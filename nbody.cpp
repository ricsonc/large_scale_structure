#include <stack>
#include <memory>
#include <cmath>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cstdio>
#include <omp.h>
#include <ctime>
#include <cassert>
#include <fstream>
#include <sstream>
#include <chrono>

#include "nbody.hpp"
#include "structures.hpp"
#include "helpers.hpp"
#include "common.hpp"

std::array<Rect, 4> split_rect(Rect rect){
        Real minx = rect.pos0.x;
        Real miny = rect.pos0.y;
        Real maxx = rect.pos1.x;
        Real maxy = rect.pos1.y;
        Real midx = (minx+maxx)/2;
        Real midy = (miny+maxy)/2;
        std::array<Rect, 4> rects;
        rects[0] = {rect.pos0,{midx,midy}};
        rects[1] = {{midx,miny},{maxx,midy}};
        rects[2] = {{minx,midy},{midx,maxy}};
        rects[3] = {{midx,midy},rect.pos1};
        return rects;
}

std::unique_ptr<Node> qtree(Rect rect, std::vector<Body*> &bptrs){
    std::unique_ptr<Node> root (new Node());
    root->rect = rect;
    if(bptrs.size() == 1){
        root->center = {bptrs[0]->p,0,1};
        return root;
    }
    std::array<Rect,4> rects = split_rect(root->rect);
    std::vector<RandomVec> centers;
    std::vector<Rect> rect_map(bptrs.size());
    for(std::size_t i = 0; i < rect_map.size(); i++){
        for(Rect &subrect: rects){
            if(subrect.contains(bptrs[i]->p)){
                rect_map[i] = subrect;
                break;
            }
        }
    }
    for(Rect &subrect: rects){
        std::vector<Body*> quadrant;
        for(std::size_t i = 0; i < bptrs.size(); i++){
            if(rect_map[i] == subrect){
                quadrant.push_back(bptrs[i]);
            }
        }
        if(quadrant.size()){
            std::unique_ptr<Node> child = qtree(subrect, quadrant);
            centers.push_back(child->center);
            root->children.push_back(std::move(child));
        }
    }
    root->center = mix_rvecs(centers);
    return root;
}

Vec pp_acceleration(Vec dp, Real plummer = 0){
    Real distancesq = dp.norm_sq();
    Real distance = std::sqrt(distancesq);
    Real mag_accel = 1/(distancesq+plummer*plummer);
    return dp*(mag_accel/distance);
}

std::vector<std::vector<Vec>> init_field(int resolution, int tiling, Real plummer = 0){
    std::vector<std::vector<Vec>> field (resolution,
                                         std::vector<Vec> (resolution, {0,0}));
    #pragma omp parallel for
    for(int i = 0; i < resolution; i++){
        for(int j = 0; j < resolution; j++){
            for(int n = -tiling; n < tiling; n++){
                for(int m = -tiling; m < tiling; m++){
                    Vec p = {n+(i+0.5)/resolution,
                             m+(j+0.5)/resolution};
                    if(p.norm_sq() < tiling){
                        field[i][j] += pp_acceleration(p, plummer);
                    }
                }
            }
        }
    }
    return field;
}

Vec NBody::accel_body_point(Body &B, Vec &P, Real mass){
    Vec dp = (P-B.p)%uargs.size;
    Vec accel;
    std::size_t n = force_field.size();
    if(dp.norm_sq() > sargs.grid_limit*uargs.size/n){
        Vec indices = dp*(n/uargs.size);
        accel = (force_field[(std::size_t)indices.x][(std::size_t)indices.y] *
                 pow(uargs.size,-2));
    } else {
        Real xu = dp.x-uargs.size;
        Real yu = dp.y-uargs.size;
        Vec close = {dp.x < -xu ? dp.x : xu, dp.y < -yu ? dp.y : yu};
        accel = pp_acceleration(close, uargs.plummer);
    }
    return accel*sargs.body_mass*mass*uargs.gconst;
}

Vec NBody::accel_body_all(Body &B){
    Vec accel = {0,0};
    std::stack<Node *> DFSS;
    DFSS.push(quadtree.get());
    while(!DFSS.empty()){
        Node *child = DFSS.top();
        DFSS.pop();
        Real distance = periodic_dist(B.p, child->center.vec, uargs.size);
        if(child->center.var == 0 ||
           distance/std::sqrt(child->center.var) > sargs.QTR){
            if(distance > 0){
                Vec a = accel_body_point(B, child->center.vec, child->center.weight);
                accel += accel_body_point(B, child->center.vec,
                                          child->center.weight);
            }
        } else {
            for(std::unique_ptr<Node> &neighbor: child->children){
                DFSS.push(neighbor.get());
            }
        }
    }
    return accel;
}

std::vector<Vec> NBody::accel_all_all(){
    size_t n = bodies.size();
    static std::vector<Vec> accs (n);
    #pragma omp parallel for
    for(std::size_t i = 0; i < n; i++){
        accs[i] = accel_body_all(bodies[i]);
    }
    return accs;
}

std::vector<Body> initialbodies(int n, Real size,
                                Real displacement_ratio, Real max_vel){
    std::vector<Body> bodies (n);
    const int side = sqrt(n);
    const Real radius = size/side;
    for(int i = 0; i < side; i++){
        for(int j = 0; j < side; j++){
            Vec d_vector = rand_vec();
            Vec displacement = d_vector*displacement_ratio;
            Vec vel_displacement = d_vector*max_vel;
            Vec p = {i+displacement.x+.5, j+displacement.y+.5};
            bodies[i*side+j] = {p*radius, vel_displacement};
        }
    }
    return bodies;
}

NBody::NBody(std::string filename,
             const Real density,
             const Real size,
             const Real plummer,
             const Real gravity,
             const Real hubble,
             const Real damping,
             const Real simtime,
             const Real timestep,
             const Real QTR,
             const int resolution,
             const int tilings,
             const int grid_limit,
             const int num_bodies,
             const std::size_t drawsize,
             const int draw_freq,
             const Real displacement,
             const Real max_velocity)
{
    const Real initsize = size*exp(-hubble*simtime); //m
    const Real mass = pow(size, 3)*density/num_bodies; //kg
    this->uargs = {initsize, hubble, plummer, gravity};
    this->sargs = {QTR, mass, simtime, timestep, grid_limit, damping};
    this->ioargs = {filename, drawsize, 0, draw_freq};
    this->bodies = initialbodies(num_bodies, initsize,
                                 displacement, max_velocity);
    this->force_field = init_field(resolution, tilings, plummer/uargs.size);
    border_wrap();
    build_qtree();
    this->accs = accel_all_all();
}

void NBody::metric_expansion(){
    Real ratio = 1+uargs.hubble*sargs.timestep;
    uargs.size *= ratio;
    for(Body &body: bodies){
        body.p = body.p*ratio;
    }
}

void NBody::border_wrap(){
    for(Body &body: bodies){
        body.p %= uargs.size;
    }
}

void NBody::build_qtree(){
    std::vector<Body*> bptrs;
    for(Body &body: bodies){
        bptrs.push_back(&body);
    }
    quadtree = qtree({{0,0},{uargs.size,uargs.size}}, bptrs);
}

void NBody::leapfrog(){
    std::vector<Vec> new_accs = accel_all_all();
    int j = 0;
    for(Body &body: bodies){
        body.v += (accs[j]+new_accs[j])*sargs.timestep*-.5;
        body.v = body.v * (1 - sargs.damping*sargs.timestep);
        j++;
    }
    accs = new_accs;
    int i = 0;
    for(Body &body: bodies){
        body.p += body.v*sargs.timestep+accs[i]*pow(sargs.timestep,2)*.5;
        i++;
    }
}

void NBody::simulate(bool verbose){
    for(int i = 0; i < sargs.simtime/sargs.timestep; i++){
        if(verbose){
            std::printf("iteration: %d\n", ioargs.frame_num);
        }
        auto t0 = std::chrono::system_clock::now();
        build_qtree();
        auto t1 = std::chrono::system_clock::now();
        leapfrog();
        auto t2 = std::chrono::system_clock::now();
        border_wrap();
        metric_expansion();
        auto t3 = std::chrono::system_clock::now();
        if(i%ioargs.frequency == 0){
            draw();
            data_dump();
        }
        auto t4 = std::chrono::system_clock::now();
        if(verbose){
            std::printf("\t quadtree construction  : %ld\n", (t1-t0).count()/1000000);
            std::printf("\t barnes hut algorithm   : %ld\n", (t2-t1).count()/1000000);
            std::printf("\t metric expansion       : %ld\n", (t3-t2).count()/1000000);
            std::printf("\t writing to file        : %ld\n", (t4-t3).count()/1000000);
        }
        ioargs.frame_num++;
    }
}

void NBody::draw(){
    std::size_t dsize = ioargs.drawsize;
    if(dsize == 0){
        return;
    }
    std::vector<std::vector<bool>> pixelarr (dsize,
                                             std::vector<bool> (dsize, true));
    for(Body &body: bodies){
        size_t x = body.p.x/uargs.size*dsize;
        size_t y = body.p.y/uargs.size*dsize;
        pixelarr[x][y] = false;
    }
    to_image(pixelarr,
             ioargs.filename+"/"+std::to_string(ioargs.frame_num)+".ppm");
}

void NBody::data_dump(){
    std::string filename = (ioargs.filename+"/"+
                            std::to_string(ioargs.frame_num)+".dump");
    std::FILE *f = std::fopen(filename.c_str(), "w");
    for(Body &body: bodies){
        std::fprintf(f, "%.12e %.12e %.12e %.12e\n",
                     body.p.x, body.p.y, body.v.x, body.v.y);
    }
    std::fclose(f);
}
