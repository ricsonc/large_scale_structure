#include <stack>
#include <memory>
#include <cmath>

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

std::unique_ptr<Node> qtree(Rect rect, std::vector<Body> &bodies){
    std::unique_ptr<Node> root (new Node());
    root->rect = rect;
    if(bodies.size() == 1){
        root->center = {bodies[0].p,0,1};
        return root;
    }
    std::array<Rect,4> rects = split_rect(root->rect);
    std::vector<RandomVec> centers;
    for(auto &rect: rects){
        std::vector<Body>quadrant;
        for(Body &body: bodies){
            if(rect.contains(body.p)){
                quadrant.push_back(body);
            }
        }
        if(quadrant.size()){
            std::unique_ptr<Node> child = qtree(rect, quadrant);
            centers.push_back(child->center);
            root->children.push_back(std::move(child));
        }
    }
    root->center = mix_rvecs(centers);
    return root; 
}

Vec pp_acceleration(Vec dp, Real plummer = 0){
    Real distancesq = dp.norm_sq();
    if(distancesq == 0){
        return {0,0};
    }
    Real distance = std::sqrt(distancesq);
    Real mag_accel = 1/(distancesq+plummer*plummer);
    return dp*(mag_accel/distance);
}

std::vector<std::vector<Vec>> init_field(int resolution, int tiling){
    std::vector<std::vector<Vec>> field (resolution, std::vector<Vec> (resolution, {0,0}));
    for(int i = 0; i < resolution; i++){
        for(int j = 0; j < resolution; j++){
            for(int n = -tiling; n <= tiling; n++){
                for(int m = -tiling; m <= tiling; m++){
                    Vec p = {(i+0.5)/resolution+n,
                             (j+0.5)/resolution+m};
                    field[i][j] += pp_acceleration(p);
                }
            }
        }
    }
    return field;
}

Vec NBody::accel_body_point(Body B, Vec P, Real mass){
    Vec dp = (P-B.p)%uargs.size;
    Vec accel;
    size_t n = this->force_field.size();
    if(dp.norm_sq() > sargs.grid_limit*uargs.size/n){
        Vec indices = dp*(force_field.size()/uargs.size);
        accel = this->force_field[(int)indices.x][(int)indices.y]*pow(uargs.size,-2);
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
        if(child->center.var == 0 || distance/std::sqrt(child->center.var) > sargs.QTR){
            accel += accel_body_point(B, child->center.vec, child->center.weight);
        } else {
            for(std::unique_ptr<Node> &neighbor: child->children){
                DFSS.push(neighbor.get());
            }
        }
    }
    return accel;
}

std::vector<Vec> NBody::accel_all_all(){ 
    static std::vector<Vec> accs (this->bodies.size());
    int i = 0;
    for(Body &body: this->bodies){
        accs[i] = accel_body_all(body);
        i++;
    }
    return accs;
}

std::vector<Body> initialbodies(int n, Real size, Real displacement_ratio, Real max_vel){
    std::vector<Body> bodies (n);
    const int side = sqrt(n);
    const Real radius = size/side;
    for(int i = 0; i < side; i++){
        for(int j = 0; j < side; j++){
            Vec d_vector = rand_vec();
            Vec displacement = d_vector*displacement_ratio;
            Vec vel_displacement = d_vector*max_vel;
            Vec p = {i+displacement.x, j+displacement.y};
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
             const Real simtime,
             const Real timestep,
             const Real QTR,
             const int resolution,
             const int tilings,
             const int grid_limit,
             const int num_bodies,
             const std::size_t drawsize,
             const Real displacement,
             const Real max_velocity) 
{
    const Real initsize = size*exp(-hubble*simtime); //m
    const Real mass = pow(size, 3)*density/num_bodies; //kg
    this->uargs = {initsize, hubble, plummer, gravity};
    this->sargs = {QTR, mass, simtime, timestep, grid_limit};
    this->ioargs = {filename, drawsize, 0};
    this->bodies = initialbodies(num_bodies, initsize, displacement, max_velocity);
    this->force_field = init_field(resolution, tilings);
}

void NBody::metric_expansion(){
    Real ratio = 1+uargs.hubble*sargs.timestep;
    uargs.size *= ratio;
    for(Body &body: bodies){
        body.p = body.p*ratio;
    }
}

void NBody::border_wrap(){
    for(Body &body: this->bodies){
        body.p %= uargs.size;
    }
}

void NBody::build_qtree(){
    quadtree = qtree({{0,0},{uargs.size,uargs.size}}, this->bodies);
}

void NBody::leapfrog(){
    std::vector<Vec> accs = accel_all_all();
    int i = 0;
    for(Body &body: this->bodies){
        body.p += body.v*sargs.timestep;
        body.v += accs[i]*sargs.timestep;
        i++;
    }
}

void NBody::simulate(bool verbose){
    for(int i = 0; i < sargs.simtime/sargs.timestep; i++){
        if(verbose){
            printf("frame: %d, time: %.3e\n", ioargs.frame_num, sargs.timestep*i);
        }
        build_qtree();
        leapfrog();
        border_wrap();
        metric_expansion();
        draw();
        ioargs.frame_num++;
    }
}

void NBody::draw(){
    std::size_t dsize = ioargs.drawsize;
    if(dsize == 0){
        return;
    }
    std::vector<std::vector<bool>> pixelarr (dsize, std::vector<bool> (dsize, true));
    for(Body &body: this->bodies){
        size_t x = body.p.x/uargs.size*dsize;
        size_t y = body.p.y/uargs.size*dsize;
        pixelarr[x][y] = false;
    }
    to_image(pixelarr, ioargs.filename+"/"+std::to_string(ioargs.frame_num)+".ppm");
}