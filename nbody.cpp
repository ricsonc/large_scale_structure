#include "nbody.hpp"
#include "helpers.hpp"

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
    if(bodies.size() == 1){
        root->center = {bodies[0].p,0,1};
        return root;
    }
    std::array<Rect,4> rects = split_rect(root->rect);
    std::vector<RandomVec> centers;
    for(auto rect: rects){
        std::vector<Body>quadrant;
        for(auto body: bodies){
            if(rect.contains(body.p)){
                quadrant.push_back(body);
            }
        }
        if(quadrant.size()){
            std::unique_ptr<Node> child = qtree(rect, quadrant);
            centers.push_back(child->center);
            root->children.push_back(child);
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
    Vec dp = P-B.p;
    Vec accel;
    if(dp.norm_sq() > this->sargs.grid_limit*this->uargs.size/this->force_field.size()){
        Vec indices = dp*(force_field.size()/this->uargs.size);
        accel = this->force_field[indices.x][indices.y]*pow(this->uargs.size,-2);
    } else {
        accel = pp_acceleration(dp, this->uargs.plummer);
    }
    return accel*this->sargs.body_mass*mass*this->uargs.gconst;
}

Vec NBody::accel_body_all(Body &B){
    Vec accel = {0,0};
    std::stack<Node *> DFSS;
    DFSS.push(this->quadtree.get());
    while(!DFSS.empty()){
        Node *child = DFSS.top();
        DFSS.pop();
        Real distance = periodic_dist(B.p, child->center.vec, this->uargs.size);
        if(child->center.var == 0 || distance/std::sqrt(child->center.var) > this->sargs.QTR){
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

std::vector<Body> initialbodies(int n, Real size, Real mass, Real displacement_ratio, Real max_vel){
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
             const int drawsize,
             const Real displacement,
             const Real max_velocity) 
{
    const Real initsize = size*exp(-hubble*simtime); //m
    const Real mass = pow(size, 3)*density/num_bodies; //kg
    universe_args uargs = {initsize, hubble, plummer, gravity};
    simulation_args sargs = {QTR, mass, simtime, timestep, grid_limit, filename, drawsize};
    std::vector<Body> bodies = initialbodies(num_bodies, initsize, mass, displacement, max_velocity);
    std::vector<std::vector<Vec>> forcefield = init_field(resolution, tilings);
}

void NBody::metric_expansion(){
    Real ratio = 1+this->uargs.hubble*this->sargs.timestep;
    this->uargs.size *= ratio;
    for(auto body: bodies){
        body.p = body.p*ratio;
    }
}

void NBody::border_wrap(){
    for(auto body: this->bodies){
        body.p.x = fmod(body.p.x, this->uargs.size);
        body.p.y = fmod(body.p.y, this->uargs.size);
    }
}

void NBody::build_qtree(){
    this->quadtree = qtree({{0,0},{this->uargs.size,this->uargs.size}}, this->bodies);
}

void NBody::leapfrog(){
    std::vector<Vec> accs = accel_all_all();
    int i = 0;
    for(Body &body: this->bodies){
        body.p += body.v*this->sargs.timestep;
        body.v += accs[i]*this->sargs.timestep;
        i++;
    }
}

void NBody::step(){
    build_qtree();
    leapfrog();
    border_wrap();
    metric_expansion();
    draw();
}

void NBody::simulate(bool verbose){
    for(int i = 0; i < this->sargs.simtime/this->sargs.timestep; i++){
        if(verbose){
            printf("%f\n", this->sargs.timestep*i);
        }
        step(this->sargs.timestep, {this.simname+to_string(i)});
    }
}

void NBody::draw(std::string filename){
    Real dsize = this->dargs.size;
    static bool init = true;
    static std::vector<bool> pixelarr;
    if(init){
        pixelarr = new bool[dsize][dsize];
        init = false;
    }
    memset(parr, false, dsize*dsize*sizeof(bool));
    for(auto body: this->bodies){
        int x = body.p.x/uargs.size*dsize;
        int y = body.p.x/uargs.size*dsize;
        to_image[x%dsize][y%dsize] = true;
    }
    to_image(parr, dsize, filename);
}