#include "helpers.hpp"
#include "structures.hpp"

std::array<Rect, 4> split_rect(Rect rect){
        Real minx = rect.pos0.x;
        Real miny = rect.pos0.y;
        Real maxx = rect.pos1.x;
        Real maxy = rect.pos1.y;
        Real midx = (minx+maxx)/2;
        Real midy = (miny+maxy)/2;

        return {{rect.pos0,{midx,midy}},
                {{midx,miny},{maxx,midy}},
                {{minx,midy},{midx,maxy}},
                {{midx,midy},rect.pos1}};
}

void qtree(Node root, std::vector<Body> &bodies, std::multiset<Node, Node> &QT){
    if(bodies.size() == 1){
        root.body = bodies[0];
        root.center = {root.body.p,0,1};
        return;
    }
    std::array<Rect,4> rects = split_rect(root.rect);
    std::vector<Rvec> centers;
    for(auto rect: rects){
        std::vector<Body>quadrant;
        for(auto body: bodies){
            if rect.contains(body){
                quadrant.push_back(body);
            }
        }
        if(quadrant.size()){
            Node child = {.rect = rect};
            build_qtree(child, quadrant, QT);
            centers.push_back(child.center);
            QT.insert(pair<Node, Node>(root, child));
        }
    }
    root.center = mix_rvecs(centers);
    return;
}

Vec pp_acceleration(dp, mass = 1, G = 1, plummer = 0){
    Real distancesq = dp.norm_sq()
    Real distance = std::sqrt(distancesq);
    Real mag_accel = G*mass/(distancesq+plummer*plummer);
    return dp*(mag_accel/distance);
}

std::vector<std::vector<Vec>> init_field(int resolution, int tiling){
    <std::vector<std::vector<Vec>> field (resolution, std::vector<Vec> (resolution, {0,0}));
    for(int i = 0; i < resolution; i++){
        for(int j = 0; j < resolution; j++){
            for(int n = -tiling, n <= tiling, n++){
                for(int m = -tiling, m <= tiling, m++){
                    Vec p = {(i+0.5)/resolution+n,
                             (j+0.5)/resolution+m};
                    field[i][j] += pp_acceleration(p);
                }
            }
        }
    }
    return field;
}

Vec NBody::accel_body_body(Body B0, Body B1){
    Vec dp = B1.p-B0.p;
    if(dp.norm_sq() > this->grid_limit*this->uargs.size/this->force_field.size()){
        Vec indices = dp*(force_field.size()/this->uargs.size);
        return this->force_field[indices.x][indices.y];
    } else {
        return pp_acceleration(dp, this->sargs.mass, this->uargs.gconst, this->uargs.plummer);
    }
}

//rewrite
Vec NBody::accel_body_all(Body B){
    Vec accel = {0,0};
    //rewrite with std::stack instead of tree_stack
    std::stack<Node> BFSstack;
    BFSstack.push()
    static Tree_stack BFSS = {0, Tree[uargs.nbody]};
    BFSS.add(quadtree);
    while(BFSS.n){
        tree = BFSS.pop();
        Real lattice_dist = lat_dist(bodies[i].p, tree.center.var, uargs.size);
        if(tree.node.body != NULL && tree.node.body->id == id){
            continue;
        } else if (tree.node.center.var == 0 || lattice_dist/sqrt(tree.center.var) > QTR){
            accel.P(lattice_accel(bodies[i].p, tree.node.center, lattice_dist));
        } else {
            for(int i = 0; i < tree.num_child; i++){
                BFSS.add(tree.trees[i]);
            }
        }
    }
    return accel;
}

std::vector<Vec> NBody::accel_all_all(){ 
    static std::vector<Vec> accs (this->bodies.size());
    int i = 0;
    for(auto body: this->bodies){
        accs[i] = body_accel(body);
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
            d_vector = rand_vec();
            displacement = d_vector*displacement_ratio;
            vel_displacement = d_vector*max_vel;
            Vec p = {i+displacement.x, j+displacement.y};
            bodies[i*side+j] = {index, mass, p*radius, vel_displacement};
        }
    }
    return bodies;
}

NBody::NBody(CReal density = 1E-26, //kg*m^-3
             CReal size = 5E+23, //m
             CReal plummer = 5E+21, //m
             CReal gravity = 6.67E-11, //m^3*kg^-1*s^-2
             CReal hubble = 2.25E-18, //s^-1
             CReal simtime = 5E+17, //s
             CReal timestep = 5E+14, //s
             CReal QTR = 3,
             const int resolution = 1024,
             const int tilings = 15,
             const int grid_limit = 8,
             const int num_bodies = 4096,
             const int drawsize = 1024,
             CReal displacement = 0.2,
             CReal max_velocity = 1E+5, //m*s^-1
             string filename = "")
{
    CReal initsize = size*exp(-hubble*simtime) //m
    CReal mass = pow(size, 3)*density/bodies; //kg
    universe_args uargs = {initsize, hubble, plummer, gravity};
    simulation_args sargs = {QTR, lattice, mass, simtime, timestep, gridlimit, filename, drawsize};
    std::vector<Body> bodies = initialbodies(num_bodies, initsize, mass, displacement, max_velocity);
    std::vector<std::vector<Vec>> forcefield = init_field(resolution, tiling);
    return {.bodies = bodies, .uargs = uargs, .sargs = sargs, .force_field = forcefield};
}

void NBody::metric_expansion(){
    Real ratio = 1+hubble*this->simulation_args.timestep;
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
    std::multimap<Node, Node> QT;
    return qtree({{0,0},{this->uargs.size,this->uargs.size}}, this->bodies, QT);
}

void NBody::leapfrog(){
    std::vector<Vec> accs = accelerations();
    int i = 0;
    for(auto body: this->bodies){
        body.p += body.v*this->simulation_args.timestep;
        body.v += accs[i]*this->simulation_args.timestep;
        i++;
    }
}

void NBody::step(){
    build_qtree();
    leapfrog();
    border_wrap();
    metric_expansion();
    if(this->dargs != NULL){
        draw(this->dargs);
    }
}

void NBody::simulate(verbose = True){
    for(int i = 0; i < this->simtime/this->timestep; i++){
        if(verbose){
            printf("%f\n", dt*i);
        }
        step(this->timestep, {simname+to_string(i)});
    }
}

void NBody::draw(string filename){
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