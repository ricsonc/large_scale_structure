#ifndef NBODY
#define NBODY

struct NBody {
        std::vector<Body> bodies;
        UniverseArgs uargs;
        SimulationArgs sargs;
        Tree *quadtree;
        
        void init_field();
        void metric_expansion();
        void border_wrap();
        void build_qtree();
        void leapfrog();
        //what do these two do?
        Vec accel(Vec p1, Vec p2, Real mass, Real distance);
        Vec lattice_accel(Vec pos, RVec center, Real lattice_dist);
        Vec body_accel()
        std::vector<Vec> accelerations();
        void step();
        void simulate(verbose = True);
        void draw();
};

#endif