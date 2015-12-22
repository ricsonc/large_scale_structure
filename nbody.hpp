#ifndef NBODY
#define NBODY

class NBody {
    public:
        std::vector<Body> bodies;
        UniverseArgs uargs;
        SimulationArgs sargs;
        void simulate(verbose = True);
        void draw();
    private:
        Tree *quadtree;
        void init_field();
        void metric_expansion();
        void border_wrap();
        void build_qtree();
        void leapfrog();
        void step();
        Vec accel_body(Body B);
        std::vector<Vec> accelerations();
};

#endif