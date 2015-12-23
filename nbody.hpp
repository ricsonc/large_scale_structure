#ifndef NBODY
#define NBODY

class NBody {
    private:
        Tree quadtree;
        std::vector<std::vector<Vec>> force_field;
        void init_field();
        void metric_expansion();
        void border_wrap();
        void build_qtree();
        void leapfrog();
        void step();
        Vec accel_body_point(Body B, Vec P, Real mass);
        Vec accel_body_all(Body B);
        std::vector<Vec> accel_all_all();
    public:
        NBody(CReal density = 1E-26, //kg*m^-3
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
              string filename = "");
        std::vector<Body> bodies;
        universe_args uargs;
        simulation_args sargs;
        void simulate(verbose = True);
        void draw();
};

#endif