#ifndef NBODY
#define NBODY

#include "structures.hpp"
#include <string>
#include "definitions.hpp"

class NBody {
    private:
        std::unique_ptr<Node> quadtree;
        std::vector<std::vector<Vec>> force_field;
        void metric_expansion();
        void border_wrap();
        void build_qtree();
        void leapfrog();
        Vec accel_body_point(Body B, Vec P, Real mass);
        Vec accel_body_all(Body &B);
        std::vector<Vec> accel_all_all();
    public:
        NBody(std::string filename = "",
              const Real density = 1E-26, //kg*m^-3
              const Real size = 5E+23, //m
              const Real plummer = 5E+21, //m
              const Real gravity = 6.67E-11, //m^3*kg^-1*s^-2
              const Real hubble = 2.25E-18, //s^-1
              const Real simtime = 5E+17, //s
              const Real timestep = 5E+15, //s
              const Real QTR = 3,
              const int resolution = 1024,
              const int tilings = 15,
              const int grid_limit = 8,
              const int num_bodies = 4096,
              const std::size_t drawsize = 1024,
              const Real displacement = 0.2,
              const Real max_velocity = 1E+5); //m*s^-1
        std::vector<Body> bodies;
        UArgs uargs;
        SArgs sargs;
        IOArgs ioargs;
        void simulate(bool verbose = true);
        void draw();
};

#endif