#ifndef NBODY
#define NBODY

#include <string>
#include <memory>

#include "common.hpp"
#include "structures.hpp"

class NBody {
    private:
        std::unique_ptr<Node> quadtree;
        std::vector<Vec> accs;
        std::vector<std::vector<Vec>> force_field;
        void metric_expansion();
        void border_wrap();
        void build_qtree();
        void leapfrog();
        Vec accel_body_point(Body &B, Vec &P, Real mass);
        Vec accel_body_all(Body &B);
        std::vector<Vec> accel_all_all();
        Real kinetic();
    public:
        NBody(std::string filename = "",
              const Real density = 1E-27, //kg*m^-3
              const Real size = 1E+25, //m
              const Real plummer = 1E+23, //m
              const Real gravity = 6.67E-11, //m^3*kg^-1*s^-2
              const Real hubble = 2.25E-18, //s^-1
              const Real damping = 5E-18,
              const Real simtime = 5E+17, //s
              const Real timestep = 1E+14, //s
              const Real QTR = 3,
              const int resolution = 2048,
              const int tilings = 31,
              const int grid_limit = 16,
              const int num_bodies = 65536,
              const std::size_t drawsize = 2048,
              const int draw_freq = 1,
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