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
    void offset(Vec);
    public:
    NBody(std::string filename,
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
          const Real max_velocity);
    std::vector<Body> bodies;
    UArgs uargs;
    SArgs sargs;
    IOArgs ioargs;
    void simulate(bool verbose = true);
    void draw();
    void data_dump();
};

#endif
