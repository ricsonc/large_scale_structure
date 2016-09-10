#ifndef STRUCTURES
#define STRUCTURES

#include <string>
#include <vector>
#include <memory>

#include "common.hpp"
#include "structures.hpp"

struct Vec {
    Real x;
    Real y;
    void operator +=(Vec v);
    Vec operator +(Vec v);
    Vec operator *(Real scalar);
    Vec operator -(Vec v);
    void operator %=(Real scalar);
    Vec operator %(Real scalar);
    bool operator ==(Vec v);
    Real norm_sq();
};
typedef struct Vec Vec;

struct RandomVec {
    Vec vec;
    Real var;
    Real weight;
};
typedef struct RVec RVec;

struct Rect {
    Vec pos0;
    Vec pos1;
    bool contains(Vec pos);
    bool operator ==(Rect r);
};
typedef struct Rect Rect;

struct Body {
    Vec p;
    Vec v;
};
typedef struct Body Body;

struct Node {
    Rect rect;
    RandomVec center;
    std::vector<std::unique_ptr<Node>> children;
};
typedef struct Node Node;

struct universe_args {
    Real size;
    Real hubble;
    Real plummer;
    Real gconst;
};
typedef struct universe_args UArgs;

struct simulation_args {
    Real QTR;
    Real body_mass;
    Real simtime;
    Real timestep;
    int grid_limit;
    Real damping;
};
typedef struct simulation_args SArgs;

struct input_output_args{
    std::string filename;
    std::size_t drawsize;
    int frame_num;
    int frequency;
};
typedef struct input_output_args IOArgs;

#endif
