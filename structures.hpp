#ifndef STRUCTURES
#define STRUCTURES

#include "definitions.hpp"
#include "structures.hpp"
#include <string>

struct Vec {
        Real x;
        Real y;
        void operator +(Vec v);
        Vec operator *(Real scalar);
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
};
typedef struct Rect Rect;

struct Body {
        Vec p;
        Vec v;
};
typedef struct Body Body;

struct Node {
        Rect rect;
        RVec center;
        Body &body;
};
typedef struct Node Node;

struct Tree {
        Node node;
        std::vector<Tree> children;
};
typedef struct Tree Tree;

struct universe_args {
    Real size;
    Real hubble;
    Real plummer;
    Real gconst;
};
typedef struct UArgs UArgs;

struct simulation_args {
    Real QTR;
    Real lattice;
    int num_bodies;
    int body_mass;
    Real simtime;
    Real timestep;
    int grid_limit;
    string filename;
    int drawsize;
};
typedef struct simulation_args simulation_args

#endif