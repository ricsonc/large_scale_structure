#ifndef STRUCTURES
#define STRUCTURES

#include "definitions.hpp"
#include "structures.hpp"
#include <string>

struct Vec {
        Real x;
        Real y;
        //functions
        vec_plus_equal(Vec vec);
        Vec vec_x_real(Real real);
        Real norm_sq(void);
};
typedef struct Vec Vec;

struct RVec {
        Vec vec;
        Real var;
        Real weight;
};
typedef struct RVec RVec;

struct Rect {
        Vec pos0;
        Vec pos1;
        //functions
        bool contains(Vec pos);
};
typedef struct Rect Rect;

struct Body {
        int id;
        Real mass;
        Vec p;
        Vec v;
};
typedef struct Body Body;

struct Node {
        Rect rect;
        RVec center;
        Body * body;
};
typedef struct Node Node;

struct Tree {
        Node node;
        Tree * trees;
        int num_child;
};
typedef struct Tree Tree;

struct Tree_stack {
        int n;
        Tree ** trees;
        //functions
        void add(Tree * tree);
        Tree * pop(void);
};
typedef struct Tree_stack Tree_stack;

struct UArgs {
        Real size;
        Real QTR;
        Real lattice;
        Real hubble;
        Real plummer;
        Real gconst;
        int num_bodies;
};
typedef struct UArgs UArgs;

struct Dargs {
        string filename;
        int size;
};
typedef struct Dargs Dargs;

#endif