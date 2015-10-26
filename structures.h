////

struct Vec {
        Real x;
        Real y;
};
typedef struct Vec Vec;

Vec::vec_plus_equal(Vec vec);
Vec Vec::vec_x_real(Real real);
Real Vec::norm_sq(void);

////

struct Tree_stack {
        int n;
        Tree [] trees;
};
typedef struct Tree_stack Tree_stack

void Tree_stack::add(Tree * tree);
Tree * Tree_stack::pop(void);

////

struct RVec {
        Vec vec;
        REAL var;
        REAL weight;
};
typedef struct RVec RVec

////

struct UArgs = {
        Real size;
        Real QTR;
        Real lattice;
        Real hubble;
        Real plummer;
        Real gconst;
        int num_bodies;
};
typedef struct UArgs UArgs;

////

struct Dargs {
        string filename;
        int size;
};
typedef struct Dargs Dargs;

////

struct Rect {
        Vec pos0;
        Vec pos1;
};
typedef struct Rect Rect;

bool Rect::contains(Vec pos);

////

struct Node {
        Rect rect;
        RVec center;
        Body * body;
};
typedef struct Node Node;

////

struct Tree {
        Node node;
        Tree [] trees;
        int num_child;
};
typedef struct Tree Tree;

////

struct Body {
        int id;
        Real mass;
        Vec p;
        Vec v;
};
typedef struct Body Body;
