#include "structures.hpp"
#include "definitions.hpp"

////

Vec::vec_plus_equal(Vec vec){
        x += vec.x;
        y += vec.y;
}

Vec Vec::vec_x_real(Real real){
        return {x*real, y*real};
}

Real Vec::norm_sq(void){
        return x*x+y*y;
}

////

void Tree_stack::add(Tree * tree){
        trees[n] = tree;
        n++;
}

Tree * Tree_stack::pop(void){
        n--;
        return trees[n];
}

////

bool Rect::contains(Vec pos){
        return (pos.x >= pos0.x &&
                pos.y >= pos0.y &&
                pos.x < pos1.x &&
                pos.y < pos1.y);
}
