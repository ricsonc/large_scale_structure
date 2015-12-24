#include "structures.hpp"

Vec Vec::operator +=(Vec v){
    this->x += v.x;
    this->y += v.y;
}

Vec Vec::operator *(Real scalar){
    return {this->x*scalar, this->y*scalar};
}

Vec Vec::operator -(Vec v){
    return {v.x-this->x, v.y-this->y};
}

Real Vec::norm_sq(){
    return pow(this->x,2)+pow(this->y,2);
}

bool Rect::contains(Vec pos){
    return (pos.x >= this->pos0.x &&
            pos.y >= this->pos0.y &&
            pos.x < this->pos1.x &&
            pos.y < this->pos1.y);
}
