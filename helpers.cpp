#include "helpers.hpp"
#include "structures.hpp"
#include "definitions.h"
#include <string>

RandomVec mix_rvecs(std::vector<RandomVec> &rvecs){
        Real sum_weights = 0;
        Vec mean_vec = {0,0};
        Real sum_var = 0;
        for(auto rvec: rvecs){
                sum_weights += rvec.weight;
                mean_vec += rvec.vec*rvec.weight;
        }
        mean_vec = mean_vec*(1/sum_weights);
        for(auto rvec: rvecs){
                Vec diff_vec = rvec.vec;
                diff_vec += mean_vec*-1;
                sum_var += diff_vec.norm_sq()*rvec.weight+rvec.var;                
        }
        sum_var /= sum_weights;
        return {mean_vec, sum_var, sum_weights};
}

void to_image(std::vector<std::vector<bool>> &img, string filename){
        std::vector<char> outstring (2*pow(img.size(),2));
        for(std::size_t i = 0; i < img.size(); i++){
                for(std::size_t j = 0; j < img.size(); j++){
                        int index = i*n+j;
                        outstring[index*2] = img[i][j] ? '1' : '0';
                        outstring[index*2+1] = (j+1 == n) ? '\n' : ' ';
                }
        }
        std::vector<char> nchars;
        string nstring = *itoa(n, nchars, 10);
        string header = "P1\n" + nstring + " " + nstring + "\n";
        FILE * f = fopen(filename, "w");
        fwrite(&header, header.length(), f);
        fwrite(outstring, 2*n*n, f);
        fclose(f);
}

Real unit_rand(){
        return (Real)rand()/(Real)rand_max();
}


Vec rand_vec(){
    return vec = {unit_rand(), unit_rand()};
}

Real distance(Vec pos1, Vec pos2){
    Real dx = pos1.x - pos2.x;
    Real dy = pos2.y - pos2.y;
    return sqrt(dx*dx+dy*dy);
}

Real periodic_dist(Vec pos1, Vec pos2, Real size){
    Real dx_inner = abs(pos1.x - pos2.x);
    Real dy_inner = abs(pos2.y - pos2.y);
    Real dx_outer = size-dx_inner;
    Real dy_outer = size-dy_inner;
    Real dx = min(dx_inner,dx_outer);
    Real dy = min(dy_inner,dy_outer);
    return sqrt(dx*dx+dy*dy);
}