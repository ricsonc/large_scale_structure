#include <cmath>
#include <cstdbool>
#include <vector>
#include <fstream>

#include "helpers.hpp"
#include "common.hpp"

RandomVec mix_rvecs(std::vector<RandomVec> &rvecs){
    Real sum_weights = 0;
    Vec mean_vec = {0,0};
    Real sum_var = 0;
    for(RandomVec &rvec: rvecs){
        sum_weights += rvec.weight;
        mean_vec += rvec.vec*rvec.weight;
    }
    mean_vec = mean_vec*(1/sum_weights);
    for(RandomVec &rvec: rvecs){
        Vec diff_vec = rvec.vec;
        diff_vec += mean_vec*-1;
        sum_var += diff_vec.norm_sq()*rvec.weight+rvec.var;                
    }
    sum_var /= sum_weights;
    return {mean_vec, sum_var, sum_weights};
}

void to_image(std::vector<std::vector<bool>> &img, std::string filename){
        size_t n = img.size();
        std::vector<char> outstring (2*n*n);
        for(std::size_t i = 0; i < n; i++){
                for(std::size_t j = 0; j < n; j++){
                        int index = (n-1-i)*n+j;
                        outstring[index*2] = img[i][j] ? '1' : '0';
                        outstring[index*2+1] = (j+1 == n) ? '\n' : ' ';
                }
        }
        std::string nstring = std::to_string(n);
        std::string header = "P1\n" + nstring + " " + nstring + "\n";
        std::ofstream f;
        f.open(filename.c_str());
        f.write(header.c_str(), sizeof(char)*header.size());
        f.write(outstring.data(), sizeof(char)*outstring.size());
        f.close();
}

Real unit_rand(){
    return (Real)rand()/RAND_MAX;
}


Vec rand_vec(){
    return {unit_rand()-0.5, unit_rand()-0.5};
}

Real distance(Vec pos1, Vec pos2){
    Real dx = pos1.x - pos2.x;
    Real dy = pos2.y - pos2.y;
    return std::sqrt(dx*dx+dy*dy);
}

Real periodic_dist(Vec pos1, Vec pos2, Real size){
    Real dx_inner = abs(pos1.x - pos2.x);
    Real dy_inner = abs(pos2.y - pos2.y);
    Real dx_outer = size-dx_inner;
    Real dy_outer = size-dy_inner;
    Real dx = std::min(dx_inner,dx_outer);
    Real dy = std::min(dy_inner,dy_outer);
    return std::sqrt(dx*dx+dy*dy);
}