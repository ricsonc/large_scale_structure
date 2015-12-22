#include "helpers.hpp"
#include "structures.h"
#include "definitions.hpp"
#include <string>

RVec mix_rvecs(RVec * rvecs, int n){
        //initialize variables
        Real sum_weights = 0;
        Vec mean_vec = {0,0};
        Real sum_var = 0;
        //compute the mixture
        for(int i = 0; i < n; i++){
                sum_weights += rvecs[i].weight;
                mean_vec.P(rvecs[i].vec.M(rvecs[i].weight));
        }
        mean_vec = mean_vec.M(1/sum_weights);
        for(int i = 0; i < n; i++){
                Vec diff_vec = rvecs[i].vec;
                diff_vec.P(mean_vec.M(-1));
                sum_var += diff_vec.norm_sq()*rvecs[i].weight+rvecs[i].var;                
        }
        sum_var /= sum_weights;
        return {mean_vec, sum_var, sum_weights};
}

void to_image(bool ** img, int n, string filename){
        //convert bool array to char representation
        char * outstring = new char[2*n*n];
        for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                        int index = i*n+j;
                        outstring[index*2] = img[i][j] ? '1' : '0';
                        outstring[index*2+1] = (j+1 == n) ? '\n' : ' ';
                }
        }
        //write header and array to file
        char * nchars;
        string nstring = *itoa(n, nchars, 10);
        string header = "P1\n" + nstring + " " + nstring + "\n";
        FILE * f = fopen(filename, "w");
        fwrite(&header, header.length(), f);
        fwrite(outstring, 2*n*n, f);
        fclose(f);
}

Real unit_rand(void){
        return (Real)rand()/(Real)rand_max();
}


Vec rand_vec(void){
        return vec = {unit_rand(), unit_rand()};
}

Real distance(Vec pos1, Vec pos2){
        Real dx = pos1.x - pos2.x;
        Real dy = pos2.y - pos2.y;
        return sqrt(dx*dx+dy*dy);
}

Real lat_dist(Vec pos1, Vec pos2, Real lattice){
        Real dx_inner = abs(pos1.x - pos2.x);
        Real dy_inner = abs(pos2.y - pos2.y);
        Real dx_outer = lattice-dx_inner;
        Real dy_outer = lattice-dy_inner;
        Real dx = min(dx_inner,dx_outer);
        Real dy = min(dy_inner,dy_outer);
        return sqrt(dx*dx+dy*dy);
}