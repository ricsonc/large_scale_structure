#include <stdio.h>
#include <math.h> 
#include <stdlib.h> 
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>

#include "nbody.hpp"

using namespace std;

Body * initialbodies(int n, Real size, Real mass, Real displ_ratio, Real max_vel){
        Body * bodies = new Body[n];
        const int side = sqrt(n);
        const Real radius = size/side;
        for(int i = 0; i < side; i++){
                for(int j = 0; j < side; j++){
                        //generate displacements
                        int index = i*side+j;
                        d_vector = rand_vec();
                        displacement = d_vector.M(displ_ratio);
                        vel_displacement = d_vector.M(max_vel);
                        //make body
                        Vec p = {i+displacement.x, j+displacement.y};
                        bodies[index] = {index, mass, p.M(radius), vel_displacement};
                }
        }
        return bodies;
}

NBody default_universe(Real simtime){
        //natural parameters
        CReal density = 1E-26; //kg*m^-3
        CReal size 5E+23; //m
        CReal plummer 5E+21; //m
        CReal gravity 6.67E-11; //m^3*kg^-1*s^-2
        CReal hubble 2.25E-18; //s^-1
        //simulation parameters
        const int num_bodies 4096;
        CReal QTR = 3;
        CReal lattice = 2;
        // initial conditions
        CReal displacement = 0.2;
        CReal max_velocity = 1E5; //m*s^-1
        //computed
        CReal initsize = size*exp(-hubble*simtime)
        CReal mass = pow(size, 3)*density/bodies;
        //created
        UArgs uargs = {initsize, QTR, lattice, hubble, plummer, gravity, num_bodies};
        Body * bodies = initialbodies(num_bodies, initsize, mass, displacement, max_velocity);
        return {.bodies = bodies, .uargs = uargs};
}

void default_sim(string simname, int disp_size){
        CReal simtime = 5E+17;
        CReal timestep = 5E+14;
        NBody universe = default_universe(Real simtime);
        universe.simulate(simtime, timestep, simname, disp_size);
}

int main(int argc, char **argv)
{
        srand(CSeed ? 0 : time(NULL));
        default_sim("filament_finder", 1024);
}