#include "nbody.hpp"
#include "helpers.hpp"
#include "structures.hpp"
#include "definitions.hpp"

std::vector<Body> initialbodies(int n, Real size, Real mass, Real displacement_ratio, Real max_vel){
    std::vector<Body> bodies (n);
    const int side = sqrt(n);
    const Real radius = size/side;
    for(int i = 0; i < side; i++){
        for(int j = 0; j < side; j++){
            d_vector = rand_vec();
            displacement = d_vector*displacement_ratio;
            vel_displacement = d_vector*max_vel;
            Vec p = {i+displacement.x, j+displacement.y};
            bodies[i*side+j] = {index, mass, p*radius, vel_displacement};
        }
    }
    return bodies;
}

NBody default_universe(string filename){
        
    CReal density = 1E-26; //kg*m^-3
    CReal size 5E+23; //m
    CReal plummer 5E+21; //m
    CReal gravity 6.67E-11; //m^3*kg^-1*s^-2
    CReal hubble 2.25E-18; //s^-1
    
    CReal simtime = 5E+17; //s
    CReal timestep = 5E+14; //s
    
    const int num_bodies 4096;
    CReal QTR = 3;
    CReal lattice = 2;
    const int drawsize = 1024;
    
    CReal displacement = 0.2;
    CReal max_velocity = 1E5; //m*s^-1
    
    CReal initsize = size*exp(-hubble*simtime) //m
    CReal mass = pow(size, 3)*density/bodies; //kg
    
    UniverseArgs uargs = {initsize, hubble, plummer, gravity};
    SimulationArgs sargs = {QTR, lattice, mass, simtime, timestep, filename, drawsize};
    std::vector<Body> bodies = initialbodies(num_bodies, initsize, mass, displacement, max_velocity);
    return {.bodies = bodies, .uargs = uargs, .sargs = sargs};
}

int main(int argc, char **argv)
{
    srand(CSeed ? 0 : time(NULL));
    NBody universe = default_universe("filament_simulation");
    universe.simulate();
}

