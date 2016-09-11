#include <ctime>
#include <omp.h>
#include <getopt.h>

#include <assert.h>

#include "nbody.hpp"
#include "helpers.hpp"
#include "structures.hpp"

#define CONST_SEED true

int main(int argc, char **argv)
{
    int c;

    std::string filename = ".";
    Real density = 1E-26; //kg*m^-3
    Real size = 2E+24; //m
    Real plummer = 2E+21; //m
    Real gravity = 6.67E-11; //m^3*kg^-1*s^-2
    Real hubble = 2.25E-18; //s^-1
    Real damping = 7E-18;
    Real simtime = 5E+17; //s
    Real timestep = 1E+14; //s
    Real QTR = 3;
    int resolution = 1<<12;
    int tilings = (1<<5)-1;
    int grid_limit = 1<<6;
    int num_bodies = 1<<20;
    std::size_t drawsize = 1<<11;
    int draw_freq = 1;
    Real displacement = 0.2;
    Real max_velocity = 5E+4; //m*s^-1
    int verbosity = 1;

    while(1){
        static struct option long_options [] = {
            {"filename", required_argument, 0, 1},
            {"density", required_argument, 0, 2},
            {"size", required_argument, 0, 3},
            {"plummer", required_argument, 0, 4},
            {"gravity", required_argument, 0, 5},
            {"hubble", required_argument, 0, 6},
            {"damping", required_argument, 0, 7},
            {"simtime", required_argument, 0, 8},
            {"timestep", required_argument, 0, 9},
            {"quadtree_ratio", required_argument, 0, 10},
            {"mesh_resolution", required_argument, 0, 11},
            {"periodic_tilings", required_argument, 0, 12},
            {"mesh_limits", required_argument, 0, 13},
            {"bodies", required_argument, 0, 14},
            {"draw_size", required_argument, 0, 15},
            {"draw_freq", required_argument, 0, 16},
            {"displacement", required_argument, 0, 17},
            {"init_velocity", required_argument, 0, 18},
            {"verbose", required_argument, 0, 19}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1){
            break;
        }
        switch(c)
        {
        case 1:
            filename = optarg;
            break;
        case 2:
            size = atof(optarg);
            break;
        case 3:
            plummer = atof(optarg);
            break;
        case 4:
            gravity = atof(optarg);
            break;
        case 5:
            hubble = atof(optarg);
            break;
        case 6:
            damping = atof(optarg);
            break;
        case 7:
            assert(optarg);
            simtime = atof(optarg);
            break;
        case 8:
            timestep = atof(optarg);
            break;
        case 9:
            QTR = atof(optarg);
            break;
        case 10:
            resolution = atof(optarg);
            break;
        case 11:
            tilings = atoi(optarg);
            break;
        case 12:
            grid_limit = atoi(optarg);
            break;
        case 13:
            num_bodies = atoi(optarg);
            break;
        case 14:
            drawsize = atoi(optarg);
            break;
        case 15:
            draw_freq = atoi(optarg);
            break;
        case 16:
            displacement = atof(optarg);
            break;
        case 17:
            max_velocity = atof(optarg);
            break;
        case 18:
            verbosity = atoi(optarg);
            break;
        default:
            abort();
        }
    }
    if (verbosity){
        printf("starting simulation\n");
    }
    srand(CONST_SEED ? 0 : time(NULL));
    NBody universe = NBody(filename,
                           density,
                           size,
                           plummer,
                           gravity,
                           hubble,
                           damping,
                           simtime,
                           timestep,
                           QTR,
                           resolution,
                           tilings,
                           grid_limit,
                           num_bodies,
                           drawsize,
                           draw_freq,
                           displacement,
                           max_velocity);
    universe.simulate();
    return 0;
}
