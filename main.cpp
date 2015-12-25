#include <ctime>

#include "nbody.hpp"
#include "helpers.hpp"
#include "structures.hpp"

#define CSeed true

int main()
{
    srand(CSeed ? 0 : time(NULL));
    NBody universe = NBody("filament_simulation");
    universe.simulate();
}