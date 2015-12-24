#include "nbody.hpp"
#include "helpers.hpp"
#include "structures.hpp"
#include "definitions.hpp"

int main()
{
    srand(CSeed ? 0 : time(NULL));
    NBody universe = NBody(filename = "filament_simulation");
    universe.simulate();
}

