#ifndef HELPERS
#define HELPERS

#include <string>
#include "structures.hpp"

RandomVec mix_rvecs(std::vector<RandomVec> rvecs);
void to_image(std::vector<std::vector<bool>> img, int n, std::string filename);
Vec rand_vec();
Real distance(Vec pos1, Vec pos2);
Real lat_dist(Vec pos1, Vec pos2, Real lattice);

#endif