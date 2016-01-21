#ifndef HELPERS
#define HELPERS

#include <cstdlib> 
#include <vector>
#include <string>
#include <fstream>

#include "structures.hpp"
#include "common.hpp"

RandomVec mix_rvecs(std::vector<RandomVec> &rvecs);
void to_image(std::vector<std::vector<bool>> &img, std::string filename);
Vec rand_vec();
Real distance(Vec pos1, Vec pos2);
Real periodic_dist(Vec pos1, Vec pos2, Real size);

#endif