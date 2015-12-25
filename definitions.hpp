#ifndef DEFINITIONS
#define DEFINITIONS

#define CSeed true
#define _GNU_SOURCE

typedef double Real;

#include <cstdio>
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <string>
#include <vector>
#include <memory>
#include <stack>
#include <cstdbool>
#include <fstream>
#include <fenv.h>

feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW  | FE_UNDERFLOW);

#endif