all:
	g++ -fopenmp -Wall -W -Wmain -Werror -Wextra -std=c++11 -pedantic-errors -pedantic -Ofast -fexpensive-optimizations -march=native *.cpp *.hpp -o large_scale_simulation
	g++ -fopenmp -Wall -W -Wmain -Werror -Wextra -std=c++11 -pedantic-errors -pedantic -O0 *.cpp *.hpp -o large_scale_simulation_debug

