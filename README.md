# large_scale_structure
===

The universe started with nearly uniformly distributed matter which clumped together due to small perturbations. This simulation attempts to numerically simulate the formation of galactic superclusters across 15 billion years.

This simulation uses a barnes-hut tree code to achieve high performance. In addition, a hybrid mesh method is used to deal with periodic boundary conditions in a computationally efficient method. Finally OMP is used to significantly speed up the simulation on multiple processors.

The numerical integration used is a symplectic leapfrog method. In order to avoid stability issues, plummer softening is used. A damping factor is also added to the simulation to preserve the conservation of energy.

===

Options: preface these with --

filename : string
the directory to which the ouput will be saved
defaults to current directory

density : double
the density of the universe
defaults to 1E-26 kg*m^-2

size : double
the size (diameter) of the universe in m
defaults to 2E+24 m

plummer : double
a smoothing parameter used to prevent singularities
defaults to 5E+21 m

gravity: double
the gravitational constant
defaults to 6.67E-11 m^3*kg^-1*s^-2

hubble: double
the hubble constant used to calculate metric expansion
defaults to 2.25E-18 s^-1

damping: double
a stability parameter used to maintain constant energy
defaults to 1E-17 s^-1

simtime: double
amount of simulation time to run for
defaults to 5E+17 s

timestep: double
amount of time per timestep
defaults to 1E+14 s

quadtree_ratio: double
a parameter which controls the degree of accuracy
defaults to 3.0

mesh_resolution : int
a parameter which controls the degree of accuracy
defaults to 4096

periodic_tilings : int
a parameter which controls the degree of accuracy
defaults to 31

mesh_limits : int
a parameter which determines whether the mesh is used for each body
defaults to 64

bodies : int
the number of gravitational bodies to simulate
defaults to 2^20 = 1048576

draw_size : int
dimensions of the ouput image in pixels
defaults to 4096

displacement : double
degree of uniformity in the initial conditions
defaults to 0.2

init_velocity : double
initial random velocity given to bodies
defaults to 5E+4

verbose : int
indicates the degree of verbosity
defaults to 1