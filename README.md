#[Video demonstration](https://drive.google.com/file/d/0B-ynXs-zChWiWFBLWXFpTTU2LVU/view?usp=sharing)

===

#Overview

The universe started with nearly uniformly distributed matter which clumped together due to small perturbations. This simulation attempts to numerically simulate the formation of galactic superclusters across 15 billion years. The simulation encompasses approximately one part in one million of the observable universe.

This simulation uses a barnes-hut tree code to achieve high performance. In addition, a hybrid mesh method is used to deal with periodic boundary conditions in a computationally efficient manner. Finally OMP is used to significantly speed up the simulation on multiple processors.

The numerical integration used is a symplectic leapfrog method. In order to avoid stability issues, plummer softening is used. A damping factor is also added to the simulation to preserve the conservation of energy.

===

#command line options :

filename : string <br />
the directory to which the ouput will be saved <br />
defaults to current directory

density : double <br />
the density of the universe <br />
defaults to 1E-26 kg*m^-2

size : double <br />
the size (diameter) of the universe in m <br />
defaults to 2E+24 m

plummer : double <br />
a smoothing parameter used to prevent singularities <br />
defaults to 5E+21 m

gravity: double <br />
the gravitational constant <br />
defaults to 6.67E-11 m^3*kg^-1*s^-2

hubble: double <br />
the hubble constant used to calculate metric expansion <br />
defaults to 2.25E-18 s^-1

damping: double <br />
a stability parameter used to maintain constant energy <br />
defaults to 1E-17 s^-1

simtime: double <br />
amount of simulation time to run for <br />
defaults to 5E+17 s

timestep: double <br />
amount of time per timestep <br />
defaults to 1E+14 s

quadtree_ratio: double <br />
a parameter which controls the degree of accuracy <br />
defaults to 3.0

mesh_resolution : int <br />
a parameter which controls the degree of accuracy <br />
defaults to 4096

periodic_tilings : int <br />
a parameter which controls the degree of accuracy <br />
defaults to 31

mesh_limits : int <br />
a parameter which determines whether the mesh is used for each body <br />
defaults to 64

bodies : int <br />
the number of gravitational bodies to simulate <br />
defaults to 2^20 = 1048576

draw_size : int <br />
dimensions of the ouput image in pixels <br />
defaults to 4096

displacement : double <br />
degree of uniformity in the initial conditions <br />
defaults to 0.2

init_velocity : double <br />
initial random velocity given to bodies <br />
defaults to 5E+4

verbose : int <br />
indicates the degree of verbosity <br />
defaults to 1
