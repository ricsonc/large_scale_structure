# large_scale_structure
===

The universe started with nearly uniformly distributed matter which clumped together due to small perturbations. This simulation attempts to numerically simulate the formation of galactic superclusters across 15 billion years.

This simulation uses a barnes-hut tree code to achieve high performance. In addition, a hybrid mesh method is used to deal with periodic boundary conditions in a computationally efficient method. Finally OMP is used to significantly speed up the simulation on multiple processors.

The numerical integration used is a symplectic leapfrog method. In order to avoid stability issues, plummer softening is used. A damping factor is also added to the simulation to preserve the conservation of energy.
