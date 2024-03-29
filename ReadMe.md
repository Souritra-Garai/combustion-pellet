# Combustion in Pellets packed with Core-Shell Structured Intermetallic Particles

[![DOI](https://zenodo.org/badge/381698482.svg)](https://zenodo.org/badge/latestdoi/381698482)

This program solves the PDEs governing flame propagation in a pellet packed with energetic intermetallic particles. The particles of concern are nickel-coated aluminium particles. The reaction kinetics for the combustion of these particles is assumed to be dominated by intermetallic diffusion. The pellet-scale heat equation is coupled to particle-level species diffusion equations and solved using an innovative quasi-implicit scheme. The physico-chemical process is represented in the following figure.

More details are available in the article -

Garai, S., & Sundaram, D. S., A Quasi-Implicit and Coupled Multiscale Scheme for Simulating Combustion of Pellets of Core-Shell Structured Intermetallic Particles. Combustion and Flame 250 (2023) 112650. DOI:[10.1016/j.combustflame.2023.112650](https://doi.org/10.1016/j.combustflame.2023.112650)

![Flame Propagation](https://github.com/Souritra-Garai/combustion-pellet/blob/master/images/Flame%20Propagation.jpg)

## Compiling the code

CMake is used to compile the present code. After cloning the repository, open a terminal in the directory where the code is present
```
mkdir build bin
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..
cmake --build .
```
You can chose your preferred compilers in place of `gcc` and `g++`.

## Running the code

Majority of thermo-physical data that the code requires as input is present in the `data` directory. Further, organisation of data in that directory is present in its own readme file.

After compilation, the `bin` directory will contain all the executables. To run the pellet flame propagation code, simply execute
```
mkdir solutions
bin/PelletFlamePropagationEXE
```
from a directory where the `data` and `solutions` directories is accessible. The temperature profile will be saved in a directory inside the `solutions` directory.

To see the available command line options like particle volume fractions of pellet etc., run
```
bin/PelletFlamePropagationEXE -help
```
The output should be generated as follows
```
Flame Propagation in a pellet packed with energetic intermetallic core-shell particles.

USAGE: Pellet-Flame-Propagation [OPTIONS]

OPTIONS:

-h, -help, --help, --usage   Display usage instructions.

-ignL ARG                    Set Initial Ignition Length Fraction of Pellet to
                             ARG.

-ignT ARG                    Set Initial Pellet Ignition Temperature to ARG K.

-phi ARG                     Set Particle Volume Fraction of Pellet to ARG.

EXAMPLES:

Pellet-Flame-Propagation -phi 0.68

Developed by Souritra Garai, 2021-23.
```
