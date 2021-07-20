/**
 * @file Core-Shell-Diffusion_Example.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Program to test functions of CoreShellDiffusion class
 * @version 0.1
 * @date 2021-07-16
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>

#include "thermo-physical-properties/Substance.hpp"
#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"
#include "pde-problems/Core-Shell-Diffusion.hpp"



CoreShellCombustionParticle<long double> Ni_clad_Al_particle(
    Al, Ni, NiAl, 2.56E-6, 102.191E3
);

CoreShellDiffusion<long double> Ni_clad_Al_particle_diffusion(
    Ni_clad_Al_particle, 1001, 32.5E-6, 30.5E-6
);

int main(int argc, char const *argv[])
{
    Ni_clad_Al_particle_diffusion.printProperties(std::cout);
    return 0;
}
