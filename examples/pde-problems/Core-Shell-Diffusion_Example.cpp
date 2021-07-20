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

Substance<long double> Al(2700, 897, 26.98, 239);
Substance<long double> Ni(8902, 440, 58.69, 90.7);
Substance<long double> NiAl(5900, 717, 85.675, 115, -118.4);

long double core_radius = 32.5E-6;
long double overall_radius = 39.5E-6;
long double mass = calcMassCoreShellParticle<long double>(Al, Ni, overall_radius, core_radius);

CoreShellCombustionParticle<long double> Ni_clad_Al_particle(
    Al, Ni, NiAl, 
    overall_radius,
    core_radius,
    mass
);

CoreShellDiffusion<long double> Ni_clad_Al_particle_diffusion(
    Ni_clad_Al_particle, 1001
);

int main(int argc, char const *argv[])
{
    Ni_clad_Al_particle_diffusion.printProperties(std::cout);
    Ni_clad_Al_particle.printProperties(std::cout);
    return 0;
}
