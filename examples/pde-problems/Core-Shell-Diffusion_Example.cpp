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

Substance<float> Al(2700, 897, 26.98, 239);
Substance<float> Ni(8902, 440, 58.69, 90.7);
Substance<float> NiAl(5900, 717, 85.675, 115, -118.4);

float core_radius = 32.5E-6;
float overall_radius = 39.5E-6;

int main(int argc, char const *argv[])
{
    CoreShellDiffusion<float>::setUpCoreShellCombustionParticle(
        Al, Ni, NiAl,
        overall_radius, core_radius
    );

    CoreShellDiffusion<float>::setGridSize(100001);
    CoreShellDiffusion<float>::setTimeStep(0.001);

    CoreShellDiffusion<float> Ni_clad_Al_particle;

    Ni_clad_Al_particle.printProperties(std::cout);

    return 0;
}
