/**
 * @file Core-Shell-Combustion-Particle_Example.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-07-09
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>

#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"

Substance<float> Al(2700, 897, 26.98, 239);
Substance<float> Ni(8908, 440, 58.69, 90.7);
Substance<float> NiAl(5900, 717, 85.675, 115, -118.4);

float core_radius = 32.5E-6;
float overall_radius = 39.5E-6;

int main(int argc, char const *argv[])
{
    CoreShellCombustionParticle<float>::setUpCoreShellCombustionParticle(
        Al, Ni, NiAl,
        overall_radius, core_radius
    );

    CoreShellCombustionParticle<float> Ni_clad_Al_particle;

    std::cout << "Enthalpy\t:\t" << Ni_clad_Al_particle.getEnthalpy(373) << "\tJ/kg" << std::endl << std::endl;

    Ni_clad_Al_particle.printProperties(std::cout);

    return 0;
}