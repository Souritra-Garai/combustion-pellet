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

#include "substances/Aluminium.hpp"
#include "substances/Nickel.hpp"
#include "substances/NickelAluminide.hpp"

#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"

long double core_radius = 32.5E-6;
long double overall_radius = 39.5E-6;

int main(int argc, char const *argv[])
{
    CoreShellCombustionParticle<long double>::setUpCoreShellCombustionParticle(
        Aluminium, Nickel, NickelAluminide,
        overall_radius, core_radius
    );

    CoreShellCombustionParticle<long double> Ni_clad_Al_particle;

    std::cout << "Enthalpy\t:\t" << Ni_clad_Al_particle.getEnthalpy(373) << "\tJ/kg" << std::endl << std::endl;

    Ni_clad_Al_particle.printProperties(std::cout);

    return 0;
}