/**
 * @file Core-Shell-Particle_Example.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-07-09
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>

#include "species/Aluminium.hpp"
#include "species/Nickel.hpp"
#include "species/NickelAluminide.hpp"

#include "thermo-physical-properties/Core-Shell-Particle.hpp"

long double core_radius = 32.5E-6;
long double overall_radius = 39.5E-6;

int main(int argc, char const *argv[])
{
    CoreShellParticle<long double>::setUpCoreShellParticle(
        Aluminium, Nickel, NickelAluminide,
        overall_radius, core_radius
    );

    CoreShellParticle<long double> Ni_clad_Al_particle;

    std::cout << "Enthalpy\t:\t" << Ni_clad_Al_particle.getEnthalpy(373) << "\tJ/kg" << std::endl << std::endl;

    Ni_clad_Al_particle.printProperties(std::cout);

    return 0;
}