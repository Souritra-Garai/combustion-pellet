/**
 * @file Core-Shell-Particle_Example.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief 
 * @version 0.1
 * @date 2021-07-09
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>

#include "thermo-physical-properties/Core-Shell-Particle.hpp"

int main(int argc, char const *argv[])
{
    CoreShellParticle<float> Ni_clad_Al_particle;

	std::cout << "Enthalpy\t:\t" << Ni_clad_Al_particle.getEnthalpy(298.15) << "\n";
	std::cout << "Density\t\t:\t" << Ni_clad_Al_particle.getDensity(298.15) << "\n";

    return 0;
}