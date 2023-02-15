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
    CoreShellParticle Ni_clad_Al_particle;

	real_t temperature;

	std::cout << "Density :\t" << Ni_clad_Al_particle.getDensity(298.15) << "\tkg/m3" << std::endl;
	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;

	std::cout << "Heat Capacity :\t" << Ni_clad_Al_particle.getHeatCapacity(temperature) << "\tJ/kg-K" << std::endl;
	std::cout << "Enthalpy :\t" <<  Ni_clad_Al_particle.getEnthalpy(temperature) << "\tJ/kg" << std::endl;

	std::cout << "Thermal Conductivity :\t" <<  Ni_clad_Al_particle.getThermalConductivity(temperature) << "\tW/m-K" << std::endl;
	
    return 0;
}