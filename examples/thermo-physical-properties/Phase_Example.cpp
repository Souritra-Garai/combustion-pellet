#include <iostream>

#include "thermo-physical-properties/Phase.hpp"
#include "utilities/Read-Data.hpp"

const real_t Phase::sharpness_coefficient = 0.1;

int main(int argc, char const *argv[])
{
	Phase nickel_solid_1 = readPhaseData("data/species/nickel/solid-1");

	real_t temperature;
	
	std::cout << "Density :\t" << nickel_solid_1.getDensity(298.15) << "\tkg/m3" << std::endl;
	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;

	std::cout << "Heat Capacity :\t" << nickel_solid_1.getHeatCapacity(temperature) << "\tJ/mol-K" << std::endl;
	std::cout << "Enthalpy :\t" <<  nickel_solid_1.getStandardEnthalpy(temperature) << "\tJ/mol" << std::endl;

	std::cout << "Thermal Conductivity :\t" <<  nickel_solid_1.getThermalConductivity(temperature) << "\tW/m-K" << std::endl;
	
	return 0;
}