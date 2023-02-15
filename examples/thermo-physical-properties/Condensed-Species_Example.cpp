#include <iostream>
#include "thermo-physical-properties/Condensed-Species.hpp"
#include "utilities/Read-Data.hpp"

const real_t Phase::sharpness_coefficient = 0.1;

CondensedSpecies nickel = readCondensedSpeciesData("data/species/nickel");

int main(int argc, char const *argv[])
{
	real_t temperature;

	std::cout << "Density :\t" << nickel.getDensity(298.15) << "\tkg/m3" << std::endl;
	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;

	std::cout << "Heat Capacity :\t" << nickel.getHeatCapacity(temperature) << "\tJ/kg-K" << std::endl;
	std::cout << "Enthalpy :\t" <<  nickel.getEnthalpy(temperature) << "\tJ/kg" << std::endl;

	std::cout << "Thermal Conductivity :\t" <<  nickel.getThermalConductivity(temperature) << "\tW/m-K" << std::endl;
	
	return 0;
}