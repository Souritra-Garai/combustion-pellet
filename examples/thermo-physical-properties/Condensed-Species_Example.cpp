#include <iostream>
#include "thermo-physical-properties/Condensed-Species.hpp"
#include "utilities/Read-Data.hpp"

template<typename real_t> const real_t Phase<real_t>::sharpness_coefficient = 0.1;

CondensedSpecies<float> nickel = readCondensedSpeciesData<float>("data/species/nickel");

int main(int argc, char const *argv[])
{
	float temperature;

	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;
	
	std::cout << "Density : " << nickel.getThermalConductivity(temperature) << "\n";
	std::cout << "Enthalpy : " << nickel.getEnthalpy(temperature) << "\tC : " << nickel.getHeatCapacity(temperature) << "\n";
	std::cout << "Conductivity : " << nickel.getThermalConductivity(temperature) << "\n";
	
	return 0;
}