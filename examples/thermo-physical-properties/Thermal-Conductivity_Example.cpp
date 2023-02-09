#include <iostream>

#include "thermo-physical-properties/Thermal-Conductivity.hpp"
#include "utilities/Read-Data.hpp"

int main(int argc, char const *argv[])
{
	ThermalConductivityQuadraticPolynomial<long double> thermal_conductivity_Ni = readThermalConductivityData<long double>("data/species/aluminium/solid");

	long double temperature;
	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;

	std::cout << "Heat Capacity : " << thermal_conductivity_Ni.getThermalConductivity(temperature) << std::endl;
	
	return 0;
}
