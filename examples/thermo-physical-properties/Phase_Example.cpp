#include <iostream>

#include "thermo-physical-properties/Phase.hpp"
#include "utilities/Read-Data.hpp"

template<typename real_t> const real_t Phase<real_t>::sharpness_coefficient = 0.1;

int main(int argc, char const *argv[])
{
	Phase<long double> nickel_solid_1 = readPhaseData<long double>("data/species/nickel/solid-1");

	long double temperature;
	
	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;

	std::cout << "Heat Capacity : " << nickel_solid_1.getHeatCapacity(temperature) << std::endl;
	std::cout << "Enthalpy : " <<  nickel_solid_1.getStandardEnthalpy(temperature) << std::endl;

	std::cout << "Density : " << nickel_solid_1.getDensity(temperature) << std::endl;
	
	return 0;
}