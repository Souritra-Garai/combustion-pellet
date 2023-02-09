#include <iostream>
#include "thermo-physical-properties/Enthalpy.hpp"
#include "utilities/Read-Data.hpp"

int main(int argc, char const *argv[])
{
	Enthalpy<long double> enthalpy_Ni = readEnthalpyData<long double>("data/species/nickel/solid-1");

	long double temperature;

	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;

	std::cout << "Heat Capacity : " << enthalpy_Ni.getHeatCapacity(temperature) << std::endl;
	std::cout << "Enthalpy : " <<  enthalpy_Ni.getStandardEnthalpy(temperature) << std::endl;
	
	return 0;
}
