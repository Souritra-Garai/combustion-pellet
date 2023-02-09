#include <iostream>

#include "thermo-physical-properties/Ideal-Gas.hpp"
#include "utilities/Read-Data.hpp"

int main(int argc, char const *argv[])
{
	IdealGas<long double> argon = readIdealGasData<long double>("data/species/argon");

	long double temperature;
	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;

	std::cout << "Heat Capacity : " << argon.getCp(temperature) << std::endl;
	std::cout << "Enthalpy : " <<  argon.getEnthalpy(temperature) << std::endl;

	std::cout << "Molar Mass : " << argon.getMolarMass() << std::endl;
	
	return 0;
}