#include <iostream>

#include "thermo-physical-properties/Ideal-Gas.hpp"
#include "utilities/Read-Data.hpp"

int main(int argc, char const *argv[])
{
	IdealGas argon = readIdealGasData("data/species/argon");

	std::cout << "Molar Mass :\t" << argon.molar_mass << "\tkg/mol" << std::endl;

	long double temperature;
	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;

	std::cout << "Heat Capacity :\t" << argon.getCp(temperature) << "\tJ/kg-K" << std::endl;
	std::cout << "Enthalpy :\t" <<  argon.getEnthalpy(temperature) << "\tJ/kg" << std::endl;
	std::cout << "Thermal Conductivity :\t" << argon.getThermalConductivity(temperature) << "\tW/m-K" << std::endl;
	
	return 0;
}