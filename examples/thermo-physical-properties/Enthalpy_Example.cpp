#include <iostream>
#include "math/Shomate-Expression.hpp"
#include "utilities/Read-Data.hpp"

int main(int argc, char const *argv[])
{
	ShomateExpression enthalpy_Ni = readShomateExpressionCoefficients("data/species/nickel/solid-1");

	long double temperature;

	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;

	temperature = ShomateExpression::normalizeInput(temperature);

	std::cout << "Heat Capacity :\t" << enthalpy_Ni.evaluateExpression(temperature) << "\tJ/mol-K" << std::endl;
	std::cout << "Enthalpy :\t" <<  enthalpy_Ni.evaluateExpressionIntegral(temperature) << "\tkJ/mol-K" << std::endl;
	
	return 0;
}
