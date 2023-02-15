#include <iostream>

#include "math/Quadratic-Expression.hpp"
#include "utilities/Read-Data.hpp"

int main(int argc, char const *argv[])
{
	QuadraticExpression thermal_conductivity_Ni = readQuadraticExpressionCoefficients("data/species/aluminium/solid");

	long double temperature;
	std::cout << "Enter Temperature (K) : ";
	std::cin >> temperature;

	std::cout << "Thermal Conductivity :\t" << thermal_conductivity_Ni.evaluateExpression(temperature) << "\tW/m-K" << std::endl;
	
	return 0;
}
