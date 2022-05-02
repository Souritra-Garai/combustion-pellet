#include <iostream>

#include "species/Aluminium.hpp"

template<typename real_t> real_t Phase<real_t>::sharpness_coefficient = 100;

int main(int argc, char const *argv[])
{
    for (double temperature = 298; temperature <= 2500; temperature += 0.1)
    {
       	std::cout     << temperature << ',' << Aluminium.getMolarMass()		* Aluminium.getHeatCapacity(temperature)         << ',' << Aluminium.getMolarMass()			* Aluminium.getEnthalpy(temperature)          << ',' << Aluminium.getThermalConductivity(temperature)         << std::endl;
	}

    return 0;
}
