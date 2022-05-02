#include <iostream>

#include <thermo-physical-properties/Phase.hpp>

Enthalpy<double> internal_energy_solid_Al(
    28.08920,
   -5.414849,
    8.560423,
    3.427370,
   -0.277375,
   -9.147187
);

ThermalConductivityQuadraticPolynomial<double> thermal_conductivity_solid_Al(248.0, -0.067, 0.0);

Phase<double> solid_Al(
    2700.0,
    internal_energy_solid_Al,
    thermal_conductivity_solid_Al,
    0,
    933,
    100
);

int main(int argc, char const *argv[])
{
    std::cout << "Temperature (K)," << "Heat Capacity (J / mol-K)," << "Standard Enthalpy (J / mol)," << "Thermal Conductivity (W / m - K)" << std::endl;

    for (double temperature = 273; temperature <= 2500; temperature += 1)
    {
        std::cout << temperature << ',' << solid_Al.getHeatCapacity(temperature) << ',' << solid_Al.getStandardEnthalpy(temperature) << ',' << solid_Al.getThermalConductivity(temperature) << std::endl;
    }

    return 0;
}
