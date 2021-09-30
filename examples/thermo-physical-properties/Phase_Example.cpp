#include <fstream>

#include <thermo-physical-properties/Phase.hpp>

#include <utilities/File_Generator.hpp>

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
    FileGenerator file_generator;

    std::ofstream my_file = file_generator.getCSVFile("Thermo_Physical_Properties");
    
    my_file << "Temperature (K)," << "Heat Capacity (J / mol-K)," << "Standard Enthalpy (J / mol)," << "Thermal Conductivity (W / m - K)" << std::endl;

    for (double temperature = 1; temperature <= 2000; temperature += 1)
    {
        my_file << temperature << ',' << solid_Al.getHeatCapacity(temperature) << ',' << solid_Al.getStandardEnthalpy(temperature) << ',' << solid_Al.getThermalConductivity(temperature) << std::endl;
    }

    my_file.close();

    return 0;
}
