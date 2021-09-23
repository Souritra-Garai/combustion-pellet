/**
 * @file Substance_Example.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Example to test functions of Substance class
 * @version 0.1
 * @date 2021-07-06
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <fstream>

#include "utilities/File_Generator.hpp"

#include "substances/Aluminium.hpp"
#include "substances/Nickel.hpp"
#include "substances/NickelAluminide.hpp"
#include "substances/Argon.hpp"

int main(int argc, char const *argv[])
{
    FileGenerator file_generator;

    std::ofstream file_Ar   = file_generator.getCSVFile("Thermo_Physical_Properties_Ar");
    std::ofstream file_Al   = file_generator.getCSVFile("Thermo_Physical_Properties_Al");
    std::ofstream file_Ni   = file_generator.getCSVFile("Thermo_Physical_Properties_Ni");
    std::ofstream file_NiAl = file_generator.getCSVFile("Thermo_Physical_Properties_NiAl");
    
    file_Ar   << "Temperature (K)," << "Heat Capacity (J / kg-K)," << "Internal Energy (J / kg)," << "Thermal Conductivity (W / m - K)" << std::endl;
    file_Al   << "Temperature (K)," << "Heat Capacity (J / kg-K)," << "Internal Energy (J / kg)," << "Thermal Conductivity (W / m - K)" << std::endl;
    file_Ni   << "Temperature (K)," << "Heat Capacity (J / kg-K)," << "Internal Energy (J / kg)," << "Thermal Conductivity (W / m - K)" << std::endl;
    file_NiAl << "Temperature (K)," << "Heat Capacity (J / kg-K)," << "Internal Energy (J / kg)," << "Thermal Conductivity (W / m - K)" << std::endl;

    for (double temperature = 300; temperature <= 2000; temperature += 0.1)
    {
        file_Ar     << temperature << ',' << Argon.getHeatCapacity(temperature)             << ',' << Argon.getInternalEnergy(temperature)              << ',' << Argon.getThermalConductivity(temperature)             << std::endl;

        file_Al     << temperature << ',' << Aluminium.getHeatCapacity(temperature)         << ',' << Aluminium.getInternalEnergy(temperature)          << ',' << Aluminium.getThermalConductivity(temperature)         << std::endl;

        file_Ni     << temperature << ',' << Nickel.getHeatCapacity(temperature)            << ',' << Nickel.getInternalEnergy(temperature)             << ',' << Nickel.getThermalConductivity(temperature)            << std::endl;

        file_NiAl   << temperature << ',' << NickelAluminide.getHeatCapacity(temperature)   << ',' << NickelAluminide.getInternalEnergy(temperature)    << ',' << NickelAluminide.getThermalConductivity(temperature)   << std::endl;
    }

    file_Ar.close();
    file_Al.close();
    file_Ni.close();
    file_NiAl.close();
    
    return 0;
}
