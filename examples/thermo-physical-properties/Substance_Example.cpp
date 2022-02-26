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
    
    file_Ar   << "Temperature (K)," << "Heat Capacity (J / mol-K)," << "Enthalpy (J / mol)," << "Thermal Conductivity (W / m - K)" << std::endl;
    file_Al   << "Temperature (K)," << "Heat Capacity (J / mol-K)," << "Enthalpy (J / mol)," << "Thermal Conductivity (W / m - K)" << std::endl;
    file_Ni   << "Temperature (K)," << "Heat Capacity (J / mol-K)," << "Enthalpy (J / mol)," << "Thermal Conductivity (W / m - K)" << std::endl;
    file_NiAl << "Temperature (K)," << "Heat Capacity (J / mol-K)," << "Enthalpy (J / mol)," << "Thermal Conductivity (W / m - K)" << std::endl;

    for (double temperature = 298; temperature <= 2500; temperature += 0.1)
    {
        file_Ar     << temperature << ',' << Argon.getMolarMass()			* Argon.getCp(temperature)			             << ',' << Argon.getMolarMass()				* Argon.getEnthalpy(temperature)              << ',' << Argon.getThermalConductivity(temperature)             << std::endl;

        file_Al     << temperature << ',' << Aluminium.getMolarMass()		* Aluminium.getHeatCapacity(temperature)         << ',' << Aluminium.getMolarMass()			* Aluminium.getEnthalpy(temperature)          << ',' << Aluminium.getThermalConductivity(temperature)         << std::endl;

        file_Ni     << temperature << ',' << Nickel.getMolarMass()			* Nickel.getHeatCapacity(temperature)            << ',' << Nickel.getMolarMass()			* Nickel.getEnthalpy(temperature)             << ',' << Nickel.getThermalConductivity(temperature)            << std::endl;

        file_NiAl   << temperature << ',' << NickelAluminide.getMolarMass() * NickelAluminide.getHeatCapacity(temperature)   << ',' << NickelAluminide.getMolarMass() 	* NickelAluminide.getEnthalpy(temperature)    << ',' << NickelAluminide.getThermalConductivity(temperature)   << std::endl;
    }

    file_Ar.close();
    file_Al.close();
    file_Ni.close();
    file_NiAl.close();
    
    return 0;
}
