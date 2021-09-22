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

int main(int argc, char const *argv[])
{
    FileGenerator file_generator;

    std::ofstream my_file = file_generator.getCSVFile("Thermo_Physical_Properties");
    
    my_file << "Temperature (K)," << "Heat Capacity (J / kg-K)," << "Internal Energy (J / kg)," << "Thermal Conductivity (W / m - K)" << std::endl;

    for (double temperature = 920; temperature <= 940; temperature += 0.01)
    {
        my_file << temperature << ',' << Aluminium.getHeatCapacity(temperature) << ',' << Aluminium.getInternalEnergy(temperature) << ',' << Aluminium.getThermalConductivity(temperature) << std::endl;
    }

    my_file.close();
    
    return 0;
}
