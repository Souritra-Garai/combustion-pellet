/**
 * @file Adiabatic_Combustion_Temperature.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-09-30
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>
#include <fstream>

#include "substances/Aluminium.hpp"
#include "substances/Nickel.hpp"
#include "substances/NickelAluminide.hpp"

#include "utilities/File_Generator.hpp"

#define MAX_ITER 1000U
#define MAX_TEMP 6E3	// K

int main(int argc, char const *argv[])
{
	FileGenerator file_generator;

    std::ofstream file = file_generator.getCSVFile("Adiabatic_Combustion_Temperature");
	
	file << "Initial Temperature (K)," << "Adiabatic Combustion Temperature (K)" << std::endl;
	
	for (double initial_temperature = 298; initial_temperature <= 2000; initial_temperature+=1)
	{
		double initial_energy = // J
			Aluminium.getInternalEnergy(initial_temperature) * Aluminium.getMolarMass() +
			Nickel.getInternalEnergy(initial_temperature) * Nickel.getMolarMass();

		// std::cout << "Intial Energy : " << initial_energy << " J\n";

		double temperature_upper = MAX_TEMP;
		double temperature_lower = initial_temperature;

		double temperature_middle = (temperature_upper + temperature_lower) / 2.0;
		double final_energy = NickelAluminide.getInternalEnergy(temperature_middle) * NickelAluminide.getMolarMass();

		unsigned int i = 0;

		while (abs(final_energy - initial_energy) > 0.1 && i < MAX_ITER)
		{
			if (final_energy > initial_energy) temperature_upper = temperature_middle;

			else temperature_lower = temperature_middle;

			temperature_middle = (temperature_upper + temperature_lower) / 2.0;
			final_energy = NickelAluminide.getInternalEnergy(temperature_middle) * NickelAluminide.getMolarMass();
		
			// std::cout << "Iteration # " << ++i << "\t\tTemperature : " << temperature_middle << " K\t\tEnergy : " << final_energy << " J\n";
		}

		file << initial_temperature << ',' << temperature_middle << std::endl;
	}

	file.close();

	return 0;
}
