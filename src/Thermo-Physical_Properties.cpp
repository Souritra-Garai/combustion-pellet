#include <iostream>
#include <fstream>

#include "utilities/File_Generator.hpp"
#include "utilities/Program_Options.hpp"

#include "species/Aluminium.hpp"
#include "species/Nickel.hpp"
#include "species/NickelAluminide.hpp"
#include "species/Argon.hpp"

template<typename real_t> real_t Phase<real_t>::sharpness_coefficient = 0.1;

long double temperature_lower_bound = 273.0;
long double temperature_upper_bound = 2500.0;

long double delta_T = 0.1;

void parseProgramOptions(int, char const *[]);

int main(int argc, char const *argv[])
{
	parseProgramOptions(argc, argv);
	
    FileGenerator file_generator;

    std::ofstream file_Ar   = file_generator.getCSVFile("Thermo_Physical_Properties_Ar");
    std::ofstream file_Al   = file_generator.getCSVFile("Thermo_Physical_Properties_Al");
    std::ofstream file_Ni   = file_generator.getCSVFile("Thermo_Physical_Properties_Ni");
    std::ofstream file_NiAl = file_generator.getCSVFile("Thermo_Physical_Properties_NiAl");
    
    file_Ar   << "Temperature (K)," << "Heat Capacity (J / mol-K)," << "Enthalpy (J / mol)," << "Thermal Conductivity (W / m - K)" << std::endl;
    file_Al   << "Temperature (K)," << "Heat Capacity (J / mol-K)," << "Enthalpy (J / mol)," << "Thermal Conductivity (W / m - K)" << std::endl;
    file_Ni   << "Temperature (K)," << "Heat Capacity (J / mol-K)," << "Enthalpy (J / mol)," << "Thermal Conductivity (W / m - K)" << std::endl;
    file_NiAl << "Temperature (K)," << "Heat Capacity (J / mol-K)," << "Enthalpy (J / mol)," << "Thermal Conductivity (W / m - K)" << std::endl;

    for (double temperature = temperature_lower_bound; temperature <= temperature_upper_bound; temperature += delta_T)
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

void parseProgramOptions(int argc, char const *argv[])
{
	ez::ezOptionParser opt;

	opt.overview	= "Print variation of thermo-physical properties with temperature for various species in csv files.";
	opt.syntax		= "Thermo-Physical_Properties [OPTIONS]";
	opt.example		= "Thermo-Physical_Properties --sharp 1.0 --Tlower 273.0 --Tupper 2500.0\n\n";
	opt.footer		= "Developed by Souritra Garai, 2021-22.\n";

	setHelpOption(opt);
	setSharpnessCoefficientOption(opt);
	setTemperatureUpperBoundOption(opt);
	setTemperatureLowerBoundOption(opt);
	setTemperatureStepOption(opt);

	opt.parse(argc, argv);

	displayHelpOption(opt);

	Phase<long double>::sharpness_coefficient = getSharpnessCoefficientOption(opt, Phase<long double>::sharpness_coefficient);
	
	temperature_lower_bound = getTemperatureLowerBoundOption(opt, temperature_lower_bound);
	temperature_upper_bound = getTemperatureUpperBoundOption(opt, temperature_upper_bound);

	delta_T = getTemperatureStepOption(opt, delta_T);
}