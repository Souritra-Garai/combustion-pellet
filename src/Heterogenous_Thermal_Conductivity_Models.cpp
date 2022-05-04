#include <fstream>
#include <iostream>

#include "utilities/File_Generator.hpp"
#include "utilities/Program_Options.hpp"

#include "species/Aluminium.hpp"
#include "species/Nickel.hpp"
#include "species/NickelAluminide.hpp"
#include "species/Argon.hpp"

#include "thermo-physical-properties/Core-Shell-Particle.hpp"
#include "thermo-physical-properties/Thermal_Conductivity_Pellet.hpp"

long double temperature = 298.15;

long double core_radius = 32.5E-6;
long double overall_radius = 39.5E-6;

void parseProgramOptions(int argc, char const *argv[]);

int main(int argc, char const *argv[])
{
	parseProgramOptions(argc, argv);

    CoreShellParticle<long double>::setUpCoreShellParticle(
        Aluminium, Nickel, NickelAluminide,
        overall_radius, core_radius
    );

    CoreShellParticle<long double> Ni_clad_Al_particle;

    FileGenerator file_generator;

    std::ofstream file = file_generator.getCSVFile("thermal_conductivity_pellet");

    file << "Particle Volume Fraction,ME1,EMT,MEB,CC,ME2" << std::endl;

    long double particle_thermal_conductivity = Ni_clad_Al_particle.getThermalConductivity(temperature);
    long double fluid_thermal_conductivity = Argon.getThermalConductivity(temperature);

    std::cout << "Particle Thermal Conductivity : " << particle_thermal_conductivity << std::endl;
    std::cout << "Fluid Thermal Conductivity : " << fluid_thermal_conductivity << std::endl;

    for (long double particle_volume_fraction = 0.001; particle_volume_fraction < 1.0; particle_volume_fraction += 0.001)
    {
        file << particle_volume_fraction << ',';
        file << getThermalConductivityME1(1.0 - particle_volume_fraction, fluid_thermal_conductivity, particle_thermal_conductivity) << ',';
        file << getThermalConductivityEMT(1.0 - particle_volume_fraction, fluid_thermal_conductivity, particle_thermal_conductivity) << ',';
        file << getThermalConductivityMEB(1.0 - particle_volume_fraction, fluid_thermal_conductivity, particle_thermal_conductivity) << ',';
        file << getThermalConductivityCC( 1.0 - particle_volume_fraction, fluid_thermal_conductivity, particle_thermal_conductivity) << ',';
        file << getThermalConductivityME2(1.0 - particle_volume_fraction, fluid_thermal_conductivity, particle_thermal_conductivity) << std::endl;
    }

    file.close();
    return 0;
}

void parseProgramOptions(int argc, char const *argv[])
{
	ez::ezOptionParser opt;

	opt.overview	= "Thermal Conductivity of a heterogenous solid-fluid mixture.";
	opt.syntax		= "Heterogenous_Thermal_Conductivity_Models [OPTIONS]";
	opt.example		= "Heterogenous_Thermal_Conductivity_Models --sharp 1.0 --T 1500.0\n\n";
	opt.footer		= "Developed by Souritra Garai, 2021-22.\n";

	setHelpOption(opt);
	setSharpnessCoefficientOption(opt);
	setTemperatureOption(opt);
	setCoreRadiusOption(opt);
	setOverallRadiusOption(opt);

	opt.parse(argc, argv);

	displayHelpOption(opt);

	Phase<long double>::sharpness_coefficient = getSharpnessCoefficientOption(opt, Phase<long double>::sharpness_coefficient);
	
	temperature = getTemperatureOption(opt, temperature);

	core_radius = getCoreRadiusOption(opt, core_radius);
	overall_radius = getOverallRadiusOption(opt, overall_radius);
}