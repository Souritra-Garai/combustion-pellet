#include <iostream>
#include <string.h>
#include <omp.h>

#include <boost/program_options.hpp>

#include "thermo-physical-properties/Substance.hpp"
#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.hpp"

#include "pde-problems/Core-Shell-Diffusion.hpp"
#include "pde-problems/Pellet-Flame-Propagation.hpp"

#include "utilities/File_Generator.hpp"
#include "utilities/Keyboard_Interrupt.hpp"

#include "substances/Argon.hpp"
#include "substances/Aluminium.hpp"
#include "substances/Nickel.hpp"
#include "substances/NickelAluminide.hpp"

#define MAX_ITER 1E8

long double core_radius = 32.5E-6;
long double overall_radius = 39.5E-6;

long double pellet_length = 6.35E-3;
long double pellet_diameter = 6.35E-3;

ArrheniusDiffusivityModel<long double> Alawieh_diffusivity(2.56E-6, 102.191E3);
ArrheniusDiffusivityModel<long double> Du_diffusivity(9.54E-8, 26E3);

long double phi = 0.68;
ArrheniusDiffusivityModel<long double> * diffusivity_model;

void printState(size_t iteration_number, PelletFlamePropagation<long double> &pellet);
int setCombustionConfiguration(int argc, char const *argv[]);

size_t data_capture_interval = round(0.005 / 0.000001);

int main(int argc, char const *argv[])
{
    CoreShellDiffusion<long double>::setUpCoreShellCombustionParticle(
        Aluminium, Nickel, NickelAluminide,
        overall_radius, core_radius
    );

    CoreShellDiffusion<long double>::setGridSize(1001);
    // CoreShellDiffusion<long double>::setTimeStep(0.0001);

    PelletFlamePropagation<long double>::setPelletDimensions(pellet_length, pellet_diameter);
    PelletFlamePropagation<long double>::setAmbientHeatLossParameters(0, 0, 0);
    PelletFlamePropagation<long double>::setTemperatureParameters(933, 298);
    PelletFlamePropagation<long double>::setDegassingFluid(Argon);

    PelletFlamePropagation<long double>::setGridSize(1001);
    PelletFlamePropagation<long double>::setTimeStep(0.000001);
    PelletFlamePropagation<long double>::setInfinitesimalChangeTemperature(0.001);
    PelletFlamePropagation<long double>::setInitialIgnitionParameters(1500, 0.1 * pellet_length);

	PelletFlamePropagation<long double>::setImplicitnessSourceTerm(1);
	PelletFlamePropagation<long double>::setImplicitnessDiffusionTerm(1);

	switch (setCombustionConfiguration(argc, argv))
	{
		case 0:	return 0;

		case 1: return 1;
		
		default:	;
	}

	PelletFlamePropagation<long double> combustion_pellet(phi);
    combustion_pellet.setDiffusivityModel(*diffusivity_model);
    combustion_pellet.initializePellet();

    FileGenerator file_generator;

    std::ofstream temperature_file = file_generator.getCSVFile("temperature");
    std::ofstream config_file = file_generator.getTXTFile("combustion_config");

	combustion_pellet.printConfiguration(config_file);
    combustion_pellet.printProperties(config_file);
    config_file.close();

    std::cout << std::endl;

    printState(0, combustion_pellet);
    
	size_t __iter = 1;
    setUpKeyboardInterrupt();
    
    try
    {
		combustion_pellet.printGridPoints(temperature_file, ',');

        while (!combustion_pellet.isCombustionComplete() && __iter <= MAX_ITER)
        {
            combustion_pellet.setUpEquations();
            combustion_pellet.solveEquations();

            if (__iter % data_capture_interval == 0)
			{
            	combustion_pellet.printTemperatureProfile(temperature_file, ',');
				std::cout << "Iteration # " << __iter << std::endl;
			}

            __iter++;
        }
    }

    catch (InterruptException& e)
    {
        std::cout << "\nCaught signal " << e.S << std::endl;

        printState(__iter, combustion_pellet);

        combustion_pellet.printTemperatureProfile(temperature_file, ',');
        temperature_file.close();

        return 1;
    }

    printState(__iter-1, combustion_pellet);
    
    combustion_pellet.printTemperatureProfile(temperature_file, ',');
    temperature_file.close();

    return 0;
}

void printState(size_t iteration_number, PelletFlamePropagation<long double> &pellet)
{
    std::cout << "Iteration # " << iteration_number << "\t";
    
    pellet.printTemperatureProfile(std::cout);

    std::cout << std::endl;
}

int setCombustionConfiguration(int argc, const char * argv[])
{
	try
	{
		boost::program_options::options_description options_description_obj("Usage: Pellet-Flame-Propagation [options]\nOptions");

		options_description_obj.add_options()
			("help",	"produce help message")
			("Du",		"set diffusivity parameters to Du et al's model (by default Alawieh et al's model is used)")
			("M",		boost::program_options::value<size_t>(),	"set number of grid points for pellet pde-solver to arg")
			("Dt",		boost::program_options::value<long double>(),	"set time step in seconds to arg")
			("phi",		boost::program_options::value<long double>(),	"set particle volume fractions to arg")
			("DT",		boost::program_options::value<long double>(),	"set infinitesimal change in temperature to arg, default is 0.1")
			("gamma",	boost::program_options::value<long double>(),	"set implicitness of source term to arg")
			("theta",	boost::program_options::value<long double>(),	"set implicitness of diffusion term to arg");

		boost::program_options::variables_map variables_map_obj;
		boost::program_options::store(
			boost::program_options::parse_command_line(argc, argv, options_description_obj),
			variables_map_obj
		);
		boost::program_options::notify(variables_map_obj);

		if (variables_map_obj.count("help"))
		{
			std::cout << options_description_obj << std::endl;
			return 0;
		}

		if (variables_map_obj.count("Du"))
		{
			diffusivity_model = &Du_diffusivity;		
			std::cout << "Using Du et al's Diffusivity parameters." << std::endl;
		}
		else
		{
			diffusivity_model = &Alawieh_diffusivity;
			std::cout << "Using Alawieh et al's Diffusivity parameters." << std::endl;
		}

		if (variables_map_obj.count("M"))
		{
			PelletFlamePropagation<long double>::setGridSize(variables_map_obj["M"].as<size_t>());
		}

		if (variables_map_obj.count("Dt"))
		{
			PelletFlamePropagation<long double>::setTimeStep(variables_map_obj["Dt"].as<long double>());
			data_capture_interval = round(0.005 / variables_map_obj["Dt"].as<long double>());
		}

		if (variables_map_obj.count("DT"))
		{
			PelletFlamePropagation<long double>::setInfinitesimalChangeTemperature(variables_map_obj["DT"].as<long double>());
		}

		if (variables_map_obj.count("phi"))
		{
			phi = variables_map_obj["phi"].as<long double>();

			if (phi > 0.0 && phi <= 1.0) 
			{
				PelletFlamePropagation<long double> combustion_pellet(phi);
				std::cout << "Particle volume fractions set to " << phi << std::endl;
			}
			else
			{
				std::cerr << "Particle volume fractions should be in the interval (0, 1]. Given value : " << phi << std::endl;
				return 1;
			}
		}

		if (variables_map_obj.count("gamma"))
		{
			PelletFlamePropagation<long double>::setImplicitnessSourceTerm(variables_map_obj["gamma"].as<long double>());
		}

		if (variables_map_obj.count("theta"))
		{
			PelletFlamePropagation<long double>::setImplicitnessDiffusionTerm(variables_map_obj["theta"].as<long double>());
		}
	}

	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	return 2;
}