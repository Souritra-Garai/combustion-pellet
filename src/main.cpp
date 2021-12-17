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

#define MAX_ITER 100000

double core_radius = 32.5E-6;
double overall_radius = 39.5E-6;

double pellet_length = 6.35E-3;
double pellet_diameter = 6.35E-3;

ArrheniusDiffusivityModel<double> Alawieh_diffusivity(2.56E-6, 102.191E3);
ArrheniusDiffusivityModel<double> Du_diffusivity(9.54E-8, 26E3);

double phi = 0.7;
ArrheniusDiffusivityModel<double> * diffusivity_model;

void printState(size_t iteration_number, PelletFlamePropagation<double> &pellet);
int setCombustionConfiguration(int argc, char const *argv[]);

size_t data_capture_interval = round(0.005 / 0.0001);

int main(int argc, char const *argv[])
{
    CoreShellDiffusion<double>::setUpCoreShellCombustionParticle(
        Aluminium, Nickel, NickelAluminide,
        overall_radius, core_radius
    );

    CoreShellDiffusion<double>::setGridSize(1001);
    // CoreShellDiffusion<double>::setTimeStep(0.0001);

    PelletFlamePropagation<double>::setPelletDimensions(pellet_length, pellet_diameter);
    PelletFlamePropagation<double>::setAmbientHeatLossParameters(0, 0);
    PelletFlamePropagation<double>::setTemperatureParameters(933, 298);
    PelletFlamePropagation<double>::setDegassingFluid(Argon);

    PelletFlamePropagation<double>::setGridSize(101);
    PelletFlamePropagation<double>::setTimeStep(0.0001);
    PelletFlamePropagation<double>::setInfinitesimalChangeTemperature(0.1);
    PelletFlamePropagation<double>::setInitialIgnitionParameters(1500, 0.1 * pellet_length);

	// PelletFlamePropagation<double>::setImplicitnessSourceTerm(0.5);
	// PelletFlamePropagation<double>::setImplicitnessDiffusionTerm(0.5);

	switch (setCombustionConfiguration(argc, argv))
	{
		case 0:	return 0;

		case 1: return 1;
		
		default:	;
	}

	PelletFlamePropagation<double> combustion_pellet(phi);
    combustion_pellet.setDiffusivityModel(*diffusivity_model);
    combustion_pellet.initializePellet();

    FileGenerator file_generator;

    std::ofstream temperature_file = file_generator.getCSVFile("temperature");
    std::ofstream config_file = file_generator.getTXTFile("combustion_config");

	combustion_pellet.printConfiguration(config_file);
    combustion_pellet.printProperties(config_file);
    config_file.close();

    std::cout << std::endl;

    size_t __iter = 1;

    printState(0, combustion_pellet);

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
				std::cout << "Iteration # " << __iter+1 << std::endl;
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

void printState(size_t iteration_number, PelletFlamePropagation<double> &pellet)
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
			("Dt",		boost::program_options::value<double>(),	"set time step in seconds to arg")
			("phi",		boost::program_options::value<double>(),	"set particle volume fractions to arg")
			("DT",		boost::program_options::value<double>(),	"set infinitesimal change in temperature to arg, default is 0.1")
			("gamma",	boost::program_options::value<double>(),	"set implicitness of source term to arg")
			("theta",	boost::program_options::value<double>(),	"set implicitness of diffusion term to arg");

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
			PelletFlamePropagation<double>::setGridSize(variables_map_obj["M"].as<size_t>());
		}

		if (variables_map_obj.count("Dt"))
		{
			PelletFlamePropagation<double>::setTimeStep(variables_map_obj["Dt"].as<double>());
			data_capture_interval = round(0.005 / variables_map_obj["Dt"].as<double>());
		}

		if (variables_map_obj.count("DT"))
		{
			PelletFlamePropagation<double>::setInfinitesimalChangeTemperature(variables_map_obj["DT"].as<double>());
		}

		if (variables_map_obj.count("phi"))
		{
			phi = variables_map_obj["phi"].as<double>();

			if (phi > 0.0 && phi <= 1.0) 
			{
				PelletFlamePropagation<double> combustion_pellet(phi);
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
			PelletFlamePropagation<double>::setImplicitnessSourceTerm(variables_map_obj["gamma"].as<double>());
		}

		if (variables_map_obj.count("theta"))
		{
			PelletFlamePropagation<double>::setImplicitnessDiffusionTerm(variables_map_obj["theta"].as<double>());
		}
	}

	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	return 2;
}