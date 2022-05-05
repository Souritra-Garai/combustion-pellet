#include <iostream>

#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.hpp"

#include "pde-problems/Core-Shell-Diffusion.hpp"
#include "pde-problems/Pellet-Flame-Propagation.hpp"

#include "utilities/File_Generator.hpp"
#include "utilities/Program_Options.hpp"
#include "utilities/Keyboard_Interrupt.hpp"

#include "species/Argon.hpp"
#include "species/Aluminium.hpp"
#include "species/Nickel.hpp"
#include "species/NickelAluminide.hpp"

#define MAX_ITER 1E6

long double delta_t = 1E-6;

size_t num_grid_points_pellet = 1001;
size_t num_grid_points_particle = 1001;

long double delta_T = 0.001;

long double core_radius = 32.5E-6;
long double overall_radius = 39.5E-6;

long double pellet_length = 6.35E-3;
long double pellet_diameter = 6.35E-3;

long double diffusion_term_implicitness = 0.5;
long double source_term_implicitness = 0.5;

ArrheniusDiffusivityModel<long double> diffusivity_model(2.56E-6, 102.191E3);

long double phi = 0.7;

void parseProgramOptions(int argc, char const *argv[]);

int main(int argc, char const *argv[])
{
	parseProgramOptions(argc, argv);

    CoreShellDiffusion<long double>::setUpCoreShellParticle(
        Aluminium, Nickel, NickelAluminide,
        overall_radius, core_radius
    );

    CoreShellDiffusion<long double>::setGridSize(num_grid_points_particle);
    // CoreShellDiffusion<long double>::setTimeStep(0.0001);

    PelletFlamePropagation<long double>::setPelletDimensions(pellet_length, pellet_diameter);
    PelletFlamePropagation<long double>::setAmbientHeatLossParameters(0, 0, 0);
    PelletFlamePropagation<long double>::setTemperatureParameters(933, 298);
    PelletFlamePropagation<long double>::setDegassingFluid(Argon);

    PelletFlamePropagation<long double>::setGridSize(num_grid_points_pellet);
    PelletFlamePropagation<long double>::setTimeStep(delta_t);
    PelletFlamePropagation<long double>::setInfinitesimalChangeTemperature(delta_T);
    PelletFlamePropagation<long double>::setInitialIgnitionParameters(1500, 0.1 * pellet_length);

	PelletFlamePropagation<long double>::setImplicitnessSourceTerm(source_term_implicitness);
	PelletFlamePropagation<long double>::setImplicitnessDiffusionTerm(diffusion_term_implicitness);

	PelletFlamePropagation<long double> combustion_pellet(phi);
    combustion_pellet.setDiffusivityModel(diffusivity_model);
    combustion_pellet.initializePellet();

    FileGenerator file_generator;

    std::ofstream temperature_file = file_generator.getCSVFile("temperature");
    std::ofstream config_file = file_generator.getTXTFile("combustion_config");

	combustion_pellet.printConfiguration(config_file);
    combustion_pellet.printProperties(config_file);
    config_file.close();

    std::cout << "Initialized Pellet. Starting iterations.\nPress Ctrl+C to stop...\n\n";

    setUpKeyboardInterrupt();
    
    try
    {
		combustion_pellet.printGridPoints(temperature_file, ',');

		size_t step = 0.001 / delta_t;

		bool combustion_not_complete = true;

		for (size_t i = 0; i < MAX_ITER && combustion_not_complete;)
		{
			size_t i_step = i + step;

			for (;i < i_step && combustion_not_complete; i++)
			{
				combustion_pellet.setUpEquations();
				combustion_pellet.solveEquations();

				combustion_not_complete = !combustion_pellet.isCombustionComplete();
			}

			combustion_pellet.printTemperatureProfile(temperature_file, ',');
			std::cout << "Iterations Completed : " << i << "\n";

		}

    	std::cout << "\nPellet combustion complete.\n";
    }

    catch (InterruptException& e)
    {
        std::cout << "\nCaught signal " << e.S << std::endl;

        combustion_pellet.printTemperatureProfile(temperature_file, ',');
        temperature_file.close();

        return 1;
    }

    combustion_pellet.printTemperatureProfile(temperature_file, ',');
    temperature_file.close();

	CoreShellDiffusion<long double>::deallocateRadiusArray();

    return 0;
}

void parseProgramOptions(int argc, char const *argv[])
{
	ez::ezOptionParser opt;

	opt.overview	= "Flame Propagation in a pellet packed with energetic intermetallic core-shell particles.";
	opt.syntax		= "Pellet-Flame-Propagation [OPTIONS]";
	opt.example		= "Pellet-Flame-Propagation --phi 0.68\n\n";
	opt.footer		= "Developed by Souritra Garai, 2021-22.\n";

	setHelpOption(opt);
	setSharpnessCoefficientOption(opt);
	setParticleGridSizeOption(opt);
	setPelletGridSizeOption(opt);
	setTimeStepOption(opt);
	setTemperatureStepOption(opt);
	setDiffusivityModelOption(opt);
	setCoreRadiusOption(opt);
	setOverallRadiusOption(opt);
	setLengthOption(opt);
	setDiameterOption(opt);
	setGammaOption(opt);
	setKappaOption(opt);
	setPhiOption(opt);

	opt.parse(argc, argv);

	displayHelpOption(opt);

	Phase<long double>::sharpness_coefficient = getSharpnessCoefficientOption(opt, Phase<long double>::sharpness_coefficient);
	
	delta_T = getTemperatureStepOption(opt, delta_T);

	num_grid_points_particle = getParticleGridSizeOption(opt, num_grid_points_particle);
	num_grid_points_pellet   = getPelletGridSizeOption(opt, num_grid_points_pellet);

	delta_t = getTimeStepOption(opt, delta_t);

	getDiffusivityModelOption(opt, diffusivity_model);

	core_radius = getCoreRadiusOption(opt, core_radius);
	overall_radius = getOverallRadiusOption(opt, overall_radius);

	pellet_length = getLengthOption(opt, pellet_length);
	pellet_diameter = getDiameterOption(opt, pellet_diameter);

	diffusion_term_implicitness = getKappaOption(opt, diffusion_term_implicitness);
	source_term_implicitness = getGammaOption(opt, source_term_implicitness);

	phi = getPhiOption(opt, phi);
}