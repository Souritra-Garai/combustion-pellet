#include <iostream>

#include "species/Aluminium.hpp"
#include "species/Nickel.hpp"
#include "species/NickelAluminide.hpp"

#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.hpp"
#include "pde-problems/Core-Shell-Diffusion.hpp"

#include "utilities/Keyboard_Interrupt.hpp"
#include "utilities/Program_Options.hpp"
#include "utilities/File_Generator.hpp"

#define MAX_ITER 1E8

size_t n = 1001;
long double delta_t = 1E-6;

long double core_radius = 32.5E-6;
long double overall_radius = 39.5E-6;

long double temperature = 1900;

ArrheniusDiffusivityModel<long double> Alawieh_diffusivity(2.56E-6, 102.191E3);
ArrheniusDiffusivityModel<long double> Du_diffusivity(9.54E-8, 26E3);
ArrheniusDiffusivityModel<long double> *diffusivity_model = &Alawieh_diffusivity;

void parseProgramOptions(int argc, char const *argv[]);

int main(int argc, char const *argv[])
{
	parseProgramOptions(argc, argv);
	
    CoreShellDiffusion<long double>::setUpCoreShellParticle(
        Aluminium, Nickel, NickelAluminide,
        overall_radius, core_radius
    );

    CoreShellDiffusion<long double>::setGridSize(n);
    CoreShellDiffusion<long double>::setTimeStep(delta_t);

    CoreShellDiffusion<long double> Ni_clad_Al_particle;

    FileGenerator file_generator;

    std::ofstream config_file = file_generator.getTXTFile("diffusion_config");
    std::ofstream conc_A_file = file_generator.getCSVFile("concentration_A");
    std::ofstream conc_B_file = file_generator.getCSVFile("concentration_B");

    Ni_clad_Al_particle.printGridPoints(conc_A_file, ',');
    Ni_clad_Al_particle.printGridPoints(conc_B_file, ',');

    Ni_clad_Al_particle.printProperties(config_file);
    config_file.close();

    long double diffusivity = diffusivity_model->getDiffusivity(temperature);

    setUpKeyboardInterrupt();

	long double simulation_time = 0.0;

    try
    {
		size_t step = 0.001 / delta_t;

		bool combustion_complete = false;

		for (size_t i = 0; i < MAX_ITER && !combustion_complete;)
		{
			size_t i_step = i + step;

			for (; i < i_step && !combustion_complete; i++)
			{
				Ni_clad_Al_particle.setUpEquations(diffusivity);
            	Ni_clad_Al_particle.solveEquations();

				simulation_time += delta_t;

				combustion_complete = Ni_clad_Al_particle.isCombustionComplete();           
			}

			Ni_clad_Al_particle.printConcentrationProfileA(conc_A_file, ',', simulation_time);
            Ni_clad_Al_particle.printConcentrationProfileB(conc_B_file, ',', simulation_time);

			std::cout << "Iterations Completed : " << i << "\n";
		}
    }

    catch (InterruptException& e)
    {
        std::cout << "\nCaught signal " << e.S << std::endl;

		Ni_clad_Al_particle.printConcentrationProfileA(conc_A_file, ',', simulation_time);
		Ni_clad_Al_particle.printConcentrationProfileB(conc_B_file, ',', simulation_time);

        conc_A_file.close();
        conc_B_file.close();

        return 1;
    }

    conc_A_file.close();
    conc_B_file.close();

	CoreShellDiffusion<long double>::deallocateRadiusArray();
    
    return 0;
}

void parseProgramOptions(int argc, char const *argv[])
{
	ez::ezOptionParser opt;

	opt.overview	= "Diffusion of Al and Ni in a Core-Shell particle.";
	opt.syntax		= "Core-Shell-Diffusion [OPTIONS]";
	opt.example		= "Core-Shell-Diffusion --sharp 1.0 --T 1500.0\n\n";
	opt.footer		= "Developed by Souritra Garai, 2021-22.\n";

	setHelpOption(opt);
	setSharpnessCoefficientOption(opt);
	setParticleGridSizeOption(opt);
	setTimeStepOption(opt);
	setTemperatureOption(opt);

	opt.parse(argc, argv);

	displayHelpOption(opt);

	Phase<long double>::sharpness_coefficient = getSharpnessCoefficientOption(opt, Phase<long double>::sharpness_coefficient);
	
	temperature = getTemperatureOption(opt, temperature);

	n = getParticleGridSizeOption(n);

	delta_t = getTimeStepOption(delta_t);
}