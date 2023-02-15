#include <iostream>

#include "pde-problems/Core-Shell-Diffusion.hpp"

#include "utilities/Keyboard-Interrupt.hpp"
#include "utilities/Program-Options.hpp"
#include "utilities/File-Generator.hpp"

#define MAX_ITER 1E8

long double temperature = 1900;

void parseProgramOptions(int argc, char const *argv[]);

int main(int argc, char const *argv[])
{
	parseProgramOptions(argc, argv);

    CoreShellDiffusion::setUpRadiusArray();

    CoreShellDiffusion Ni_clad_Al_particle;

    FileGenerator file_generator;

    std::ofstream conc_A_file = file_generator.getCSVFile("concentration_A");
    std::ofstream conc_B_file = file_generator.getCSVFile("concentration_B");

    Ni_clad_Al_particle.printGridPoints(conc_A_file, ',');
    Ni_clad_Al_particle.printGridPoints(conc_B_file, ',');

    setUpKeyboardInterrupt();

	long double simulation_time = 0.0;

    try
    {
		size_t step = 0.001 / CoreShellDiffusion::delta_t;

		bool combustion_complete = false;

		for (size_t i = 0; i < MAX_ITER && !combustion_complete;)
		{
			size_t i_step = i + step;

			for (; i < i_step && !combustion_complete; i++)
			{
				Ni_clad_Al_particle.setUpEquations(temperature);
            	Ni_clad_Al_particle.solveEquations();

				simulation_time += CoreShellDiffusion::delta_t;

				combustion_complete = Ni_clad_Al_particle.isCombustionComplete();           
			}

			Ni_clad_Al_particle.printConcentrationProfileA(conc_A_file, ',', simulation_time);
            Ni_clad_Al_particle.printConcentrationProfileB(conc_B_file, ',', simulation_time);

			std::cout << "Iterations Completed : " << i << "\n";
		}
    }

    catch (InterruptException& e)
    {
        std::cout << "\nCaught signal " << e.S << "\nQuitting..." << std::endl;

		Ni_clad_Al_particle.printConcentrationProfileA(conc_A_file, ',', simulation_time);
		Ni_clad_Al_particle.printConcentrationProfileB(conc_B_file, ',', simulation_time);

        conc_A_file.close();
        conc_B_file.close();

        return 1;
    }

    conc_A_file.close();
    conc_B_file.close();

	CoreShellDiffusion::deallocateRadiusArray();
    
    return 0;
}

void parseProgramOptions(int argc, char const *argv[])
{
	ez::ezOptionParser opt;

	opt.overview	= "Diffusion of Al and Ni in a Core-Shell particle.";
	opt.syntax		= "Core-Shell-Diffusion [OPTIONS]";
	opt.example		= "Core-Shell-Diffusion --T 1500.0\n\n";
	opt.footer		= "Developed by Souritra Garai, 2021-23.\n";

	setHelpOption(opt);
	setTemperatureOption(opt);
	
	opt.parse(argc, argv);

	displayHelpOption(opt);
	
	temperature = getTemperatureOption(opt, temperature);
}