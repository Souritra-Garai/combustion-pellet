#include <iostream>
#include <chrono>

#include "pde-problems/Core-Shell-Diffusion.hpp"
#include "pde-problems/Pellet-Flame-Propagation.hpp"

#include "utilities/File-Generator.hpp"
#include "utilities/Program-Options.hpp"
#include "utilities/Keyboard-Interrupt.hpp"

#define MAX_ITER 1E6

long double phi = 0.7;

long double initial_ignition_temperature = 1500;
long double initial_ignition_length_fraction = 0.1;

void parseProgramOptions(int argc, char const *argv[]);

int main(int argc, char const *argv[])
{
	parseProgramOptions(argc, argv);

    CoreShellDiffusion<long double>::setUpRadiusArray();

	PelletFlamePropagation<long double> combustion_pellet(phi);

    combustion_pellet.initializePellet(
		initial_ignition_temperature,
		initial_ignition_length_fraction
	);

    FileGenerator file_generator;

    std::ofstream temperature_file = file_generator.getCSVFile("temperature");

    std::cout << "Initialized Pellet. Starting iterations.\nPress Ctrl+C to stop...\n\n";

    setUpKeyboardInterrupt();

	size_t i = 0;
    
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    try
    {
		combustion_pellet.printGridPoints(temperature_file, ',');

		size_t step = 0.001 / PelletFlamePropagation<long double>::delta_t;

		bool combustion_not_complete = true;

		for (; i < MAX_ITER && combustion_not_complete;)
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

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

	std::ofstream time_file = file_generator.getTXTFile("runtime");

	time_file << "Time difference\t= " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << " [ms]" << std::endl;
	time_file << "Time per iteration\t= " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / i << " [ms]" << std::endl;
    time_file.close();

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
	opt.example		= "Pellet-Flame-Propagation -phi 0.68\n\n";
	opt.footer		= "Developed by Souritra Garai, 2021-23.\n";

	setHelpOption(opt);
	setPhiOption(opt);
	setIgnitionTemperatureOption(opt);
	setIgnitionLengthOption(opt);

	opt.parse(argc, argv);

	displayHelpOption(opt);

	phi = getPhiOption(opt, phi);

	initial_ignition_temperature = getIgnitionTemperatureOption(opt, initial_ignition_temperature);
	initial_ignition_length_fraction = getIgnitionLengthOption(opt, initial_ignition_length_fraction);
}