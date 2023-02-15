#include <iostream>

#include "pde-problems/Core-Shell-Diffusion.hpp"
#include "pde-problems/Pellet-Flame-Propagation.hpp"

#include "utilities/File-Generator.hpp"
#include "utilities/Keyboard-Interrupt.hpp"

#define MAX_ITER 1E8

void printState(size_t iteration_number, PelletFlamePropagation &pellet);

int main(int argc, char const *argv[])
{
    CoreShellDiffusion::setUpRadiusArray();

    PelletFlamePropagation combustion_pellet(0.8);
   
    combustion_pellet.initializePellet(1900, 0.1);

    FileGenerator file_generator;

    std::ofstream temperature_file = file_generator.getCSVFile("temperature");
	combustion_pellet.printGridPoints(temperature_file, ',');
    
    size_t __iter = 1;

    printState(0, combustion_pellet);

    setUpKeyboardInterrupt();

    try
    {
        while (!combustion_pellet.isCombustionComplete() && __iter <= MAX_ITER)
        {
            combustion_pellet.setUpEquations();
            combustion_pellet.solveEquations();

            if (__iter % 1000 == 0)
            {
                printState(__iter, combustion_pellet);
	            combustion_pellet.printTemperatureProfile(temperature_file, ',');
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

	CoreShellDiffusion::deallocateRadiusArray();

    return 0;
}

void printState(size_t iteration_number, PelletFlamePropagation &pellet)
{
	std::cout << "Iteration # " << iteration_number << "\t";

	pellet.printTemperatureProfile(std::cout);

	std::cout << std::endl;
}