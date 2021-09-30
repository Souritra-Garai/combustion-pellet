#include <iostream>
#include <omp.h>

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

#define MAX_ITER 3000

double core_radius = 32.5E-6;
double overall_radius = 39.5E-6;

double pellet_length = 6.35E-3;
double pellet_diameter = 6.35E-3;

ArrheniusDiffusivityModel<double> Alawieh_diffusivity(2.56E-6, 102.191E3);
ArrheniusDiffusivityModel<double> Du_diffusivity(9.54E-8, 26E3);

void printState(size_t iteration_number, PelletFlamePropagation<double> &pellet);

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

    // omp_set_num_threads(1);

    PelletFlamePropagation<double> combustion_pellet(0.5);
    combustion_pellet.setDiffusivityModel(Du_diffusivity);

    combustion_pellet.initializePellet();

    FileGenerator file_generator;

    std::ofstream temperature_file = file_generator.getCSVFile("temperature");
    std::ofstream config_file = file_generator.getTXTFile("combustion_config");

    combustion_pellet.printProperties(config_file);
    config_file.close();

    std::cout << std::endl;

    size_t __iter = 1;

    printState(0, combustion_pellet);

    setUpKeyboardInterrupt();
    
    try
    {
        while (!combustion_pellet.isCombustionComplete() && __iter <= MAX_ITER)
        {
            combustion_pellet.setUpEquations();
            combustion_pellet.solveEquations();

            combustion_pellet.printTemperatureProfile(temperature_file, ',');

            if (__iter % 30 == 0)
            {
                std::cout << "Iteration # " << __iter << std::endl;
                // printState(__iter, combustion_pellet);
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