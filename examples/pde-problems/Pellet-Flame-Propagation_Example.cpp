#include <iostream>
#include <omp.h>

#include "thermo-physical-properties/Substance.hpp"
#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.hpp"

#include "pde-problems/Core-Shell-Diffusion.hpp"
#include "pde-problems/Pellet-Flame-Propagation.hpp"

#include "utilities/File_Generator.hpp"
#include "utilities/Keyboard_Interrupt.hpp"

#define MAX_ITER 3000

Substance<float> Al(2700, 1060, 26.98E-3, 220);
Substance<float> Ni(8902, 440, 58.69E-3, 66);
Substance<float> NiAl(5900, 717, 85.67E-3, 115, -118.4E3 / 85.67E-3);

float core_radius = 32.5E-6;
float overall_radius = 39.5E-6;

Substance<float> Ar(0.5, 520, 39.95E-3, 0.3, 0);

float pellet_length = 6.35E-3;
float pellet_diameter = 6.35E-3;

ArrheniusDiffusivityModel<float> Alawieh_diffusivity(2.56E-6, 102.191E3);
ArrheniusDiffusivityModel<float> Du_diffusivity(9.54E-8, 26E3);

void printState(size_t iteration_number, PelletFlamePropagation<float> &pellet);

int main(int argc, char const *argv[])
{
    CoreShellDiffusion<float>::setUpCoreShellCombustionParticle(
        Al, Ni, NiAl,
        overall_radius, core_radius
    );

    CoreShellDiffusion<float>::setGridSize(1001);
    // CoreShellDiffusion<float>::setTimeStep(0.0001);

    PelletFlamePropagation<float>::setPelletDimensions(pellet_length, pellet_diameter);
    PelletFlamePropagation<float>::setAmbientHeatLossParameters(0, 0);
    PelletFlamePropagation<float>::setTemperatureParameters(933, 298);
    PelletFlamePropagation<float>::setDegassingFluid(Ar);

    PelletFlamePropagation<float>::setGridSize(101);
    PelletFlamePropagation<float>::setTimeStep(0.0001);
    PelletFlamePropagation<float>::setInfinitesimalChangeTemperature(0.1);
    PelletFlamePropagation<float>::setInitialIgnitionParameters(1911, 0.1 * pellet_length);

    // omp_set_num_threads(1);

    PelletFlamePropagation<float> combustion_pellet(0.5);
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

void printState(size_t iteration_number, PelletFlamePropagation<float> &pellet)
{
    std::cout << "Iteration # " << iteration_number << "\t";
    
    pellet.printTemperatureProfile(std::cout);

    std::cout << std::endl;
}