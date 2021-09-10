/**
 * @file Core-Shell-Diffusion_Example.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Program to test functions of CoreShellDiffusion class
 * @version 0.1
 * @date 2021-07-16
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>
#include <omp.h>

#include "thermo-physical-properties/Substance.hpp"
#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.hpp"
#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"
#include "pde-problems/Core-Shell-Diffusion.hpp"

#include "utilities/Keyboard_Interrupt.hpp"
#include "utilities/File_Generator.hpp"

#define MAX_ITER 5000
#define Dt 0.0001

Substance<float> Al(2700, 897, 26.98E-3, 239);
Substance<float> Ni(8902, 440, 58.69E-3, 90.7);
Substance<float> NiAl(5900, 717, 85.67E-3, 115, -118.4E3 / 85.67E-3);

float core_radius = 32.5E-6;
float overall_radius = 39.5E-6;

ArrheniusDiffusivityModel<float> Alawieh_diffusivity(2.56E-6, 102.191E3);
ArrheniusDiffusivityModel<float> Du_diffusivity(9.54E-8, 26E3);

void printState(size_t iteration_number, CoreShellDiffusion<float> &particle);

int main(int argc, char const *argv[])
{
    CoreShellDiffusion<float>::setUpCoreShellCombustionParticle(
        Al, Ni, NiAl,
        overall_radius, core_radius
    );

    CoreShellDiffusion<float>::setGridSize(1001);
    CoreShellDiffusion<float>::setTimeStep(Dt);

    CoreShellDiffusion<float> Ni_clad_Al_particle;

    FileGenerator file_generator;

    std::ofstream config_file = file_generator.getTXTFile("diffusion_config");
    std::ofstream conc_A_file = file_generator.getCSVFile("concentration_A");
    std::ofstream conc_B_file = file_generator.getCSVFile("concentration_B");

    Ni_clad_Al_particle.printGridPoints(conc_A_file, ',');
    Ni_clad_Al_particle.printGridPoints(conc_B_file, ',');

    Ni_clad_Al_particle.printProperties(config_file);
    config_file.close();

    size_t __iter = 1;

    float temperature = 1400;
    float diffusivity = Alawieh_diffusivity.getDiffusivity(temperature);

    setUpKeyboardInterrupt();
    
    try
    {
        while (!Ni_clad_Al_particle.isCombustionComplete() && __iter <= MAX_ITER)
        {
            Ni_clad_Al_particle.setUpEquations(diffusivity);
            Ni_clad_Al_particle.solveEquations();
            
            Ni_clad_Al_particle.printConcentrationProfileA(conc_A_file, ',', Dt * __iter);
            Ni_clad_Al_particle.printConcentrationProfileB(conc_B_file, ',', Dt * __iter);

            if (__iter % 100 == 0)
            {
                std::cout << "Iteration # " << __iter << std::endl;
            }

            __iter++;
        }
    }

    catch (InterruptException& e)
    {
        std::cout << "\nCaught signal " << e.S << std::endl;

        conc_A_file.close();
        conc_B_file.close();

        return 1;
    }

    conc_A_file.close();
    conc_B_file.close();
    
    return 0;
}

void printState(size_t iteration_number, CoreShellDiffusion<float> &particle)
{
    float Y_Al = particle.getMassFractionsCoreMaterial();
    float Y_Ni = particle.getMassFractionsShellMaterial();
    float Y_NiAl = particle.getMassFractionsProductMaterial();

    std::cout << "Iteration # " << iteration_number;

    std::cout << "\tAl\t:\t" << Y_Al;
    std::cout << "\tNi\t:\t" << Y_Ni;
    std::cout << "\tNiAl\t:\t" << Y_NiAl;
    std::cout << "\tSum\t:\t" << Y_Al + Y_Ni + Y_NiAl;

    std::cout << std::endl;
}