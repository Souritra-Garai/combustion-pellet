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

#define MAX_ITER 10000

Substance<float> Al(2700, 897, 26.98E-3, 239);
Substance<float> Ni(8902, 440, 58.69E-3, 90.7);
Substance<float> NiAl(5900, 717, 85.67E-3, 115, -118.4E3 / 85.67E-3);

float core_radius = 32.5E-6;
float overall_radius = 39.5E-6;

void printState(size_t iteration_number, CoreShellDiffusion<float> &particle);

int main(int argc, char const *argv[])
{
    CoreShellDiffusion<float>::setUpCoreShellCombustionParticle(
        Al, Ni, NiAl,
        overall_radius, core_radius
    );

    CoreShellDiffusion<float>::setGridSize(101);
    CoreShellDiffusion<float>::setTimeStep(0.0001);

    ArrheniusDiffusivityModel<float> diffusivity_model(2.56E-6, 102.191E3);

    CoreShellDiffusion<float> Ni_clad_Al_particle, Ni_clad_Al_particle_copy;

    FileGenerator file_generator;

    std::ofstream config_file = file_generator.getTXTFile("diffusion_config");
    std::ofstream conc_A_file = file_generator.getCSVFile("concentration_A");
    std::ofstream conc_B_file = file_generator.getCSVFile("concentration_B");
    std::ofstream state_file = file_generator.getCSVFile("mass_fractions");

    Ni_clad_Al_particle.printProperties(std::cout);
    Ni_clad_Al_particle.printProperties(config_file);
    config_file.close();

    size_t __iter = 1;

    printState(0, Ni_clad_Al_particle);
    state_file << "Al, Ni, NiAl, Sum\n";

    float temperature = 1600;
    float diffusivity = diffusivity_model.getDiffusivity(temperature);

    setUpKeyboardInterrupt();
    
    try
    {
        while (!Ni_clad_Al_particle_copy.isCombustionComplete() && __iter <= MAX_ITER)
        {
            Ni_clad_Al_particle.setUpEquations(diffusivity);
            Ni_clad_Al_particle.solveEquations();

            Ni_clad_Al_particle.copyTo(Ni_clad_Al_particle_copy);

            Ni_clad_Al_particle_copy.setUpEquations(diffusivity);
            Ni_clad_Al_particle_copy.solveEquations();

            float Y_Al = Ni_clad_Al_particle_copy.getMassFractionsCoreMaterial();
            float Y_Ni = Ni_clad_Al_particle_copy.getMassFractionsShellMaterial();
            float Y_NiAl = Ni_clad_Al_particle_copy.getMassFractionsProductMaterial();

            state_file << Y_Al << ',' << Y_Ni << ',' << Y_NiAl << ',' << Y_Al + Y_Ni + Y_NiAl << std::endl;

            Ni_clad_Al_particle_copy.printConcentrationProfileA(conc_A_file, ',');
            Ni_clad_Al_particle_copy.printConcentrationProfileB(conc_B_file, ',');

            if (__iter % 20 == 0)
            {
                printState(__iter, Ni_clad_Al_particle_copy);
            }

            __iter++;
        }
    }

    catch (InterruptException& e)
    {
        std::cout << "\nCaught signal " << e.S << std::endl;

        printState(__iter, Ni_clad_Al_particle_copy);

        float Y_Al = Ni_clad_Al_particle_copy.getMassFractionsCoreMaterial();
        float Y_Ni = Ni_clad_Al_particle_copy.getMassFractionsShellMaterial();
        float Y_NiAl = Ni_clad_Al_particle_copy.getMassFractionsProductMaterial();

        state_file << Y_Al << ',' << Y_Ni << ',' << Y_NiAl << ',' << Y_Al + Y_Ni + Y_NiAl << std::endl;

        Ni_clad_Al_particle_copy.printConcentrationProfileA(conc_A_file, ',');
        Ni_clad_Al_particle_copy.printConcentrationProfileB(conc_B_file, ',');

        conc_A_file.close();
        conc_B_file.close();
        state_file.close();

        return 1;
    }

    printState(__iter, Ni_clad_Al_particle);

    float Y_Al = Ni_clad_Al_particle_copy.getMassFractionsCoreMaterial();
    float Y_Ni = Ni_clad_Al_particle_copy.getMassFractionsShellMaterial();
    float Y_NiAl = Ni_clad_Al_particle_copy.getMassFractionsProductMaterial();

    state_file << Y_Al << ',' << Y_Ni << ',' << Y_NiAl << ',' << Y_Al + Y_Ni + Y_NiAl << std::endl;

    Ni_clad_Al_particle_copy.printConcentrationProfileA(conc_A_file, ',');
    Ni_clad_Al_particle_copy.printConcentrationProfileB(conc_B_file, ',');

    conc_A_file.close();
    conc_B_file.close();
    state_file.close();
    
    return 0;
}

void printState(size_t iteration_number, CoreShellDiffusion<float> &particle)
{
    std::cout << "Iteration # " << iteration_number;
    float Y_Al = particle.getMassFractionsCoreMaterial();
    float Y_Ni = particle.getMassFractionsShellMaterial();
    float Y_NiAl = particle.getMassFractionsProductMaterial();

    std::cout << "\tAl : " << Y_Al;
    std::cout << "\tNi : " << Y_Ni;
    std::cout << "\tNiAl : " << Y_NiAl;
    std::cout << "\tSum : " << Y_Al + Y_Ni + Y_NiAl;

    float m_Al = particle.getDiffusionMassA();
    float m_Ni = particle.getDiffusionMassB();

    std::cout << "\tAl Mass : " << m_Al;
    std::cout << "\tNi Mass : " << m_Ni;
    std::cout << "\tSum : " << m_Al + m_Ni;

    std::cout << std::endl;
}