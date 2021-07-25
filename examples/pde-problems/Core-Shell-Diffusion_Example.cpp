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

#include "thermo-physical-properties/Substance.hpp"
#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"
#include "pde-problems/Core-Shell-Diffusion.hpp"

#include "utilities/Keyboard_Interrupt.hpp"

#define MAX_ITER 10000

Substance<float> Al(2700, 897, 26.98E-3, 239);
Substance<float> Ni(8902, 440, 58.69E-3, 90.7);
Substance<float> NiAl(5900, 717, 85.67E-3, 115, -118.4E3 / 85.67E-3);

float core_radius = 32.5E-6;
float overall_radius = 39.5E-6;

int main(int argc, char const *argv[])
{
    CoreShellDiffusion<float>::setUpCoreShellCombustionParticle(
        Al, Ni, NiAl,
        overall_radius, core_radius
    );

    CoreShellDiffusion<float>::setGridSize(1001);
    CoreShellDiffusion<float>::setTimeStep(0.000001);

    CoreShellDiffusion<float>::setDiffusivityParameters(2.56E-6, 102.191E3);

    CoreShellDiffusion<float> Ni_clad_Al_particle;

    Ni_clad_Al_particle.printProperties(std::cout);

    size_t __iter = 1;

    std::cout << "Iteration # " << 0;
    std::cout << "\tAl : " << Ni_clad_Al_particle.getMassFractionsCoreMaterial();
    std::cout << "\tNi : " << Ni_clad_Al_particle.getMassFractionsShellMaterial();
    std::cout << "\tNiAl : " << Ni_clad_Al_particle.getMassFractionsProductMaterial();
    std::cout << "\tMass : " << Ni_clad_Al_particle.numCalcMass();
    std::cout << std::endl;

    float temperature = 1600;

    setUpKeyboardInterrupt();

    try
    {
        while (!Ni_clad_Al_particle.isCombustionComplete() && __iter <= MAX_ITER)
        {
            Ni_clad_Al_particle.setUpEquations(temperature);
            Ni_clad_Al_particle.solveEquations();

            if (__iter % 100 == 0) {
                std::cout << "Iteration # " << __iter;
                std::cout << "\tAl : " << Ni_clad_Al_particle.getMassFractionsCoreMaterial();
                std::cout << "\tNi : " << Ni_clad_Al_particle.getMassFractionsShellMaterial();
                std::cout << "\tNiAl : " << Ni_clad_Al_particle.getMassFractionsProductMaterial();
                std::cout << "\tMass : " << Ni_clad_Al_particle.numCalcCoreMass();
                std::cout << "\tMass : " << Ni_clad_Al_particle.numCalcShellMass();
                std::cout << std::endl;

                // Ni_clad_Al_particle.printConcentrationProfiles(std::cout);
            }

            __iter++;
        }

        std::cout << "Iteration # " << __iter;
        std::cout << "\tAl : " << Ni_clad_Al_particle.getMassFractionsCoreMaterial();
        std::cout << "\tNi : " << Ni_clad_Al_particle.getMassFractionsShellMaterial();
        std::cout << "\tNiAl : " << Ni_clad_Al_particle.getMassFractionsProductMaterial();
        std::cout << "\tMass : " << Ni_clad_Al_particle.numCalcCoreMass();
        std::cout << "\tMass : " << Ni_clad_Al_particle.numCalcShellMass();
        std::cout << std::endl;

        // Ni_clad_Al_particle.printConcentrationProfiles(std::cout);
    }

    catch (InterruptException& e)
    {
        std::cout << "\nCaught signal " << e.S << std::endl;

        std::cout << "Iteration # " << __iter;
        std::cout << "\tAl : " << Ni_clad_Al_particle.getMassFractionsCoreMaterial();
        std::cout << "\tNi : " << Ni_clad_Al_particle.getMassFractionsShellMaterial();
        std::cout << "\tNiAl : " << Ni_clad_Al_particle.getMassFractionsProductMaterial();
        std::cout << "\tMass : " << Ni_clad_Al_particle.numCalcCoreMass();
        std::cout << "\tMass : " << Ni_clad_Al_particle.numCalcShellMass();
        std::cout << std::endl;

        return 1;
    }
    
    return 0;
}
