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
#include <signal.h>

#include "thermo-physical-properties/Substance.hpp"
#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"
#include "pde-problems/Core-Shell-Diffusion.hpp"

#define MAX_ITER 10000

class InterruptException : public std::exception
{
public:
  InterruptException(int s) : S(s) {}
  int S;
};

void sig_to_exception(int s)
{
  throw InterruptException(s);
}

Substance<float> Al(2700, 897, 26.98, 239);
Substance<float> Ni(8902, 440, 58.69, 90.7);
Substance<float> NiAl(5900, 717, 85.675, 115, -118.4);

float core_radius = 32.5E-6;
float overall_radius = 39.5E-6;

int main(int argc, char const *argv[])
{
    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = sig_to_exception;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, NULL);

    CoreShellDiffusion<float>::setUpCoreShellCombustionParticle(
        Al, Ni, NiAl,
        overall_radius, core_radius
    );

    CoreShellDiffusion<float>::setGridSize(101);
    CoreShellDiffusion<float>::setTimeStep(0.00001);

    CoreShellDiffusion<float>::setDiffusivityParameters(2.56E-6, 102.191E3);

    CoreShellDiffusion<float> Ni_clad_Al_particle;

    Ni_clad_Al_particle.printProperties(std::cout);

    size_t __iter = 1;

    float temperature = 1600;

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
                std::cout << std::endl;

                // Ni_clad_Al_particle.printConcentrationProfiles(std::cout);
            }

            __iter++;
        }

        std::cout << "Iteration # " << __iter;
        std::cout << "\tAl : " << Ni_clad_Al_particle.getMassFractionsCoreMaterial();
        std::cout << "\tNi : " << Ni_clad_Al_particle.getMassFractionsShellMaterial();
        std::cout << "\tNiAl : " << Ni_clad_Al_particle.getMassFractionsProductMaterial();
        std::cout << std::endl;

        Ni_clad_Al_particle.printConcentrationProfiles(std::cout);
    }

    catch (InterruptException& e)
    {
        std::cout << "\nCaught signal " << e.S << std::endl;
        return 1;
    }
    
    return 0;
}
