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

#include "pde-problems/Core-Shell-Diffusion.hpp"

#define MAX_ITER 1E6

void printState(long double time, CoreShellDiffusion<long double> &particle);

int main(int argc, char const *argv[])
{
	
    CoreShellDiffusion<long double>::setUpRadiusArray();

    CoreShellDiffusion<long double> Ni_clad_Al_particle;

    long double temperature = 1900;

    size_t step = 0.001 / CoreShellDiffusion<long double>::delta_t;

	long double time = 0.0;

	bool combustion_complete = false;

	printState(time, Ni_clad_Al_particle);

	for (size_t i = 0; i < MAX_ITER && !combustion_complete;)
	{
		size_t i_step = i + step;

		for (; i < i_step && !combustion_complete; i++)
		{
			Ni_clad_Al_particle.setUpEquations(temperature);
			Ni_clad_Al_particle.solveEquations();

			time += CoreShellDiffusion<long double>::delta_t;

			combustion_complete = Ni_clad_Al_particle.isCombustionComplete();           
		}

		std::cout << "Iterations Completed : " << i++ << "\n";

		printState(time, Ni_clad_Al_particle);
	}

	CoreShellDiffusion<long double>::deallocateRadiusArray();
	
    return 0;
}

void printState(long double time, CoreShellDiffusion<long double> &particle)
{
    long double Y_Al = particle.getMassFractionsCoreMaterial();
    long double Y_Ni = particle.getMassFractionsShellMaterial();
    long double Y_NiAl = particle.getMassFractionsProductMaterial();

	std::cout << "Time : " << time << "\ts";

    std::cout << "\tAl\t:\t" << Y_Al;
    std::cout << "\tNi\t:\t" << Y_Ni;
    std::cout << "\tNiAl\t:\t" << Y_NiAl;
    std::cout << "\tSum\t:\t" << Y_Al + Y_Ni + Y_NiAl;

    std::cout << std::endl;
}