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

void printState(real_t time, CoreShellDiffusion &particle);

int main(int argc, char const *argv[])
{
	
    CoreShellDiffusion::setUpRadiusArray();

    CoreShellDiffusion Ni_clad_Al_particle;

    real_t temperature = 1900;

    size_t step = 0.001 / CoreShellDiffusion::delta_t;

	real_t time = 0.0;

	bool combustion_complete = false;

	printState(time, Ni_clad_Al_particle);

	for (size_t i = 0; i < MAX_ITER && !combustion_complete;)
	{
		size_t i_step = i + step;

		for (; i < i_step && !combustion_complete; i++)
		{
			Ni_clad_Al_particle.setUpEquations(temperature);
			Ni_clad_Al_particle.solveEquations();

			time += CoreShellDiffusion::delta_t;

			combustion_complete = Ni_clad_Al_particle.isCombustionComplete();           
		}

		std::cout << "Iterations Completed : " << i++ << "\n";

		printState(time, Ni_clad_Al_particle);
	}

	CoreShellDiffusion::deallocateRadiusArray();
	
    return 0;
}

void printState(real_t time, CoreShellDiffusion &particle)
{
    real_t m_Al = particle.getAtomMassA();
    real_t m_Ni = particle.getAtomMassB();
    real_t m = m_Al + m_Ni;

	static const real_t initial_mass = CoreShellParticle::mass;

	std::cout << "Time : " << time << "\ts";

    std::cout << "\tAl:\t" << m_Al;
    std::cout << "\tNi:\t" << m_Ni;
    std::cout << "\tMass:\t" << m;
    std::cout << "\tDiff%:\t" << (m-initial_mass)*100/initial_mass;

    std::cout << std::endl;
}