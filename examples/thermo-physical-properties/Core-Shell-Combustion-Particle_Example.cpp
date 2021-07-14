/**
 * @file Core-Shell-Combustion-Particle_Example.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-07-09
 * 
 * @copyright Copyright (c) 2021
 * 
 */


#include <iostream>

#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"

Substance<float> A(10, 10, 10, 10);
Substance<float> B(10, 10, 10, 10);
Substance<float> AB(10, 10, 10, 10, -120);

CoreShellCombustionParticle<float> AB_particle(A, B, AB, 10, 10);

int main(int argc, char const *argv[])
{
    std::cout << "Enthalpy\t:\t" << AB_particle.getEnthalpy(350) << "\tJ / kg" << std::endl;
    AB_particle.printProperties(std::cout);
    return 0;
}
