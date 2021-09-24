#include <iostream>

#include "thermo-physical-properties/Packed-Pellet.hpp"
#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"

#include "pde-problems/Core-Shell-Diffusion.hpp"

#include "substances/Argon.hpp"
#include "substances/Aluminium.hpp"
#include "substances/Nickel.hpp"
#include "substances/NickelAluminide.hpp"

double core_radius = 32.5E-6;
double overall_radius = 39.5E-6;

double pellet_length = 6.35E-3;
double pellet_diameter = 6.35E-3;

int main(int argc, char const *argv[])
{
    // CoreShellCombustionParticle<double>::setUpCoreShellCombustionParticle(
    //     Aluminium, Nickel, NickelAluminide,
    //     overall_radius, core_radius
    // );

    CoreShellDiffusion<double>::setUpCoreShellCombustionParticle(
        Aluminium, Nickel, NickelAluminide,
        overall_radius, core_radius
    );

    CoreShellDiffusion<double> particle;

    PackedPellet<double>::setPelletDimensions(pellet_length, pellet_diameter);
    PackedPellet<double>::setAmbientHeatLossParameters(19.68, 0.25);
    PackedPellet<double>::setTemperatureParameters(933.0, 298.0);
    PackedPellet<double>::setDegassingFluid(Argon);

    PackedPellet<double> pellet(0.5);
    
    pellet.printProperties(std::cout);

    std::cout << "\n\nEnthalpy\t:\t" << pellet.getInternalEnergy(&particle, 298.15) << "\tJ/kg" << std::endl;  
    return 0;
}