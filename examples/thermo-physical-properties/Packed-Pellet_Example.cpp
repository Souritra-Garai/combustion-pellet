#include <iostream>

#include "thermo-physical-properties/Substance.hpp"
#include "thermo-physical-properties/Packed-Pellet.hpp"
#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"

#include "pde-problems/Core-Shell-Diffusion.hpp"

Substance<float> Al(2700, 897, 26.98, 239);
Substance<float> Ni(8908, 440, 58.69, 90.7);
Substance<float> NiAl(5900, 717, 85.675, 115, -118.4);

float core_radius = 32.5E-6;
float overall_radius = 39.5E-6;

Substance<float> Ar(0.5, 520, 39.95, 0.3, 0);

float pellet_length = 6.35E-3;
float pellet_diameter = 6.35E-3;

int main(int argc, char const *argv[])
{
    CoreShellCombustionParticle<float>::setUpCoreShellCombustionParticle(
        Al, Ni, NiAl,
        overall_radius, core_radius
    );

    // CoreShellDiffusion<float>::setUpCoreShellCombustionParticle(
    //     Al, Ni, NiAl,
    //     overall_radius, core_radius
    // );

    CoreShellDiffusion<float> particle;

    PackedPellet<float>::setPelletDimensions(pellet_length, pellet_diameter);
    PackedPellet<float>::setAmbientHeatLossProperties(19.68, 0.25, 298);
    PackedPellet<float>::setDegassingFluid(Ar);

    PackedPellet<float> pellet(0.5);
    
    pellet.printProperties(std::cout);

    std::cout << "\n\nEnthalpy\t:\t" << pellet.getEnthalpy(&particle, 373) << "\tJ/kg" << std::endl;   

    return 0;
}