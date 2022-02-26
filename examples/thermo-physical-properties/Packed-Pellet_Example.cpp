#include <iostream>

#include "thermo-physical-properties/Packed-Pellet.hpp"
#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"

#include "substances/Argon.hpp"
#include "substances/Aluminium.hpp"
#include "substances/Nickel.hpp"
#include "substances/NickelAluminide.hpp"

long double core_radius = 32.5E-6;
long double overall_radius = 39.5E-6;

long double pellet_length = 6.35E-3;
long double pellet_diameter = 6.35E-3;

int main(int argc, char const *argv[])
{
    CoreShellCombustionParticle<long double>::setUpCoreShellCombustionParticle(
        Aluminium, Nickel, NickelAluminide,
        overall_radius, core_radius
    );

    PackedPellet<long double>::setPelletDimensions(pellet_length, pellet_diameter);
    PackedPellet<long double>::setAmbientHeatLossParameters(19.68, 0, 0.25);
    PackedPellet<long double>::setTemperatureParameters(933.0, 298.0);
    PackedPellet<long double>::setDegassingFluid(Argon);

    PackedPellet<long double> pellet(0.5);
    
    pellet.printProperties(std::cout);

    return 0;
}