#include <iostream>

#include "species/Aluminium.hpp"
#include "species/Nickel.hpp"
#include "species/NickelAluminide.hpp"
#include "species/Argon.hpp"

#include "thermo-physical-properties/Core-Shell-Particle.hpp"
#include "thermo-physical-properties/Thermal_Conductivity_Pellet.hpp"

long double core_radius = 32.5E-6;
long double overall_radius = 39.5E-6;

int main(int argc, char const *argv[])
{
    CoreShellParticle<long double>::setUpCoreShellParticle(
        Aluminium, Nickel, NickelAluminide,
        overall_radius, core_radius
    );

    CoreShellParticle<long double> Ni_clad_Al_particle;

    long double particle_thermal_conductivity = Ni_clad_Al_particle.getThermalConductivity(298);
    long double fluid_thermal_conductivity = Argon.getThermalConductivity(298);

    std::cout << "Particle Thermal Conductivity : " << particle_thermal_conductivity << std::endl;
    std::cout << "Fluid Thermal Conductivity : " << fluid_thermal_conductivity << std::endl;

    for (long double particle_volume_fraction = 0.001; particle_volume_fraction < 1.0; particle_volume_fraction += 0.001)
    {
        std::cout << particle_volume_fraction << ',';
        std::cout << getThermalConductivityME1(1.0 - particle_volume_fraction, fluid_thermal_conductivity, particle_thermal_conductivity) << ',';
        std::cout << getThermalConductivityEMT(1.0 - particle_volume_fraction, fluid_thermal_conductivity, particle_thermal_conductivity) << ',';
        std::cout << getThermalConductivityMEB(1.0 - particle_volume_fraction, fluid_thermal_conductivity, particle_thermal_conductivity) << ',';
        std::cout << getThermalConductivityCC( 1.0 - particle_volume_fraction, fluid_thermal_conductivity, particle_thermal_conductivity) << ',';
        std::cout << getThermalConductivityME2(1.0 - particle_volume_fraction, fluid_thermal_conductivity, particle_thermal_conductivity) << std::endl;
    }
	
    return 0;
}