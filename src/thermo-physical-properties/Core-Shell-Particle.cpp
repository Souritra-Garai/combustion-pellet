#include "thermo-physical-properties/Core-Shell-Particle.hpp"
#include "utilities/Read-Data.hpp"

#include <cmath> // PI and pow()

// Order of memory initialization of static variables is crucial
// as variables initialized later are dependent on variables initialized first

const real_t Phase::sharpness_coefficient = readScalarData<real_t>("data/core-shell-particle/", "sharpness-coefficient.txt");

ArrheniusDiffusivityModel CoreShellParticle::diffusivity_model = readArrheniusDiffusivityModelParameters("data/core-shell-particle/diffusivity-parameters");

CondensedSpecies CoreShellParticle::core_species = readCondensedSpeciesData("data/core-shell-particle/core-species");
CondensedSpecies CoreShellParticle::shell_species = readCondensedSpeciesData("data/core-shell-particle/shell-species");
CondensedSpecies CoreShellParticle::product_species = readCondensedSpeciesData("data/core-shell-particle/product-species");

const real_t CoreShellParticle::overall_radius = readScalarData<real_t>("data/core-shell-particle/", "overall-radius.txt");
const real_t CoreShellParticle::core_radius    = readScalarData<real_t>("data/core-shell-particle/", "core-radius.txt");

// Returns volume of core of core-shell particle in m^3
real_t calcCoreVolume()
{
    return 4.0 * M_PI * std::pow(CoreShellParticle::core_radius, 3) / 3.0;
}

// Returns volume of shell of core-shell particle in m^3
real_t calcShellVolume()
{
    return 4.0 * M_PI * (std::pow(CoreShellParticle::overall_radius, 3) - std::pow(CoreShellParticle::core_radius, 3)) / 3.0;
}

// Returns mass of core of core-shell particle in kg
real_t calcCoreMass()
{
    return CoreShellParticle::core_species.getDensity(298.15) * calcCoreVolume();
}

// Returns mass of species of core-shell particle in kg
real_t calcShellMass()
{
    return CoreShellParticle::shell_species.getDensity(298.15) * calcShellVolume();
}

// Returns mass of core-shell particle in kg
real_t calcParticleMass()
{
    return calcCoreMass() + calcShellMass();
}

const real_t CoreShellParticle::mass = calcParticleMass();

CoreShellParticle::CoreShellParticle()
{
    // Calculate and set mass fractions of the 
    // core, shell and product material
    _mass_fraction_core_material    = calcCoreMass()  / mass;
    _mass_fraction_shell_material   = calcShellMass() / mass;
    _mass_fraction_product_material = 0.0;
}
