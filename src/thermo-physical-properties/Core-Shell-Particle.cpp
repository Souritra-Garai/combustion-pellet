/**
 * @file Core-Shell-Combustion-Particle.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief 
 * @version 0.1
 * @date 2021-07-08
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "thermo-physical-properties/Core-Shell-Particle.hpp"
#include "utilities/Read-Data.hpp"

#include <cmath>

const real_t Phase::sharpness_coefficient = readScalarData<real_t>("data/core-shell-particle/", "sharpness-coefficient.txt");

ArrheniusDiffusivityModel CoreShellParticle::diffusivity_model = readArrheniusDiffusivityModelParameters("data/core-shell-particle/diffusivity-parameters");

CondensedSpecies CoreShellParticle::core_species = readCondensedSpeciesData("data/core-shell-particle/core-species");
CondensedSpecies CoreShellParticle::shell_species = readCondensedSpeciesData("data/core-shell-particle/shell-species");
CondensedSpecies CoreShellParticle::product_species = readCondensedSpeciesData("data/core-shell-particle/product-species");

const real_t CoreShellParticle::overall_radius = readScalarData<real_t>("data/core-shell-particle/", "overall-radius.txt");
const real_t CoreShellParticle::core_radius    = readScalarData<real_t>("data/core-shell-particle/", "core-radius.txt");

real_t calcCoreVolume()
{
    return 4.0 * M_PI * std::pow(CoreShellParticle::core_radius, 3) / 3.0;
}

real_t calcShellVolume()
{
    return 4.0 * M_PI * (std::pow(CoreShellParticle::overall_radius, 3) - std::pow(CoreShellParticle::core_radius, 3)) / 3.0;
}

real_t calcCoreMass()
{
    return CoreShellParticle::core_species.getDensity(298.15) * calcCoreVolume();
}

real_t calcShellMass()
{
    return CoreShellParticle::shell_species.getDensity(298.15) * calcShellVolume();
}

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
