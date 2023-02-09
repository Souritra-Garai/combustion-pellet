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

// Required pow() function and value of \f$ \pi \f$
#include <math.h>

template<typename real_t> const real_t Phase<real_t>::sharpness_coefficient = readScalarData<real_t>("data/core-shell-particle/", "sharpness-coefficient.txt");

// Instantiate the static private members of the CoreShellParticle class
// Creating the Substance members
template<typename real_t> CondensedSpecies<real_t> CoreShellParticle<real_t>::_core_species = readCondensedSpeciesData<real_t>("data/core-shell-particle/core-species");
template<typename real_t> CondensedSpecies<real_t> CoreShellParticle<real_t>::_shell_species = readCondensedSpeciesData<real_t>("data/core-shell-particle/shell-species");
template<typename real_t> CondensedSpecies<real_t> CoreShellParticle<real_t>::_product_species = readCondensedSpeciesData<real_t>("data/core-shell-particle/product-species");
// Creating the real_t members
template<typename real_t> const real_t CoreShellParticle<real_t>::_overall_radius = readScalarData<real_t>("data/core-shell-particle/", "overall-radius.txt");
template<typename real_t> const real_t CoreShellParticle<real_t>::_core_radius    = readScalarData<real_t>("data/core-shell-particle/", "core-radius.txt");
template<typename real_t> const real_t CoreShellParticle<real_t>::_mass = CoreShellParticle<real_t>::calcParticleMass();

template<typename real_t> ArrheniusDiffusivityModel<real_t> CoreShellParticle<real_t>::_diffusivity_model = readArrheniusDiffusivityModelParameters<real_t>("data/core-shell-particle/diffusivity-parameters");

// Calculate and return volume of sphere 
// with static member _core_radius as the radius
template<typename real_t>
real_t CoreShellParticle<real_t>::calcCoreVolume()
{
    // Volume of sphere = \f$ \frac{4}{3} \pi r^3
    return 4.0 * M_PI * pow(_core_radius, 3) / 3.0;
}

// Calculate and return volume fo a shell
// with static members _overall_radius and _core_radius
// as the outer and inner radii respectively
template<typename real_t>
real_t CoreShellParticle<real_t>::calcShellVolume()
{
    // Volume of spherical shell = \f$ \frac{4}{3} \pi \left( r_{outer}^3 - r_{inner}^3 \right) \f$
    return 4.0 * M_PI * (pow(_overall_radius, 3) - pow(_core_radius, 3)) / 3.0;
}

// Calculate and return mass of uniform sphere
// with density of core material
template<typename real_t>
real_t CoreShellParticle<real_t>::calcCoreMass()
{
    // Mass of core = Density of core material \f$ \times \f$ Volume of core
    return _core_species.getDensity(298.15) * calcCoreVolume();
}

// Calculate and return mass of spherical shell
// with uniform density equal to that of shell material
template<typename real_t>
real_t CoreShellParticle<real_t>::calcShellMass()
{
    // Mass of shell = Density of shell material \f$ \times \f$ Volume of shell
    return _shell_species.getDensity(298.15) * calcShellVolume();
}

// Sum and return the masses of core part and the shell part
template<typename real_t>
real_t CoreShellParticle<real_t>::calcParticleMass()
{
    // Mass of Core-Shell Particle = Mass of core + Mass of shell
    return calcCoreMass() + calcShellMass();
}

// Create an instance of CoreShellParticle
// and initialize the mass fraction variables
template<typename real_t>
CoreShellParticle<real_t>::CoreShellParticle()
{
    // Calculate and set mass fractions of the 
    // core, shell and product material
    _mass_fraction_core_material    = calcCoreMass()  / _mass;
    _mass_fraction_shell_material   = calcShellMass() / _mass;
    _mass_fraction_product_material = 0.0;
}

template class Phase<float>;
template class Phase<double>;
template class Phase<long double>;

template class CoreShellParticle<float>;
template class CoreShellParticle<double>;
template class CoreShellParticle<long double>;
