/**
 * @file Packed-Pellet.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Implementation details for the class PackedPellet
 * @version 0.1
 * @date 2021-07-31
 * 
 * @copyright Copyright (c) 2021
 * 
 */

// For definition of PackedPellet class
#include "thermo-physical-properties/Packed-Pellet.hpp"

#include "utilities/Read-Data.hpp"

/*****************************************************************************************************/
// Instantiation of static members
// Only one copy of these variables are shared across all classes

// Length of pellet in \f$ m \f$
template<typename real_t> const real_t PackedPellet<real_t>::length = readScalarData<real_t>("data/pellet", "length.txt");
// Diameter of pellet in \f$ m \f$
template<typename real_t> const real_t PackedPellet<real_t>::diameter = readScalarData<real_t>("data/pellet", "diameter.txt");

// Ambient heat loss parameters
// Convective heat transfer coefficient
template<typename real_t> const real_t PackedPellet<real_t>::convective_heat_transfer_coefficient_curved_surface	= readScalarData<real_t>("data/pellet", "convective-heat-transfer-coefficient-curved-surface.txt");
template<typename real_t> const real_t PackedPellet<real_t>::convective_heat_transfer_coefficient_flat_surface		= readScalarData<real_t>("data/pellet", "convective-heat-transfer-coefficient-flat-surface.txt");
// Emissivity
template<typename real_t> const real_t PackedPellet<real_t>::radiative_emissivity = readScalarData<real_t>("data/pellet", "radiative-emissivity.txt");

// Ambient temperature in K
template<typename real_t> const real_t PackedPellet<real_t>::ambient_temperature = readScalarData<real_t>("data/pellet", "ambient-temperature.txt");

// Substance filling the voids in the packed pellet
template<typename real_t> IdealGas<real_t> PackedPellet<real_t>::interstitial_gas = readIdealGasData<real_t>("data/pellet/interstitial-gas");

/*****************************************************************************************************/
// Definitions of Static Member Functions

template<typename real_t>
real_t PackedPellet<real_t>::calcOverallParticleDensity(real_t particle_volume_fractions)
{
    // Get the density of a defualt Core-Shell Combustion Particle instance
    // This assumes the static members of the CoreShellParticle class
    // have been initialized

    return particle_volume_fractions * CoreShellParticle<real_t>().getDensity(298.15);
}

/*****************************************************************************************************/
// Constructors and destructors

template<typename real_t>
PackedPellet<real_t>::PackedPellet(
    real_t particle_volume_fractions
) : // Mem Initialization list for initialising the constant members
    overall_particle_density(calcOverallParticleDensity(particle_volume_fractions)),
	interstitial_volume_fractions(1.0 - particle_volume_fractions)
{ 
    // Nothing to do
    ; 
}

/*****************************************************************************************************/
template class PackedPellet<float>;
template class PackedPellet<double>;
template class PackedPellet<long double>;