#include "thermo-physical-properties/Packed-Pellet.hpp"

#include "utilities/Read-Data.hpp"

const real_t PackedPellet::length = readScalarData<real_t>("data/pellet", "length.txt");
const real_t PackedPellet::diameter = readScalarData<real_t>("data/pellet", "diameter.txt");

const real_t PackedPellet::convective_heat_transfer_coefficient_curved_surface	= readScalarData<real_t>("data/pellet", "convective-heat-transfer-coefficient-curved-surface.txt");
const real_t PackedPellet::convective_heat_transfer_coefficient_flat_surface	= readScalarData<real_t>("data/pellet", "convective-heat-transfer-coefficient-flat-surface.txt");
const real_t PackedPellet::radiative_emissivity = readScalarData<real_t>("data/pellet", "radiative-emissivity.txt");

const real_t PackedPellet::ambient_pressure = readScalarData<real_t>("data/pellet", "ambient-pressure.txt");
const real_t PackedPellet::ambient_temperature = readScalarData<real_t>("data/pellet", "ambient-temperature.txt");

IdealGas PackedPellet::interstitial_gas = readIdealGasData("data/pellet/interstitial-gas");

// Input temperature T in K
// Returns density in kg/m^3
real_t calcOverallParticleDensity(real_t particle_volume_fractions, real_t temperature = 298.15)
{
	return particle_volume_fractions * CoreShellParticle().getDensity(temperature);
}

PackedPellet::PackedPellet(
	real_t particle_volume_fractions
) : overall_particle_density(calcOverallParticleDensity(particle_volume_fractions, ambient_temperature)),
	interstitial_volume_fractions(1.0 - particle_volume_fractions)
{ ; }
