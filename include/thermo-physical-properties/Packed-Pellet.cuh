#ifndef __PACKED_PELLET__
#define __PACKED_PELLET__

#include <ostream>

#include "thermo-physical-properties/Ideal_Gas.cuh"
#include "thermo-physical-properties/Species.cuh"
#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "thermo-physical-properties/Thermal_Conductivity_Pellet.cuh"

namespace PackedPellet
{
	__device__ double length;
	__device__ double diameter;

	__device__ double convective_heat_transfer_coefficient_curved_surface;
	__device__ double convective_heat_transfer_coefficient_flat_surface;
	__device__ double radiative_emissivity;

	__device__ double ambient_temperature;
	__device__ double ignition_temperature;

	__device__ IdealGas *degassing_fluid;

	__device__ void setPelletDimensions(double length, double diameter)
	{
		PackedPellet::length = length;
		PackedPellet::diameter = diameter;
	}

	__device__ void setHeatLossParameters(
		double convective_heat_transfer_coefficient_curved_surface,
		double convective_heat_transfer_coefficient_flat_surface,
		double radiative_emissivity
	) {
		PackedPellet::convective_heat_transfer_coefficient_curved_surface = convective_heat_transfer_coefficient_curved_surface;
		PackedPellet::convective_heat_transfer_coefficient_flat_surface = convective_heat_transfer_coefficient_flat_surface;
		PackedPellet::radiative_emissivity = radiative_emissivity;
	}

	__device__ void setTemperatureParameters(
		double ambient_temperature,
		double ignition_temperature
	) {
		PackedPellet::ambient_temperature = ambient_temperature;
		PackedPellet::ignition_temperature = ignition_temperature;
	}

	__device__ void setDegassingFluid(
		IdealGas *degassing_fluid
	) {
		PackedPellet::degassing_fluid = degassing_fluid;
	}

	class Pellet
	{
		protected:
	
			double _degassing_fluid_volume_fractions;
			double _overall_particle_density;

		public:

			__device__ Pellet(double particle_volume_fractions)
			{
				_degassing_fluid_volume_fractions = 1.0 - particle_volume_fractions;
				
				_overall_particle_density = particle_volume_fractions * CoreShellParticle::Particle().getDensity(ambient_temperature);
			}

			__device__ __forceinline__ double getThermalConductivity(CoreShellParticle::Particle *particle, double temperature)
			{
				return PelletThermalConductivity::getThermalConductivityMEB(
					_degassing_fluid_volume_fractions,
					degassing_fluid->getThermalConductivity(temperature),
					particle->getThermalConductivity(temperature)
				);
			}
	};
} // namespace PackedPellet

#endif