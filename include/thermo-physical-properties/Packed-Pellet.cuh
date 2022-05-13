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

	__host__ void printConfiguration(std::ostream &output_stream, double phi)
	{
		double length;
		cudaMemcpyFromSymbol(&length, PackedPellet::length, sizeof(double));
		double diameter;
		cudaMemcpyFromSymbol(&diameter, PackedPellet::diameter, sizeof(double));

		double convective_heat_transfer_coefficient_curved_surface;
		cudaMemcpyFromSymbol(
			&convective_heat_transfer_coefficient_curved_surface,
			PackedPellet::convective_heat_transfer_coefficient_curved_surface,
			sizeof(double)
		);
		
		double convective_heat_transfer_coefficient_flat_surface;
		cudaMemcpyFromSymbol(
			&convective_heat_transfer_coefficient_flat_surface,
			PackedPellet::convective_heat_transfer_coefficient_flat_surface,
			sizeof(double)
		);

		double radiative_emissivity;
		cudaMemcpyFromSymbol(
			&radiative_emissivity,
			PackedPellet::radiative_emissivity,
			sizeof(double)
		);

		double ambient_temperature;
		cudaMemcpyFromSymbol(
			&ambient_temperature,
			PackedPellet::ambient_temperature,
			sizeof(double)
		);

		double ignition_temperature;
		cudaMemcpyFromSymbol(
			&ignition_temperature,
			PackedPellet::ignition_temperature,
			sizeof(double)
		);

		output_stream << "Pellet Properties\n\n";

		output_stream << "Particle Volume Fractions : " << phi << '\n';

		output_stream << "Length : " << length << " m\nDiameter : " << diameter << " m\n";

		output_stream << "Convective Heat Trasfer Coefficients for Heat Loss\n";
		output_stream << "\tCurved Surface : " << convective_heat_transfer_coefficient_curved_surface << " W / m2-K\n";
		output_stream << "\tFlat Surface : " << convective_heat_transfer_coefficient_flat_surface << " W / m2-K\n";

		output_stream << "Radiative Emissivity : " << radiative_emissivity << "\n";

		output_stream << "Ambient Temperature : " << ambient_temperature << " K\n";
		output_stream << "Ignition Temperature : " << ignition_temperature << " K\n\n";
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