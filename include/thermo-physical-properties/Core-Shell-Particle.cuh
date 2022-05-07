#ifndef __CORE_SHELL_PARTICLE__
#define __CORE_SHELL_PARTICLE__

#include "thermo-physical-properties/Species.cuh"

#include <ostream>

namespace CoreShellParticle
{
	__device__ Species *core_material;
	__device__ Species *shell_material;
	__device__ Species *product_material;

	__device__ double overall_radius;
	__device__ double core_radius;

	__device__ double mass;
	__device__ double core_mass;
	__device__ double shell_mass;

	__device__ void setParticleMass(
		double core_radius,
		double overall_radius)
	{
		CoreShellParticle::core_radius = core_radius;
		CoreShellParticle::overall_radius = overall_radius;

		core_mass = (4.0 / 3.0) * 3.14 * pow(core_radius, 3) * core_material->getDensity(298.15);
		shell_mass = (4.0 / 3.0) * 3.14 * (pow(overall_radius, 3) - pow(core_radius, 3)) * shell_material->getDensity(298.15);

		mass = core_mass + shell_mass;
	}

	__device__ void initialize(
		Species *core_material,
		Species *shell_material,
		Species *product_material,
		double core_radius,
		double overall_radius
	)
	{
		CoreShellParticle::core_material = core_material;
		CoreShellParticle::shell_material = shell_material;
		CoreShellParticle::product_material = product_material;

		setParticleMass(core_radius, overall_radius);
	}

	__host__ void printConfiguration(std::ostream &output_buffer)
	{
		double overall_radius;
		double core_radius;

		double mass;
		double core_mass;
		double shell_mass;

		cudaMemcpyFromSymbol(&overall_radius, CoreShellParticle::overall_radius, sizeof(double));
		cudaMemcpyFromSymbol(&core_radius, CoreShellParticle::core_radius, sizeof(double));

		cudaMemcpyFromSymbol(&mass, CoreShellParticle::mass, sizeof(double));
		cudaMemcpyFromSymbol(&core_mass, CoreShellParticle::core_mass, sizeof(double));
		cudaMemcpyFromSymbol(&shell_mass, CoreShellParticle::shell_mass, sizeof(double));
		
		output_buffer << "Particle Configuration\n\n";

		output_buffer << "Overall Radius :\t" << overall_radius << "\tm\n";
		output_buffer << "Core Radius :\t" << core_radius << "\tm\n\n";

		output_buffer << "Overall Mass :\t" << mass << "\tkg\n";
		output_buffer << "Core Mass :\t" << core_mass << "\tkg\n";
		output_buffer << "Shell Mass :\t" << shell_mass << "\tkg\n\n\n";
	}

	class Particle
	{

		protected:

			double _mass_fraction_core_material;
			double _mass_fraction_shell_material;
			double _mass_fraction_product_material;

		public:

			__device__ Particle()
			{
				_mass_fraction_core_material = core_mass / mass;
				_mass_fraction_shell_material = shell_mass / mass;

				_mass_fraction_product_material = 0.0;
			}

			__device__ __forceinline__ double getMassFractionsCoreMaterial() { return _mass_fraction_core_material; }
			__device__ __forceinline__ double getMassFractionsShellMaterial() { return _mass_fraction_shell_material; }
			__device__ __forceinline__ double getMassFractionsProductMaterial() { return _mass_fraction_product_material; }

			__device__ __forceinline__ double getDensity(double temperature)
			{
				double volume_fraction_core_material = getMassFractionsCoreMaterial() / core_material->getDensity(temperature);
				double volume_fraction_shell_material = getMassFractionsShellMaterial() / shell_material->getDensity(temperature);
				double volume_fraction_product_material = getMassFractionsProductMaterial() / product_material->getDensity(temperature);

				double sum = volume_fraction_core_material + volume_fraction_shell_material + volume_fraction_product_material;

				volume_fraction_core_material /= sum;
				volume_fraction_shell_material /= sum;
				volume_fraction_product_material /= sum;

				return 1.0 / (volume_fraction_core_material / core_material->getDensity(temperature) +
							volume_fraction_shell_material / shell_material->getDensity(temperature) +
							volume_fraction_product_material / product_material->getDensity(temperature));
			}

			__device__ __forceinline__ double getThermalConductivity(double temperature)
			{
				double volume_fraction_core_material = getMassFractionsCoreMaterial() / core_material->getDensity(temperature);
				double volume_fraction_shell_material = getMassFractionsShellMaterial() / shell_material->getDensity(temperature);
				double volume_fraction_product_material = getMassFractionsProductMaterial() / product_material->getDensity(temperature);

				double sum = volume_fraction_core_material + volume_fraction_shell_material + volume_fraction_product_material;

				volume_fraction_core_material /= sum;
				volume_fraction_shell_material /= sum;
				volume_fraction_product_material /= sum;

				return volume_fraction_core_material * core_material->getThermalConductivity(temperature) +
					volume_fraction_shell_material * shell_material->getThermalConductivity(temperature) +
					volume_fraction_product_material * product_material->getThermalConductivity(temperature);
			}

			__device__ __forceinline__ double getHeatCapacity(double temperature)
			{
				return 
				getMassFractionsCoreMaterial()    * core_material->getHeatCapacity(temperature) +
				getMassFractionsShellMaterial()   * shell_material->getHeatCapacity(temperature) +
				getMassFractionsProductMaterial() * product_material->getHeatCapacity(temperature); 
			}

			__device__ __forceinline__ double getEnthalpy(double temperature)
			{
				return 
				getMassFractionsCoreMaterial()    * core_material->getEnthalpy(temperature) +
				getMassFractionsShellMaterial()   * shell_material->getEnthalpy(temperature) +
				getMassFractionsProductMaterial() * product_material->getEnthalpy(temperature);
			}

			__device__ __forceinline__ void setMassFractions(
				double mass_fraction_core_material,
				double mass_fraction_shell_material,
				double mass_fraction_product_material
			) {
				_mass_fraction_core_material = mass_fraction_core_material;
				_mass_fraction_shell_material = mass_fraction_shell_material;
				_mass_fraction_product_material = mass_fraction_product_material;
			}

			__device__ __forceinline__ bool isReactionComplete(double tolerance = 0.001)
			{
				if (_mass_fraction_core_material < tolerance || _mass_fraction_shell_material < tolerance) return true;

				else return false;
			}
	};
}

#endif