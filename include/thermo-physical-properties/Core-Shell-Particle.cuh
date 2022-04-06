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

	Species *core_material_host_ptr;
	Species *shell_material_host_ptr;
	Species *product_material_host_ptr;

	__global__ void setParticleMass(
		double core_radius,
		double overall_radius)
	{
		CoreShellParticle::core_radius = core_radius;
		CoreShellParticle::overall_radius = overall_radius;

		core_mass = (4.0 / 3.0) * M_PI * pow(core_radius, 3) * core_material->getDensity(298.15);
		shell_mass = (4.0 / 3.0) * M_PI * (pow(overall_radius, 3) - pow(core_radius, 3)) * shell_material->getDensity(298.15);

		mass = core_mass + shell_mass;
	}

	__host__ void initialize(
		Species *pointer_to_core_material,
		Species *pointer_to_shell_material,
		Species *pointer_to_product_material,
		double core_radius,
		double overall_radius
	)
	{
		cudaMalloc(&core_material_host_ptr, sizeof(Species));
		cudaMalloc(&shell_material_host_ptr, sizeof(Species));
		cudaMalloc(&product_material_host_ptr, sizeof(Species));

		cudaMemcpy(core_material_host_ptr, pointer_to_core_material, sizeof(Species), cudaMemcpyHostToDevice);
		cudaMemcpy(shell_material_host_ptr, pointer_to_shell_material, sizeof(Species), cudaMemcpyHostToDevice);
		cudaMemcpy(product_material_host_ptr, pointer_to_product_material, sizeof(Species), cudaMemcpyHostToDevice);

		cudaMemcpyToSymbol(core_material, &core_material_host_ptr, sizeof(Species *));
		cudaMemcpyToSymbol(shell_material, &shell_material_host_ptr, sizeof(Species *));
		cudaMemcpyToSymbol(product_material, &product_material_host_ptr, sizeof(Species *));

		setParticleMass<<<1,1>>>(core_radius, overall_radius);
	}

	__host__ void deallocate()
	{
		cudaFree(core_material_host_ptr);
		cudaFree(shell_material_host_ptr);
		cudaFree(product_material_host_ptr);
	}

	class Particle
	{

		protected:

			double _mass_fraction_core_material;
			double _mass_fraction_shell_material;
			double _mass_fraction_product_material;

		public:

			__device__ void initialize()
			{
				_mass_fraction_core_material = core_mass / mass;
				_mass_fraction_shell_material = shell_mass / mass;

				_mass_fraction_product_material = 0.0;
			}

			__device__ __host__ __forceinline__ double getMassFractionsCoreMaterial() { return _mass_fraction_core_material; }
			__device__ __host__ __forceinline__ double getMassFractionsShellMaterial() { return _mass_fraction_shell_material; }
			__device__ __host__ __forceinline__ double getMassFractionsProductMaterial() { return _mass_fraction_product_material; }

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

			__device__ __host__ __forceinline__ void setMassFractions(
				double mass_fraction_core_material,
				double mass_fraction_shell_material,
				double mass_fraction_product_material
			) {
				_mass_fraction_core_material = mass_fraction_core_material;
				_mass_fraction_shell_material = mass_fraction_shell_material;
				_mass_fraction_product_material = mass_fraction_product_material;
			}
	};
}

#endif