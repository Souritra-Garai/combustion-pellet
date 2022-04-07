#ifndef __SPECIES__
#define __SPECIES__

#include <thermo-physical-properties/Phase.cuh>

class Species
{
	private :

		unsigned int _num_phases;

		Phase *_phases;

		double _molar_mass;

	public :

		__host__ __device__ Species(
			unsigned int number_of_phases,
			Phase *array_of_phases,
			double molar_mass
		) {
			_num_phases = number_of_phases;

			_phases = array_of_phases;

			_molar_mass = molar_mass;
		}

		__host__ __device__ __forceinline__ double getMolarMass()
		{
			return _molar_mass;
		}
		
		__host__ __device__ __forceinline__ double getDensity(double temperature)
		{
			double density = 0;

			#pragma unroll
			
				for (unsigned int i = 0; i < _num_phases; i++)
				{
					density += _phases[i].getDensity(temperature);
				}

			return density;
		}
		
		__host__ __device__ __forceinline__ double getMolarDensity(double temperature)
		{
			return getDensity(temperature) / getMolarMass();
		}

		__host__ __device__ __forceinline__ double getHeatCapacity(double temperature)
		{
			double heat_capacity = 0;

			#pragma unroll

				for (unsigned int i = 0; i < _num_phases; i++)
				{
					heat_capacity += _phases[i].getHeatCapacity(temperature);
				}

			return heat_capacity / getMolarMass();
		}

		__host__ __device__ __forceinline__ double getEnthalpy(double temperature)
		{
			double enthalpy = 0;

			#pragma unroll

				for (unsigned int i = 0; i < _num_phases; i++)
				{
					enthalpy += _phases[i].getEnthalpy(temperature);
				}

			return enthalpy / getMolarMass();
		}

		__host__ __device__ __forceinline__ double getThermalConductivity(double temperature)
		{ 
			double heat_conductivity = 0;

			#pragma unroll

				for (unsigned int i = 0; i < _num_phases; i++)
				{
					heat_conductivity += _phases[i].getThermalConductivity(temperature);
				}

			return heat_conductivity;
		}
};

#endif