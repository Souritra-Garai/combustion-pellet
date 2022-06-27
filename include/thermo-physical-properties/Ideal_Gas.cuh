#ifndef __IDEAL_GAS__
#define __IDEAL_GAS__

#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Thermal_Conductivity.cuh"

#define PRESSURE	101325.0

class IdealGas
{
	private :

		double _molar_mass;

		Enthalpy _enthalpy;
		ThermalConductivity _thermal_conductivity;

	public :

		__host__ __device__ IdealGas(
			double molar_mass,
			Enthalpy species_enthalpy,
			ThermalConductivity species_thermal_conductivity
		) :	_enthalpy(species_enthalpy),
			_thermal_conductivity(species_thermal_conductivity)
		{
			_molar_mass = molar_mass;
		}

		__host__ __device__ __forceinline__ double getMolarMass()
		{
			return _molar_mass;
		}

		__host__ __device__ __forceinline__ double getMolarDensity(double temperature)
        { 
            return PRESSURE / (8.314 * temperature);
        }

		__host__ __device__ __forceinline__ double getDensity(double temperature)
		{
			return getMolarMass() * getMolarDensity(temperature);
		}

		__host__ __device__ __forceinline__ double getEnthalpy(double temperature)
		{
			return _enthalpy.getEnthalpy(temperature) / getMolarMass();
		}

		__host__ __device__ __forceinline__ double getHeatCapacity(double temperature)
		{
			return _enthalpy.getHeatCapacity(temperature) / getMolarMass();
		}

		__host__ __device__ __forceinline__ double getThermalConductivity(double temperature)
		{
			return _thermal_conductivity.getThermalConductivity(temperature);
		}
};

#endif