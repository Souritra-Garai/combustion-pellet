/**
 * @file IdealGas.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __IDEAL_GAS__
#define __IDEAL_GAS__

// Required for printing to file / screen
#include <ostream>

#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Thermal_Conductivity.cuh"

class IdealGas
{
	private :

		double _molar_mass;

		Enthalpy _enthalpy;
		ThermalConductivity _thermal_conductivity;

	public :

		__host__ __device__ void initialize(
			double molar_mass,
			Enthalpy enthalpy,
			ThermalConductivity thermal_conductivity
		) {
			_molar_mass = molar_mass;
			
			_enthalpy = enthalpy;
			_thermal_conductivity = thermal_conductivity;
		}

		__device__ __forceinline__ double getMolarMass()
		{
			return _molar_mass;
		}

		__device__ __forceinline__ double getMolarDensity(double temperature)
        { 
            return 1.01325E5 / (8.314 * temperature);
        }

		__device__ __forceinline__ double getDensity(double temperature)
		{
			return getMolarMass() * getMolarDensity(temperature);
		}

		__device__ __forceinline__ double getEnthalpy(double temperature)
		{
			return _enthalpy.getEnthalpy(temperature) / getMolarMass();
		}

		__device__ __forceinline__ double getHeatCapacity(double temperature)
		{
			return _enthalpy.getHeatCapacity(temperature) / getMolarMass();
		}

		__device__ __forceinline__ double getThermalConductivity(double temperature)
		{
			return _thermal_conductivity.getThermalConductivity(temperature);
		}
};

#endif