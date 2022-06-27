#ifndef __VACCUM__
#define __VACCUM__

#include "thermo-physical-properties/Thermal_Conductivity.cuh"
#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Ideal_Gas.cuh"

__device__ IdealGas *vaccum;

__device__ void loadVaccum()
{
	Enthalpy enthalpy_v(
		1E-10,
		0,
		0,
		0,
		0,
		0
	);

	ThermalConductivity thermal_conductivity_v(1E-10, 0, 0);

	vaccum = new IdealGas(
		1E-10,
		enthalpy_v,
		thermal_conductivity_v
	);
}

__device__ void unloadVaccum()
{
	delete vaccum;
}

#endif