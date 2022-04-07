#ifndef __ARGON__
#define __ARGON__

#include "thermo-physical-properties/Thermal_Conductivity.cuh"
#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Ideal_Gas.cuh"

__device__ IdealGas *argon;

__device__ void loadArgon()
{
	Enthalpy enthalpy_Ar(
		20.78600,
		2.825911E-7,
		-1.464191E-7,
		1.092131E-8,
		-3.661371E-8,
		-6.197350
	);

	ThermalConductivity thermal_conductivity_Ar(1.49E-3, 5.98E-5, -1.92E-8);

	argon = new IdealGas(
		39.948E-3,
		enthalpy_Ar,
		thermal_conductivity_Ar
	);
}

__device__ void unloadArgon()
{
	delete argon;
}

#endif