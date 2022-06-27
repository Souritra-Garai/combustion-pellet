#ifndef __HELIUM__
#define __HELIUM__

#include "thermo-physical-properties/Thermal_Conductivity.cuh"
#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Ideal_Gas.cuh"

__device__ IdealGas *helium;

__device__ void loadHelium()
{
	Enthalpy enthalpy_He(
		20.78603,
		4.850638E-10,
		-1.582916E-10,
		1.525102E-11,
		3.196347E-11,
		-6.197350
	);

	ThermalConductivity thermal_conductivity_He(0.0545, 3.58E-4, -5.08E-8);

	helium = new IdealGas(
		4.002602E-3,
		enthalpy_He,
		thermal_conductivity_He
	);
}

__device__ void unloadHelium()
{
	delete helium;
}

#endif