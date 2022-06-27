#ifndef __NITROGEN__
#define __NITROGEN__

#include "thermo-physical-properties/Thermal_Conductivity.cuh"
#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Ideal_Gas.cuh"

__device__ IdealGas *nitrogen;

__device__ void loadNitrogen()
{
	Enthalpy enthalpy_N(
		19.50583,
		19.88705,
		-8.598535,
		1.369784,
		0.527601,
		-4.935202
	);

	ThermalConductivity thermal_conductivity_N(8.38E-3, 6.37E-5, -6.85E-9);

	nitrogen = new IdealGas(
		28.0134E-3,
		enthalpy_N,
		thermal_conductivity_N
	);
}

__device__ void unloadNitrogen()
{
	delete nitrogen;
}

#endif