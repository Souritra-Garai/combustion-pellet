#ifndef __ARGON__
#define __ARGON__

#include "thermo-physical-properties/Thermal_Conductivity.cuh"
#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Ideal_Gas.cuh"

IdealGas argon;

__host__ void loadArgon()
{
	Enthalpy enthalpy_Ar;
	ThermalConductivity thermal_conductivity_Ar;

	enthalpy_Ar.assignCoefficients(
		20.78600,
		2.825911E-7,
		-1.464191E-7,
		1.092131E-8,
		-3.661371E-8,
		-6.197350
	);

	thermal_conductivity_Ar.assignCoefficients(1.49E-3, 5.98E-5, -1.92E-8);

	argon.initialize(
		39.948E-3,
		enthalpy_Ar,
		thermal_conductivity_Ar
	);
}

#endif