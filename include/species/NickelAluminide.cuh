#ifndef __NICKEL_ALUMINIDE__
#define __NICKEL_ALUMINIDE__

#include "thermo-physical-properties/Thermal_Conductivity.cuh"
#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Phase.cuh"

#include "thermo-physical-properties/Species.cuh"

__device__ Species *nickel_aluminide;
__device__ Phase *phases_NiAl;

__device__ void loadNickelAlumnide(double sharpness_coefficient = 100.0)
{
	Enthalpy enthalpy_solid_NiAl(
		41.86,
		13.81,
		0.0,
		0.0,
		0.0,
		-131.5
	);

	Enthalpy enthalpy_liquid_NiAl(
		71.16,
		0.0,
		0.0,
		0.0,
		0.0,
		-99.54
	);

	ThermalConductivity thermal_conductivity_solid_NiAl(72.8, 0.0156, -1.35E-5);
	ThermalConductivity thermal_conductivity_liquid_NiAl(23.9, 0.0, 0.0);

	Phase solid_NiAl(
		5650.0,
		enthalpy_solid_NiAl,
		thermal_conductivity_solid_NiAl,
		-INFINITY,
		1912.0,
		sharpness_coefficient
	);

	Phase liquid_NiAl(
		4798.3,
		enthalpy_liquid_NiAl,
		thermal_conductivity_liquid_NiAl,
		1912.0,
		INFINITY,
		sharpness_coefficient
	);

	phases_NiAl = new Phase[2]{solid_NiAl, liquid_NiAl};

	nickel_aluminide = new Species(2, phases_NiAl, 85.6748386E-3);
}

__device__ void unloadNickelAluminide()
{
	delete nickel_aluminide;
	delete [] phases_NiAl;
}

#endif