#ifndef __ALUMINIUM__
#define __ALUMINIUM__

#include "thermo-physical-properties/Thermal_Conductivity.cuh"
#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Phase.cuh"

#include "thermo-physical-properties/Species.cuh"

__device__ Species *aluminium;
__device__ Phase *phases_Al;

__device__ void loadAluminium(double sharpness_coefficient = 100.0)
{
	Enthalpy enthalpy_solid_Al(
		28.08920,
		-5.414849,
		8.560423,
		3.427370,
		-0.277375,
		-9.147187
	);

	Enthalpy enthalpy_liquid_Al(
		31.75104,
		3.935826E-8,
		-1.786515E-8,
		2.694171E-9,
		5.480037E-9,
		-0.945684
	);

	ThermalConductivity thermal_conductivity_solid_Al(248.0, -0.067, 0.0);
	ThermalConductivity thermal_conductivity_liquid_Al(33.9, 7.89E-2, -2.099E-5);

	Phase solid_Al(
		2700.0,
		enthalpy_solid_Al,
		thermal_conductivity_solid_Al,
		-INFINITY,
		933.47,
		sharpness_coefficient
	);

	Phase liquid_Al(
		2375.0,
		enthalpy_liquid_Al,
		thermal_conductivity_liquid_Al,
		933.47,
		INFINITY,
		sharpness_coefficient
	);

	phases_Al = new Phase[2]{solid_Al, liquid_Al};

	aluminium = new Species(2, phases_Al, 26.9815386E-3);
}

__device__ void unloadAluminium()
{
	delete aluminium;
	delete [] phases_Al;
}

#endif