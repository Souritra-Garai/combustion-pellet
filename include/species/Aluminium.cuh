#ifndef __ALUMINIUM__
#define __ALUMINIUM__

#include "thermo-physical-properties/Thermal_Conductivity.cuh"
#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Phase.cuh"

#include "thermo-physical-properties/Species.cuh"

Species aluminium;

__host__ void loadAluminium(double sharpness_coefficient = 100.0)
{
	Enthalpy enthalpy_solid_Al, enthalpy_liquid_Al;
	ThermalConductivity thermal_conductivity_solid_Al, thermal_conductivity_liquid_Al;

	Phase phases_Al[2];

	enthalpy_solid_Al.assignCoefficients(
		28.08920,
		-5.414849,
		8.560423,
		3.427370,
		-0.277375,
		-9.147187
	);

	enthalpy_liquid_Al.assignCoefficients(
		31.75104,
		3.935826E-8,
		-1.786515E-8,
		2.694171E-9,
		5.480037E-9,
		-0.945684
	);

	thermal_conductivity_solid_Al.assignCoefficients(248.0, -0.067, 0.0);
	thermal_conductivity_liquid_Al.assignCoefficients(33.9, 7.89E-2, -2.099E-5);

	phases_Al[0].initialize(
		2700.0,
		enthalpy_solid_Al,
		thermal_conductivity_solid_Al,
		273.0,
		933.47,
		sharpness_coefficient
	);

	phases_Al[1].initialize(
		2375.0,
		enthalpy_liquid_Al,
		thermal_conductivity_liquid_Al,
		933.47,
		INFINITY,
		sharpness_coefficient
	);

	aluminium.initialize(2, phases_Al, 26.9815386E-3);
}

#endif