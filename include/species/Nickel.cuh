#ifndef __NICKEL__
#define __NICKEL__

#include "thermo-physical-properties/Thermal_Conductivity.cuh"
#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Phase.cuh"

#include "thermo-physical-properties/Species.cuh"

Species nickel;

__host__ void loadNickel(double sharpness_coefficient = 100.0)
{
	Enthalpy enthalpy_solid_Ni_1, enthalpy_solid_Ni_2, enthalpy_solid_Ni_3, enthalpy_liquid_Ni;
	
	ThermalConductivity thermal_conductivity_solid_Ni_1, thermal_conductivity_solid_Ni_2, thermal_conductivity_liquid_Ni;

	Phase phases_Ni[4];

	enthalpy_solid_Ni_1.assignCoefficients(
		13.69160,
		82.46509,
		-174.9548,
		161.6011,
		-0.092417,
		-6.833644
	);

	enthalpy_solid_Ni_2.assignCoefficients(
		1248.045,
		-1257.510,
		0,
		0,
		-165.1266,
		-788.8263
	);

	enthalpy_solid_Ni_3.assignCoefficients(
		16.49839,
		18.74913,
		-6.639841,
		1.717278,
		1.872051,
		-0.467675
	);

	enthalpy_liquid_Ni.assignCoefficients(
		38.91103,
		0.0,
		0.0,
		0.0,
		0.0,
		-2.722630
	);

	thermal_conductivity_solid_Ni_1.assignCoefficients(107, -0.096, 3.2E-5);
	thermal_conductivity_solid_Ni_2.assignCoefficients(59.5, -7.67E-3, 1.7E-5);
	thermal_conductivity_liquid_Ni.assignCoefficients(17.95, 2.097E-2, 0.0);

	phases_Ni[0].initialize(
		8902.0,
		enthalpy_solid_Ni_1,
		thermal_conductivity_solid_Ni_1,
		273.0,
		600.0,
		sharpness_coefficient
	);

	phases_Ni[1].initialize(
		8902.0,
		enthalpy_solid_Ni_2,
		thermal_conductivity_solid_Ni_2,
		600.0,
		700.0,
		sharpness_coefficient
	);

	phases_Ni[2].initialize(
		8902.0,
		enthalpy_solid_Ni_3,
		thermal_conductivity_solid_Ni_2,
		700.0,
		1728.0,
		sharpness_coefficient
	);

	phases_Ni[3].initialize(
		7810.0,
		enthalpy_liquid_Ni,
		thermal_conductivity_liquid_Ni,
		1728.0,
		INFINITY,
		sharpness_coefficient
	);

	nickel.initialize(4, phases_Ni, 58.6934E-3);
}

#endif