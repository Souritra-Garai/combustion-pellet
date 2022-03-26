#ifndef __NICKEL_ALUMINIDE__
#define __NICKEL_ALUMINIDE__

#include "thermo-physical-properties/Thermal_Conductivity.cuh"
#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Phase.cuh"

#include "thermo-physical-properties/Species.cuh"

Species nickel_aluminide;

__host__ void loadNickelAlumnide(double sharpness_coefficient = 100.0)
{
	Enthalpy enthalpy_solid_NiAl, enthalpy_liquid_NiAl;

	ThermalConductivity thermal_conductivity_solid_NiAl, thermal_conductivity_liquid_NiAl;

	Phase phases_NiAl[2];

	enthalpy_solid_NiAl.assignCoefficients(
		41.86,
		13.81,
		0.0,
		0.0,
		0.0,
		-131.5
	);

	enthalpy_liquid_NiAl.assignCoefficients(
		71.16,
		0.0,
		0.0,
		0.0,
		0.0,
		-99.54
	);

	thermal_conductivity_solid_NiAl.assignCoefficients(72.8, 0.0156, -1.35E-5);
	thermal_conductivity_liquid_NiAl.assignCoefficients(23.9, 0.0, 0.0);

	phases_NiAl[0].initialize(
		5650.0,
		enthalpy_solid_NiAl,
		thermal_conductivity_solid_NiAl,
		273.0,
		1912.0,
		sharpness_coefficient
	);

	phases_NiAl[1].initialize(
		4798.3,
		enthalpy_liquid_NiAl,
		thermal_conductivity_liquid_NiAl,
		1912.0,
		INFINITY,
		sharpness_coefficient
	);

	nickel_aluminide.initialize(2, phases_NiAl, 85.6748386E-3);
}

#endif