#ifndef __THERMAL_CONDUCTIVITY_PELLET__
#define __THERMAL_CONDUCTIVITY_PELLET__

#include "math/Data-Type.hpp"

real_t getThermalConductivityEMT(
	real_t porosity,
	real_t thermal_conductivity_fluid,
	real_t thermal_conductivity_particle
);

real_t getThermalConductivityCC(
	real_t porosity,
	real_t thermal_conductivity_fluid,
	real_t thermal_conductivity_particle
);

real_t getThermalConductivityME1(
	real_t porosity,
	real_t thermal_conductivity_fluid,
	real_t thermal_conductivity_particle
);

real_t getThermalConductivityME2(
	real_t porosity,
	real_t thermal_conductivity_fluid,
	real_t thermal_conductivity_particle
);

real_t getThermalConductivityMEB(
	real_t porosity,
	real_t thermal_conductivity_fluid,
	real_t thermal_conductivity_particle
);

real_t getThermalConductivityCCB(
	real_t porosity,
	real_t thermal_conductivity_fluid,
	real_t thermal_conductivity_particle
);

#endif