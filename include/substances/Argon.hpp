#ifndef __ARGON__
#define __ARGON__

#include "thermo-physical-properties/Thermal_Conductivity.hpp"
#include "thermo-physical-properties/Enthalpy.hpp"
#include "thermo-physical-properties/IdealGas.hpp"

#define SHARPNESS_COEFFICIENT_AR 100.0

Enthalpy<long double> enthalpy_Ar(
    20.78600,
    2.825911E-7,
   -1.464191E-7,
    1.092131E-8,
   -3.661371E-8,
   -6.197350
);

ThermalConductivityQuadraticPolynomial<long double> thermal_conductivity_Ar(1.49E-3, 5.98E-5, -1.92E-8);

IdealGas<long double> Argon(
	39.948E-3,
	(long double) 5.0/3.0,
	enthalpy_Ar,
	thermal_conductivity_Ar
);

#endif