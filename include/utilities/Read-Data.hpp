#ifndef __READ_DATA__
#define __READ_DATA__

#include <thermo-physical-properties/Enthalpy.hpp>
#include <thermo-physical-properties/Thermal-Conductivity.hpp>

#include <thermo-physical-properties/Ideal-Gas.hpp>

#include <thermo-physical-properties/Phase.hpp>
#include <thermo-physical-properties/Condensed-Species.hpp>

#include <thermo-physical-properties/Arrhenius-Diffusivity-Model.hpp>

template<typename real_t>
real_t readScalarData(const char *directory_name, const char *var_name);

template<typename real_t>
Enthalpy<real_t> readEnthalpyData(const char *directory_name);

template<typename real_t>
ThermalConductivityQuadraticPolynomial<real_t> readThermalConductivityData(const char *directory_name);

template<typename real_t>
IdealGas<real_t> readIdealGasData(const char *directory_name);


template<typename real_t>
Phase<real_t> readPhaseData(const char *directory_name);

template<typename real_t>
CondensedSpecies<real_t> readCondensedSpeciesData(const char *directory_name);

template<typename real_t>
ArrheniusDiffusivityModel<real_t> readArrheniusDiffusivityModelParameters(const char *directory_name);

#endif