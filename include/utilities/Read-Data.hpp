#ifndef __READ_DATA__
#define __READ_DATA__

#include "math/Shomate-Expression.hpp"
#include "math/Quadratic-Expression.hpp"

#include "thermo-physical-properties/Ideal-Gas.hpp"

#include "thermo-physical-properties/Phase.hpp"
#include "thermo-physical-properties/Condensed-Species.hpp"

#include "thermo-physical-properties/Arrhenius-Diffusivity-Model.hpp"

template<typename data_t>
data_t readScalarData(const char *directory_name, const char *variable_file_name);

ShomateExpression readShomateExpressionCoefficients(const char *directory_name);

QuadraticExpression readQuadraticExpressionCoefficients(const char *directory_name);

IdealGas readIdealGasData(const char *directory_name);

Phase readPhaseData(const char *directory_name);

CondensedSpecies readCondensedSpeciesData(const char *directory_name);

ArrheniusDiffusivityModel readArrheniusDiffusivityModelParameters(const char *directory_name);

#endif