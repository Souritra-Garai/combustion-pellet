#ifndef __SIGMOID_FUNCTION__
#define __SIGMOID_FUNCTION__

#include <cmath>

#include "math/Data-Type.hpp"

inline real_t getSigmoid(real_t x, real_t origin = 0, real_t scale = 0)
{
	// return 1 / (1 + std::exp( - scale * (x - origin)));
	return 0.5 * (1. + std::tanh(0.5 * scale * (x - origin)));
}

inline real_t getSigmoidDerivative(real_t x, real_t origin = 0, real_t scale = 0)
{
	real_t sigma = getSigmoid(x, origin, scale);

	return scale * sigma * (1 - sigma);
}

#endif