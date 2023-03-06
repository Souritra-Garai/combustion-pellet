#ifndef __SHOMATE_EXPRESSION__
#define __SHOMATE_EXPRESSION__

#include <cmath>

#include "math/Data-Type.hpp"

class ShomateExpression
{
	private:

		const real_t _A;
		const real_t _B;
		const real_t _C;
		const real_t _D;
		const real_t _E;
		const real_t _F;

		const real_t _integral_B;
		const real_t _integral_C;
		const real_t _integral_D;

	public:
		
		ShomateExpression(
			real_t A = 0,
			real_t B = 0,
			real_t C = 0,
			real_t D = 0,
			real_t E = 0,
			real_t F = 0
		) :	_A(A),
			_B(B),
			_C(C),
			_D(D),
			_E(E),
			_F(F),
			_integral_B(B / 2.),
			_integral_C(C / 3.),
			_integral_D(D / 4.)
		{
			;
		}

		// Normalizes T(K) for input
		// T (in K) -> T / 1000
		inline static real_t normalizeInput(real_t x)
		{
			return x * 1E-3;
		}

		// Expects input as T(K)/1000
		// Returns output in J/mol.-K
		inline real_t evaluateExpression(real_t x) const
		{
			return 
				_A +
				_B * x +
				_C * std::pow(x, 2) +
				_D * std::pow(x, 3) +
				_E * std::pow(x,-2)
			;
		}

		// Expects input as T(K)/1000
		// Returns output in kJ/mol.
		inline real_t evaluateExpressionIntegral(real_t x) const
		{
			return
				_A * x +
				_integral_B * std::pow(x, 2) +
				_integral_C * std::pow(x, 3) +
				_integral_D * std::pow(x, 4) -
				_E / x +
				_F
			;
		}

		// Converts to SI units
		// kJ/mol. -> J/mol.
		inline static real_t normalizeIntegralOutput(real_t h)
		{
			return h * 1E3;
		}
};

#endif