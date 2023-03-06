#ifndef __LINEAR_EXPRESSION__
#define __LINEAR_EXPRESSION__

#include "math/Data-Type.hpp"

class LinearExpression
{
	public:

		real_t a_0;
		real_t a_1;

		LinearExpression(real_t a0 = 0., real_t a1 = 0.) : a_0(a0), a_1(a1) {;}

		inline real_t evaluateExpression(real_t x = 0.) const
		{
			return a_0 + a_1 * x;
		}

		inline LinearExpression operator + () const
		{
			return LinearExpression(a_0, a_1);
		}

		inline LinearExpression operator + (LinearExpression z) const
		{
			return LinearExpression(a_0 + z.a_0, a_1 + z.a_1);
		}
		
		inline LinearExpression operator - () const
		{
			return LinearExpression(-a_0, -a_1);
		}

		inline LinearExpression operator - (LinearExpression z) const
		{
			return this->operator+(-z);
		}

		inline LinearExpression operator + (real_t scalar) const
		{
			return LinearExpression(a_0 + scalar, a_1 + scalar);
		}

		inline LinearExpression operator - (real_t scalar) const
		{
			return this->operator+(-scalar);
		}

		inline LinearExpression operator * (real_t scalar) const
		{
			return LinearExpression(a_0 * scalar, a_1 * scalar);
		}

		inline LinearExpression operator / (real_t scalar) const
		{
			return this->operator*(1. / scalar);
		}

		inline void operator += (LinearExpression z)
		{
			a_0 += z.a_0;
			a_1 += z.a_1;
		}

		inline void operator -= (LinearExpression z)
		{
			this->operator+=(-z);
		}

		inline void operator += (real_t scalar)
		{
			a_0 += scalar;
			a_1 += scalar;
		}

		inline void operator -= (real_t scalar)
		{
			this->operator+=(-scalar);
		}

		inline void operator *= (real_t scalar)
		{
			a_0 *= scalar;
			a_1 *= scalar;
		}

		inline void operator /= (real_t scalar)
		{
			this->operator*=(1. / scalar);
		}
};

#endif