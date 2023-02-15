#ifndef __QUADRATIC_EXPRESSION__
#define __QUADRATIC_EXPRESSION__

#include "math/Data-Type.hpp"

class QuadraticExpression
{
    private :

        const real_t _a_0;
        const real_t _a_1;
        const real_t _a_2;

    public :

        QuadraticExpression(
            real_t a_0,
            real_t a_1,
            real_t a_2
        ) :	_a_0(a_0), _a_1(a_1), _a_2(a_2)
		{
            ;
        }

        inline real_t evaluateExpression(real_t x) const
        {
            return _a_0 + _a_1 * x + _a_2 * x * x;
        }
};

#endif