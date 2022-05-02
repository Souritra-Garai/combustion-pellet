/**
 * @file Enthalpy.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief 
 * @version 0.1
 * @date 2021-09-15
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __ENTHALPY__
#define __ENTHALPY__

#include <math.h>

template <typename real_t>
class Enthalpy
{
    private:
    
        const real_t _A;
        const real_t _B;
        const real_t _C;
        const real_t _D;
        const real_t _E;
        const real_t _F;

    public:
        
        Enthalpy(
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
			_F(F)
        {
			;
        }

        inline real_t getStandardEnthalpy(real_t temperature)
        {
            temperature /= 1000;

            // A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H
            return 1000.0 * (
                _A * temperature +
                _B * pow(temperature, 2) / 2 +
                _C * pow(temperature, 3) / 3 +
                _D * pow(temperature, 4) / 4 -
                _E / temperature +
                _F
            );
        }

        inline real_t getHeatCapacity(real_t temperature)
        {
            temperature /= 1000;
            
            // A + B*t + C*t2 + D*t3 + E/t2
            return 
                _A +
                _B * temperature +
                _C * pow(temperature, 2) +
                _D * pow(temperature, 3) +
                _E / pow(temperature, 2)
            ;
        }
};

#endif