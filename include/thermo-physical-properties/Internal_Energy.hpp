/**
 * @file Internal_Energy.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-09-15
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __INTERNAL_ENERGY__
#define __INTERNAL_ENERGY__

#include <math.h>

template <typename real_t>
class InternalEnergy
{
    private:
    
        real_t _A;
        real_t _B;
        real_t _C;
        real_t _D;
        real_t _E;
        real_t _F;

    public:
        
        InternalEnergy(
            real_t A = 0,
            real_t B = 0,
            real_t C = 0,
            real_t D = 0,
            real_t E = 0,
            real_t F = 0
        )
        {
            _A = A;
            _B = B;
            _C = C;
            _D = D;
            _E = E;
            _F = F;
        }

        real_t getInternalEnergy(real_t temperature)
        {
            temperature /= 1000;

            // A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H
            return
                _A * temperature +
                _B * pow(temperature, 2) / 2 +
                _C * pow(temperature, 3) / 3 +
                _D * pow(temperature, 4) / 4 -
                _E / temperature +
                _F
            ;
        }

        real_t getHeatCapacity(real_t temperature)
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