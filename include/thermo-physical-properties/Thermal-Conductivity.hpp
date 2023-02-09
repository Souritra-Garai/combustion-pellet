/**
 * @file Thermal-Conductivity.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief 
 * @version 2.1
 * @date 2021-09-17
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __THERMAL_CONDUCTIVITY__
#define __THERMAL_CONDUCTIVITY__

template <typename real_t>
class ThermalConductivityQuadraticPolynomial
{
    private :

        const real_t _a_0;
        const real_t _a_1;
        const real_t _a_2;

    public :

        ThermalConductivityQuadraticPolynomial(
            real_t a_0,
            real_t a_1,
            real_t a_2
        ) :	_a_0(a_0), _a_1(a_1), _a_2(a_2)
		{
            ;
        }

        inline real_t getThermalConductivity(real_t temperature)
        {
            return _a_0 + _a_1 * temperature + _a_2 * temperature * temperature;
        }
};

#endif