/**
 * @file Thermal_Conductivity.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
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

        real_t _a_0;
        real_t _a_1;
        real_t _a_2;

    public :

        ThermalConductivityQuadraticPolynomial(
            real_t a_0,
            real_t a_1,
            real_t a_2
        ) {
            _a_0 = a_0;
            _a_1 = a_1;
            _a_2 = a_2;
        }

        real_t getThermalConductivity(real_t temperature)
        {
            return _a_0 + _a_1 * temperature + _a_2 * temperature * temperature;
        }
};

template <typename real_t>
class ThermalConductivityLagrangeInterpolation
{
    private:

        unsigned int _num_data_points;

        real_t * _temperature_array;
        real_t * _thermal_conductivity_array;

    public:

        ThermalConductivityLagrangeInterpolation(
            unsigned int num_data_points,
            real_t temperature_array[],
            real_t thermal_conductivity_array[]
        ) {
            _num_data_points = num_data_points;

            _temperature_array          = new real_t[num_data_points];
            _thermal_conductivity_array = new real_t[num_data_points];
            
            for (unsigned int i = 0; i < _num_data_points; i++)
            {
                _temperature_array[i] = temperature_array[i];

                _thermal_conductivity_array[i] = thermal_conductivity_array[i];
            }
        }

        ~ThermalConductivityLagrangeInterpolation()
        {
            delete [] _temperature_array;
            delete [] _thermal_conductivity_array;
        }

        real_t getThermalConductivity(real_t temperature)
        {
            real_t thermal_conductivity = 0;

            for (unsigned int i = 0; i < _num_data_points; i++)
            {
                real_t term = _thermal_conductivity_array[i];

                for (unsigned int j = 0; j < _num_data_points; j++)
                {
                    if (i != j)
                    {
                        term *= (temperature - _temperature_array[j]) / (_temperature_array[i] - _temperature_array[j]);
                    }
                }

                thermal_conductivity += term;
            }

            return thermal_conductivity;
        }
};

#endif