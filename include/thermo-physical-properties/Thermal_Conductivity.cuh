#ifndef __THERMAL_CONDUCTIVITY__
#define __THERMAL_CONDUCTIVITY__

class ThermalConductivity
{
    private :

        double _a_0;
        double _a_1;
        double _a_2;

    public :

        __host__ __device__ void assignCoefficients(
			double a_0,
			double a_1,
			double a_2
		) {
			_a_0 = a_0;
			_a_1 = a_1;
			_a_2 = a_2;
		}

        __device__ double getThermalConductivity(double temperature)
        {
            return _a_0 + _a_1 * temperature + _a_2 * temperature * temperature;
        }
};

#endif