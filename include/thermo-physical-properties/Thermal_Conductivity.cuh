#ifndef __THERMAL_CONDUCTIVITY__
#define __THERMAL_CONDUCTIVITY__

class ThermalConductivity
{
    private :

        double _a_0;
        double _a_1;
        double _a_2;

    public :

        __host__ __device__ ThermalConductivity(
			double a_0 = 0.0,
			double a_1 = 0.0,
			double a_2 = 0.0
		) {
			_a_0 = a_0;
			_a_1 = a_1;
			_a_2 = a_2;
		}

        __device__ __host__ __forceinline__ double getThermalConductivity(double temperature)
        {
            return _a_0 + _a_1 * temperature + _a_2 * temperature * temperature;
        }
};

#endif