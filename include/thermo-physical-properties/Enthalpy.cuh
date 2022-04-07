#ifndef __ENTHALPY__
#define __ENTHALPY__

class Enthalpy
{
    private:
    
        double _A;
        double _B;
        double _C;
        double _D;
        double _E;
        double _F;

    public:
        
        __host__ __device__ Enthalpy(
            double A = 0.0,
            double B = 0.0,
            double C = 0.0,
            double D = 0.0,
            double E = 0.0,
            double F = 0.0
        ) {
            _A = A;
			_B = B;
			_C = C;
			_D = D;
			_E = E;
			_F = F;
        }

		__host__ __device__ __forceinline__ double getEnthalpy(double temperature)
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

		__host__ __device__ __forceinline__ double getHeatCapacity(double temperature)
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