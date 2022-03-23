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

// #include <math.h>

class Enthalpy
{
    private:
    
        double *_coefficients;

    public:
        
        __host__ Enthalpy(
            double A = 0,
            double B = 0,
            double C = 0,
            double D = 0,
            double E = 0,
            double F = 0
        )
        {
            cudaMalloc(&_coefficients, 6U * sizeof(double));

			double host_coefficients[] = {A, B, C, D, E, F};
			cudaMemcpy(_coefficients, host_coefficients, 6U * sizeof(double), cudaMemcpyHostToDevice);
        }

		__host__ ~Enthalpy()
		{
			cudaFree(_coefficients);
		}

		__host__ double* getCoefficients()
		{
			return _coefficients;
		}
};

__device__ double getStandardEnthalpy(double temperature, double coefficients[6])
{
	temperature /= 1000;

	// A*t + B*t2/2 + C*t3/3 + D*t4/4 − E/t + F − H
	return 1000.0 * (
		coefficients[0] * temperature +
		coefficients[1] * pow(temperature, 2) / 2 +
		coefficients[2] * pow(temperature, 3) / 3 +
		coefficients[3] * pow(temperature, 4) / 4 -
		coefficients[4] / temperature +
		coefficients[5]
	);
}

__device__ double getHeatCapacity(double temperature, double coefficients[6])
{
	temperature /= 1000;
	
	// A + B*t + C*t2 + D*t3 + E/t2
	return 
		coefficients[0] +
		coefficients[1] * temperature +
		coefficients[2] * pow(temperature, 2) +
		coefficients[3] * pow(temperature, 3) +
		coefficients[4] / pow(temperature, 2)
	;
}


#endif