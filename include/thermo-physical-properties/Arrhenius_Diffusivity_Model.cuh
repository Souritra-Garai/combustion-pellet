#ifndef __ARRHENIUS_DIFFUSIVITY_MODEL__
#define __ARRHENIUS_DIFFUSIVITY_MODEL__

class ArrheniusDiffusivityModel
{
    private:

		double _pre_exponential_factor_low;
		double _activation_energy_low;

        double _pre_exponential_factor;
        double _activation_energy;

    public:
       
        __host__ __device__ ArrheniusDiffusivityModel(
            double pre_exponential_factor,
            double activation_energy
        ) {
            setParameters(pre_exponential_factor, activation_energy);
			setParametersLow(0, 0);
        }

        __host__ __device__ ArrheniusDiffusivityModel() {
            setParameters(0, 0);
			setParametersLow(0, 0);
        }

        __host__ __device__ void setParameters(
            double pre_exponential_factor,
            double activation_energy
        ) {
            _pre_exponential_factor = pre_exponential_factor;
            _activation_energy = activation_energy;
        }

		__host__ __device__ void setParametersLow(
			double pre_exponential_factor,
            double activation_energy
        ) {
            _pre_exponential_factor_low = pre_exponential_factor;
            _activation_energy_low = activation_energy;
        }
		
        inline __host__ __device__ double getDiffusivity(double temperature) {
			
			if (temperature >= 933.0)

            	return _pre_exponential_factor * exp(- _activation_energy / (8.314 * temperature));

			else

				return _pre_exponential_factor_low * exp(- _activation_energy_low / (8.314 * temperature));
        }
};

#endif