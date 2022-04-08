#ifndef __ARRHENIUS_DIFFUSIVITY_MODEL__
#define __ARRHENIUS_DIFFUSIVITY_MODEL__

class ArrheniusDiffusivityModel
{
    private:

        double _pre_exponential_factor;
        double _activation_energy;

    public:
       
        __host__ __device__ ArrheniusDiffusivityModel(
            double pre_exponential_factor,
            double activation_energy
        ) {
            setParameters(pre_exponential_factor, activation_energy);
        }

        __host__ __device__ ArrheniusDiffusivityModel() {
            setParameters(0, 0);
        }

        __host__ __device__ void setParameters(
            double pre_exponential_factor,
            double activation_energy
        ) {
            _pre_exponential_factor = pre_exponential_factor;
            _activation_energy = activation_energy;
        }
		
        __host__ __device__ double getDiffusivity(double temperature) {
            return _pre_exponential_factor * exp(- _activation_energy / (8.314 * temperature));
        }

        __host__ __device__ double getPreExponentialFactor() {
            return _pre_exponential_factor;
        }
		
        __host__ __device__ double getActivationEnergy() {
            return _activation_energy;
        }
};

#endif