#ifndef __ARRHENIUS_DIFFUSIVITY_MODEL__
#define __ARRHENIUS_DIFFUSIVITY_MODEL__

#include <cmath>

#include <math/Data-Type.hpp>

class ArrheniusDiffusivityModel
{
	private:
		
		const real_t _critical_temperature;
		const real_t _pre_exponential_factor_low;
		const real_t _pre_exponential_factor_high;
		const real_t _activation_energy_low;
		const real_t _activation_energy_high;

	public:
		
		// Inputs 
		// temperature in K
		// pre-exponential factor in m^2/s
		// activation energy in J/mol.
		ArrheniusDiffusivityModel(
			real_t critical_temperature,
			real_t pre_exponential_factor_low,
			real_t pre_exponential_factor_high,
			real_t activation_energy_low,
			real_t activation_energy_high
		) : _critical_temperature(critical_temperature),
			_pre_exponential_factor_low(pre_exponential_factor_low),
			_pre_exponential_factor_high(pre_exponential_factor_high),
			_activation_energy_low(activation_energy_low / 8.314),
			_activation_energy_high(activation_energy_high / 8.314)
		{ ; }

		// Requires input T in K
		// Returns output in m^2/s
		inline real_t getDiffusivity(real_t temperature) const
		{
			if (temperature > _critical_temperature)

				return _pre_exponential_factor_high * std::exp(- _activation_energy_high / temperature);

			else

				return _pre_exponential_factor_low * std::exp(- _activation_energy_low / temperature);
		}
};

#endif