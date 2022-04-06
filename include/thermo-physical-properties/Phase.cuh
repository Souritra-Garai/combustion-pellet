#ifndef __PHASE__
#define __PHASE__

#include "thermo-physical-properties/Enthalpy.cuh"
#include "thermo-physical-properties/Thermal_Conductivity.cuh"

class Phase
{
	private:

		double _density;

		double _temperature_lower_bound;
		double _temperature_upper_bound;

		double _sharpness_coefficient;

		Enthalpy _enthalpy;

		ThermalConductivity _thermal_conductivity;

		__device__ __host__ __forceinline__ static double _getSigmoid(double x, double origin = 0, double scale = 0)
		{
			return 1 / (1 + exp( - scale * (x - origin)));
		}

		__device__ __host__ __forceinline__ static double _getSigmoidDerivative(double x, double origin = 0, double scale = 0)
		{
			double sigma = _getSigmoid(x, origin, scale);

			return scale * sigma * (1 - sigma);
		}

    public:

		__host__ void initialize(
			double density,
			Enthalpy enthalpy,
			ThermalConductivity thermal_conductivity,
			double temperature_lower_bound = 0,
			double temeprature_upper_bound = INFINITY,
			double sharpness_coefficient = 1
		) {
			_density = density;

			_enthalpy = enthalpy;
			_thermal_conductivity = thermal_conductivity;

			_temperature_lower_bound = temperature_lower_bound;
			_temperature_upper_bound = temeprature_upper_bound;

			_sharpness_coefficient = sharpness_coefficient;
		}

		__device__ __host__ __forceinline__ double getDensity(double temperature) 
		{
			return _density * (
				_getSigmoid(temperature, _temperature_lower_bound, _sharpness_coefficient) -
				_getSigmoid(temperature, _temperature_upper_bound, _sharpness_coefficient)
			);
		}

		__device__ __host__ __forceinline__ double getEnthalpy(double temperature)
		{
			return _enthalpy.getEnthalpy(temperature) * (
				_getSigmoid(temperature, _temperature_lower_bound, _sharpness_coefficient) -
				_getSigmoid(temperature, _temperature_upper_bound, _sharpness_coefficient)
			);
		}

		__device__ __host__ __forceinline__ double getHeatCapacity(double temperature)
		{
			return _enthalpy.getHeatCapacity(temperature) * (
				_getSigmoid(temperature, _temperature_lower_bound, _sharpness_coefficient) -
				_getSigmoid(temperature, _temperature_upper_bound, _sharpness_coefficient)
			) + _enthalpy.getEnthalpy(temperature) * (
				_getSigmoidDerivative(temperature, _temperature_lower_bound, _sharpness_coefficient) -
				_getSigmoidDerivative(temperature, _temperature_upper_bound, _sharpness_coefficient)
			);
		}

		__device__ __host__ __forceinline__ double getThermalConductivity(double temperature)
		{
			return _thermal_conductivity.getThermalConductivity(temperature) * (
				_getSigmoid(temperature, _temperature_lower_bound, _sharpness_coefficient) -
				_getSigmoid(temperature, _temperature_upper_bound, _sharpness_coefficient)
			);
		}
};

#endif