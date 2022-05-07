#ifndef __PELLET_FLAME_PROPAGATION__
#define __PELLET_FLAME_PROPAGATION__

#include <ostream>

#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.cuh"
#include "thermo-physical-properties/Packed-Pellet.cuh"
#include "pde-problems/Core-Shell-Diffusion.cuh"
#include "lu-solver/LU_Solver.cuh"

#define STEFAN_BOLTZMANN_CONSTANT 5.670374419E-8 // W / m2 - K4

namespace PelletFlamePropagation
{
	__device__ double kappa;
	__device__ double gamma;

	__device__ size_t n;
	__device__ double delta_x;

	__device__ double delta_T;

	__device__ double initial_ignition_length;
	__device__ double initial_ignition_temperature;

	__device__ void setImplicitness(
		double transient_term_implicitness,
		double diffusion_term_implicitness
	) {
		gamma = transient_term_implicitness;
		kappa = diffusion_term_implicitness;
	}

	__device__ void setNumGridPoints(size_t N)
	{
		n = N;
		delta_x = PackedPellet::length / ((double) N - 1.0);
	}

	__device__ void setInfinitesimalChangeInTemperature(double Delta_T)
	{
		delta_T = Delta_T;
	}

	__device__ void setInitialConditions(
		double ignition_temperature,
		double ignition_length
	) {
		initial_ignition_length = ignition_length;
		initial_ignition_temperature = ignition_temperature;
	}

	__host__ void printConfiguration(std::ostream &output_stream)
	{
		;
	}

	__host__ void printGridPoints(std::ostream &output_stream, char delimiter=',')
	{
		double length;
		cudaMemcpyFromSymbol(&length, PackedPellet::length, sizeof(double));

		size_t n;
		cudaMemcpyFromSymbol(&n, PelletFlamePropagation::n, sizeof(size_t));

		output_stream << NAN << delimiter;

		for (size_t i = 0; i < n-1; i++) output_stream << length * (double) i / ((double) n - 1.0) << delimiter;
		
		output_stream << length << "\n";
	}

	__host__ void printTemperatureArray(std::ostream &output_stream, double *temperature_array, double time, char delimiter='\t')
	{
		output_stream << time << delimiter;

		size_t n;
		cudaMemcpyFromSymbol(&n, PelletFlamePropagation::n, sizeof(size_t));
		
		for (size_t i = 0; i < n-1; i++) output_stream << temperature_array[i] << delimiter;

		output_stream << temperature_array[n-1] << '\n';
	}

	class FlamePropagation : public PackedPellet::Pellet
	{
		private:

			double *_temperature_array;
			
			LUSolver _solver;

			double *_thermal_conductivity_array;

			double *_prev_step_particle_enthalpy_array;

			CoreShellDIffusion::Diffusion **_diffusion_problems_array[3];
			// 0 - actual particles
			// 1 - particles for constant temperature evolution
			// 2 - particles for elevated temperature evolution

			ArrheniusDiffusivityModel *_diffusivity_model;

			__device__ __forceinline__ double getXCoordinate(size_t index)
			{
				return (double) index * PackedPellet::length / ((double) n - 1.0);
			}

			__device__ __forceinline__ bool inReactionZone(size_t index)
			{
				// printf("%d\t%b\t%b\n", index);
				return 
					_temperature_array[index] >= PackedPellet::ignition_temperature &&
					!_diffusion_problems_array[0][index]->isReactionComplete();
			}

			__device__ __forceinline__ double getParticleEnthalpyTemperatureDerivate(size_t index)
			{
				return
				(	_diffusion_problems_array[2][index]->getEnthalpy(_temperature_array[index]) -
					_diffusion_problems_array[1][index]->getEnthalpy(_temperature_array[index])
				) / delta_T;
			}

			__device__ __forceinline__ double getParticleEnthalpyTimeDerivative(size_t index)
			{
				return 
				(	_diffusion_problems_array[1][index]->getEnthalpy(_temperature_array[index]) -
					_diffusion_problems_array[0][index]->getEnthalpy(_temperature_array[index])
				) / CoreShellDIffusion::delta_t;
			}

			__device__ __forceinline__ double getParticleEnthalpyExplicitDerivative(size_t index)
			{
				return
				(	_diffusion_problems_array[0][index]->getEnthalpy(_temperature_array[index]) -
					_prev_step_particle_enthalpy_array[index]
				) / CoreShellDIffusion::delta_t;
			}
			
			__device__ __forceinline__ double2 getTransientTerm(size_t index)
			{
				double2 expr = {0, 0};

				if (index > 0 && index < n-1)
				{
					expr.x = _diffusion_problems_array[0][index]->getHeatCapacity(_temperature_array[index]) / CoreShellDIffusion::delta_t;

					expr.y = 0.0;

					if (inReactionZone(index))
					{
						expr.x += gamma * getParticleEnthalpyTemperatureDerivate(index) / CoreShellDIffusion::delta_t;

						expr.y = (1.0 - gamma) * getParticleEnthalpyExplicitDerivative(index) + gamma * getParticleEnthalpyTimeDerivative(index);
					}
				}

				return expr;
			}

			__device__ __forceinline__ double2 getHeatLossTerm(size_t index)
			{
				double2 expr;

				double h = (index == 0) || (index == n-1) ?
					PackedPellet::convective_heat_transfer_coefficient_flat_surface :
					PackedPellet::convective_heat_transfer_coefficient_curved_surface;

				expr.x = gamma * (h + 4.0 * PackedPellet::radiative_emissivity * STEFAN_BOLTZMANN_CONSTANT * pow(_temperature_array[index], 3));

				expr.y = 	h * (_temperature_array[index] - PackedPellet::ambient_temperature) + 
							PackedPellet::radiative_emissivity * STEFAN_BOLTZMANN_CONSTANT * (pow(_temperature_array[index], 4) - pow(PackedPellet::ambient_temperature, 4));
				
				return expr;
			}

			__device__ __forceinline__ void setUpBoundaryConditionX0()
			{
				double2 expr = getHeatLossTerm(0);

				double lambda_by_delta_x = 0.5 * (_thermal_conductivity_array[0] + _thermal_conductivity_array[1]) / delta_x;

				_solver.setEquationFirstRow(
					lambda_by_delta_x + expr.x,
					- lambda_by_delta_x,
					expr.x * _temperature_array[0] - expr.y
				);
			}

			__device__ __forceinline__ void setUpBoundaryConditionXN()
			{
				double2 expr = getHeatLossTerm(n-1);

				double lambda_by_delta_x = 0.5 * (_thermal_conductivity_array[n-2] + _thermal_conductivity_array[n-1]) / delta_x;

				_solver.setEquationLastRow(
					- lambda_by_delta_x,
					lambda_by_delta_x - expr.x,
					expr.y - expr.x * _temperature_array[n-1]
				);
			}

		public:

			__device__ __forceinline__ void collectEnthalpy(size_t index)
			{
				if (index > 0 && index < n-1)
				{
					_prev_step_particle_enthalpy_array[index] = _diffusion_problems_array[0][index]->getEnthalpy(_temperature_array[index]);
				}
			}

			__device__ __forceinline__ void setUpEquation(size_t index)
			{
				if (index == 0) setUpBoundaryConditionX0();
				
				else if (index < n-1)
				{
					double2 alpha = getTransientTerm(index);
					double2 beta = getHeatLossTerm(index);

					double lambda_by_delta_x_sqr_forward  = 0.5 * (_thermal_conductivity_array[index] + _thermal_conductivity_array[index+1]) / pow(delta_x, 2);
					double lambda_by_delta_x_sqr_backward = 0.5 * (_thermal_conductivity_array[index-1] + _thermal_conductivity_array[index]) / pow(delta_x, 2);

					_solver.setEquation(
						index,
						- kappa * lambda_by_delta_x_sqr_backward,
						_overall_particle_density * alpha.x - 4.0 * beta.x / PackedPellet::diameter + kappa * (lambda_by_delta_x_sqr_backward + lambda_by_delta_x_sqr_forward) +
						_degassing_fluid_volume_fractions * PackedPellet::degassing_fluid->getDensity(_temperature_array[index]) * PackedPellet::degassing_fluid->getHeatCapacity(_temperature_array[index]) / CoreShellDIffusion::delta_t,
						- kappa * lambda_by_delta_x_sqr_forward,
						- _overall_particle_density * (alpha.y - alpha.x * _temperature_array[index]) + 4.0 * (beta.y - beta.x * _temperature_array[index]) / PackedPellet::diameter +
						_degassing_fluid_volume_fractions * PackedPellet::degassing_fluid->getDensity(_temperature_array[index]) * PackedPellet::degassing_fluid->getHeatCapacity(_temperature_array[index]) * _temperature_array[index] / CoreShellDIffusion::delta_t +
						(1.0 - kappa) * (
							lambda_by_delta_x_sqr_forward  * (_temperature_array[index + 1] - _temperature_array[index]) -
							lambda_by_delta_x_sqr_backward * (_temperature_array[index] - _temperature_array[index - 1])
						)
					);
				}

				else if (index == n-1) setUpBoundaryConditionXN();
			}

			__device__ __forceinline__ void setParticleDiffusivity(size_t particle_index, size_t pellet_position_index)
			{
				if (particle_index < 3 && pellet_position_index > 0 && pellet_position_index < n-1)
				{
					if (inReactionZone(pellet_position_index))
					{
						_diffusion_problems_array[particle_index][pellet_position_index]->setDiffusivity(
							_diffusivity_model->getDiffusivity(
								_temperature_array[pellet_position_index] + (particle_index==2) * delta_T
							)
						);
					}
				}
			}

			__device__ __forceinline__ void setParticleEquation(size_t particle_index, size_t pellet_position_index, size_t eqn_index)
			{
				if (particle_index < 3 && pellet_position_index > 0 && pellet_position_index < n-1)
				{
					if (inReactionZone(pellet_position_index))
					{
						_diffusion_problems_array[particle_index][pellet_position_index]->setUpEquations(
							eqn_index,
							_diffusion_problems_array[0][pellet_position_index]
						);
					}
				}
			}

			__device__ __forceinline__ void solveParticleEquation(size_t particle_index, size_t pellet_position_index, size_t solver_index)
			{
				if (particle_index < 3 && pellet_position_index > 0 && pellet_position_index < n-1)
				{
					if (inReactionZone(pellet_position_index))
					{
						_diffusion_problems_array[particle_index][pellet_position_index]->solveEquations(solver_index);
					}
				}
			}

			__device__ __forceinline__ void updateParticleMassFraction(size_t particle_index, size_t pellet_position_index)
			{
				if (particle_index < 3 && pellet_position_index > 0 && pellet_position_index < n-1)
				{
					if (inReactionZone(pellet_position_index))
					{
						_diffusion_problems_array[particle_index][pellet_position_index]->updateMassFractions();
					}
				}
			}
			
			__device__ __forceinline__ void solveEquation()
			{
				_solver.getSolution(_temperature_array);
			}

			__device__ __forceinline__ void updateThermalConductivity(size_t index)
			{
				if (index == 0)
				{
					_thermal_conductivity_array[0] = getThermalConductivity(_diffusion_problems_array[0][1], _temperature_array[0]);
				}

				else if (index > 0 && index < n-1)
				{
					_thermal_conductivity_array[index] = getThermalConductivity(_diffusion_problems_array[0][index], _temperature_array[index]);
				}

				else if (index == n-1)
				{
					_thermal_conductivity_array[n-1] = getThermalConductivity(_diffusion_problems_array[0][n-2], _temperature_array[n-1]);
				}
			}

			__device__ FlamePropagation(
				double particle_volume_fractions
			) : PackedPellet::Pellet(particle_volume_fractions),
				_solver(n)
			{
				_thermal_conductivity_array = new double[n];
				_prev_step_particle_enthalpy_array = new double[n];

				_diffusion_problems_array[0] = new CoreShellDIffusion::Diffusion*[n];
				_diffusion_problems_array[1] = new CoreShellDIffusion::Diffusion*[n];
				_diffusion_problems_array[2] = new CoreShellDIffusion::Diffusion*[n];
			}

			__device__ ~FlamePropagation()
			{
				delete [] _diffusion_problems_array[2];
				delete [] _diffusion_problems_array[1];
				delete [] _diffusion_problems_array[0];

				delete [] _prev_step_particle_enthalpy_array;
				delete [] _thermal_conductivity_array;
			}

			__device__ void setArrayAddresses(
				double *temperature_array,
				double *concentration_array_A,
				double *concentration_array_B
			) {
				_temperature_array = temperature_array;

				for (size_t i = 1; i < n-1; i++)
				{
					_diffusion_problems_array[0][i]->setArrayAddresses(
						concentration_array_A + 3 * i * n,
						concentration_array_B + 3 * i * n
					);

					_diffusion_problems_array[1][i]->setArrayAddresses(
						concentration_array_A + (3 * i + 1) * n,
						concentration_array_B + (3 * i + 1) * n
					);

					_diffusion_problems_array[2][i]->setArrayAddresses(
						concentration_array_A + (3 * i + 2) * n,
						concentration_array_B + (3 * i + 2) * n
					);
				}
			}

			__device__ void setDiffusivityModel(ArrheniusDiffusivityModel *diffusivity_model)
			{
				_diffusivity_model = diffusivity_model;
			}
	
			__device__ void initializePellet(size_t pellet_position_index, size_t particle_position_index)
			{
				if (pellet_position_index == 0 && particle_position_index == 0) 
				
					_temperature_array[0] = initial_ignition_temperature;

				else if (pellet_position_index > 0 && pellet_position_index < n-1)
				{
					if (particle_position_index == 0) 

						_temperature_array[pellet_position_index] = 
							getXCoordinate(pellet_position_index) <= initial_ignition_length ?
							initial_ignition_temperature : PackedPellet::ambient_temperature;

					_diffusion_problems_array[0][pellet_position_index]->setInitialState(particle_position_index);
					_diffusion_problems_array[1][pellet_position_index]->setInitialState(particle_position_index);
					_diffusion_problems_array[2][pellet_position_index]->setInitialState(particle_position_index);
				}

				else if (pellet_position_index == n-1)

					_temperature_array[n-1] = PackedPellet::ambient_temperature;
			}

			__device__ __forceinline__ void allocateParticleMemory(size_t pellet_position_index)
			{
				if (pellet_position_index > 0 && pellet_position_index < n-1)
				{
					_diffusion_problems_array[0][pellet_position_index] = new CoreShellDIffusion::Diffusion();
					_diffusion_problems_array[1][pellet_position_index] = new CoreShellDIffusion::Diffusion();
					_diffusion_problems_array[2][pellet_position_index] = new CoreShellDIffusion::Diffusion();
				}
			}

			__device__ __forceinline__ void deallocateParticleMemory(size_t pellet_position_index)
			{
				if (pellet_position_index > 0 && pellet_position_index < n-1)
				{
					delete _diffusion_problems_array[0][pellet_position_index];
					delete _diffusion_problems_array[1][pellet_position_index];
					delete _diffusion_problems_array[2][pellet_position_index];
				}
			}

			__device__ __forceinline__ bool isCombustionComplete()
			{
				bool flag = true;
			
				for (size_t i = 1; i < n-1 && flag; i++)
				{
					flag = !inReactionZone(i);
				}
				
				return flag;
			}
	};
} // namespace PelletFlamePropagation

#endif