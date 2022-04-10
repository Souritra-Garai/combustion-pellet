#ifndef __CORE_SHELL_DIFFUSION__
#define __CORE_SHELL_DIFFUSION__

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* a, double b) { return b; }
#endif

#include <ostream>

#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "lusolver/LU_Solver.cuh"

namespace CoreShellDIffusion
{
	__device__ double time;
	__device__ double delta_t;
	__device__ double _coefficient_2;

	__device__ size_t n;
	__device__ double delta_r;
	__device__ double *radial_coordinate_sqr;

	__device__ void setGridSize(size_t num_points)
	{
		n = num_points;

		delta_r = CoreShellParticle::overall_radius / ((double) n - 1.0);

		radial_coordinate_sqr = new double[n];

		for (int i = 0; i < n; i++)
		{
			radial_coordinate_sqr[i] = pow(CoreShellParticle::overall_radius * (double) i / ((double) num_points - 1.0), 2);
		}
	}

	__device__  void deallocate()
	{
		delete [] radial_coordinate_sqr;
	}

	__device__ void setTimeStep(double time_step)
	{
		time = 0.0;
		delta_t = time_step;
		_coefficient_2 = 1.0 / delta_t;
	}

	__host__ void printConfiguration(std::ostream &output_stream)
	{
		double delta_t, delta_r;
		size_t n;

		cudaMemcpyFromSymbol(&delta_t, CoreShellDIffusion::delta_t, sizeof(double));
		cudaMemcpyFromSymbol(&delta_r, CoreShellDIffusion::delta_r, sizeof(double));
		cudaMemcpyFromSymbol(&n, CoreShellDIffusion::n, sizeof(size_t));

		output_stream << "Time Step\t:\t" << delta_t << "\ts" << std::endl;

		output_stream << "Number of Grid Points\t:\t" << n << std::endl;

		output_stream << "Grid Size\t:\t" << delta_r << "\tm" << std::endl;
	}

	__host__ void printGridPoints(std::ostream &output_stream, char delimiter=',')
	{
		double radius;
		cudaMemcpyFromSymbol(&radius, CoreShellParticle::overall_radius, sizeof(double));

		size_t n;
		cudaMemcpyFromSymbol(&n, CoreShellDIffusion::n, sizeof(size_t));

		output_stream << NAN << delimiter;

		for (size_t i = 0; i < n-1; i++) output_stream << radius * (double) i / ((double) n - 1.0) << delimiter;
		output_stream << "\n";
	}

	__host__ void printConcentrationArray(std::ostream &output_stream, double *concentration_array, double time = 0.0, char delimiter = '\t')
	{
		output_stream << time << delimiter;

		size_t n;
		cudaMemcpyFromSymbol(&n, CoreShellDIffusion::n, sizeof(size_t));
		
		for (size_t i = 0; i < n-1; i++) output_stream << concentration_array[i] << delimiter;

		output_stream << concentration_array[n-1] << '\n';
	}

	class Diffusion : public CoreShellParticle::Particle
	{
		private :

			double *_concentration_array_A;
			double *_concentration_array_B;

			LUSolver _solver_A;
			LUSolver _solver_B;

			double _coefficient_1;

		public :

			__device__ Diffusion() : CoreShellParticle::Particle(), _solver_A(n), _solver_B(n)
			{
				;
			}

			__device__ void setArrayAddresses(
				double *concentration_array_A,
				double *concentration_array_B
			) {
				_concentration_array_A = concentration_array_A;
				_concentration_array_B = concentration_array_B;
			}

			__device__ __forceinline__ double getRxnConcA(size_t index)
			{
				return max(_concentration_array_A[index] - _concentration_array_B[index], 0.0);
			}

			__device__ __forceinline__ double getRxnConcB(size_t index)
			{
				return max(_concentration_array_B[index] - _concentration_array_A[index], 0.0);
			}

			__device__ __forceinline__ double getRxnConcAB(size_t index)
			{
				return min(_concentration_array_A[index], _concentration_array_B[index]);
			}

			__device__ __forceinline__ double getAtmConcA(size_t index)
			{
				return _concentration_array_A[index];
			}

			__device__ __forceinline__ double getAtmConcB(size_t index)
			{
				return _concentration_array_B[index];
			}

			__device__ __forceinline__ void addAtmMolesA(size_t index, double *sum_A)
			{
				double mol_A = getAtmConcA(index) * radial_coordinate_sqr[index];
				
				atomicAdd(sum_A, mol_A);
			}

			__device__ __forceinline__ void setInitialState(size_t index)
			{
				if (radial_coordinate_sqr[index] <= pow(CoreShellParticle::core_radius, 2))
				{
					_concentration_array_A[index] = CoreShellParticle::core_material->getMolarDensity(298.15);
					_concentration_array_B[index] = 0.0;
				}
				else
				{
					_concentration_array_A[index] = 0.0;
					_concentration_array_B[index] = CoreShellParticle::shell_material->getMolarDensity(298.15);
				}
			}

			__device__ double getAtmMassA()
			{
				double moles = 0.5 * getAtmConcA(n-1) * radial_coordinate_sqr[n-1];

				for (size_t i = 1; i < n-1; i++)	moles += getAtmConcA(i) * radial_coordinate_sqr[i];

				return 4.0 * M_PI * delta_r * CoreShellParticle::core_material->getMolarMass() * moles;
			}

			__device__ double getAtmMassB()
			{
				double moles = 0.5 * getAtmConcA(n-1) * radial_coordinate_sqr[n-1];

				for (size_t i = 1; i < n-1; i++)	moles += getAtmConcB(i) * radial_coordinate_sqr[i];

				return 4.0 * M_PI * delta_r * CoreShellParticle::shell_material->getMolarMass() * moles;
			}

			__device__ __forceinline__ void setDiffusivity(double diffusivity)
			{
				_coefficient_1 = - 0.5 * diffusivity / pow(delta_r, 2);
			}

			__device__ __forceinline__ void setUpEquations(size_t index, Diffusion *diffusion_problem)
			{
				if (index == 0)
				{
					_solver_A.setEquationFirstRow(1, -1, 0);
					_solver_B.setEquationFirstRow(1, -1, 0);
				}
				
				else if (index < n-1)
				{
					double coefficient_3 = _coefficient_1 * radial_coordinate_sqr[index+1] / radial_coordinate_sqr[index];
					double coefficient_4 = _coefficient_2 - _coefficient_1 - coefficient_3;
					double coefficient_5 = _coefficient_2 + _coefficient_1 + coefficient_3;

					_solver_A.setEquation(
						index,
						_coefficient_1,
						coefficient_4,
						coefficient_3,
						- _coefficient_1 * diffusion_problem->_concentration_array_A[index-1] +
						coefficient_5 * diffusion_problem->_concentration_array_A[index] +
						- coefficient_3 * diffusion_problem->_concentration_array_A[index+1]
					);

					_solver_B.setEquation(
						index,
						_coefficient_1,
						coefficient_4,
						coefficient_3,
						- _coefficient_1 * diffusion_problem->_concentration_array_B[index-1] +
						coefficient_5 * diffusion_problem->_concentration_array_B[index] +
						- coefficient_3 * diffusion_problem->_concentration_array_B[index+1]
					);
				}

				else if (index == n-1)
				{
					_solver_A.setEquationLastRow(-1, 1, 0);
					_solver_B.setEquationLastRow(-1, 1, 0);
				}
			}

			__device__ __forceinline__ void setUpEquations(size_t index)
			{
				setUpEquations(index, this);
			}
			
			__device__ __forceinline__ void solveEquations(size_t index)
			{
				if (index == 0)			_solver_A.getSolution(_concentration_array_A);
				else if (index == 1)	_solver_B.getSolution(_concentration_array_B);
			}

			__device__ __forceinline__ void updateMassFractions()
			{
				double Y_A  = 0.5 * getRxnConcA(n-1)  * radial_coordinate_sqr[n-1];
				double Y_B  = 0.5 * getRxnConcA(n-1)  * radial_coordinate_sqr[n-1];
				double Y_AB = 0.5 * getRxnConcAB(n-1) * radial_coordinate_sqr[n-1];

				for (size_t i = 1; i < n-1; i++)
				{
					Y_A  += getRxnConcA(i)  * radial_coordinate_sqr[i];
					Y_B  += getRxnConcB(i)  * radial_coordinate_sqr[i];
					Y_AB += getRxnConcAB(i) * radial_coordinate_sqr[i];
				}

				Y_A  *= CoreShellParticle::core_material->getMolarMass();
				Y_B  *= CoreShellParticle::shell_material->getMolarMass();
				Y_AB *= CoreShellParticle::product_material->getMolarMass();

				double sum = Y_A + Y_B + Y_AB;
				setMassFractions(Y_A / sum, Y_B / sum, Y_AB / sum);
			}
	};
}

#endif