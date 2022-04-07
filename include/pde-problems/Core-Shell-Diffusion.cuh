#ifndef __CORE_SHELL_DIFFUSION__
#define __CORE_SHELL_DIFFUSION__

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* a, double b) { return b; }
#endif

#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "lusolver/LU_Solver.cuh"

namespace CoreShellDIffusion
{
	__device__ double delta_t;

	__device__ size_t n;
	__device__ double delta_r;
	__device__ double *radial_coordinate_sqr;

	__device__ double coefficient_2;

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
		delta_t = time_step;
		coefficient_2 = 1.0 / delta_t;
	}

	// static void printConfiguration(std::ostream &output_stream);

	class Diffusion : public CoreShellParticle::Particle
	{
		private :

			double *_concentration_array_A;
			double *_concentration_array_B;

			LUSolver _solver_A;
			LUSolver _solver_B;

			double coefficient_1;

		public :

			__device__ Diffusion() : CoreShellParticle::Particle(), _solver_A(n), _solver_B(n)
			{
				_concentration_array_A = new double[n];
				_concentration_array_B = new double[n];
			}

			__device__ ~Diffusion()
			{
				delete [] _concentration_array_B;
				delete [] _concentration_array_A;
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

			__device__ __forceinline__ void setCoefficient_1(double diffusivity)
			{
				coefficient_1 = - 0.5 * diffusivity / pow(delta_r, 2);
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
					double coefficient_3 = coefficient_1 * radial_coordinate_sqr[index+1] / radial_coordinate_sqr[index];
					double coefficient_4 = coefficient_2 - coefficient_1 - coefficient_3;
					double coefficient_5 = coefficient_2 + coefficient_1 + coefficient_3;

					_solver_A.setEquation(
						index,
						coefficient_1,
						coefficient_4,
						coefficient_3,
						- coefficient_1 * diffusion_problem->_concentration_array_A[index-1] +
						coefficient_5 * diffusion_problem->_concentration_array_A[index] +
						- coefficient_3 * diffusion_problem->_concentration_array_A[index+1]
					);

					_solver_B.setEquation(
						index,
						coefficient_1,
						coefficient_4,
						coefficient_3,
						- coefficient_1 * diffusion_problem->_concentration_array_B[index-1] +
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

			// void copyFrom(Diffusion &diffusion_problem);
			// void copyTo(Diffusion &diffusion_problem);

			// void printConcentrationProfileA(std::ostream &output_stream, char delimiter = '\t', double curr_time = 0);
			// void printConcentrationProfileB(std::ostream &output_stream, char delimiter = '\t', double curr_time = 0);

			// void printGridPoints(std::ostream &output_stream, char delimiter = '\t');
	};
}

#endif