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

	__device__ __forceinline__ double getRadialCoordinate(size_t index)
	{
		return CoreShellParticle::overall_radius * (double) index / ((double) n - 1.0);
	}

	__global__ void setGridSize(size_t num_points)
	{
		n = num_points;

		delta_r = CoreShellParticle::overall_radius / ((double) n - 1.0);

		radial_coordinate_sqr = new double[n];

		for (int i = 0; i < n; i++)
		{
			radial_coordinate_sqr[i] = pow(CoreShellParticle::overall_radius * (double) i / ((double) num_points - 1.0), 2);
		}
	}

	__global__  void deallocate()
	{
		delete [] radial_coordinate_sqr;
	}

	__global__ void setTimeStep(double time_step)
	{
		delta_t = time_step;
	}

	class Diffusion : public CoreShellParticle::Particle
	{
		private :

			double *_concentration_array_A;
			double *_concentration_array_B;

			LUSolver _solver_A;
			LUSolver _solver_B;

		public :

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

			// static void printConfiguration(std::ostream &output_stream);

			__device__ void initializeParticleFromDevice()
			{
				initialize();

				printf("Line 105\n");

				_concentration_array_A = new double[n];
				_concentration_array_B = new double[n];

				printf("Line 110\n");

				_solver_A.allocateMemoryFromDevice(n);
				_solver_B.allocateMemoryFromDevice(n);

				printf("Line 115\n");

				for (size_t i = 0; i < n; i++)
				{
					if (i < n)
					{
						if (radial_coordinate_sqr[i] < pow(CoreShellParticle::core_radius,2))
						{
							_concentration_array_A[i] = CoreShellParticle::core_material->getMolarDensity(298.15);
							_concentration_array_B[i] = 0;
						}

						else
						{
							_concentration_array_A[i] = 0;
							_concentration_array_B[i] = CoreShellParticle::shell_material->getMolarDensity(298.15);
						}
					}
				}

				printf("Line 135\n");
			}

			__device__ void deallocateMemoryFromDevice()
			{
				delete [] _concentration_array_A;
				delete [] _concentration_array_B;

				_solver_A.deallocateMemoryFromDevice();
				_solver_B.deallocateMemoryFromDevice();
			}

			// void setUpEquations(double diffusivity);
			// void setUpEquations(double diffusivity, Diffusion &diffusion_problem);

			// void solveEquations();

			// double getDiffusionMassA();

			// double getDiffusionMassB();

			// void copyFrom(Diffusion &diffusion_problem);
			// void copyTo(Diffusion &diffusion_problem);

			// void printConcentrationProfileA(std::ostream &output_stream, char delimiter = '\t', double curr_time = 0);
			// void printConcentrationProfileB(std::ostream &output_stream, char delimiter = '\t', double curr_time = 0);

			// void printGridPoints(std::ostream &output_stream, char delimiter = '\t');
	};

	__global__ void addShellMassA(CoreShellDIffusion::Diffusion *diffusion_problem, double *sum)
	{
		size_t i = blockDim.x + blockIdx.x + threadIdx.x;

		if (i < n-1) 
		{
			double s = diffusion_problem->getAtmConcA(i+1) * radial_coordinate_sqr[i+1];
			atomicAdd(sum, s);
		}
	}

	// __host__ double getMassA(CoreShellDIffusion::Diffusion *diffusion_problem)
	// {
	// 	double sum = 0.5 * diffusion_problem->getAtmConcA(n-1) * radial_coordinate_sqr[n-1];

	// 	double *dev_sum;
	// 	cudaMalloc(&dev_sum, sizeof(double));
	// 	cudaMemcpy(dev_sum, &sum, sizeof(double), cudaMemcpyHostToDevice);
		
	// 	addShellMassA<<<n+255/256, 256>>>(diffusion_problem, dev_sum);

	// 	cudaMemcpy(&sum, dev_sum, sizeof(double), cudaMemcpyDeviceToHost);
	// 	cudaFree(dev_sum);

	// 	return sum;
	// }
}

#endif