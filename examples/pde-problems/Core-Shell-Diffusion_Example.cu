#include <iostream>

#include "species/Aluminium.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "pde-problems/Core-Shell-Diffusion.cuh"

// #include "utilities/Keyboard_Interrupt.hpp"
// #include "utilities/File_Generator.hpp"

#define MAX_ITER 5000
#define Dt 0.0001

long double core_radius = 32.5E-6;
long double overall_radius = 39.5E-6;

// ArrheniusDiffusivityModel<long double> Alawieh_diffusivity(2.56E-6, 102.191E3);
// ArrheniusDiffusivityModel<long double> Du_diffusivity(9.54E-8, 26E3);

// void printState(size_t iteration_number, CoreShellDiffusion<long double> &particle);

__global__ void initialize(CoreShellDIffusion::Diffusion *diffusion_problem)
{
	diffusion_problem->initializeParticleFromDevice();
}

__global__ void deinitialize(CoreShellDIffusion::Diffusion *diffusion_problem)
{
	diffusion_problem->deallocateMemoryFromDevice();
}

int main(int argc, char const *argv[])
{
	loadAluminium();
	loadNickel();
	loadNickelAlumnide();

	CoreShellParticle::initialize(
		&aluminium, &nickel, &nickel_aluminide,
		core_radius, overall_radius
	);

	CoreShellDIffusion::setGridSize<<<1,1>>>(1001);
	CoreShellDIffusion::setTimeStep<<<1,1>>>(Dt);

	CoreShellDIffusion::Diffusion *diffusion_problem;
	cudaMalloc(&diffusion_problem, sizeof(CoreShellDIffusion::Diffusion));

	initialize<<<1,1>>>(diffusion_problem);
	deinitialize<<<1,1>>>(diffusion_problem);

	aluminium.deallocateMemory();
	nickel.deallocateMemory();
	nickel_aluminide.deallocateMemory();


	CoreShellDIffusion::deallocate<<<1,1>>>();
	CoreShellParticle::deallocate();

    // FileGenerator file_generator;

    // std::ofstream config_file = file_generator.getTXTFile("diffusion_config");
    // std::ofstream conc_A_file = file_generator.getCSVFile("concentration_A");
    // std::ofstream conc_B_file = file_generator.getCSVFile("concentration_B");

    // Ni_clad_Al_particle.printGridPoints(conc_A_file, ',');
    // Ni_clad_Al_particle.printGridPoints(conc_B_file, ',');

    // Ni_clad_Al_particle.printProperties(config_file);
    // config_file.close();

    // size_t __iter = 1;

    // long double temperature = 1400;
    // long double diffusivity = Alawieh_diffusivity.getDiffusivity(temperature);

    // setUpKeyboardInterrupt();
    
    // try
    // {
    //     while (!Ni_clad_Al_particle.isCombustionComplete() && __iter <= MAX_ITER)
    //     {
    //         Ni_clad_Al_particle.setUpEquations(diffusivity);
    //         Ni_clad_Al_particle.solveEquations();
            
    //         Ni_clad_Al_particle.printConcentrationProfileA(conc_A_file, ',', Dt * __iter);
    //         Ni_clad_Al_particle.printConcentrationProfileB(conc_B_file, ',', Dt * __iter);

    //         if (__iter % 100 == 0)
    //         {
    //             std::cout << "Iteration # " << __iter << std::endl;
    //         }

    //         __iter++;
    //     }
    // }

    // catch (InterruptException& e)
    // {
    //     std::cout << "\nCaught signal " << e.S << std::endl;

    //     conc_A_file.close();
    //     conc_B_file.close();

    //     return 1;
    // }

    // conc_A_file.close();
    // conc_B_file.close();
    
    return 0;
}

// void printState(size_t iteration_number, CoreShellDiffusion<long double> &particle)
// {
//     long double Y_Al = particle.getMassFractionsCoreMaterial();
//     long double Y_Ni = particle.getMassFractionsShellMaterial();
//     long double Y_NiAl = particle.getMassFractionsProductMaterial();

//     std::cout << "Iteration # " << iteration_number;

//     std::cout << "\tAl\t:\t" << Y_Al;
//     std::cout << "\tNi\t:\t" << Y_Ni;
//     std::cout << "\tNiAl\t:\t" << Y_NiAl;
//     std::cout << "\tSum\t:\t" << Y_Al + Y_Ni + Y_NiAl;

//     std::cout << std::endl;
// }