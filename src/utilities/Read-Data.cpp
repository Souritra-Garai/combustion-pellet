#include <utilities/Read-Data.hpp>

#include <iostream>
#include <cstring>
#include <fstream>
#include <dirent.h>

void openFile(std::ifstream &file, const char * directory_name, const char *file_name)
{
	char file_name_concat[_POSIX_PATH_MAX];
	std::strcpy(file_name_concat, directory_name);
	std::strcat(file_name_concat, "/");
	std::strcat(file_name_concat, file_name);

	file.open(file_name_concat, std::ios::in);

	if (file.fail())
	{
		std::cerr << "[ERROR] Could not open file - " << file_name_concat << std::endl;
		std::abort();
	}

	file.seekg(0, std::ios::beg);
}

template<typename data_t>
data_t readScalarData(const char *directory_name, const char *file_name)
{
	data_t scalar_val;

	std::ifstream scalar_file;
	openFile(scalar_file, directory_name, file_name);
	
	scalar_file >> scalar_val;
	scalar_file.close();
	
	return scalar_val;
}

template float readScalarData<float>(const char *directory_name, const char *file_name);
template size_t readScalarData<size_t>(const char *directory_name, const char *file_name);
template double readScalarData<double>(const char *directory_name, const char *file_name);
template long double readScalarData<long double>(const char *directory_name, const char *file_name);

ShomateExpression readShomateExpressionCoefficients(const char *directory_name)
{
	std::ifstream enthalpy_file;
	openFile(enthalpy_file, directory_name, "enthalpy.txt");

	real_t A, B, C, D, E, F;

	enthalpy_file >> A >> B >> C >> D >> E >> F;

	enthalpy_file.close();

	return ShomateExpression(A, B, C, D, E, F);
}

QuadraticExpression readQuadraticExpressionCoefficients(const char *directory_name)
{
	std::ifstream thermal_conductivity_file;
	openFile(thermal_conductivity_file, directory_name, "thermal-conductivity.txt");

	real_t a_0, a_1, a_2;

	thermal_conductivity_file >> a_0 >> a_1 >> a_2;

	thermal_conductivity_file.close();

	return QuadraticExpression(a_0, a_1, a_2);
}

IdealGas readIdealGasData(const char *directory_name)
{
	real_t molar_mass = readScalarData<real_t>(directory_name, "molar-mass.txt");
	real_t gamma = readScalarData<real_t>(directory_name, "gamma.txt");

	ShomateExpression	enthalpy = readShomateExpressionCoefficients(directory_name);
	QuadraticExpression thermal_conductivity = readQuadraticExpressionCoefficients(directory_name);
	
	return IdealGas(
		molar_mass,
		gamma,
		enthalpy,
		thermal_conductivity
	);
}

Phase readPhaseData(const char *directory_name)
{
	real_t density = readScalarData<real_t>(directory_name, "density.txt");
	
	std::ifstream temperature_bounds_file;
	openFile(temperature_bounds_file, directory_name, "temperature-bounds.txt");
	
	real_t temperature_lower_bound, temperature_upper_bound;
	temperature_bounds_file >> temperature_lower_bound >> temperature_upper_bound;
	temperature_bounds_file.close();

	ShomateExpression	enthalpy = readShomateExpressionCoefficients(directory_name);
	QuadraticExpression	thermal_conductivity = readQuadraticExpressionCoefficients(directory_name);
	
	return Phase(
		density,
		enthalpy,
		thermal_conductivity,
		temperature_lower_bound,
		temperature_upper_bound
	);
}

CondensedSpecies readCondensedSpeciesData(const char *directory_name)
{
	char phase_directory_name[10][_POSIX_PATH_MAX];
	unsigned int n = 0;

	DIR* directory = opendir(directory_name);
	dirent* dirent_p;

	while ((dirent_p = readdir(directory)) != NULL && n < 10)
	{
		if (!strcmp(dirent_p->d_name, ".") || !strcmp(dirent_p->d_name, "..")) continue;
		if (dirent_p->d_type != DT_DIR) continue;

		std::strcpy(phase_directory_name[n], directory_name);
		std::strcat(phase_directory_name[n], "/");
		std::strcat(phase_directory_name[n], dirent_p->d_name);
		n++;
	}

	closedir(directory);

	Phase *phases = (Phase*) malloc(n * sizeof(Phase));

	for (unsigned int i = 0; i < n; i++)
	{
		Phase phase = readPhaseData(phase_directory_name[i]);
		memcpy(phases + i, &phase, sizeof(Phase));
	}

	real_t molar_mass = readScalarData<real_t>(directory_name, "molar-mass.txt");

	CondensedSpecies condensed_species(n, phases, molar_mass);

	free(phases);

	return condensed_species;
}

ArrheniusDiffusivityModel readArrheniusDiffusivityModelParameters(const char *directory_name)
{
	real_t t_C, pre_exp_low, pre_exp_high, act_eng_low, act_eng_high;
	
	t_C = readScalarData<real_t>(directory_name, "critical-temperature.txt");

	std::ifstream param_file;
	
	openFile(param_file, directory_name, "pre-exponential-factor.txt");
	param_file >> pre_exp_low >> pre_exp_high;
	param_file.close();

	openFile(param_file, directory_name, "activation-energy.txt");
	param_file >> act_eng_low >> act_eng_high;
	param_file.close();

	return ArrheniusDiffusivityModel(
		t_C,
		pre_exp_low,
		pre_exp_high,
		act_eng_low,
		act_eng_high
	);
}
