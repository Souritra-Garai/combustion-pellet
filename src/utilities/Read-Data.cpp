#include <utilities/Read-Data.hpp>

#include <iostream>
#include <string.h>
#include <fstream>
#include <dirent.h>

void openFile(std::ifstream &file, const char * directory_name, const char *file_name)
{
	char file_name_concat[_POSIX_PATH_MAX];
	strcpy(file_name_concat, directory_name);
	strcat(file_name_concat, "/");
	strcat(file_name_concat, file_name);

	file.open(file_name_concat, std::ios::in);

	if (file.fail())
	{
		std::cerr << "[ERROR] Could not open file - " << file_name_concat << std::endl;
		std::abort();
	}

	file.seekg(0, std::ios::beg);
}

template<typename real_t>
real_t readScalarData(const char *directory_name, const char *file_name)
{
	real_t scalar_val;

	std::ifstream scalar_file;
	openFile(scalar_file, directory_name, file_name);
	
	scalar_file >> scalar_val;
	scalar_file.close();
	
	return scalar_val;
}

template float readScalarData(const char *, const char *);
template size_t readScalarData(const char *, const char *);
template double readScalarData(const char *, const char *);
template long double readScalarData(const char *, const char *);

template<typename real_t>
Enthalpy<real_t> readEnthalpyData(const char *directory_name)
{
	std::ifstream enthalpy_file;
	openFile(enthalpy_file, directory_name, "enthalpy.txt");

	real_t A, B, C, D, E, F;

	enthalpy_file >> A >> B >> C >> D >> E >> F;

	enthalpy_file.close();

	return Enthalpy<real_t>(A, B, C, D, E, F);
}

template Enthalpy<float> readEnthalpyData(const char *);
template Enthalpy<double> readEnthalpyData(const char *);
template Enthalpy<long double> readEnthalpyData(const char *);

template<typename real_t>
ThermalConductivityQuadraticPolynomial<real_t> readThermalConductivityData(const char *directory_name)
{
	std::ifstream thermal_conductivity_file;
	openFile(thermal_conductivity_file, directory_name, "thermal-conductivity.txt");

	real_t a_0, a_1, a_2;

	thermal_conductivity_file >> a_0 >> a_1 >> a_2;

	thermal_conductivity_file.close();

	return ThermalConductivityQuadraticPolynomial<real_t>(a_0, a_1, a_2);
}

template ThermalConductivityQuadraticPolynomial<float> readThermalConductivityData(const char*);
template ThermalConductivityQuadraticPolynomial<double> readThermalConductivityData(const char*);
template ThermalConductivityQuadraticPolynomial<long double> readThermalConductivityData(const char*);

template<typename real_t>
IdealGas<real_t> readIdealGasData(const char *directory_name)
{
	real_t molar_mass = readScalarData<real_t>(directory_name, "molar-mass.txt");
	real_t gamma = readScalarData<real_t>(directory_name, "gamma.txt");

	Enthalpy<real_t> enthalpy = readEnthalpyData<real_t>(directory_name);
	ThermalConductivityQuadraticPolynomial<real_t> thermal_conductivity = readThermalConductivityData<real_t>(directory_name);
	
	return IdealGas<real_t>(
		molar_mass,
		gamma,
		enthalpy,
		thermal_conductivity
	);
}

template IdealGas<float> readIdealGasData(const char *);
template IdealGas<double> readIdealGasData(const char *);
template IdealGas<long double> readIdealGasData(const char *);

template<typename real_t>
Phase<real_t> readPhaseData(const char *directory_name)
{
	real_t density = readScalarData<real_t>(directory_name, "density.txt");
	
	std::ifstream temperature_bounds_file;
	openFile(temperature_bounds_file, directory_name, "temperature-bounds.txt");
	
	real_t temperature_lower_bound, temperature_upper_bound;
	temperature_bounds_file >> temperature_lower_bound >> temperature_upper_bound;
	temperature_bounds_file.close();

	Enthalpy<real_t> enthalpy = readEnthalpyData<real_t>(directory_name);
	ThermalConductivityQuadraticPolynomial<real_t> thermal_conductivity = readThermalConductivityData<real_t>(directory_name);
	
	return Phase<real_t>(
		density,
		enthalpy,
		thermal_conductivity,
		temperature_lower_bound,
		temperature_upper_bound
	);
}

template Phase<float> readPhaseData(const char *);
template Phase<double> readPhaseData(const char *);
template Phase<long double> readPhaseData(const char *);

template<typename real_t>
CondensedSpecies<real_t> readCondensedSpeciesData(const char *directory_name)
{
	char phase_directory_name[10][_POSIX_PATH_MAX];
	unsigned int n = 0;

	DIR* directory = opendir(directory_name);
	dirent* dirent_p;

	while ((dirent_p = readdir(directory)) != NULL && n < 10)
	{
		if (!strcmp(dirent_p->d_name, ".") || !strcmp(dirent_p->d_name, "..")) continue;
		if (dirent_p->d_type != DT_DIR) continue;

		strcpy(phase_directory_name[n], directory_name);
		strcat(phase_directory_name[n], "/");
		strcat(phase_directory_name[n], dirent_p->d_name);
		n++;
	}

	closedir(directory);

	Phase<real_t> *phases = (Phase<real_t>*) malloc(n * sizeof(Phase<real_t>));

	for (unsigned int i = 0; i < n; i++)
	{
		Phase<real_t> phase = readPhaseData<real_t>(phase_directory_name[i]);
		memcpy(phases + i, &phase, sizeof(Phase<real_t>));
	}

	real_t molar_mass = readScalarData<real_t>(directory_name, "molar-mass.txt");

	CondensedSpecies<real_t> condensed_species(n, phases, molar_mass);

	free(phases);

	return condensed_species;
}

template CondensedSpecies<float> readCondensedSpeciesData(const char *);
template CondensedSpecies<double> readCondensedSpeciesData(const char *);
template CondensedSpecies<long double> readCondensedSpeciesData(const char *);

template<typename real_t>
ArrheniusDiffusivityModel<real_t> readArrheniusDiffusivityModelParameters(const char *directory_name)
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

	return ArrheniusDiffusivityModel<real_t>(
		t_C,
		pre_exp_low,
		pre_exp_high,
		act_eng_low,
		act_eng_high
	);
}

template ArrheniusDiffusivityModel<float> readArrheniusDiffusivityModelParameters(const char *);
template ArrheniusDiffusivityModel<double> readArrheniusDiffusivityModelParameters(const char *);
template ArrheniusDiffusivityModel<long double> readArrheniusDiffusivityModelParameters(const char *);

