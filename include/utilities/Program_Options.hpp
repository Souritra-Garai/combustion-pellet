#ifndef __PROGRAM_OPTIONS__
#define __PROGRAM_OPTIONS__

#include <iostream>

#include "utilities/ezOptionParser.hpp"
#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.hpp"

void setHelpOption(ez::ezOptionParser &opt)
{
	opt.add(
		"",
		0,
		0,
		0,
		"Display usage instructions.",
		"-h",
		"-help",
		"--help",
		"--usage"
	);
}

void displayHelpOption(ez::ezOptionParser &opt)
{
	if (opt.isSet("-h"))
	{
		std::string usage;
		opt.getUsage(usage);
		std::cout << usage;

		exit(0);
	}
}

void setSharpnessCoefficientOption(ez::ezOptionParser &opt)
{
	opt.add(
		"0.1",
		0,
		1,
		0,
		"Set sharpness coefficient to ARG.",
		"--sharp"
	);
}

double getSharpnessCoefficientOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("--sharp"))
	{
		opt.get("--sharp")->getDouble(default_value);
		std::cout << "Sharpness coefficient is set to " << default_value << "\n";
	}

	return default_value;
}

void setTemperatureUpperBoundOption(ez::ezOptionParser &opt)
{
	opt.add(
		"2500.0",
		0,
		1,
		0,
		"Set temperature range upper bound to ARG kelvin.",
		"--Tupper"
	);
}

double getTemperatureUpperBoundOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("--Tupper"))
	{
		opt.get("--Tupper")->getDouble(default_value);
		std::cout << "Temperature range upper bound is set to " << default_value << "\n";
	}

	return default_value;
}

void setTemperatureLowerBoundOption(ez::ezOptionParser &opt)
{
	opt.add(
		"273.0",
		0,
		1,
		0,
		"Set temperature range lower bound to ARG kelvin.",
		"--Tlower"
	);
}

double getTemperatureLowerBoundOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("--Tlower"))
	{
		opt.get("--Tlower")->getDouble(default_value);
		std::cout << "Temperature range lower bound is set to " << default_value << "\n";
	}

	return default_value;
}

void setTemperatureOption(ez::ezOptionParser &opt)
{
	opt.add(
		"1500.0",
		0,
		1,
		0,
		"Set temperature to ARG.",
		"--T"
	);
}

double getTemperatureOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("--T"))
	{
		opt.get("--T")->getDouble(default_value);
		std::cout << "Temperature is set to " << default_value << " K\n";
	}

	return default_value;
}

void setTemperatureStepOption(ez::ezOptionParser &opt)
{
	opt.add(
		".001",
		0,
		1,
		0,
		"Set temperature step to ARG kelvin.",
		"--DT"
	);
}

double getTemperatureStepOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("--DT"))
	{
		opt.get("--DT")->getDouble(default_value);
		std::cout << "Temperature step is set to " << default_value << " kelvin\n";
	}

	return default_value;
}

void setParticleGridSizeOption(ez::ezOptionParser &opt)
{
	opt.add(
		"1001",
		0,
		1,
		0,
		"Set number of grid points in particle to ARG.",
		"--N"
	);
}

int getParticleGridSizeOption(ez::ezOptionParser &opt, int default_value)
{
	if (opt.isSet("--N"))
	{
		opt.get("--N")->getInt(default_value);
		std::cout << "Particle grid size is set to " << default_value << ".\n";
	}

	return default_value;
}

void setTimeStepOption(ez::ezOptionParser &opt)
{
	opt.add(
		"0.000001",
		0,
		1,
		0,
		"Set the time step to ARG seconds.",
		"--Dt"
	);
}

double getTimeStepOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("--Dt"))
	{
		opt.get("--Dt")->getDouble(default_value);
		std::cout << "Time step is set to " << default_value << " seconds.\n";
	}

	return default_value;
}

void setDiffusivityModelOption(ez::ezOptionParser &opt)
{
	opt.add(
		"",
		0,
		0,
		0,
		"Set the diffusivity model to Du. Default is Alawieh.",
		"--Du"
	);
}

void getDiffusivityModelOption(ez::ezOptionParser &opt, ArrheniusDiffusivityModel<long double> &diffusivity_model)
{
	if (opt.isSet("--Du"))
	{
		diffusivity_model.setParameters(9.54E-8, 26E3);
		std::cout << "Diffusivity parameters are set to Du's model.\n";
	}
}

void setCoreRadiusOption(ez::ezOptionParser &opt)
{
	opt.add(
		"32.5E-6",
		0,
		1,
		0,
		"Set the core radius of Core-Shell Particle to ARG metres.",
		"--r_C"
	);
}

double getCoreRadiusOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("--r_C"))
	{
		opt.get("--r_C")->getDouble(default_value);
		std::cout << "Core radius of Core-Shell Particle is set to " << default_value << " metres.\n";
	}

	return default_value;
}

void setOverallRadiusOption(ez::ezOptionParser &opt)
{
	opt.add(
		"39.5E-6",
		0,
		1,
		0,
		"Set the overall radius of Core-Shell Particle to ARG metres.",
		"--r_P"
	);
}

double getOverallRadiusOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("--r_P"))
	{
		opt.get("--r_P")->getDouble(default_value);
		std::cout << "Overall radius of Core-Shell Particle is set to " << default_value << " metres.\n";
	}

	return default_value;
}

#endif