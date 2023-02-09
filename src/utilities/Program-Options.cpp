#include "utilities/Program-Options.hpp"
#include "utilities/Read-Data.hpp"
#include <iostream>

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

void setTemperatureOption(ez::ezOptionParser &opt)
{
	opt.add(
		"1500.0",
		0,
		1,
		0,
		"Set Temperature to ARG K.",
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

void setPhiOption(ez::ezOptionParser &opt)
{
	opt.add(
		"0.7",
		0,
		1,
		0,
		"Set Particle Volume Fraction of Pellet to ARG.",
		"-phi"
	);
}

double getPhiOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("-phi"))
	{
		opt.get("-phi")->getDouble(default_value);

		if (default_value <= 0 || default_value >= 1)
		{
			std::cerr << "Particle Volume Fraction should be in (0,1). Given " << default_value << "\n";
			std::abort();
		}
		
		std::cout << "Particle Volume Fraction of Pellet is set to " << default_value << "\n";
	}

	return default_value;
}

void setIgnitionTemperatureOption(ez::ezOptionParser &opt)
{
	opt.add(
		"1500.0",
		0,
		1,
		0,
		"Set Initial Pellet Ignition Temperature to ARG K.",
		"-ignT"
	);
}

double getIgnitionTemperatureOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("-ignT"))
	{
		opt.get("-ignT")->getDouble(default_value);

		double critical_temperature = readScalarData<double>("data/core-shell-particle/diffusivity-parameters", "critical-temperature.txt");

		if (default_value < critical_temperature)
		{
			std::cerr << "Ignition Temperature should be greater than " << critical_temperature << " K. Given " << default_value << " K\n";
			std::abort();
		}

		std::cout << "Initial Pellet Ignition Temperature is set to " << default_value << " K\n";
	}

	return default_value;
}

void setIgnitionLengthOption(ez::ezOptionParser &opt)
{
	opt.add(
		"0.1",
		0,
		1,
		0,
		"Set Initial Ignition Length Fraction of Pellet to ARG.",
		"-ignL"
	);
}

double getIgnitionLengthOption(ez::ezOptionParser &opt, double default_value)
{
	if (opt.isSet("-ignL"))
	{
		opt.get("-ignL")->getDouble(default_value);

		if (default_value <= 0 || default_value >= 1)
		{
			std::cerr << "Pellet Length Fraction should be in (0,1). Given " << default_value << "\n";
			std::abort();
		}
		
		std::cout << "Ignition Length Fraction of Pellet is set to " << default_value << "\n";
	}

	return default_value;
}