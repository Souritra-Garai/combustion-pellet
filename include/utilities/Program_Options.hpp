#ifndef __PROGRAM_OPTIONS__
#define __PROGRAM_OPTIONS__

#include <iostream>

#include "utilities/ezOptionParser.hpp"

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
		"Set temperature range upper bound to ARG.",
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
		"Set temperature range lower bound to ARG.",
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

#endif