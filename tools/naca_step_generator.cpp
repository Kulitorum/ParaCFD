#include "core/geometry/naca_step.h"

#include <cstdio>
#include <cstdlib>
#include <string>

using paracfd::core::NacaStepOptions;
using paracfd::core::NacaCoordinateOptions;
using paracfd::core::NacaGeometryConvention;
using paracfd::core::naca_geometry_convention_name;
using paracfd::core::write_naca_4digit_coordinates;
using paracfd::core::write_naca_4digit_step;

namespace
{
	void usage()
	{
		std::fprintf(stderr,
			"usage: naca_step_generator --naca DIGITS --output FILE [--chord METRES] "
			"[--span METRES] [--points-per-side N] [--xfoil-dat FILE] "
			"[--geometry closed|original|nasa-tmr-0012|nrel-4415]\n"
			"example: naca_step_generator --naca 2412 --chord 1 --span 2 "
			"--output Test-Data/NACA2412.step\n");
	}
}

int main(int argc, char** argv)
{
	NacaStepOptions options;
	std::string xfoil_dat_path;
	for (int index = 1; index < argc; ++index)
	{
		const std::string argument = argv[index];
		if (argument == "--naca" && index + 1 < argc) options.four_digit_code = argv[++index];
		else if (argument == "--output" && index + 1 < argc) options.output_path = argv[++index];
		else if (argument == "--chord" && index + 1 < argc) options.chord = std::atof(argv[++index]);
		else if (argument == "--span" && index + 1 < argc) options.span = std::atof(argv[++index]);
		else if (argument == "--points-per-side" && index + 1 < argc) options.points_per_side = std::atoi(argv[++index]);
		else if (argument == "--xfoil-dat" && index + 1 < argc) xfoil_dat_path = argv[++index];
		else if (argument == "--geometry" && index + 1 < argc)
		{
			const std::string value = argv[++index];
			if (value == "closed") options.geometry = NacaGeometryConvention::closed_sharp;
			else if (value == "original") options.geometry = NacaGeometryConvention::original_finite_trailing_edge;
			else if (value == "nasa-tmr-0012") options.geometry = NacaGeometryConvention::nasa_tmr_naca0012;
			else if (value == "nrel-4415") options.geometry = NacaGeometryConvention::nrel_measured_naca4415;
			else { usage(); return 2; }
		}
		else if (argument == "--help" || argument == "-h") { usage(); return 0; }
		else { usage(); return 2; }
	}
	if (options.output_path.empty()) { usage(); return 2; }
	std::string error;
	if (!write_naca_4digit_step(options, &error))
	{
		std::fprintf(stderr, "naca_step_generator: %s\n", error.c_str());
		return 1;
	}
	if (!xfoil_dat_path.empty())
	{
		NacaCoordinateOptions coordinates;
		coordinates.four_digit_code = options.four_digit_code;
		coordinates.geometry = options.geometry;
		coordinates.chord = options.chord;
		coordinates.points_per_side = options.points_per_side;
		coordinates.output_path = xfoil_dat_path;
		if (!write_naca_4digit_coordinates(coordinates, &error))
		{
			std::fprintf(stderr, "naca_step_generator: %s\n", error.c_str());
			return 1;
		}
	}
	std::printf("Wrote NACA %s (%s): chord %.9g m, span %.9g m, %d samples/side -> %s\n",
		options.four_digit_code.c_str(), naca_geometry_convention_name(options.geometry), options.chord, options.span,
		options.points_per_side, options.output_path.c_str());
	if (!xfoil_dat_path.empty())
		std::printf("Wrote matching section coordinates -> %s\n",
			xfoil_dat_path.c_str());
	return 0;
}
