// naca_step.h -- deterministic OpenCascade STEP generator for validation wings.
//
// The public interface is deliberately OpenCascade-free. Dimensions are SI metres;
// the implementation writes a millimetre-based STEP file so the normal ParaCFD STEP
// importer returns the requested dimensions after its established mm -> m conversion.
#pragma once

#include <string>

namespace paracfd::core
{
	struct NacaStepOptions
	{
		std::string four_digit_code = "2412";
		double chord = 1.0;
		double span = 2.0;
		int points_per_side = 129;
		std::string output_path;
	};

	// Generate a smooth, closed four-digit NACA section and extrude it along +Y.
	// The leading edge is at X=0, the trailing edge at X=chord, and the span is
	// centred on Y=0. The closed end faces make a clean finite wing benchmark.
	bool write_naca_4digit_step(const NacaStepOptions& options,
		std::string* error = nullptr);
}
