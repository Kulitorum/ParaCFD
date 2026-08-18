// naca_step.h -- deterministic OpenCascade STEP generator for validation wings.
//
// The public interface is deliberately OpenCascade-free. Dimensions are SI metres;
// the implementation writes a millimetre-based STEP file so the normal ParaCFD STEP
// importer returns the requested dimensions after its established mm -> m conversion.
#pragma once

#include <string>

namespace paracfd::core
{
	enum class NacaGeometryConvention
	{
		// ParaCFD's original analytical section: the fourth thickness
		// coefficient is changed to -0.1036 so the trailing edge closes.
		closed_sharp,
		// The original analytical four-digit thickness distribution from
		// NACA Report 460. The -0.1015 coefficient leaves a finite gap, which
		// this STEP generator closes with a straight trailing-edge face. The
		// Report 460 tunnel models instead faired that finite ordinate to a
		// rounded zero-ordinate edge, so this is not an exact test-model replica.
		original_finite_trailing_edge,
		// Exact normalized geometry prescribed by NASA's Turbulence Modeling
		// Resource NACA 0012 validation case (not the generic NACA formula).
		nasa_tmr_naca0012,
		// Measured ordinates of the 18-inch Ohio State/NREL NACA 4415 test
		// article (Table A1 of the 7x10-ft wind-tunnel report).
		nrel_measured_naca4415
	};

	struct NacaStepOptions
	{
		std::string four_digit_code = "2412";
		NacaGeometryConvention geometry = NacaGeometryConvention::closed_sharp;
		double chord = 1.0;
		double span = 2.0;
		int points_per_side = 129;
		std::string output_path;
	};

	struct NacaCoordinateOptions
	{
		std::string four_digit_code = "2412";
		NacaGeometryConvention geometry = NacaGeometryConvention::closed_sharp;
		double chord = 1.0;
		int points_per_side = 129;
		std::string output_path;
	};

	// Analytic closed primitive used by the imported-CAD/embedded-boundary
	// regression suite.  Keeping the public type OCCT-free lets tests exercise the
	// exact STEP importer and CAD-edge provenance path used by production geometry.
	struct ClosedCylinderStepOptions
	{
		double centre_x = 0.0;
		double centre_y = 0.0;
		double base_z = 0.0;
		double radius = 0.75;
		double height = 5.6;
		std::string output_path;
	};

	// Generate a smooth, closed four-digit NACA section and extrude it along +Y.
	// The leading edge is at X=0, the trailing edge at X=chord, and the span is
	// centred on Y=0. The closed end faces make a clean finite wing benchmark.
	bool write_naca_4digit_step(const NacaStepOptions& options,
		std::string* error = nullptr);

	// Write the exact analytical section used by the STEP generator in the
	// conventional XFOIL/Selig order (trailing edge -> upper -> nose -> lower
	// -> trailing edge). This lets the validation oracle analyse precisely the
	// geometry sent through OpenCascade rather than XFOIL's finite-gap built-in.
	bool write_naca_4digit_coordinates(const NacaCoordinateOptions& options,
		std::string* error = nullptr);

	bool write_closed_cylinder_step(const ClosedCylinderStepOptions& options,
		std::string* error = nullptr);

	const char* naca_geometry_convention_name(NacaGeometryConvention convention);
}
