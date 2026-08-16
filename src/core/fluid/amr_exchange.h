#pragma once

#include "core/fluid/amr_fields.h"

#include <array>

namespace paracfd::core
{
	void restrict_fine_pressure_to_coarse(const AmrHierarchy& hierarchy, AmrHostFields& fields, int fine_level);
	void prolong_coarse_pressure_to_fine(const AmrHierarchy& hierarchy, AmrHostFields& fields, int fine_level);

	struct FluxMatchResult
	{
		double original_coarse_flux = 0.0;
		double fine_flux_sum = 0.0;
		double corrected_coarse_velocity = 0.0;
		double reflux_correction = 0.0;
	};

	// One coarse face is tiled by four fine faces. Areas include EB aperture fractions.
	FluxMatchResult match_two_to_one_flux(double coarse_velocity, double coarse_open_area,
		const std::array<double, 4>& fine_velocity, const std::array<double, 4>& fine_open_area);
}
