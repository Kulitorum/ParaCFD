#pragma once

#include "core/geometry/embedded_boundary.h"
#include "core/paraglider_config.h"

#include <limits>
#include <vector>

namespace paracfd::core
{
	struct CompositeAmrPressureSystem;
	struct EbPressureState
	{
		std::vector<double> regular;   // one value per Cartesian cell; ignored for split cells
		std::vector<double> irregular; // one value per explicit fluid fragment
		double pressure(FragmentRef ref) const;
	};

	struct TriangleAeroResult
	{
		double p_plus = std::numeric_limits<double>::quiet_NaN();
		double p_minus = std::numeric_limits<double>::quiet_NaN();
		double delta_p = std::numeric_limits<double>::quiet_NaN();
		double cp_plus = std::numeric_limits<double>::quiet_NaN();
		double cp_minus = std::numeric_limits<double>::quiet_NaN();
		double delta_cp = std::numeric_limits<double>::quiet_NaN();
		Vec3d pressure_force{};
		double represented_area = 0.0;
	};

	struct AerodynamicLoads
	{
		std::vector<TriangleAeroResult> triangles;
		Vec3d pressure_force{}; // pressure only; skin friction is not implemented
		Vec3d pressure_moment{};
		double pressure_drag = 0.0, pressure_side = 0.0, pressure_lift = 0.0;
		bool force_coefficients_valid = false;
		double cd_pressure = std::numeric_limits<double>::quiet_NaN();
		double cs_pressure = std::numeric_limits<double>::quiet_NaN();
		double cl_pressure = std::numeric_limits<double>::quiet_NaN();
	};

	AerodynamicLoads compute_pressure_loads(const EmbeddedBoundary& eb, const EbPressureState& pressure,
		std::size_t source_triangle_count, const FreestreamConfig& freestream,
		const AeroReferenceConfig& reference, double pressure_reference = 0.0);
	AerodynamicLoads compute_pressure_loads(const CompositeAmrPressureSystem& system,
		const std::vector<double>& pressure, std::size_t source_triangle_count,
		const FreestreamConfig& freestream, const AeroReferenceConfig& reference,
		double pressure_reference = 0.0);
}
