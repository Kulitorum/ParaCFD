#pragma once

#include "core/geometry/embedded_boundary.h"
#include "core/paraglider_config.h"

#include <limits>
#include <vector>

namespace paracfd::core
{
	struct CompositeAmrPressureSystem;
	struct SmoothFabricWallPatchLoad;
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
		Vec3d viscous_force{};
		double represented_area = 0.0;
	};

	struct AerodynamicLoads
	{
		std::vector<TriangleAeroResult> triangles;
		Vec3d pressure_force{}; // pressure component, kept separate from wall shear
		Vec3d pressure_moment{};
		double pressure_drag = 0.0, pressure_side = 0.0, pressure_lift = 0.0;
		Vec3d viscous_force{}, viscous_moment{};
		double viscous_drag = 0.0, viscous_side = 0.0, viscous_lift = 0.0;
		Vec3d total_force{}, total_moment{};
		double drag = 0.0, side = 0.0, lift = 0.0;
		bool viscous_loads_valid = false;
		bool force_coefficients_valid = false;
		double cd_pressure = std::numeric_limits<double>::quiet_NaN();
		double cs_pressure = std::numeric_limits<double>::quiet_NaN();
		double cl_pressure = std::numeric_limits<double>::quiet_NaN();
		double cd_viscous = std::numeric_limits<double>::quiet_NaN();
		double cs_viscous = std::numeric_limits<double>::quiet_NaN();
		double cl_viscous = std::numeric_limits<double>::quiet_NaN();
		double cd = std::numeric_limits<double>::quiet_NaN();
		double cs = std::numeric_limits<double>::quiet_NaN();
		double cl = std::numeric_limits<double>::quiet_NaN();
	};

	AerodynamicLoads compute_pressure_loads(const EmbeddedBoundary& eb, const EbPressureState& pressure,
		std::size_t source_triangle_count, const FreestreamConfig& freestream,
		const AeroReferenceConfig& reference, double pressure_reference = 0.0);
	AerodynamicLoads compute_pressure_loads(const CompositeAmrPressureSystem& system,
		const std::vector<double>& pressure, std::size_t source_triangle_count,
		const FreestreamConfig& freestream, const AeroReferenceConfig& reference,
		double pressure_reference = 0.0);
	void accumulate_viscous_loads(AerodynamicLoads& loads,
		const std::vector<SmoothFabricWallPatchLoad>& patches,
		const FreestreamConfig& freestream, const AeroReferenceConfig& reference);
	// Convert integrated loads from the retained +Y half to their full-wing mirror pair.
	// Per-triangle fields remain the actually simulated half and are not duplicated.
	void reconstruct_y_symmetric_integrated_loads(AerodynamicLoads& loads,
		double symmetry_plane_y, double moment_origin_y);
}
