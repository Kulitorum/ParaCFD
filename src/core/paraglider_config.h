#pragma once

#include "core/geometry/model_placement.h"
#include "core/geometry/triangle_bvh.h"

#include <string>

namespace paracfd::core
{
	struct FreestreamConfig
	{
		double speed = 10.0;      // m/s, +X internally
		double rho = 1.225;       // kg/m^3
		double nu = 1.5e-5;      // m^2/s
	};

	struct DomainConfig
	{
		double upstream_margin = 3.0;   // m from placed bbox
		double downstream_margin = 8.0; // m
		double lateral_margin = 3.0;    // m per side
		double vertical_margin = 3.0;   // m per side; there is no ground plane
		bool half_wing_symmetry = false; // input mesh is +Y half; Y-min is its mirror plane
	};

	struct AmrConfig
	{
		double base_cell_size = 0.25; // m
		int max_levels = 3;           // ratio exactly 2 between adjacent levels
		int topology_refinement_levels = 0; // extra surface levels with a one-parent-brick EB halo
		int brick_size = 32;          // interior cells per axis
		int ghost_cells = 1;
		double wing_refinement_distance = 1.0;    // m
		double surface_refinement_distance = 0.35; // m
		double wake_length = 8.0;                 // m, +X from wing bbox
		double wake_radius = 2.0;                 // m around wing bbox y/z centre
		double min_volume_fraction = 0.25;         // conservative same-side merge threshold
		// Reporting threshold only: positive-area apertures below this h^2 fraction are
		// counted in diagnostics but remain in the conservative pressure/flux graph.
		double min_aperture_area_fraction = 1e-4;
	};

	struct ExternalSolverConfig
	{
		double cfl = 0.7;
		double smagorinsky_cs = 0.10;
		// The production FP32 cut-cell operator reaches a measured residual floor of
		// O(1e-4) on the full PlanB solid. Tighter requests stall rather than improve
		// the corrected flux field; use an FP64 validation build when tighter solves
		// are required.
		double projection_tolerance = 5e-4;
		int projection_max_iterations = 600;
	};

	struct AeroReferenceConfig
	{
		double area = 0.0;   // m^2; <=0 means coefficients deliberately unavailable
		double length = 0.0; // m; <=0 means moment coefficients unavailable
		Vec3d moment_origin{};
	};

	struct ParagliderConfig
	{
		// The sole geometry input: one closed aerodynamic OCCT solid. Construction
		// sheets, ribs, baffles, and open vent geometry are not accepted here.
		std::string step_path;
		double tessellation_deflection_mm = 2.0;
		ModelPlacement placement;
		FreestreamConfig freestream;
		DomainConfig domain;
		AmrConfig amr;
		ExternalSolverConfig solver;
		AeroReferenceConfig reference;
	};

	bool load_paraglider_config(const std::string& path, ParagliderConfig& out, std::string* error = nullptr);
	bool save_paraglider_config(const std::string& path, const ParagliderConfig& config, std::string* error = nullptr);
	Aabb3d automatic_flow_domain(const TriMesh& placed_mesh, const DomainConfig& margins);
}
