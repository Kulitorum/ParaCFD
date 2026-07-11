// sim_setup.h — build a runnable M2 open-channel simulation for the G1 viewer from a
// JSON config. Qt-free. Produces a ChannelFluidCore (inlet/Orlanski BCs + optional
// procedural cylinder obstacle) plus display metadata (domain extents, reference speed)
// the GUI needs for camera framing and colour-map ranges.
//
// It also exposes a SimRecipe (the immutable ingredients: grid, BCs, params, config
// obstacle, init state) + make_core(), so the GUI can REBUILD the core on the worker
// thread with a different obstacle mask (e.g. when the loaded STEP model is voxelized and
// injected as an obstacle) without re-parsing the config.
#pragma once

#include "core/fluid/channel_core.h"
#include "core/fluid/grid_metrics.h"
#include "core/geometry/model_placement.h"
#include "core/geometry/tri_mesh.h"

#include <memory>
#include <string>
#include <vector>

namespace windcfd::gui
{
	struct SimInfo
	{
		int nx = 0, ny = 0, nz = 0;
		double h = 0.05;               // voxel edge [m]
		double Lx = 0, Ly = 0, Lz = 0; // domain extents [m]
		double U = 1.0;                // reference inlet speed [m/s]
		double rho = 1.0;              // density [kg/m^3]
		double nu = 0.0;               // kinematic viscosity used [m^2/s]
		std::string name = "g1_viewer";
	};

	// Everything needed to (re)construct the core with a chosen obstacle mask. Copied into the
	// worker's rebuild factory so a re-inject is a pure worker-thread operation.
	struct SimRecipe
	{
		windcfd::core::MacGrid grid;               // GRADED: the DEVICE-view MacGrid (kernels); UNIFORM: null-metric grid.
		// Graded fine-core grid (null ⇒ uniform). GridMetrics is move-only, so it is held by shared_ptr:
		// the recipe is copied (into recipe_, the worker rebuild factory, scene_io) and every copy shares
		// this one instance, keeping the host+device metric arrays alive as long as any device-view MacGrid
		// (the core's g_, snapshots) references them. Host consumers use metrics->host_view().
		std::shared_ptr<windcfd::core::GridMetrics> metrics; // null ⇒ uniform grid
		windcfd::core::FineCoreSpec fine_core;               // the parsed fine-core spec (for Apply re-grid)
		bool graded = false;                                 // true ⇒ metrics is set (graded grid active)
		windcfd::core::ChannelBC bc;               // solid_mode is overridden per make_core() call
		windcfd::core::ChannelParams pr;
		std::vector<unsigned char> base_solid;   // config obstacle (procedural cylinder); all-zero if none
		int base_solid_mode = windcfd::core::SOLID_NOSLIP; // config's obstacle surface mode
		double init_u = 1.0;                     // init_uniform u0
		double init_v_blip = 0.0;                // init_uniform transverse blip (absolute m/s)
		SimInfo info;

		// Provenance so the GUI can rebuild this exact sim at a different grid (domain + h) by
		// re-running the SAME setup path with a GridOverride — see build_sim.
		std::string source_config;               // config JSON path this recipe was built from
	};

	// Optional domain+resolution override applied AFTER the config's own domain/h are read, so the
	// GUI's "Apply" can change ONLY the grid (nx,ny,nz from Lx/Ly/Lz and the voxel/cell size h) while
	// every other physics parameter (Re, obstacle, ν …) stays exactly as the config specifies. h is the
	// uniform grid spacing = fluid cell size = the voxel size the model is voxelized at (one and the same
	// in this simulator). The inlet current speed U is ALSO override-able (U > 0) so an Apply keeps the
	// "Input speed" the user has dialed in; the same knob is live via SimWorker::setInletSpeed (no reset).
	// All other physics is preserved.
	struct GridOverride
	{
		bool active = false;
		double Lx = 0.0, Ly = 0.0, Lz = 0.0; // domain extents [m]
		double h = 0.05;                     // voxel / cell size [m]
		double U = 0.0;                      // inlet (current) speed [m/s]; > 0 overrides the config U
	};

	// Derive the integer grid dimensions from a domain + voxel size, EXACTLY as the sim builder does
	// (floor with the same minimum clamps). Shared with the GUI so its live "→ nx×ny×nz" readout is
	// truthful (matches the grid Apply will actually build).
	void grid_dims_for(double Lx, double Ly, double Lz, double h, int& nx, int& ny, int& nz);

	// Construct a fresh, flow-initialised core from the recipe with `solid` as the obstacle mask
	// and `solid_mode` as its surface condition (SOLID_FREESLIP / SOLID_NOSLIP). Ready to step.
	std::unique_ptr<windcfd::core::ChannelFluidCore> make_core(const SimRecipe& r,
		const std::vector<unsigned char>& solid, int solid_mode);

	// Parse `config_path` (JSON), fill `recipe` (incl. recipe.info), and return the initial core
	// (config obstacle, config solid_mode). On any error, falls back to a small default open-
	// channel-with-cylinder setup and records the reason in `warn`. When `ov` is non-null and
	// active, its domain + voxel size replace the config's (everything else unchanged) so the GUI
	// can rebuild the same sim at a coarser/finer grid.
	std::unique_ptr<windcfd::core::ChannelFluidCore> build_sim(const std::string& config_path,
		SimRecipe& recipe, std::string& warn, const GridOverride* ov = nullptr);
}
