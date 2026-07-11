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
#include "core/geometry/model_placement.h"
#include "core/geometry/tri_mesh.h"
#include "core/sediment/seabed_engine.h"

#include <memory>
#include <string>
#include <vector>

namespace scour::gui
{
	struct SimInfo
	{
		int nx = 0, ny = 0, nz = 0;
		double h = 0.05;               // voxel edge [m]
		double Lx = 0, Ly = 0, Lz = 0; // domain extents [m]
		double U = 1.0;                // reference inlet speed [m/s]
		double rho = 1.0;              // density [kg/m^3]
		double nu = 0.0;               // kinematic viscosity used [m^2/s]
		bool cylinder = false;
		std::string name = "g1_viewer";
	};

	// Everything needed to (re)construct the core with a chosen obstacle mask. Copied into the
	// worker's rebuild factory so a re-inject is a pure worker-thread operation.
	struct SimRecipe
	{
		scour::core::MacGrid grid;
		scour::core::ChannelBC bc;               // solid_mode is overridden per make_core() call
		scour::core::ChannelParams pr;
		std::vector<unsigned char> base_solid;   // config obstacle (procedural cylinder); all-zero if none
		int base_solid_mode = scour::core::SOLID_NOSLIP; // config's obstacle surface mode
		double init_u = 1.0;                     // init_uniform u0
		double init_v_blip = 0.0;                // init_uniform transverse blip (absolute m/s)
		SimInfo info;

		// Provenance so the GUI can rebuild this exact sim at a different grid (domain + h) by
		// re-running the SAME setup path with a GridOverride — see build_sim/build_seabed_sim.
		std::string source_config;               // config/scenario JSON path this recipe was built from
		bool is_scenario = false;                // true ⇒ built by build_seabed_sim (erodible seabed)
	};

	// Optional domain+resolution override applied AFTER the config's own domain/h are read, so the
	// GUI's "Apply" can change ONLY the grid (nx,ny,nz from Lx/Ly/Lz and the voxel/cell size h) while
	// every other physics/scenario parameter (Re, obstacle, d50, MORFAC, α, Cs, ν, structure …) stays
	// exactly as the config specifies. h is the uniform grid spacing = fluid cell size = the voxel size
	// the model + bed are voxelized at (one and the same in this simulator). The inlet current speed U
	// is ALSO override-able (U > 0) so an Apply keeps the "Input speed" the user has dialed in; the
	// same knob is live via SimWorker::setInletSpeed (no reset). All other physics is preserved.
	struct GridOverride
	{
		bool active = false;
		double Lx = 0.0, Ly = 0.0, Lz = 0.0; // domain extents [m]
		double h = 0.05;                     // voxel / cell size [m]
		double U = 0.0;                      // inlet (current) speed [m/s]; > 0 overrides the config U
	};

	// Derive the integer grid dimensions from a domain + voxel size, EXACTLY as the sim builders do
	// (floor with the same minimum clamps). Shared with the GUI so its live "→ nx×ny×nz" readout is
	// truthful (matches the grid Apply will actually build).
	void grid_dims_for(double Lx, double Ly, double Lz, double h, int& nx, int& ny, int& nz);

	// Construct a fresh, flow-initialised core from the recipe with `solid` as the obstacle mask
	// and `solid_mode` as its surface condition (SOLID_FREESLIP / SOLID_NOSLIP). Ready to step.
	std::unique_ptr<scour::core::ChannelFluidCore> make_core(const SimRecipe& r,
		const std::vector<unsigned char>& solid, int solid_mode);

	// Parse `config_path` (JSON), fill `recipe` (incl. recipe.info), and return the initial core
	// (config obstacle, config solid_mode). On any error, falls back to a small default open-
	// channel-with-cylinder setup and records the reason in `warn`. When `ov` is non-null and
	// active, its domain + voxel size replace the config's (everything else unchanged) so the GUI
	// can rebuild the same sim at a coarser/finer grid.
	std::unique_ptr<scour::core::ChannelFluidCore> build_sim(const std::string& config_path,
		SimRecipe& recipe, std::string& warn, const GridOverride* ov = nullptr);

	// --- Erodible-seabed scenario (the live morphodynamic demo) ------------------------------------
	// Everything the GUI needs to run + visualise a seabed scenario: the morphodynamic engine (moved
	// to the worker), the structure mesh + its on-the-bed placement (for the display), and the initial
	// bed elevation for the height-field colour map.
	struct SeabedScenario
	{
		bool active = false;
		bool has_structure = false;
		int spinup = 200;                                      // flow-only steps before morphology
		std::unique_ptr<scour::core::SeabedMorpho> engine;     // moved into the SimWorker
		scour::core::TriMesh structure_mesh;                   // for display (empty ⇒ bare seabed)
		scour::core::ModelPlacement structure_place;           // structure seated on the bed
		double bed_z0 = 1.0;                                   // initial bed elevation [m]
		scour::core::SeabedParams params;                      // the sediment params used (for a later "Add sand")
	};

	// Build a seabed scenario from `config_path` (domain, h, sand_depth, U, d50, MORFAC, α, Cs, the
	// STEP structure path). Voxelizes the structure onto the bed, constructs the morphodynamic engine
	// and the flow core (bed + structure solid; bed-inlet mask on), fills `recipe.info` for the viewer
	// and `scen` for the worker/viz. Returns the ready-to-step core, or (on error) falls back to
	// build_sim and leaves scen.active=false with the reason in `warn`. When `ov` is non-null and
	// active, its domain + voxel size replace the config's (the structure is re-voxelized and the
	// sand reservoir re-initialised at the new h; every other scenario parameter is unchanged).
	// When `structure_place` is non-null it OVERRIDES the default centre-on-bed placement of the
	// structure (it already includes the sand-depth seating), so a GUI Apply can re-voxelize the
	// structure exactly where the user positioned it with the gizmo. Null ⇒ the config's centre-on-bed
	// placement (an initial scenario load + the physics gates are byte-identical).
	std::unique_ptr<scour::core::ChannelFluidCore> build_seabed_sim(const std::string& config_path,
		SimRecipe& recipe, SeabedScenario& scen, std::string& warn, const GridOverride* ov = nullptr,
		const scour::core::ModelPlacement* structure_place = nullptr);

	// Build an erodible-seabed scenario on an EXISTING grid from a given RIGID structure mask (the GUI's
	// "Add sand"): fills the lower `sand_depth` m of every non-rigid column with sand (rigid/structure
	// cells stay solid), constructs the morphodynamic engine + flow core, and fills `recipe`/`scen` just
	// like build_seabed_sim — but with no config/STEP (the structure comes from the live sim's solids, so
	// a fluid viewer with a voxelized obstacle can be turned into a live sediment run). `base` carries the
	// current sediment params (d50/α/…); null ⇒ nominal defaults. Structure display mesh (if any) is the
	// caller's to attach. PRE-CALIBRATION / qualitative (nominal constants), like build_seabed_sim.
	// When `build_core` is false the flow core is NOT constructed (returns nullptr): the recipe (bc/pr/
	// base_solid/info) and scen (engine, params, bed_z0) are still fully filled, so the GUI can convert a
	// running fluid viewer to a seabed run IN PLACE — reusing the existing, already-developed core via
	// ChannelFluidCore::update_solid — without resetting the flow. `recipe.base_solid` is the initial
	// sand+structure flow-solid mask, and `recipe.bc` carries the log-law inlet z0/bed_datum to apply live.
	std::unique_ptr<scour::core::ChannelFluidCore> build_seabed_from_structure(
		const scour::core::MacGrid& g, double U, double sand_depth,
		const std::vector<unsigned char>& structure, const scour::core::SeabedParams* base,
		const std::string& source_config, int spinup, double nu_fluid,
		SimRecipe& recipe, SeabedScenario& scen, std::string& warn, bool build_core = true);
}
