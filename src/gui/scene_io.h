// scene_io.h — save / restore a full simulation "scene" to a single self-contained .scn file.
//
// A .scn is EVERYTHING needed to reconstruct AND continue a run: the setup (grid, BCs, params,
// the config obstacle, provenance), the loaded STEP model (embedded triangle mesh, display-only),
// and all live field DATA at the saved step (velocities, pressure, the flow solid mask). Loading a
// .scn tears the sim down and rebuilds it through the SAME build path (sim_setup make_core), then
// injects the saved fields — so the flow resumes exactly where it was, not from t=0.
//
// Two uses of the ONE format:
//   * a user "Save Scene As" writes  <name>.scn                     (added to Recent Files)
//   * an auto-save writes            <name>.<step>.scn  every N steps (a restart point; not Recent)
// Both are full, independently loadable .scn files.
//
// File layout (little-endian, the workstation is x64):
//   magic "SCOURSCN"  (8 bytes)
//   uint32 version
//   uint64 json_len ; json_len bytes of UTF-8 metadata (all scalars/structs + blob directory flags)
//   then a sequence of TLV blobs until EOF, each:
//     uint32 key_len ; key bytes ; uint8 dtype (0=f32,1=f64,2=u8,3=u32) ; uint64 count ; count*sz data
// Big arrays (u/v/w/p, mesh positions/normals) are stored as float32 (the user-chosen precision);
// masks/indices as u8/u32. Qt-free (std streams + nlohmann json) so it is unit-testable and never
// pulls Qt into the physics libs.
#pragma once

#include "core/fluid/channel_bc.h"
#include "core/fluid/mac_grid.h"
#include "core/geometry/model_placement.h"
#include "core/geometry/tri_mesh.h"
#include "gui/sim_setup.h" // SimRecipe

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace windcfd::gui
{
	// The live DYNAMIC state at a saved step, gathered from the core on the worker thread.
	// Field arrays are host DOUBLE here (as read from the device); the writer downcasts the big ones.
	struct CheckpointState
	{
		windcfd::core::MacGrid grid;
		long long steps = 0;
		double sim_time = 0.0;
		windcfd::core::ChannelBC bc;                     // live inlet speed/profile at save time
		int solid_mode = windcfd::core::SOLID_NOSLIP;    // interior-obstacle surface mode (= bc.solid_mode)
		bool bed_inlet_mask = false;

		std::vector<double> u, v, w, p;                // MAC face velocities + cell pressure
		std::vector<unsigned char> solid;              // live flow solid mask (p_count, 1=solid)
	};
	using CheckpointStatePtr = std::shared_ptr<CheckpointState>;

	// The static scene definition: enough to reconstruct the sim + display the STEP model.
	struct SceneDefinition
	{
		SimRecipe recipe;                              // grid, bc, pr, base_solid, init, provenance, info
		bool has_mesh = false;
		windcfd::core::TriMesh mesh;                     // display STEP mesh (metres); empty ⇒ none
		windcfd::core::ModelPlacement place;             // where the mesh sits (display translate)
	};

	// A complete .scn = the static definition + the full dynamic state.
	struct SceneFile
	{
		SceneDefinition def;
		CheckpointState state;
	};

	// Lightweight header for the restart-point picker (field blobs NOT read).
	struct SceneHeader
	{
		int version = 0;
		std::string name;
		long long steps = 0;
		double sim_time = 0.0;
		int nx = 0, ny = 0, nz = 0;
		double h = 0.0;
		bool has_mesh = false;
	};

	// Write / read a full .scn. Returns false and fills `warn` on any I/O or format error.
	bool write_scene(const std::string& path, const SceneFile& sf, std::string& warn);
	bool read_scene(const std::string& path, SceneFile& sf, std::string& warn);
	bool read_scene_header(const std::string& path, SceneHeader& hdr, std::string& warn);

	// Enumerate the checkpoint siblings of a scene path: the base <name>.scn plus every
	// <name>.<step>.scn next to it, sorted by step ascending (base scene reports its own header step).
	struct CheckpointRef { std::string path; long long step = 0; };
	std::vector<CheckpointRef> list_checkpoints(const std::string& scene_path);

	// Strip a .scn / .<digits>.scn suffix from `scene_path` and return "<base>.<step>.scn".
	std::string checkpoint_path_for(const std::string& scene_path, long long step);

	// True if `path`'s file name matches the numbered-checkpoint pattern <name>.<digits>.scn
	// (as opposed to a base <name>.scn). Numbered checkpoints are NOT added to Recent Files.
	bool is_numbered_checkpoint(const std::string& path);
}
