// settle.h — the rigid-body drop/settle "preparation phase" API (PLAN G3).
//
// Scatter N scaled copies of a printed protection unit (e.g. Holcim XStone) above the bed, let
// them fall + settle under gravity with body–body + body–bed contact, and return the RESTING pose
// of each as an affine core::ModelPlacement — ready to hand straight to voxelize_mesh(). The pile
// becomes the rigid structure of the morphodynamic run.
//
// This header is deliberately free of any Jolt Physics include (exactly as step_import.h is free
// of OpenCascade): all Jolt usage is isolated in the scour_physics static lib's implementation TU,
// so libscour + every physics gate stay Jolt-free. It depends only on the OCC-free/Qt-free geometry
// contracts (TriMesh, ModelPlacement).
//
// ⚠ Settling is FULLY DECOUPLED from the CUDA fluid: it runs on the CPU and only emits transforms.
// The COLLISION geometry is a convex compound (from load_step_solids(), or a single whole-mesh hull
// as a bootstrap); the flow still voxelizes the TRUE perforated display mesh, so the animal holes
// are preserved downstream — only the settling contact uses the convex pieces.
//
// Two entry points: settle_blocks() runs the whole drop to rest in one call (headless); SettleWorld
// steps it incrementally so the GUI can ANIMATE the drop (and, later, G3.3 can re-settle a persistent
// pile as the bed erodes).
//
// Units: SI (metres, seconds, kg, m/s²). Determinism: fixed `seed` + identical inputs ⇒ identical
// output (Jolt runs single-threaded here), so an experiment (this pile) is reproducible.
#pragma once

#include "core/geometry/model_placement.h"
#include "core/geometry/tri_mesh.h"

#include <array>
#include <memory>
#include <vector>

namespace scour::core
{
	// A rigid protection unit as a UNION OF CONVEX PIECES for collision. Each piece's vertex cloud
	// becomes one Jolt ConvexHullShape; the union is the body's compound collision shape. Build it
	// from load_step_solids() (one convex piece per STEP solid), or, as a bootstrap before the
	// pre-cut multi-solid STEP exists, from a single whole-mesh entry (its convex hull). The
	// vertices are taken in the unit's own local frame (metres, before any per-instance scale/pose).
	struct ConvexShape
	{
		std::vector<TriMesh> pieces; // ≥1 near-convex collision pieces (vertices → hull point clouds)

		bool empty() const { return pieces.empty(); }
	};

	// One block to drop: a uniform scale applied to `ConvexShape` and an initial world placement
	// (M = rotation·scale, t = spawn position, metres) from which it falls. Scattered spawns come
	// from scatter_blocks(); the caller may also hand-author instances.
	struct BlockInstance
	{
		double scale = 1.0;      // uniform size multiplier (mass ∝ density·V·scale³)
		ModelPlacement pose;     // initial world transform before settling
	};

	// A rigid ground for the blocks to land on. v1 = a horizontal plane at z = bed_z (the sand-bed
	// top). The developing scour hole is applied later via SettleWorld::resettleOnBed (a full surface).
	struct GroundSpec
	{
		double bed_z = 0.0; // plane height (metres); blocks rest with min-z ≈ bed_z
	};

	// A sampled bed-elevation SURFACE (world Z-up) for the persistent re-settle (PLAN G3.3): as the
	// morphodynamic bed erodes, the flat settling ground is replaced by this surface so the blocks sink
	// into the developing scour hole. Sample (i,j) sits at world (x0 + i·h, y0 + j·h) with height
	// z[j·nx + i]. Downsample the flow's z_b as coarse as the block scale to keep the collision mesh cheap.
	struct BedField
	{
		int nx = 0, ny = 0;        // grid dimensions (≥ 2 each)
		double x0 = 0.0, y0 = 0.0; // world position of sample (0,0) [m]
		double h = 0.05;           // sample spacing [m]
		std::vector<float> z;      // nx·ny bed elevation, row-major (j·nx + i) [m]
	};

	// Settling controls. Defaults are nominal marine-concrete values; the calibrated set is the
	// orchestrator's/MH's call, never an agent's (mirrors the M6 rule for physics constants).
	struct SettleParams
	{
		double density = 2400.0;      // marine concrete, kg/m³
		double gravity = 9.81;        // m/s² downward (−z)
		double friction = 0.6;        // Coulomb μ, block–block and block–bed (interlocking driver)
		double restitution = 0.1;     // low — rocks barely bounce
		double linear_damping = 0.05; // 1/s
		double angular_damping = 0.05;// 1/s
		double dt = 1.0 / 120.0;      // fixed physics tick, s (each is subdivided into `collision_steps`)
		int max_steps = 8000;         // hard cap (~67 s of ticks) if a body never sleeps

		// Total bodies this world will EVER hold, incl. ones added later via SettleWorld::addBlock (a
		// timed release starts from an empty world, so instances.size() is 0 at construction). Jolt's
		// body-pair + contact-constraint buffers are sized from this: too small ⇒ Jolt SILENTLY DROPS
		// contacts ⇒ blocks fall through each other and the floor. Set it to the final unit count (e.g.
		// 800) BEFORE creating the world; 0 ⇒ size from instances.size(). See settle.cpp for the caps.
		int reserve_bodies = 0;

		// --- Solve quality (units dropped 5 m above the domain hit the sand at ~10–13 m/s, so both
		// tunnelling through the floor and resting overlap need active control) ------------------------
		int collision_steps = 2;      // Jolt sub-steps PER dt: subdivides each tick for a higher-quality
		                              // integrate+solve (the "add substeps" knob). Sub-tick = dt/this.
		bool continuous = true;       // CCD (EMotionQuality::LinearCast) on every block — a swept
		                              // collision test so a fast drop CANNOT tunnel through the floor or
		                              // another unit in one tick (the fix for penetrating-the-floor).
		double penetration_slop = 5.0e-3; // m: allowed resting overlap (Jolt default 0.02 m = 2 cm is
		                                  // coarse for ≤ dm units → visible interpenetration; tighten it).
		int velocity_steps = 12;      // contact-solver velocity iterations (Jolt default 10)
		int position_steps = 4;       // contact-solver POSITION iterations — pushes out residual
		                              // penetration each tick (Jolt default 2); the anti-overlap knob.

		double sleep_linear = 5.0e-3; // m/s: Jolt's point-velocity sleep threshold (folds in rotation)
		double sleep_angular = 1.0e-2;// rad/s (informational; Jolt sleeps on point velocity)
		unsigned seed = 12345u;       // reproducible scatter + settle
	};

	// Outcome of a settle run: the resting pose per instance (SAME index order as the input
	// `instances`) plus diagnostics the G3.1 gate checks (converged? all above bed? no deep overlap?).
	struct SettleResult
	{
		std::vector<ModelPlacement> placements; // resting world transform per block
		int steps = 0;                          // physics steps actually taken
		bool all_asleep = false;                // true ⇒ every body slept before max_steps
		double min_z = 0.0;                     // lowest resting body-vertex z (≈ bed_z at rest)
		double max_penetration = 0.0;           // deepest body–body overlap at rest (m); ~0 when clean
	};

	// Generate `count` blocks scattered over a square footprint of side `footprint` centred at
	// (cx, cy), each dropped from a random height in [drop_z_lo, drop_z_hi] with a random 3-D
	// orientation and a uniform scale in [scale_min, scale_max]. `shape_bbox_size` (metres) sizes the
	// spawn spacing so blocks don't start interpenetrating. The returned list is SORTED LARGEST-FIRST
	// (descending scale), so a timed release drops the big units first and the smaller ones later (an
	// armour-layer pour). Deterministic for a fixed `seed`.
	std::vector<BlockInstance> scatter_blocks(int count, double cx, double cy, double footprint,
		double drop_z_lo, double drop_z_hi, double scale_min, double scale_max,
		const std::array<float, 3>& shape_bbox_size, unsigned seed);

	// Drop `instances` of `shape` under gravity onto `ground` with body–body + body–ground contact,
	// stepping until every body sleeps or `params.max_steps` is reached. Returns the resting pose per
	// instance (+ diagnostics). One-shot convenience wrapper around SettleWorld.
	SettleResult settle_blocks(const ConvexShape& shape, const std::vector<BlockInstance>& instances,
		const GroundSpec& ground, const SettleParams& params);

	// Stateful settling world for LIVE, incremental stepping — the animated drop (the GUI steps it a
	// few substeps per frame and draws the mesh at each body's current pose) and, later (PLAN G3.3),
	// a persistent pile that re-settles as the bed erodes. Pimpl-hides ALL Jolt state so this header
	// stays Jolt-free. The ONLY caller-visible physics object; its Impl is the one place linking Jolt.
	class SettleWorld
	{
	public:
		SettleWorld(const ConvexShape& shape, const std::vector<BlockInstance>& instances,
			const GroundSpec& ground, const SettleParams& params);
		~SettleWorld();
		SettleWorld(const SettleWorld&) = delete;
		SettleWorld& operator=(const SettleWorld&) = delete;

		// Add one block to the LIVE world (a dynamic body dropped from block.pose at block.scale). Used
		// by the timed release to drop units one at a time, largest-first, so the big ones settle before
		// the smaller ones land among/on them. Appends to the block list (placements()/result() grow to
		// match). Returns false if the block's collision shape was degenerate (index still appended).
		bool addBlock(const BlockInstance& block);

		// Advance `substeps` physics steps (each params.dt). Returns true while the pile is still
		// settling; false once every body has slept OR the max-steps budget is spent. Loop on it.
		bool step(int substeps = 1);

		bool all_asleep() const;                 // every body asleep (settled)
		int steps() const;                       // physics steps taken so far
		std::vector<ModelPlacement> placements(); // current world pose per block (live while settling)
		SettleResult result();                   // placements + diagnostics (min_z, penetration, …)

		// Replace the ground with a bed-elevation SURFACE (world Z-up triangle mesh from `bed`), wake all
		// blocks, and re-settle them onto it — they sink into any scour hole that has developed. Steps
		// until asleep or the budget; returns the new resting placements + diagnostics. The persistent
		// world keeps the blocks between calls, so this continues from their CURRENT poses (PLAN G3.3).
		SettleResult resettleOnBed(const BedField& bed);

	private:
		struct Impl;
		std::unique_ptr<Impl> impl_;
	};
}
