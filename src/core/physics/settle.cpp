// settle.cpp — the ONLY translation unit that links Jolt Physics (v5.5.0), exactly as
// step_import.cpp is the only one that links OpenCascade. Implements the PLAN G3 drop/settle
// preparation phase: scatter N scaled protection units above the bed, drop them under gravity with
// body–body + body–bed contact, and read back each body's transform as a core::ModelPlacement.
//
// The stateful SettleWorld::Impl holds the live Jolt world (so the GUI can step it a few substeps
// per frame and animate the drop, and G3.3 can re-settle a persistent pile); settle_blocks() is a
// one-shot wrapper that steps it to rest. Jolt boilerplate follows the canonical HelloWorld template
// (RegisterDefaultAllocator → Factory → RegisterTypes → PhysicsSystem).
//
// We run SINGLE-THREADED (JobSystemSingleThreaded) so a fixed seed reproduces the same pile bit-for-
// bit — Jolt is deterministic for a given build + thread schedule, and single-threaded removes the
// thread-count variable entirely (the G3.1 gate checks it). Contact callbacks then fire on the
// calling thread, so the penetration tracker needs no synchronisation.
//
// Coordinate frame: the whole simulator is Z-up, so gravity is (0,0,−g) and the floor normal is +Z
// (Jolt is frame-agnostic — gravity is just a vector). Scale handling: a block's per-instance `scale`
// is baked into the CONVEX-HULL point cloud and folded into the returned placement as M = scale·R, so
// the voxelizer — which applies world(v)=M·v+t to the UNSCALED display mesh — reproduces the pose.
#include "core/physics/settle.h"

#include <Jolt/Jolt.h>

#include <Jolt/Core/Factory.h>
#include <Jolt/Core/JobSystemSingleThreaded.h>
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Body/BodyType.h>
#include <Jolt/Physics/Collision/ContactListener.h>
#include <Jolt/Geometry/IndexedTriangle.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Collision/Shape/ConvexHullShape.h>
#include <Jolt/Physics/Collision/Shape/MeshShape.h>
#include <Jolt/Physics/Collision/Shape/StaticCompoundShape.h>
#include <Jolt/Physics/PhysicsSettings.h>
#include <Jolt/Physics/PhysicsSystem.h>
#include <Jolt/RegisterTypes.h>

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <memory>
#include <random>
#include <vector>

JPH_SUPPRESS_WARNINGS

using namespace JPH;
using namespace JPH::literals; // the _r Real-literal suffix

namespace
{
	// Jolt calls Trace for warnings; route it to a no-op so the library stays silent.
	void TraceNoop(const char*, ...) {}

	// Object + broad-phase layers: the standard two-layer setup (static bed vs moving blocks).
	namespace Layers
	{
		static constexpr ObjectLayer NON_MOVING = 0;
		static constexpr ObjectLayer MOVING = 1;
		static constexpr ObjectLayer NUM_LAYERS = 2;
	}

	namespace BroadPhaseLayers
	{
		static constexpr BroadPhaseLayer NON_MOVING(0);
		static constexpr BroadPhaseLayer MOVING(1);
		static constexpr uint NUM_LAYERS(2);
	}

	class ObjectLayerPairFilterImpl : public ObjectLayerPairFilter
	{
	public:
		bool ShouldCollide(ObjectLayer o1, ObjectLayer o2) const override
		{
			switch (o1)
			{
			case Layers::NON_MOVING: return o2 == Layers::MOVING; // bed only collides with blocks
			case Layers::MOVING: return true;                     // blocks collide with everything
			default: return false;
			}
		}
	};

	class BPLayerInterfaceImpl final : public BroadPhaseLayerInterface
	{
	public:
		BPLayerInterfaceImpl()
		{
			mObjectToBroadPhase[Layers::NON_MOVING] = BroadPhaseLayers::NON_MOVING;
			mObjectToBroadPhase[Layers::MOVING] = BroadPhaseLayers::MOVING;
		}
		uint GetNumBroadPhaseLayers() const override { return BroadPhaseLayers::NUM_LAYERS; }
		BroadPhaseLayer GetBroadPhaseLayer(ObjectLayer l) const override { return mObjectToBroadPhase[l]; }
#if defined(JPH_EXTERNAL_PROFILE) || defined(JPH_PROFILE_ENABLED)
		const char* GetBroadPhaseLayerName(BroadPhaseLayer) const override { return "layer"; }
#endif
	private:
		BroadPhaseLayer mObjectToBroadPhase[Layers::NUM_LAYERS];
	};

	class ObjectVsBroadPhaseLayerFilterImpl : public ObjectVsBroadPhaseLayerFilter
	{
	public:
		bool ShouldCollide(ObjectLayer l1, BroadPhaseLayer l2) const override
		{
			switch (l1)
			{
			case Layers::NON_MOVING: return l2 == BroadPhaseLayers::MOVING;
			case Layers::MOVING: return true;
			default: return false;
			}
		}
	};

	// Records the deepest body–body/body–bed overlap seen in the CURRENT step. Callbacks fire on the
	// calling thread (single-threaded job system) ⇒ no synchronisation needed. The settle loop reads
	// `cur` right after each Update, so the final value is the resting penetration (≈0 when clean).
	class PenetrationTracker : public ContactListener
	{
	public:
		float cur = 0.0f;
		void OnContactAdded(const Body&, const Body&, const ContactManifold& m, ContactSettings&) override { cur = std::max(cur, m.mPenetrationDepth); }
		void OnContactPersisted(const Body&, const Body&, const ContactManifold& m, ContactSettings&) override { cur = std::max(cur, m.mPenetrationDepth); }
	};

	// One-time, process-lifetime Jolt global init (allocator + factory + type registration). C++11
	// guarantees the function-local static is initialised exactly once, thread-safely. We never
	// UnregisterTypes/delete the factory — process lifetime is fine for a tool, and it lets multiple
	// SettleWorlds (G3.3 re-settles) be created without re-registering.
	void ensure_jolt_initialized()
	{
		static const bool once = []() {
			RegisterDefaultAllocator();
			Trace = TraceNoop;
			Factory::sInstance = new Factory();
			RegisterTypes();
			return true;
		}();
		(void)once;
	}

	// A do-nothing member whose sole job is to run ensure_jolt_initialized() from its constructor.
	// Declared as the FIRST Jolt-touching member of SettleWorld::Impl so the allocator/factory are
	// registered BEFORE the Jolt members (TempAllocatorImpl, PhysicsSystem, …) are constructed —
	// member initialisers run before the constructor body, so a body-level init call would be too late
	// (TempAllocatorImpl would allocate through an unregistered allocator → crash).
	struct JoltInit
	{
		JoltInit() { ensure_jolt_initialized(); }
	};

	// Pure rotation Quat from a ModelPlacement's linear part (columns normalised to strip any scale).
	Quat rotation_from_placement(const scour::core::ModelPlacement& p)
	{
		Vec3 c0((float)p.m[0], (float)p.m[3], (float)p.m[6]); // column 0 = M·e_x
		Vec3 c1((float)p.m[1], (float)p.m[4], (float)p.m[7]);
		Vec3 c2((float)p.m[2], (float)p.m[5], (float)p.m[8]);
		c0 = c0.Length() > 1.0e-12f ? c0.Normalized() : Vec3(1, 0, 0);
		c1 = c1.Length() > 1.0e-12f ? c1.Normalized() : Vec3(0, 1, 0);
		c2 = c2.Length() > 1.0e-12f ? c2.Normalized() : Vec3(0, 0, 1);
		Mat44 m = Mat44::sIdentity();
		m.SetColumn3(0, c0);
		m.SetColumn3(1, c1);
		m.SetColumn3(2, c2);
		return m.GetQuaternion().Normalized();
	}

	// Build a block's collision shape: each convex piece → a ConvexHullShape (points × scale), unioned
	// into a StaticCompoundShape when there is more than one. Returns an invalid ref on failure.
	ShapeRefC build_block_shape(const scour::core::ConvexShape& shape, double scale, double density)
	{
		std::vector<ShapeRefC> hulls;
		hulls.reserve(shape.pieces.size());
		for (const scour::core::TriMesh& piece : shape.pieces)
		{
			if (piece.vertex_count() < 4) continue; // need a 3-D hull
			ConvexHullShapeSettings hs;
			hs.mPoints.reserve((uint)piece.vertex_count());
			for (std::size_t v = 0; v < piece.vertex_count(); ++v)
				hs.mPoints.push_back(Vec3((float)(piece.positions[3 * v + 0] * scale),
					(float)(piece.positions[3 * v + 1] * scale),
					(float)(piece.positions[3 * v + 2] * scale)));
			hs.mDensity = (float)density;
			hs.mMaxConvexRadius = 0.0f; // keep sharp edges (blocks are ≥ cm-scale; no rounding margin)
			ShapeSettings::ShapeResult r = hs.Create();
			if (r.IsValid()) hulls.push_back(r.Get());
		}

		if (hulls.empty()) return ShapeRefC();
		if (hulls.size() == 1) return hulls.front();

		StaticCompoundShapeSettings cs;
		for (const ShapeRefC& h : hulls) cs.AddShape(Vec3::sZero(), Quat::sIdentity(), h);
		ShapeSettings::ShapeResult r = cs.Create();
		return r.IsValid() ? r.Get() : ShapeRefC();
	}
}

namespace scour::core
{
	// ------------------------------------------------------------------------------------------
	// SettleWorld::Impl — the live Jolt world. Members declared so that physics_system (which holds
	// references to the filters + the contact listener) is destroyed BEFORE them (reverse order).
	// ------------------------------------------------------------------------------------------
	struct SettleWorld::Impl
	{
		ConvexShape shape;
		std::vector<BlockInstance> instances;
		GroundSpec ground;
		SettleParams params;

		JoltInit jolt_init; // MUST precede every Jolt member below (registers the allocator first)
		std::unique_ptr<TempAllocatorImpl> temp_allocator; // sized from the body count in the ctor body
		JobSystemSingleThreaded job_system{ cMaxPhysicsJobs }; // deterministic (no worker threads)
		BPLayerInterfaceImpl broad_phase_layer_interface;
		ObjectVsBroadPhaseLayerFilterImpl object_vs_broadphase_layer_filter;
		ObjectLayerPairFilterImpl object_vs_object_layer_filter;
		PenetrationTracker tracker;
		PhysicsSystem physics_system;

		std::vector<BodyID> ids;
		std::vector<bool> alive;
		BodyID ground_id;   // current static ground body (flat floor, then the bed-surface mesh)
		int steps = 0;

		Impl(const ConvexShape& s, const std::vector<BlockInstance>& insts, const GroundSpec& g, const SettleParams& p)
			: shape(s), instances(insts), ground(g), params(p)
		{
			// Jolt is already initialised by the jolt_init member (constructed before this body).
			// CAPACITY sizing: the buffers must hold the FINAL body count — params.reserve_bodies covers
			// bodies added later via addBlock, so a timed release (which starts with an EMPTY instances
			// list, cap ← reserve_bodies) still gets big-enough buffers. A dense pile has many broad-phase
			// pairs per body, and if cMaxBodyPairs / cMaxContactConstraints overflow Jolt SILENTLY DROPS
			// contacts (Trace is a no-op) ⇒ interpenetration + floor tunnelling. Size pairs generously
			// (×24 — an AABB overlaps many neighbours in a pile) and constraints (×12). NOTE `cap` is the
			// capacity ONLY; the bodies created NOW come from `n_inst` = instances.size() (see the loop).
			const uint n_inst = (uint)instances.size();
			const uint cap = std::max<uint>(n_inst, (uint)std::max(0, params.reserve_bodies));
			const uint cMaxBodies = std::max<uint>(1024, cap + 256);
			const uint cMaxBodyPairs = std::max<uint>(8192, cap * 24);
			const uint cMaxContactConstraints = std::max<uint>(4096, cap * 12);

			// Per-step scratch: grows with the pile (constraint/island arrays). 16 MB floor + ~24 KB per
			// body so a large release (e.g. 800–few-k units) never exhausts the arena mid-Update.
			temp_allocator = std::make_unique<TempAllocatorImpl>((uint)std::max<size_t>(16u * 1024 * 1024, (size_t)cap * 24 * 1024));

			physics_system.Init(cMaxBodies, 0, cMaxBodyPairs, cMaxContactConstraints,
				broad_phase_layer_interface, object_vs_broadphase_layer_filter, object_vs_object_layer_filter);
			physics_system.SetGravity(Vec3(0.0f, 0.0f, -(float)params.gravity));

			PhysicsSettings ps = physics_system.GetPhysicsSettings();
			ps.mPointVelocitySleepThreshold = (float)params.sleep_linear;
			ps.mTimeBeforeSleep = 0.5f;
			// Higher-quality contact solve: more velocity+position iterations and a tighter penetration
			// slop so ≤ dm units stop overlapping visibly at rest (Jolt's 0.02 m default is coarse here).
			ps.mNumVelocitySteps = (uint)std::max(1, params.velocity_steps);
			ps.mNumPositionSteps = (uint)std::max(1, params.position_steps);
			ps.mPenetrationSlop = (float)params.penetration_slop;
			physics_system.SetPhysicsSettings(ps);

			physics_system.SetContactListener(&tracker);

			BodyInterface& bi = physics_system.GetBodyInterface();

			// Static floor: a large thin box whose TOP face sits at ground.bed_z.
			const float floor_half = 0.5f;
			BoxShapeSettings floor_shape_settings(Vec3(1000.0f, 1000.0f, floor_half));
			floor_shape_settings.SetEmbedded();
			ShapeRefC floor_shape = floor_shape_settings.Create().Get();
			BodyCreationSettings floor_settings(floor_shape,
				RVec3(0.0_r, 0.0_r, (Real)(ground.bed_z - floor_half)), Quat::sIdentity(),
				EMotionType::Static, Layers::NON_MOVING);
			floor_settings.mFriction = (float)params.friction;
			if (Body* floor = bi.CreateBody(floor_settings))
			{
				bi.AddBody(floor->GetID(), EActivation::DontActivate);
				ground_id = floor->GetID();
			}

			// Dynamic blocks: create only the bodies present NOW (n_inst); the buffers already reserve
			// headroom for any added later via addBlock. An empty world (timed release) skips this loop.
			ids.assign(n_inst, BodyID());
			alive.assign(n_inst, false);
			for (uint i = 0; i < n_inst; ++i)
			{
				ShapeRefC block_shape = build_block_shape(shape, instances[i].scale, params.density);
				if (block_shape.GetPtr() == nullptr) continue; // degenerate: skip (keeps indices aligned)

				const ModelPlacement& pose = instances[i].pose;
				BodyCreationSettings bcs(block_shape,
					RVec3((Real)pose.tx, (Real)pose.ty, (Real)pose.tz), rotation_from_placement(pose),
					EMotionType::Dynamic, Layers::MOVING);
				bcs.mFriction = (float)params.friction;
				bcs.mRestitution = (float)params.restitution;
				bcs.mLinearDamping = (float)params.linear_damping;
				bcs.mAngularDamping = (float)params.angular_damping;
				bcs.mAllowSleeping = true;
				// CCD: a swept collision test so a unit dropped fast (≈10–13 m/s from 5 m up) can't tunnel
				// through the floor or another unit within one tick. Auto-engages only when fast (cheap).
				bcs.mMotionQuality = params.continuous ? EMotionQuality::LinearCast : EMotionQuality::Discrete;

				BodyID id = bi.CreateAndAddBody(bcs, EActivation::Activate);
				if (!id.IsInvalid())
				{
					ids[i] = id;
					alive[i] = true;
				}
			}

			physics_system.OptimizeBroadPhase();
			// Bodies stay added for the world's lifetime; the PhysicsSystem destructor frees them.
		}

		bool add_block(const BlockInstance& b)
		{
			BodyInterface& bi = physics_system.GetBodyInterface();
			ShapeRefC block_shape = build_block_shape(shape, b.scale, params.density);
			if (block_shape.GetPtr() == nullptr)
			{
				instances.push_back(b); ids.push_back(BodyID()); alive.push_back(false);
				return false;
			}
			BodyCreationSettings bcs(block_shape,
				RVec3((Real)b.pose.tx, (Real)b.pose.ty, (Real)b.pose.tz), rotation_from_placement(b.pose),
				EMotionType::Dynamic, Layers::MOVING);
			bcs.mFriction = (float)params.friction;
			bcs.mRestitution = (float)params.restitution;
			bcs.mLinearDamping = (float)params.linear_damping;
			bcs.mAngularDamping = (float)params.angular_damping;
			bcs.mAllowSleeping = true;
			bcs.mMotionQuality = params.continuous ? EMotionQuality::LinearCast : EMotionQuality::Discrete; // CCD (see ctor)
			BodyID id = bi.CreateAndAddBody(bcs, EActivation::Activate);
			instances.push_back(b);
			ids.push_back(id);
			alive.push_back(!id.IsInvalid());
			return !id.IsInvalid();
		}

		bool all_asleep() const
		{
			return physics_system.GetNumActiveBodies(EBodyType::RigidBody) == 0;
		}

		bool advance(int substeps)
		{
			const int csteps = std::max(1, params.collision_steps);
			for (int s = 0; s < substeps && steps < params.max_steps; ++s)
			{
				if (all_asleep()) break;
				tracker.cur = 0.0f;
				physics_system.Update((float)params.dt, csteps, temp_allocator.get(), &job_system);
				++steps;
			}
			return !all_asleep() && steps < params.max_steps;
		}

		std::vector<ModelPlacement> placements()
		{
			std::vector<ModelPlacement> out(instances.size());
			BodyInterface& bi = physics_system.GetBodyInterface();
			for (std::size_t i = 0; i < instances.size(); ++i)
			{
				ModelPlacement pl = instances[i].pose; // fallback = spawn pose if the body was skipped
				if (alive[i])
				{
					const RMat44 wt = bi.GetWorldTransform(ids[i]);
					const Vec3 t = wt.GetTranslation();
					const Vec3 rc0 = wt.GetColumn3(0), rc1 = wt.GetColumn3(1), rc2 = wt.GetColumn3(2);
					const double sc = instances[i].scale;
					pl.tx = t.GetX(); pl.ty = t.GetY(); pl.tz = t.GetZ();
					// M = sc·R (row-major); R column j = rc_j = (R0j,R1j,R2j).
					pl.m[0] = sc * rc0.GetX(); pl.m[1] = sc * rc1.GetX(); pl.m[2] = sc * rc2.GetX();
					pl.m[3] = sc * rc0.GetY(); pl.m[4] = sc * rc1.GetY(); pl.m[5] = sc * rc2.GetY();
					pl.m[6] = sc * rc0.GetZ(); pl.m[7] = sc * rc1.GetZ(); pl.m[8] = sc * rc2.GetZ();
				}
				out[i] = pl;
			}
			return out;
		}

		SettleResult result()
		{
			SettleResult r;
			r.placements = placements();
			r.steps = steps;
			r.all_asleep = all_asleep();
			r.max_penetration = tracker.cur;

			double min_z = ground.bed_z;
			bool any = false;
			for (std::size_t i = 0; i < instances.size(); ++i)
			{
				if (!alive[i]) continue;
				const ModelPlacement& pl = r.placements[i];
				for (const TriMesh& piece : shape.pieces)
					for (std::size_t v = 0; v < piece.vertex_count(); ++v)
					{
						const double vx = piece.positions[3 * v + 0], vy = piece.positions[3 * v + 1], vz = piece.positions[3 * v + 2];
						const double wz = pl.m[6] * vx + pl.m[7] * vy + pl.m[8] * vz + pl.tz;
						if (!any || wz < min_z) { min_z = wz; any = true; }
					}
			}
			r.min_z = min_z;
			return r;
		}

		SettleResult resettle(const BedField& bed)
		{
			BodyInterface& bi = physics_system.GetBodyInterface();

			// Replace the ground with a triangle-mesh SURFACE of the bed (world Z-up), so the blocks
			// re-settle onto the new bed and sink into any scour hole that has developed.
			if (bed.nx >= 2 && bed.ny >= 2 && (int)bed.z.size() >= bed.nx * bed.ny)
			{
				Array<Float3> verts;
				verts.reserve((size_t)bed.nx * bed.ny);
				for (int j = 0; j < bed.ny; ++j)
					for (int i = 0; i < bed.nx; ++i)
						verts.push_back(Float3((float)(bed.x0 + i * bed.h), (float)(bed.y0 + j * bed.h),
							bed.z[(size_t)j * bed.nx + i]));

				Array<IndexedTriangle> tris;
				tris.reserve((size_t)(bed.nx - 1) * (bed.ny - 1) * 2);
				for (int j = 0; j < bed.ny - 1; ++j)
					for (int i = 0; i < bed.nx - 1; ++i)
					{
						const uint32 a = (uint32)(j * bed.nx + i);
						const uint32 b = (uint32)(j * bed.nx + i + 1);
						const uint32 c = (uint32)((j + 1) * bed.nx + i);
						const uint32 d = (uint32)((j + 1) * bed.nx + i + 1);
						tris.push_back(IndexedTriangle(a, b, c, 0)); // up-facing (normal ≈ +Z)
						tris.push_back(IndexedTriangle(b, d, c, 0));
					}

				MeshShapeSettings mss(verts, tris);
				mss.SetEmbedded();
				ShapeSettings::ShapeResult mr = mss.Create();
				if (mr.IsValid())
				{
					if (!ground_id.IsInvalid()) { bi.RemoveBody(ground_id); bi.DestroyBody(ground_id); ground_id = BodyID(); }
					BodyCreationSettings gcs(mr.Get(), RVec3::sZero(), Quat::sIdentity(), EMotionType::Static, Layers::NON_MOVING);
					gcs.mFriction = (float)params.friction;
					if (Body* g = bi.CreateBody(gcs))
					{
						bi.AddBody(g->GetID(), EActivation::DontActivate);
						ground_id = g->GetID();
					}
				}
			}

			// Wake every block and re-settle onto the new ground (continues from their current poses).
			for (std::size_t i = 0; i < ids.size(); ++i)
				if (alive[i]) bi.ActivateBody(ids[i]);

			const int csteps = std::max(1, params.collision_steps);
			for (int local = 0; local < params.max_steps; ++local)
			{
				if (physics_system.GetNumActiveBodies(EBodyType::RigidBody) == 0) break;
				tracker.cur = 0.0f;
				physics_system.Update((float)params.dt, csteps, temp_allocator.get(), &job_system);
				++steps;
			}
			return result();
		}
	};

	SettleWorld::SettleWorld(const ConvexShape& shape, const std::vector<BlockInstance>& instances,
		const GroundSpec& ground, const SettleParams& params)
		: impl_(std::make_unique<Impl>(shape, instances, ground, params))
	{
	}

	SettleWorld::~SettleWorld() = default;

	bool SettleWorld::addBlock(const BlockInstance& block) { return impl_->add_block(block); }
	bool SettleWorld::step(int substeps) { return impl_->advance(substeps); }
	bool SettleWorld::all_asleep() const { return impl_->all_asleep(); }
	int SettleWorld::steps() const { return impl_->steps; }
	std::vector<ModelPlacement> SettleWorld::placements() { return impl_->placements(); }
	SettleResult SettleWorld::result() { return impl_->result(); }
	SettleResult SettleWorld::resettleOnBed(const BedField& bed) { return impl_->resettle(bed); }

	// ------------------------------------------------------------------------------------------
	std::vector<BlockInstance> scatter_blocks(int count, double cx, double cy, double footprint,
		double drop_z_lo, double drop_z_hi, double scale_min, double scale_max,
		const std::array<float, 3>& shape_bbox_size, unsigned seed)
	{
		std::vector<BlockInstance> out;
		if (count <= 0) return out;
		out.reserve((std::size_t)count);

		std::mt19937 rng(seed);
		std::uniform_real_distribution<double> u01(0.0, 1.0);

		// Jittered grid over the footprint so blocks do not start interpenetrating. Cell jitter is
		// clamped by the largest block extent so two neighbours can never overlap at spawn.
		const int cols = std::max(1, (int)std::ceil(std::sqrt((double)count)));
		const double cell = footprint / cols;
		const double maxdim = std::max({ shape_bbox_size[0], shape_bbox_size[1], shape_bbox_size[2] }) * scale_max;
		const double jit = std::max(0.0, cell - maxdim) * 0.5; // ± this stays separated

		for (int i = 0; i < count; ++i)
		{
			const int c = i % cols, r = i / cols;
			BlockInstance bi;
			bi.scale = scale_min + u01(rng) * (scale_max - scale_min);

			const double gx = cx - 0.5 * footprint + (c + 0.5) * cell + (u01(rng) * 2.0 - 1.0) * jit;
			const double gy = cy - 0.5 * footprint + (r + 0.5) * cell + (u01(rng) * 2.0 - 1.0) * jit;
			const double gz = drop_z_lo + u01(rng) * std::max(0.0, drop_z_hi - drop_z_lo);
			bi.pose.tx = gx;
			bi.pose.ty = gy;
			bi.pose.tz = gz;

			// Uniformly-random 3-D orientation (Shoemake) → row-major rotation matrix, so blocks
			// tumble and land on random faces/edges the way dropped rocks do.
			const double u1 = u01(rng), u2 = u01(rng), u3 = u01(rng);
			const double s1 = std::sqrt(1.0 - u1), s2 = std::sqrt(u1);
			const double twoPi = 6.283185307179586;
			const double qx = s1 * std::sin(twoPi * u2), qy = s1 * std::cos(twoPi * u2);
			const double qz = s2 * std::sin(twoPi * u3), qw = s2 * std::cos(twoPi * u3);
			bi.pose.m[0] = 1 - 2 * (qy * qy + qz * qz); bi.pose.m[1] = 2 * (qx * qy - qz * qw); bi.pose.m[2] = 2 * (qx * qz + qy * qw);
			bi.pose.m[3] = 2 * (qx * qy + qz * qw); bi.pose.m[4] = 1 - 2 * (qx * qx + qz * qz); bi.pose.m[5] = 2 * (qy * qz - qx * qw);
			bi.pose.m[6] = 2 * (qx * qz - qy * qw); bi.pose.m[7] = 2 * (qy * qz + qx * qw); bi.pose.m[8] = 1 - 2 * (qx * qx + qy * qy);

			out.push_back(bi);
		}

		// Largest-first, so a timed release drops the big units before the small ones (armour-layer pour).
		// The grid positions were assigned by original (random-scale) index, so after sorting the sizes
		// are spatially mixed — no big/small spatial gradient, just a big→small RELEASE order.
		std::sort(out.begin(), out.end(), [](const BlockInstance& a, const BlockInstance& b) { return a.scale > b.scale; });
		return out;
	}

	SettleResult settle_blocks(const ConvexShape& shape, const std::vector<BlockInstance>& instances,
		const GroundSpec& ground, const SettleParams& params)
	{
		SettleResult result;
		if (shape.empty() || instances.empty()) return result;

		SettleWorld world(shape, instances, ground, params);
		while (world.step(64)) {} // advance in 64-substep chunks until asleep / budget
		return world.result();
	}
}
