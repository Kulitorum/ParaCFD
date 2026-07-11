// test_settle.cpp — unit test / PLAN G3.1 gate for the Jolt rigid-body drop/settle prep phase
// (scour::core::settle_blocks, in the scour_physics lib). OCC-free: it builds a PROCEDURAL box
// ConvexShape, scatters a handful of blocks above a rigid bed, drops them, and asserts they settle
// to sleep, come to rest ON the bed without interpenetrating, and reproduce bit-for-bit for a fixed
// seed. Links only scour_physics (→ Jolt) + gtest; no OpenCascade, no CUDA.
#include "core/physics/settle.h"

#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

namespace
{
	// One convex box of the given half-extents, centred at its local origin. The settle solver builds
	// a convex hull from the corner positions and uses them for the resting min-z, so only the 8
	// corner vertices are needed (no indices/normals).
	scour::core::ConvexShape make_box(float hx, float hy, float hz)
	{
		scour::core::TriMesh m;
		for (int sx = -1; sx <= 1; sx += 2)
			for (int sy = -1; sy <= 1; sy += 2)
				for (int sz = -1; sz <= 1; sz += 2)
				{
					m.positions.push_back(sx * hx);
					m.positions.push_back(sy * hy);
					m.positions.push_back(sz * hz);
				}
		m.bbox_min = { { -hx, -hy, -hz } };
		m.bbox_max = { { hx, hy, hz } };
		scour::core::ConvexShape s;
		s.pieces.push_back(std::move(m));
		return s;
	}

	scour::core::SettleResult run_settle(unsigned seed)
	{
		scour::core::ConvexShape box = make_box(0.2f, 0.2f, 0.2f); // 0.4 m cube
		std::vector<scour::core::BlockInstance> insts = scour::core::scatter_blocks(
			/*count*/ 9, /*cx*/ 0.0, /*cy*/ 0.0, /*footprint*/ 2.0,
			/*drop_z_lo*/ 0.5, /*drop_z_hi*/ 1.2, /*scale_min*/ 1.0, /*scale_max*/ 1.0,
			box.pieces.front().bbox_size(), seed);

		scour::core::GroundSpec ground;
		ground.bed_z = 0.0;
		scour::core::SettleParams params; // nominal defaults
		return scour::core::settle_blocks(box, insts, ground, params);
	}
}

TEST(Settle, BlocksDropSettleAndRestOnBed)
{
	scour::core::SettleResult r = run_settle(42u);

	ASSERT_EQ(r.placements.size(), 9u);
	EXPECT_TRUE(r.all_asleep) << "did not settle within the step budget (" << r.steps << " steps)";
	EXPECT_LT(r.steps, 8000);

	// The pile rests ON the bed (z = 0): the lowest vertex is at the bed, not through it, not floating.
	EXPECT_GT(r.min_z, -0.05) << "a block fell through the bed (min_z = " << r.min_z << ")";
	EXPECT_LT(r.min_z, 0.05) << "no block reached the bed (min_z = " << r.min_z << ")";

	// No significant interpenetration at rest.
	EXPECT_LT(r.max_penetration, 0.05) << "blocks overlap at rest (" << r.max_penetration << " m)";

	// Every block centre ends up above the bed (a 0.4 m cube's centre is ≥ ~0.2 m up).
	for (const scour::core::ModelPlacement& p : r.placements)
		EXPECT_GT(p.tz, -0.05);
}

TEST(Settle, DeterministicForFixedSeed)
{
	scour::core::SettleResult a = run_settle(1234u);
	scour::core::SettleResult b = run_settle(1234u);

	ASSERT_EQ(a.placements.size(), b.placements.size());
	EXPECT_EQ(a.steps, b.steps);
	for (std::size_t i = 0; i < a.placements.size(); ++i)
	{
		EXPECT_NEAR(a.placements[i].tx, b.placements[i].tx, 1e-6);
		EXPECT_NEAR(a.placements[i].ty, b.placements[i].ty, 1e-6);
		EXPECT_NEAR(a.placements[i].tz, b.placements[i].tz, 1e-6);
		for (int k = 0; k < 9; ++k)
			EXPECT_NEAR(a.placements[i].m[k], b.placements[i].m[k], 1e-6);
	}
}

// PLAN G3.3: the persistent world re-settles onto an eroded bed — a block over a scour hole sinks
// into it, a block on the still-flat bed stays put.
TEST(Settle, BlocksSinkIntoScourHoleOnResettle)
{
	scour::core::ConvexShape box = make_box(0.2f, 0.2f, 0.2f); // 0.4 m cube

	// Four blocks in a row along x at y = 0, dropped onto a flat bed at z = 0.
	std::vector<scour::core::BlockInstance> insts;
	for (int i = 0; i < 4; ++i)
	{
		scour::core::BlockInstance bi;
		bi.scale = 1.0;
		bi.pose.tx = -0.9 + i * 0.6; // x = -0.9, -0.3, 0.3, 0.9
		bi.pose.ty = 0.0;
		bi.pose.tz = 0.6 + 0.05 * i; // small stagger so they don't spawn interpenetrating
		insts.push_back(bi);
	}

	scour::core::GroundSpec ground; ground.bed_z = 0.0;
	scour::core::SettleParams params;
	scour::core::SettleWorld world(box, insts, ground, params);
	while (world.step(64)) {}
	const std::vector<scour::core::ModelPlacement> before = world.placements();

	// Erode a deep scour hole (z = -0.6) under the middle of the row (x ∈ [-0.4, 1.0]); flat elsewhere.
	scour::core::BedField bed;
	bed.nx = 41; bed.ny = 41; bed.x0 = -2.0; bed.y0 = -2.0; bed.h = 0.1;
	bed.z.assign((size_t)bed.nx * bed.ny, 0.0f);
	for (int j = 0; j < bed.ny; ++j)
		for (int i = 0; i < bed.nx; ++i)
		{
			const double x = bed.x0 + i * bed.h, y = bed.y0 + j * bed.h;
			if (x > -0.4 && x < 1.0 && y > -0.6 && y < 0.6) bed.z[(size_t)j * bed.nx + i] = -0.6f;
		}

	const scour::core::SettleResult after = world.resettleOnBed(bed);
	ASSERT_EQ(after.placements.size(), before.size());
	EXPECT_TRUE(after.all_asleep) << "re-settle did not converge (" << after.steps << " steps)";

	// Block 2 (x ≈ 0.3) sits over the hole → it drops well below its previous rest height.
	EXPECT_LT(after.placements[2].tz, before[2].tz - 0.2)
		<< "middle block should sink into the scour hole (before=" << before[2].tz << " after=" << after.placements[2].tz << ")";
	// Block 0 (x ≈ -0.9) is on the still-flat bed → it stays at ~the same height.
	EXPECT_NEAR(after.placements[0].tz, before[0].tz, 0.15)
		<< "edge block on the flat bed should stay put (before=" << before[0].tz << " after=" << after.placements[0].tz << ")";
}

// scatter_blocks must return the units sorted LARGEST-FIRST so a timed release drops big before small.
TEST(Settle, ScatterIsSortedLargestFirst)
{
	scour::core::ConvexShape box = make_box(0.1f, 0.1f, 0.1f);
	std::vector<scour::core::BlockInstance> insts = scour::core::scatter_blocks(
		12, 0.0, 0.0, 3.0, 1.0, 2.0, 0.4, 1.5, box.pieces.front().bbox_size(), 99u);
	ASSERT_EQ(insts.size(), 12u);
	for (std::size_t i = 1; i < insts.size(); ++i)
		EXPECT_GE(insts[i - 1].scale, insts[i].scale) << "not largest-first at index " << i;
}

// addBlock grows the LIVE world (the timed-release path): start empty, drop 3 units, all settle.
TEST(Settle, AddBlockGrowsAndSettles)
{
	scour::core::ConvexShape box = make_box(0.2f, 0.2f, 0.2f);
	scour::core::GroundSpec ground; ground.bed_z = 0.0;
	scour::core::SettleParams params;
	scour::core::SettleWorld world(box, {}, ground, params); // start with no blocks

	for (int i = 0; i < 3; ++i)
	{
		scour::core::BlockInstance bi;
		bi.scale = 1.0;
		bi.pose.tx = -0.6 + i * 0.6;
		bi.pose.tz = 0.6;
		EXPECT_TRUE(world.addBlock(bi));
	}
	while (world.step(64)) {}

	const scour::core::SettleResult res = world.result();
	ASSERT_EQ(res.placements.size(), 3u);
	EXPECT_TRUE(res.all_asleep);
	EXPECT_GT(res.min_z, -0.05);
	EXPECT_LT(res.min_z, 0.10);
}

// Regression for the 800-unit pour: the GUI builds an EMPTY world and adds every unit via addBlock, so
// the Jolt pair/contact buffers must be sized from params.reserve_bodies — NOT from the (zero) initial
// instance count. A dense pile of hundreds of small blocks generates far more than the old 1024-buffer
// floor of broad-phase pairs; if the buffers overflow Jolt silently drops contacts and the blocks
// interpenetrate + sink through the floor. With the fix they settle cleanly.
TEST(Settle, DensePileDoesNotInterpenetrate)
{
	scour::core::ConvexShape box = make_box(0.05f, 0.05f, 0.05f); // 0.1 m units (small ⇒ dense pile)

	const int N = 300;
	scour::core::GroundSpec ground; ground.bed_z = 0.0;
	scour::core::SettleParams params;
	params.reserve_bodies = N; // size the buffers for the full pour

	scour::core::SettleWorld world(box, {}, ground, params);

	// Scatter N over a footprint that keeps spawns SEPARATED (like the GUI), staggered in z; they fall and
	// pack into a dense touching layer — far more than the old 1024 broad-phase-pair floor, the stressor.
	std::vector<scour::core::BlockInstance> insts = scour::core::scatter_blocks(
		N, 0.0, 0.0, /*footprint*/ 1.9, /*z_lo*/ 0.4, /*z_hi*/ 2.2, 1.0, 1.0,
		box.pieces.front().bbox_size(), 314u);
	for (const scour::core::BlockInstance& b : insts) world.addBlock(b);

	// Track the WORST overlap seen across the WHOLE drop (transient interpenetration is what the eye sees).
	double worst_pen = 0.0;
	while (world.step(64)) worst_pen = std::max(worst_pen, world.result().max_penetration);

	const scour::core::SettleResult r = world.result();
	ASSERT_EQ(r.placements.size(), (std::size_t)N);
	EXPECT_TRUE(r.all_asleep) << "dense pile did not settle (" << r.steps << " steps)";
	// No block tunnels the floor, and none overlaps a neighbour by more than a small fraction of its size
	// (0.1 m unit ⇒ well under the block extent — the overflow bug produced multi-cm sinking/overlap).
	EXPECT_GT(r.min_z, -0.03) << "a block sank through the floor (min_z = " << r.min_z << ")";
	EXPECT_LT(worst_pen, 0.03) << "blocks interpenetrated during the drop (" << worst_pen << " m)";
}

TEST(Settle, EmptyInputsReturnEmpty)
{
	scour::core::ConvexShape box = make_box(0.1f, 0.1f, 0.1f);
	scour::core::GroundSpec ground;
	scour::core::SettleParams params;

	// No instances ⇒ no placements.
	scour::core::SettleResult r = scour::core::settle_blocks(box, {}, ground, params);
	EXPECT_TRUE(r.placements.empty());

	// No shape ⇒ no placements either.
	std::vector<scour::core::BlockInstance> one(1);
	scour::core::SettleResult r2 = scour::core::settle_blocks(scour::core::ConvexShape{}, one, ground, params);
	EXPECT_TRUE(r2.placements.empty());
}
