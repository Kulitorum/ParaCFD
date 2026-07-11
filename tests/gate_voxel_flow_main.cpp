// gate_voxel_flow_main.cpp — headless integration check for the voxelized-obstacle path:
// build a procedural box TriMesh, voxelize it into a small open channel, step the M2 core a
// few hundred steps, and assert the live flow actually treats the mask as a solid barrier:
//   (a) mean cell-centred |u| INSIDE the solid cells ≈ 0 (< 1e-3·U) — no penetration;
//   (b) the peak fluid speed exceeds the inlet U — the flow diverts/accelerates around it.
// OCC-free + Qt-free (links only libscour): it uses the same voxelize_mesh + ChannelFluidCore
// path the GUI uses, so it validates injection end-to-end without a display. Prints PASS/FAIL.
#include "core/fluid/channel_core.h"
#include "core/geometry/model_placement.h"
#include "core/geometry/voxelize.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace scour::core;

namespace
{
	// Axis-aligned box [0,Lx]x[0,Ly]x[0,Lz] (metres) as 12 triangles.
	TriMesh make_box(float sx, float sy, float sz)
	{
		TriMesh m;
		const float v[8][3] = {
			{ 0, 0, 0 }, { sx, 0, 0 }, { sx, sy, 0 }, { 0, sy, 0 },
			{ 0, 0, sz }, { sx, 0, sz }, { sx, sy, sz }, { 0, sy, sz }
		};
		for (auto& p : v) { m.positions.push_back(p[0]); m.positions.push_back(p[1]); m.positions.push_back(p[2]); }
		const unsigned q[6][4] = {
			{ 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 2, 3, 7, 6 }, { 1, 2, 6, 5 }, { 0, 4, 7, 3 }
		};
		for (auto& f : q)
		{
			m.indices.push_back(f[0]); m.indices.push_back(f[1]); m.indices.push_back(f[2]);
			m.indices.push_back(f[0]); m.indices.push_back(f[2]); m.indices.push_back(f[3]);
		}
		m.bbox_min = { { 0, 0, 0 } };
		m.bbox_max = { { sx, sy, sz } };
		return m;
	}
}

int main()
{
	// --- Small open channel with a solid box obstacle ------------------------
	const double U = 0.3;
	MacGrid g; g.h = 0.02; g.nx = 50; g.ny = 25; g.nz = 15; // 1.0 x 0.5 x 0.3 m
	const double Lx = g.nx * g.h, Ly = g.ny * g.h;

	TriMesh box = make_box(0.12f, 0.12f, 0.12f); // 6 cells per side → interior solid cells exist
	ModelPlacement place = place_model_on_bed(box, Lx, Ly);

	double mesh_vol = 0.0, voxel_vol = 0.0;
	int thin = 0;
	std::vector<float> frac;
	std::vector<unsigned char> solid = voxelize_mesh(box, g, place, &mesh_vol, &voxel_vol, &frac, &thin);
	long long nsolid = 0;
	for (unsigned char c : solid) nsolid += (c != 0);
	std::printf("gate_voxel_flow: box voxelized -> %lld solid cells (mesh vol %.4g vs voxel vol %.4g m^3, thin=%d)\n",
		nsolid, mesh_vol, voxel_vol, thin);
	if (nsolid <= 0) { std::printf("gate_voxel_flow: FAIL (no solid cells)\n"); return 1; }

	ChannelBC bc;
	bc.inlet_mode = INLET_UNIFORM; bc.U_inlet = U; bc.Uc = U;
	bc.solid_mode = SOLID_FREESLIP; // shipping default (RESEARCH §3)
	bc.ymin = bc.ymax = WALL_FREESLIP;
	bc.zmin = bc.zmax = WALL_FREESLIP;

	ChannelParams pr;
	pr.rho = 1000.0; pr.nu = 5.0e-4; pr.Cs = 0.0; pr.cfl = 1.0; pr.safety = 0.9;
	pr.proj_tol = 1e-3; pr.proj_max_iter = 40; pr.fixed_dt = 0.0; pr.advect_band = 1;

	ChannelFluidCore core(g, bc, pr, solid);
	core.init_uniform(U, 0.03 * U);

	const int nsteps = 400;
	for (int s = 0; s < nsteps; ++s) core.step();

	// --- Sample the settled field (host snapshot; no direct CUDA API here) ---
	FluidSnapshot snap = core.snapshot();
	const std::vector<double>& u = snap.u;
	const std::vector<double>& v = snap.v;
	const std::vector<double>& w = snap.w;

	double sum_solid = 0.0, max_fluid = 0.0;
	long long ns = 0;
	for (int k = 0; k < g.nz; ++k)
		for (int j = 0; j < g.ny; ++j)
			for (int i = 0; i < g.nx; ++i)
			{
				double uc = 0.5 * (u[g.uidx(i, j, k)] + u[g.uidx(i + 1, j, k)]);
				double vc = 0.5 * (v[g.vidx(i, j, k)] + v[g.vidx(i, j + 1, k)]);
				double wc = 0.5 * (w[g.widx(i, j, k)] + w[g.widx(i, j, k + 1)]);
				double sp = std::sqrt(uc * uc + vc * vc + wc * wc);
				if (solid[g.pidx(i, j, k)]) { ++ns; sum_solid += sp; }
				else if (sp > max_fluid) max_fluid = sp;
			}
	double mean_solid = ns > 0 ? sum_solid / (double)ns : 0.0;

	const double pen_tol = 1e-3 * U;
	bool pass_pen = mean_solid < pen_tol;
	bool pass_wake = max_fluid > U; // diversion accelerates the flow past the inlet speed
	std::printf("gate_voxel_flow: after %d steps: mean|u| in solid = %.3e m/s (< %.3e) [%s]\n",
		nsteps, mean_solid, pen_tol, pass_pen ? "PASS" : "FAIL");
	std::printf("gate_voxel_flow: max fluid speed = %.4f m/s (> U = %.4f) [%s]\n",
		max_fluid, U, pass_wake ? "PASS" : "FAIL");

	bool pass = pass_pen && pass_wake;
	std::printf("gate_voxel_flow: RESULT %s\n", pass ? "PASS" : "FAIL");
	return pass ? 0 : 1;
}
