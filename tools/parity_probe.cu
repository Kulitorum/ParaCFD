// parity_probe.cu — Phase-0 migration oracle for the graded-structured-grid change.
//
// The graded-grid refactor rewrites ~300 scalar-`h` call sites into metric-aware form
// with a claimed "zero behaviour change" on a uniform grid. CLAUDE.md advertises a
// "1e-5 GPU-vs-CPU parity test per kernel" invariant, but NO test target existed — the
// CPU twins were defined and never run. This restores that invariant AND adds the
// behaviour-preservation gate, so "zero behaviour change" becomes falsifiable.
//
// For every MAC-solver kernel in the Phase-A rewrite scope (both the M1 `mac_*`/
// `poisson_*` family and the M2 masked `ch_*` channel family), on a fixed seeded input
// over a small grid with a nontrivial solid mask, this probe checks TWO things:
//
//   (1) PARITY   — the GPU kernel matches its CPU reference twin at rel. max-norm 1e-5.
//                  Catches launch / indexing / memory bugs introduced by the rewrite.
//   (2) GOLDEN   — the GPU output matches a blessed snapshot at rel. max-norm 1e-5.
//                  Catches formula changes: the "grading collapses to uniform"
//                  regression (design D4 / tasks 0.4 & 5.1). Bless once against the
//                  current uniform solver (`--bless`); every later run must still match.
//
// Usage:
//   parity_probe --bless [--golden <dir>]   # write goldens from the current solver
//   parity_probe          [--golden <dir>]  # verify parity + goldens (exit != 0 on fail)
//
// Default golden dir is "tests/golden" relative to the working directory; the ctest
// entry sets WORKING_DIRECTORY to the repo root so that path resolves.
#include "core/fluid/channel_bc.h"
#include "core/fluid/channel_ops.h"
#include "core/fluid/channel_pressure.h"
#include "core/fluid/grid_metrics.h"
#include "core/fluid/mac_grid.h"
#include "core/fluid/mac_ops.h"
#include "core/geometry/building.h"
#include "core/geometry/voxelize.h" // seal_enclosed_voids (enclosed-void flood-fill)
#include "core/windloads.h"

#include <cuda_runtime.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <vector>

using namespace windcfd::core;

// --- CUDA + device-buffer plumbing -------------------------------------------
#define CU(expr)                                                                          \
	do                                                                                    \
	{                                                                                     \
		cudaError_t _e = (expr);                                                           \
		if (_e != cudaSuccess)                                                             \
		{                                                                                 \
			std::fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(_e), __FILE__, __LINE__); \
			std::exit(2);                                                                  \
		}                                                                                 \
	} while (0)

namespace
{
	std::vector<void*> g_allocs; // tracked device allocations, freed at teardown

	double* dev(const std::vector<double>& h)
	{
		double* d = nullptr;
		CU(cudaMalloc(&d, h.size() * sizeof(double)));
		CU(cudaMemcpy(d, h.data(), h.size() * sizeof(double), cudaMemcpyHostToDevice));
		g_allocs.push_back(d);
		return d;
	}
	double* dev_zero(int n)
	{
		double* d = nullptr;
		CU(cudaMalloc(&d, (size_t)n * sizeof(double)));
		CU(cudaMemset(d, 0, (size_t)n * sizeof(double)));
		g_allocs.push_back(d);
		return d;
	}
	unsigned char* dev_u8(const std::vector<unsigned char>& h)
	{
		unsigned char* d = nullptr;
		CU(cudaMalloc(&d, h.size()));
		CU(cudaMemcpy(d, h.data(), h.size(), cudaMemcpyHostToDevice));
		g_allocs.push_back(d);
		return d;
	}
	std::vector<double> back(const double* d, int n)
	{
		std::vector<double> h(n);
		CU(cudaMemcpy(h.data(), d, (size_t)n * sizeof(double), cudaMemcpyDeviceToHost));
		return h;
	}
	void free_all()
	{
		for (void* p : g_allocs) cudaFree(p);
		g_allocs.clear();
	}

	// Deterministic per-field fill so every run (and every machine) sees identical inputs.
	std::vector<double> filled(int n, uint32_t seed, double lo, double hi)
	{
		std::mt19937 rng(seed);
		std::uniform_real_distribution<double> dist(lo, hi);
		std::vector<double> v(n);
		for (double& x : v) x = dist(rng);
		return v;
	}

	double rel_maxnorm(const std::vector<double>& a, const std::vector<double>& b)
	{
		double num = 0.0, den = 1e-30;
		size_t n = a.size() < b.size() ? a.size() : b.size();
		for (size_t i = 0; i < n; ++i)
		{
			double d = std::fabs(a[i] - b[i]);
			if (d > num) num = d;
			double m = std::fabs(b[i]);
			if (m > den) den = m;
		}
		return num / den;
	}

	// --- Reporter: parity (GPU vs CPU) + golden (GPU vs blessed snapshot) -----
	struct Reporter
	{
		std::string golden_dir = "tests/golden";
		bool bless = false;
		int total = 0, failed = 0;

		std::string path(const std::string& name) const { return golden_dir + "/" + name + ".f64"; }

		bool load_golden(const std::string& name, std::vector<double>& out) const
		{
			std::ifstream f(path(name), std::ios::binary);
			if (!f) return false;
			f.seekg(0, std::ios::end);
			std::streamoff bytes = f.tellg();
			f.seekg(0, std::ios::beg);
			out.resize((size_t)bytes / sizeof(double));
			f.read(reinterpret_cast<char*>(out.data()), bytes);
			return (bool)f;
		}
		void save_golden(const std::string& name, const std::vector<double>& v) const
		{
			std::ofstream f(path(name), std::ios::binary);
			f.write(reinterpret_cast<const char*>(v.data()), (std::streamsize)(v.size() * sizeof(double)));
		}

		// One kernel: gpu output vs cpu twin (parity) and vs golden (behaviour). force_verify=true
		// always compares to the existing golden (never blesses) — used by the uniform-metric-array
		// pass so it is checked against the null-path golden even under --bless.
		void report(const std::string& name, const std::vector<double>& gpu, const std::vector<double>& cpu, bool force_verify = false)
		{
			++total;
			double par = rel_maxnorm(gpu, cpu);
			bool ok = par <= 1e-5;
			std::string gtag;
			if (bless && !force_verify)
			{
				save_golden(name, gpu);
				gtag = "golden=BLESSED";
			}
			else
			{
				std::vector<double> gold;
				if (!load_golden(name, gold))
				{
					gtag = "golden=MISSING";
					ok = false;
				}
				else
				{
					double gd = rel_maxnorm(gpu, gold);
					bool gok = gd <= 1e-5 && gold.size() == gpu.size();
					char buf[64];
					std::snprintf(buf, sizeof buf, "golden=%.2e", gd);
					gtag = buf;
					ok = ok && gok;
				}
			}
			if (!ok) ++failed;
			std::printf("  %-22s parity=%.2e  %-16s  %s\n", name.c_str(), par, gtag.c_str(), ok ? "PASS" : "FAIL");
		}

		// GPU-vs-CPU parity ONLY (no golden) — for a genuinely graded grid, which has no scalar-h
		// reference. Proves the metric-reading kernels are self-consistent on non-uniform spacing.
		void report_parity(const std::string& name, const std::vector<double>& gpu, const std::vector<double>& cpu)
		{
			++total;
			double par = rel_maxnorm(gpu, cpu);
			bool ok = par <= 1e-5;
			if (!ok) ++failed;
			std::printf("  %-22s parity=%.2e  (graded)          %s\n", name.c_str(), par, ok ? "PASS" : "FAIL");
		}

		void check(const std::string& name, bool ok, const std::string& detail)
		{
			++total;
			if (!ok) ++failed;
			std::printf("  %-22s %-24s  %s\n", name.c_str(), detail.c_str(), ok ? "PASS" : "FAIL");
		}
	};
}

int main(int argc, char** argv)
{
	Reporter rep;
	for (int i = 1; i < argc; ++i)
	{
		if (std::strcmp(argv[i], "--bless") == 0) rep.bless = true;
		else if (std::strcmp(argv[i], "--golden") == 0 && i + 1 < argc) rep.golden_dir = argv[++i];
	}

	// --- Test grid + physics knobs (small, but with a real solid + real ghosts) ---
	MacGrid g;
	g.nx = 16; g.ny = 16; g.nz = 16; g.h = 0.05;
	const int NP = g.p_count(), NU = g.u_count(), NV = g.v_count(), NW = g.w_count();
	const int NMAX = NU > NV ? (NU > NW ? NU : NW) : (NV > NW ? NV : NW);

	const double rho = 1.225, nu = 1e-3, dt = 1e-3, Cs = 0.16, omega = 2.0 / 3.0;
	const int band = 2, sweeps = 2, dir_xmax = 1;

	// M1 closed-box BC (cavity ghosts: noslip walls + moving lid), M2 open-channel BC.
	BC mbc; // defaults: xmin/xmax noslip, y freeslip, zmin noslip, zmax movlid
	ChannelBC cbc; // defaults: uniform inlet, freeslip laterals, noslip solid surfaces
	cbc.U_inlet = 1.0; cbc.Uc = 1.0;

	// Solid mask: a 4^3 block near the centre; nearsolid = 1-cell dilation.
	std::vector<unsigned char> solid(NP, 0), nearsolid(NP, 0);
	auto is_solid_box = [](int i, int j, int k) { return i >= 6 && i < 10 && j >= 6 && j < 10 && k >= 6 && k < 10; };
	for (int k = 0; k < g.nz; ++k)
		for (int j = 0; j < g.ny; ++j)
			for (int i = 0; i < g.nx; ++i)
				if (is_solid_box(i, j, k)) solid[g.pidx(i, j, k)] = 1;
	for (int k = 0; k < g.nz; ++k)
		for (int j = 0; j < g.ny; ++j)
			for (int i = 0; i < g.nx; ++i)
			{
				bool nb = false;
				for (int dk = -1; dk <= 1 && !nb; ++dk)
					for (int dj = -1; dj <= 1 && !nb; ++dj)
						for (int di = -1; di <= 1 && !nb; ++di)
						{
							int ii = i + di, jj = j + dj, kk = k + dk;
							if (ii >= 0 && ii < g.nx && jj >= 0 && jj < g.ny && kk >= 0 && kk < g.nz && solid[g.pidx(ii, jj, kk)]) nb = true;
						}
				if (nb) nearsolid[g.pidx(i, j, k)] = 1;
			}

	// Seeded host inputs (independent per field).
	std::vector<double> u = filled(NU, 101, -1.0, 1.0);
	std::vector<double> v = filled(NV, 102, -1.0, 1.0);
	std::vector<double> w = filled(NW, 103, -1.0, 1.0);
	std::vector<double> p = filled(NP, 104, -1.0, 1.0);
	std::vector<double> rhs = filled(NP, 105, -1.0, 1.0);
	std::vector<double> nut = filled(NP, 106, 0.0, 0.05);

	unsigned char* dsolid = dev_u8(solid);
	unsigned char* dnear = dev_u8(nearsolid);

	std::printf("parity_probe: grid %dx%dx%d h=%.3f  %s  golden=\"%s\"\n",
		g.nx, g.ny, g.nz, g.h, rep.bless ? "[BLESS]" : "[verify]", rep.golden_dir.c_str());

	// ============================ M1 mac_* / poisson_* ============================
	std::printf("-- M1 mac_* / poisson_* --\n");
	{
		// mac_advect (u,v,w out)
		double* du = dev(u), *dv = dev(v), *dw = dev(w);
		double* duo = dev_zero(NU), *dvo = dev_zero(NV), *dwo = dev_zero(NW), *dscr = dev_zero(NMAX);
		mac_advect_gpu(du, dv, dw, duo, dvo, dwo, dscr, g, mbc, dt, 1);
		std::vector<double> cu(NU), cv(NV), cw(NW);
		mac_advect_cpu(u, v, w, cu, cv, cw, g, mbc, dt, 1);
		rep.report("mac_advect_u", back(duo, NU), cu);
		rep.report("mac_advect_v", back(dvo, NV), cv);
		rep.report("mac_advect_w", back(dwo, NW), cw);
	}
	{
		// smagorinsky_nut
		double* du = dev(u), *dv = dev(v), *dw = dev(w), *dnut = dev_zero(NP);
		smagorinsky_nut_gpu(du, dv, dw, dnut, g, mbc, Cs);
		std::vector<double> cn(NP);
		smagorinsky_nut_cpu(u, v, w, cn, g, mbc, Cs);
		rep.report("smagorinsky_nut", back(dnut, NP), cn);
	}
	{
		// mac_diffuse (with a real nut field)
		double* du = dev(u), *dv = dev(v), *dw = dev(w), *dnut = dev(nut);
		double* duo = dev_zero(NU), *dvo = dev_zero(NV), *dwo = dev_zero(NW);
		mac_diffuse_gpu(du, dv, dw, duo, dvo, dwo, dnut, g, mbc, dt, nu);
		std::vector<double> cu(NU), cv(NV), cw(NW);
		mac_diffuse_cpu(u, v, w, cu, cv, cw, &nut, g, mbc, dt, nu);
		rep.report("mac_diffuse_u", back(duo, NU), cu);
		rep.report("mac_diffuse_v", back(dvo, NV), cv);
		rep.report("mac_diffuse_w", back(dwo, NW), cw);
	}
	{
		// poisson_rhs
		double* du = dev(u), *dv = dev(v), *dw = dev(w), *dr = dev_zero(NP);
		poisson_rhs_gpu(du, dv, dw, dr, g, rho, dt);
		std::vector<double> cr(NP);
		poisson_rhs_cpu(u, v, w, cr, g, rho, dt);
		rep.report("poisson_rhs", back(dr, NP), cr);
	}
	{
		// poisson_apply
		double* dp = dev(p), *dap = dev_zero(NP);
		poisson_apply_gpu(dp, dap, g);
		std::vector<double> cap(NP);
		poisson_apply_cpu(p, cap, g);
		rep.report("poisson_apply", back(dap, NP), cap);
	}
	{
		// jacobi_smooth (mutates p)
		double* dp = dev(p), *dr = dev(rhs), *dscr = dev_zero(NP);
		jacobi_smooth_gpu(dp, dr, dscr, g, omega, sweeps);
		std::vector<double> cp = p;
		jacobi_smooth_cpu(cp, rhs, g, omega, sweeps);
		rep.report("jacobi_smooth", back(dp, NP), cp);
	}
	{
		// gs_band (mutates p)
		double* dp = dev(p), *dr = dev(rhs);
		gs_band_gpu(dp, dr, g, band, sweeps, true);
		std::vector<double> cp = p;
		gs_band_cpu(cp, rhs, g, band, sweeps, true);
		rep.report("gs_band", back(dp, NP), cp);
	}
	{
		// poisson_residual
		double* dp = dev(p), *dr = dev(rhs), *dres = dev_zero(NP);
		poisson_residual_gpu(dp, dr, dres, g);
		std::vector<double> cres(NP);
		poisson_residual_cpu(p, rhs, cres, g);
		rep.report("poisson_residual", back(dres, NP), cres);
	}
	{
		// restrict (fine -> coarse, factor 2)
		MacGrid gc = g; gc.nx = g.nx / 2; gc.ny = g.ny / 2; gc.nz = g.nz / 2; gc.h = g.h * 2;
		int NPC = gc.p_count();
		std::vector<double> fine = filled(NP, 107, -1.0, 1.0);
		double* dfine = dev(fine), *dcoarse = dev_zero(NPC);
		restrict_gpu(dfine, dcoarse, g, gc);
		std::vector<double> cc(NPC);
		restrict_cpu(fine, cc, g, gc);
		rep.report("restrict", back(dcoarse, NPC), cc);
	}
	{
		// prolong_add (coarse -> fine, accumulate; same fine input to both)
		MacGrid gc = g; gc.nx = g.nx / 2; gc.ny = g.ny / 2; gc.nz = g.nz / 2; gc.h = g.h * 2;
		int NPC = gc.p_count();
		std::vector<double> coarse = filled(NPC, 108, -1.0, 1.0);
		std::vector<double> fine = filled(NP, 109, -1.0, 1.0);
		double* dcoarse = dev(coarse), *dfine = dev(fine);
		prolong_add_gpu(dcoarse, dfine, gc, g);
		std::vector<double> cf = fine;
		prolong_add_cpu(coarse, cf, gc, g);
		rep.report("prolong_add", back(dfine, NP), cf);
	}
	{
		// subtract_gradient (mutates u,v,w)
		double* du = dev(u), *dv = dev(v), *dw = dev(w), *dp = dev(p);
		subtract_gradient_gpu(du, dv, dw, dp, g, rho, dt);
		std::vector<double> cu = u, cv = v, cw = w;
		subtract_gradient_cpu(cu, cv, cw, p, g, rho, dt);
		rep.report("subtract_gradient_u", back(du, NU), cu);
		rep.report("subtract_gradient_v", back(dv, NV), cv);
		rep.report("subtract_gradient_w", back(dw, NW), cw);
	}
	{
		// max_abs_divergence (scalar)
		double* du = dev(u), *dv = dev(v), *dw = dev(w), *dscr = dev_zero(NP);
		double gpu = max_abs_divergence_gpu(du, dv, dw, dscr, g);
		double cpu = max_abs_divergence_cpu(u, v, w, g);
		rep.report("max_abs_divergence", {gpu}, {cpu});
	}

	// ============================ M2 masked ch_* ==================================
	std::printf("-- M2 masked ch_* (channel) --\n");
	{
		// ch_advect (masked MacCormack)
		double* du = dev(u), *dv = dev(v), *dw = dev(w);
		double* duo = dev_zero(NU), *dvo = dev_zero(NV), *dwo = dev_zero(NW), *dscr = dev_zero(NMAX);
		ch_advect_gpu(du, dv, dw, duo, dvo, dwo, dscr, dsolid, dnear, g, cbc, dt);
		std::vector<double> cu(NU), cv(NV), cw(NW);
		ch_advect_cpu(u, v, w, cu, cv, cw, solid, nearsolid, g, cbc, dt);
		rep.report("ch_advect_u", back(duo, NU), cu);
		rep.report("ch_advect_v", back(dvo, NV), cv);
		rep.report("ch_advect_w", back(dwo, NW), cw);
	}
	{
		// ch_smagorinsky
		double* du = dev(u), *dv = dev(v), *dw = dev(w), *dnut = dev_zero(NP);
		ch_smagorinsky_gpu(du, dv, dw, dnut, dsolid, g, cbc, Cs);
		std::vector<double> cn(NP);
		ch_smagorinsky_cpu(u, v, w, cn, solid, g, cbc, Cs);
		rep.report("ch_smagorinsky", back(dnut, NP), cn);
	}
	{
		// ch_diffuse (with a real nut field)
		double* du = dev(u), *dv = dev(v), *dw = dev(w), *dnut = dev(nut);
		double* duo = dev_zero(NU), *dvo = dev_zero(NV), *dwo = dev_zero(NW);
		ch_diffuse_gpu(du, dv, dw, duo, dvo, dwo, dnut, dsolid, g, cbc, dt, nu);
		std::vector<double> cu(NU), cv(NV), cw(NW);
		ch_diffuse_cpu(u, v, w, cu, cv, cw, &nut, solid, g, cbc, dt, nu);
		rep.report("ch_diffuse_u", back(duo, NU), cu);
		rep.report("ch_diffuse_v", back(dvo, NV), cv);
		rep.report("ch_diffuse_w", back(dwo, NW), cw);
	}
	{
		// ch_poisson_rhs
		double* du = dev(u), *dv = dev(v), *dw = dev(w), *dr = dev_zero(NP);
		ch_poisson_rhs_gpu(du, dv, dw, dr, dsolid, g, rho, dt);
		std::vector<double> cr(NP);
		ch_poisson_rhs_cpu(u, v, w, cr, solid, g, rho, dt);
		rep.report("ch_poisson_rhs", back(dr, NP), cr);
	}
	{
		// ch_poisson_apply
		double* dp = dev(p), *dap = dev_zero(NP);
		ch_poisson_apply_gpu(dp, dap, dsolid, g, dir_xmax);
		std::vector<double> cap(NP);
		ch_poisson_apply_cpu(p, cap, solid, g, dir_xmax);
		rep.report("ch_poisson_apply", back(dap, NP), cap);
	}
	{
		// ch_jacobi (mutates p)
		double* dp = dev(p), *dr = dev(rhs), *dscr = dev_zero(NP);
		ch_jacobi_gpu(dp, dr, dscr, dsolid, g, dir_xmax, omega, sweeps);
		std::vector<double> cp = p;
		ch_jacobi_cpu(cp, rhs, solid, g, dir_xmax, omega, sweeps);
		rep.report("ch_jacobi", back(dp, NP), cp);
	}
	{
		// ch_subtract_gradient (mutates u,v,w)
		double* du = dev(u), *dv = dev(v), *dw = dev(w), *dp = dev(p);
		ch_subtract_gradient_gpu(du, dv, dw, dp, dsolid, g, cbc, rho, dt, dir_xmax);
		std::vector<double> cu = u, cv = v, cw = w;
		ch_subtract_gradient_cpu(cu, cv, cw, p, solid, g, cbc, rho, dt, dir_xmax);
		rep.report("ch_subtract_gradient_u", back(du, NU), cu);
		rep.report("ch_subtract_gradient_v", back(dv, NV), cv);
		rep.report("ch_subtract_gradient_w", back(dw, NW), cw);
	}
	{
		// ch_max_div (scalar)
		double* du = dev(u), *dv = dev(v), *dw = dev(w), *dscr = dev_zero(NP);
		double gpu = ch_max_div_gpu(du, dv, dw, dscr, dsolid, g);
		double cpu = ch_max_div_cpu(u, v, w, solid, g);
		rep.report("ch_max_div", {gpu}, {cpu});
	}

	// ==================== Phase-A representation validation ======================
	std::printf("-- world<->index round-trip + uniform-metric-array collapse --\n");
	{
		// (1) locate_frac round-trip on a genuinely GRADED 1-D metric (task 1.3): world → fractional
		// index → world must recover the original coordinate on non-uniform spacing.
		int N = 12; std::vector<double> xf(N + 1); xf[0] = 0.0; double wdt = 0.02;
		for (int i = 0; i < N; ++i) { xf[i + 1] = xf[i] + wdt; wdt *= 1.12; } // geometric growth
		double maxerr = 0.0;
		for (int s = 0; s <= 40; ++s)
		{
			double x = xf[0] + (xf[N] - xf[0]) * (s / 40.0) * 0.999 + 1e-6;
			double t = locate_frac(xf.data(), N + 1, x);
			int lo = (int)std::floor(t); if (lo < 0) lo = 0; if (lo > N - 1) lo = N - 1;
			double xr = xf[lo] + (t - lo) * (xf[lo + 1] - xf[lo]);
			maxerr = std::max(maxerr, std::fabs(xr - x));
		}
		bool ok = maxerr < 1e-9; ++rep.total; if (!ok) ++rep.failed;
		std::printf("  %-22s err=%.2e       %s\n", "locate_roundtrip", maxerr, ok ? "PASS" : "FAIL");
	}
	{
		// (2) Uniform metric ARRAYS (not the null-pointer fallback) must reproduce the golden. Build
		// dx=h / xc / xf on host (for the CPU twins) and device (for the GPU kernels) and re-run the
		// formula-changed operators through the metric-READING path, comparing to the null-path golden.
		std::vector<double> Hdx(g.nx, g.h), Hdy(g.ny, g.h), Hdz(g.nz, g.h);
		std::vector<double> Hxc(g.nx), Hyc(g.ny), Hzc(g.nz), Hxf(g.nx + 1), Hyf(g.ny + 1), Hzf(g.nz + 1);
		for (int i = 0; i < g.nx; ++i) Hxc[i] = (i + 0.5) * g.h;
		for (int i = 0; i <= g.nx; ++i) Hxf[i] = i * g.h;
		for (int j = 0; j < g.ny; ++j) Hyc[j] = (j + 0.5) * g.h;
		for (int j = 0; j <= g.ny; ++j) Hyf[j] = j * g.h;
		for (int k = 0; k < g.nz; ++k) Hzc[k] = (k + 0.5) * g.h;
		for (int k = 0; k <= g.nz; ++k) Hzf[k] = k * g.h;
		MacGrid gh = g, gd = g; // host-array view (CPU twins) + device-array view (GPU kernels)
		gh.dxa = Hdx.data(); gh.dya = Hdy.data(); gh.dza = Hdz.data();
		gh.xca = Hxc.data(); gh.yca = Hyc.data(); gh.zca = Hzc.data();
		gh.xfa = Hxf.data(); gh.yfa = Hyf.data(); gh.zfa = Hzf.data(); gh.hmin = g.h;
		gd.dxa = dev(Hdx); gd.dya = dev(Hdy); gd.dza = dev(Hdz);
		gd.xca = dev(Hxc); gd.yca = dev(Hyc); gd.zca = dev(Hzc);
		gd.xfa = dev(Hxf); gd.yfa = dev(Hyf); gd.zfa = dev(Hzf); gd.hmin = g.h;
		{ double* dp = dev(p), *dap = dev_zero(NP); poisson_apply_gpu(dp, dap, gd);
			std::vector<double> cap(NP); poisson_apply_cpu(p, cap, gh); rep.report("poisson_apply", back(dap, NP), cap, true); }
		{ double* du = dev(u), *dv = dev(v), *dw = dev(w), *dp = dev(p); subtract_gradient_gpu(du, dv, dw, dp, gd, rho, dt);
			std::vector<double> cu = u, cv = v, cw = w; subtract_gradient_cpu(cu, cv, cw, p, gh, rho, dt); rep.report("subtract_gradient_u", back(du, NU), cu, true); }
		{ double* du = dev(u), *dv = dev(v), *dw = dev(w); double* duo = dev_zero(NU), *dvo = dev_zero(NV), *dwo = dev_zero(NW), *dscr = dev_zero(NMAX);
			mac_advect_gpu(du, dv, dw, duo, dvo, dwo, dscr, gd, mbc, dt, 1);
			std::vector<double> cu(NU), cv(NV), cw(NW); mac_advect_cpu(u, v, w, cu, cv, cw, gh, mbc, dt, 1); rep.report("mac_advect_u", back(duo, NU), cu, true); }
		{ double* du = dev(u), *dv = dev(v), *dw = dev(w), *dnut = dev(nut); double* duo = dev_zero(NU), *dvo = dev_zero(NV), *dwo = dev_zero(NW);
			mac_diffuse_gpu(du, dv, dw, duo, dvo, dwo, dnut, gd, mbc, dt, nu);
			std::vector<double> cu(NU), cv(NV), cw(NW); mac_diffuse_cpu(u, v, w, cu, cv, cw, &nut, gh, mbc, dt, nu); rep.report("mac_diffuse_u", back(duo, NU), cu, true); }
		{ double* dp = dev(p), *dap = dev_zero(NP); ch_poisson_apply_gpu(dp, dap, dsolid, gd, dir_xmax);
			std::vector<double> cap(NP); ch_poisson_apply_cpu(p, cap, solid, gh, dir_xmax); rep.report("ch_poisson_apply", back(dap, NP), cap, true); }
		{ double* du = dev(u), *dv = dev(v), *dw = dev(w), *dnut = dev_zero(NP); ch_smagorinsky_gpu(du, dv, dw, dnut, dsolid, gd, cbc, Cs);
			std::vector<double> cn(NP); ch_smagorinsky_cpu(u, v, w, cn, solid, gh, cbc, Cs); rep.report("ch_smagorinsky", back(dnut, NP), cn, true); }
		{ double* du = dev(u), *dv = dev(v), *dw = dev(w), *dnut = dev(nut); double* duo = dev_zero(NU), *dvo = dev_zero(NV), *dwo = dev_zero(NW);
			ch_diffuse_gpu(du, dv, dw, duo, dvo, dwo, dnut, dsolid, gd, cbc, dt, nu);
			std::vector<double> cu(NU), cv(NV), cw(NW); ch_diffuse_cpu(u, v, w, cu, cv, cw, &nut, solid, gh, cbc, dt, nu); rep.report("ch_diffuse_u", back(duo, NU), cu, true); }
		{ double* du = dev(u), *dv = dev(v), *dw = dev(w), *duo = dev_zero(NU), *dvo = dev_zero(NV), *dwo = dev_zero(NW), *dscr = dev_zero(NMAX);
			ch_advect_gpu(du, dv, dw, duo, dvo, dwo, dscr, dsolid, dnear, gd, cbc, dt);
			std::vector<double> cu(NU), cv(NV), cw(NW); ch_advect_cpu(u, v, w, cu, cv, cw, solid, nearsolid, gh, cbc, dt); rep.report("ch_advect_u", back(duo, NU), cu, true); }
		{ double* du = dev(u), *dv = dev(v), *dw = dev(w), *dscr = dev_zero(NP); double gp = ch_max_div_gpu(du, dv, dw, dscr, dsolid, gd);
			double cp = ch_max_div_cpu(u, v, w, solid, gh); rep.report("ch_max_div", {gp}, {cp}, true); }
	}

	// ==================== Phase-B: fine-core generator + graded validation =======
	std::printf("-- graded-grid generator properties (task 3.1) --\n");
	{
		const double L = 1.0, a = 0.4, b = 0.6, hf = 0.025, gr = 1.15;
		std::vector<double> xf = graded_axis_faces(L, a, b, hf, gr);
		int n = (int)xf.size() - 1;
		std::vector<double> dxa(n);
		for (int i = 0; i < n; ++i) dxa[i] = xf[i + 1] - xf[i];
		bool mono = std::fabs(xf[0]) < 1e-12 && std::fabs(xf[n] - L) < 1e-9;
		for (int i = 0; i < n; ++i) if (dxa[i] <= 0.0) mono = false;
		double maxratio = 0.0;
		for (int i = 1; i < n; ++i) { double r = dxa[i] > dxa[i - 1] ? dxa[i] / dxa[i - 1] : dxa[i - 1] / dxa[i]; maxratio = std::max(maxratio, r); }
		int nfine = 0;
		for (int i = 0; i < n; ++i) if (std::fabs(dxa[i] - hf) < 1e-9) ++nfine;
		int n_fine_exp = std::max(1, (int)std::lround((b - a) / hf));
		bool coreok = nfine >= n_fine_exp; // the uniform core is ≥ n_fine cells of EXACTLY h_fine
		char d1[64]; std::snprintf(d1, sizeof d1, "%d cells, ends 0..L", n);
		char d3[64]; std::snprintf(d3, sizeof d3, "max ratio %.4f <= %.2f", maxratio, gr);
		rep.check("gen_monotonic", mono, d1);
		rep.check("gen_core_uniform", coreok, "core cells == h_fine");
		rep.check("gen_growth_bound", maxratio <= gr * 1.0001, d3);
	}
	std::printf("-- graded generator: no sub-h_fine sliver at a domain edge (CFL dt-collapse guard) --\n");
	{
		// A core box that doesn't divide the axis evenly (4.0 m core in a 5 m axis at h=0.03) used to leave a
		// ~0.005 m sliver cell at the edge after centring — collapsing h_min() and the CFL timestep (dt ∝
		// h_min), so the flow crawled though it stayed stable. The fix slides the core flush to that edge.
		// Assert the smallest cell is ≥ h_fine (no sliver) and the edge cell is exactly h_fine.
		const double L = 5.0, a = 0.0, b = 4.0, hf = 0.03, gr = 1.15;
		std::vector<double> xf = graded_axis_faces(L, a, b, hf, gr);
		int n = (int)xf.size() - 1;
		double dmin = 1e9;
		for (int i = 0; i < n; ++i) dmin = std::min(dmin, xf[i + 1] - xf[i]);
		bool no_sliver = dmin >= 0.9 * hf;                        // no sliver (the ~2% fill_gap junction cell is fine; a real sliver is 6x smaller)
		bool edge_fine = std::fabs((xf[1] - xf[0]) - hf) < 1e-9;  // the ground/edge cell is exactly h_fine
		char d[80]; std::snprintf(d, sizeof d, "min dx=%.4f (h_fine=%.3f), edge=%.4f", dmin, hf, xf[1] - xf[0]);
		rep.check("gen_no_edge_sliver", no_sliver && edge_fine, d);
	}

	std::printf("-- GPU-vs-CPU parity on a GENUINELY GRADED grid (task 5.3) --\n");
	{
		FineCoreSpec spec;
		spec.Lx = spec.Ly = spec.Lz = 1.0;
		spec.x0 = spec.y0 = spec.z0 = 0.4; spec.x1 = spec.y1 = spec.z1 = 0.6;
		spec.h_fine = 0.03; spec.growth = 1.15;
		GridMetrics gm = GridMetrics::generate(spec);
		MacGrid gd = gm.device_view(), gh = gm.host_view();
		int gNP = gh.p_count(), gNU = gh.u_count(), gNV = gh.v_count(), gNW = gh.w_count();
		int gNMAX = gNU > gNV ? (gNU > gNW ? gNU : gNW) : (gNV > gNW ? gNV : gNW);
		std::printf("   graded grid %dx%dx%d  h_fine=%.3f  h_min=%.4f  Lx=%.3f\n", gh.nx, gh.ny, gh.nz, gm.h_fine(), gm.h_min(), gh.Lx());

		// solid box near the fine core + its 1-cell dilation
		std::vector<unsigned char> gsolid(gNP, 0), gnear(gNP, 0);
		for (int k = 0; k < gh.nz; ++k) for (int j = 0; j < gh.ny; ++j) for (int i = 0; i < gh.nx; ++i)
		{
			double cx = gh.xc(i), cy = gh.yc(j), cz = gh.zc(k);
			if (cx > 0.45 && cx < 0.55 && cy > 0.45 && cy < 0.55 && cz > 0.45 && cz < 0.55) gsolid[gh.pidx(i, j, k)] = 1;
		}
		for (int k = 0; k < gh.nz; ++k) for (int j = 0; j < gh.ny; ++j) for (int i = 0; i < gh.nx; ++i)
		{
			bool nb = false;
			for (int dk = -1; dk <= 1 && !nb; ++dk) for (int dj = -1; dj <= 1 && !nb; ++dj) for (int di = -1; di <= 1 && !nb; ++di)
			{ int ii = i + di, jj = j + dj, kk = k + dk; if (ii >= 0 && ii < gh.nx && jj >= 0 && jj < gh.ny && kk >= 0 && kk < gh.nz && gsolid[gh.pidx(ii, jj, kk)]) nb = true; }
			if (nb) gnear[gh.pidx(i, j, k)] = 1;
		}
		unsigned char* gds = dev_u8(gsolid);
		unsigned char* gdn = dev_u8(gnear);
		std::vector<double> gu = filled(gNU, 201, -1, 1), gv = filled(gNV, 202, -1, 1), gw = filled(gNW, 203, -1, 1);
		std::vector<double> gp = filled(gNP, 204, -1, 1), grhs = filled(gNP, 205, -1, 1), gnut = filled(gNP, 206, 0, 0.05);

		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *duo = dev_zero(gNU), *dvo = dev_zero(gNV), *dwo = dev_zero(gNW), *ds = dev_zero(gNMAX);
			mac_advect_gpu(du, dv, dw, duo, dvo, dwo, ds, gd, mbc, dt, 1);
			std::vector<double> cu(gNU), cv(gNV), cw(gNW); mac_advect_cpu(gu, gv, gw, cu, cv, cw, gh, mbc, dt, 1); rep.report_parity("mac_advect_u", back(duo, gNU), cu); }
		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *dn = dev_zero(gNP); smagorinsky_nut_gpu(du, dv, dw, dn, gd, mbc, Cs);
			std::vector<double> cn(gNP); smagorinsky_nut_cpu(gu, gv, gw, cn, gh, mbc, Cs); rep.report_parity("smagorinsky_nut", back(dn, gNP), cn); }
		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *dn = dev(gnut), *duo = dev_zero(gNU), *dvo = dev_zero(gNV), *dwo = dev_zero(gNW);
			mac_diffuse_gpu(du, dv, dw, duo, dvo, dwo, dn, gd, mbc, dt, nu);
			std::vector<double> cu(gNU), cv(gNV), cw(gNW); mac_diffuse_cpu(gu, gv, gw, cu, cv, cw, &gnut, gh, mbc, dt, nu); rep.report_parity("mac_diffuse_u", back(duo, gNU), cu); }
		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *dr = dev_zero(gNP); poisson_rhs_gpu(du, dv, dw, dr, gd, rho, dt);
			std::vector<double> cr(gNP); poisson_rhs_cpu(gu, gv, gw, cr, gh, rho, dt); rep.report_parity("poisson_rhs", back(dr, gNP), cr); }
		{ double* dp = dev(gp), *da = dev_zero(gNP); poisson_apply_gpu(dp, da, gd);
			std::vector<double> ca(gNP); poisson_apply_cpu(gp, ca, gh); rep.report_parity("poisson_apply", back(da, gNP), ca); }
		{ double* dp = dev(gp), *dr = dev(grhs), *ds = dev_zero(gNP); jacobi_smooth_gpu(dp, dr, ds, gd, omega, sweeps);
			std::vector<double> cp = gp; jacobi_smooth_cpu(cp, grhs, gh, omega, sweeps); rep.report_parity("jacobi_smooth", back(dp, gNP), cp); }
		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *dp = dev(gp); subtract_gradient_gpu(du, dv, dw, dp, gd, rho, dt);
			std::vector<double> cu = gu, cv = gv, cw = gw; subtract_gradient_cpu(cu, cv, cw, gp, gh, rho, dt); rep.report_parity("subtract_gradient_u", back(du, gNU), cu); }
		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *duo = dev_zero(gNU), *dvo = dev_zero(gNV), *dwo = dev_zero(gNW), *ds = dev_zero(gNMAX);
			ch_advect_gpu(du, dv, dw, duo, dvo, dwo, ds, gds, gdn, gd, cbc, dt);
			std::vector<double> cu(gNU), cv(gNV), cw(gNW); ch_advect_cpu(gu, gv, gw, cu, cv, cw, gsolid, gnear, gh, cbc, dt); rep.report_parity("ch_advect_u", back(duo, gNU), cu); }
		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *dn = dev_zero(gNP); ch_smagorinsky_gpu(du, dv, dw, dn, gds, gd, cbc, Cs);
			std::vector<double> cn(gNP); ch_smagorinsky_cpu(gu, gv, gw, cn, gsolid, gh, cbc, Cs); rep.report_parity("ch_smagorinsky", back(dn, gNP), cn); }
		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *dn = dev(gnut), *duo = dev_zero(gNU), *dvo = dev_zero(gNV), *dwo = dev_zero(gNW);
			ch_diffuse_gpu(du, dv, dw, duo, dvo, dwo, dn, gds, gd, cbc, dt, nu);
			std::vector<double> cu(gNU), cv(gNV), cw(gNW); ch_diffuse_cpu(gu, gv, gw, cu, cv, cw, &gnut, gsolid, gh, cbc, dt, nu); rep.report_parity("ch_diffuse_u", back(duo, gNU), cu); }
		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *dr = dev_zero(gNP); ch_poisson_rhs_gpu(du, dv, dw, dr, gds, gd, rho, dt);
			std::vector<double> cr(gNP); ch_poisson_rhs_cpu(gu, gv, gw, cr, gsolid, gh, rho, dt); rep.report_parity("ch_poisson_rhs", back(dr, gNP), cr); }
		{ double* dp = dev(gp), *da = dev_zero(gNP); ch_poisson_apply_gpu(dp, da, gds, gd, dir_xmax);
			std::vector<double> ca(gNP); ch_poisson_apply_cpu(gp, ca, gsolid, gh, dir_xmax); rep.report_parity("ch_poisson_apply", back(da, gNP), ca); }
		{ double* dp = dev(gp), *dr = dev(grhs), *ds = dev_zero(gNP); ch_jacobi_gpu(dp, dr, ds, gds, gd, dir_xmax, omega, sweeps);
			std::vector<double> cp = gp; ch_jacobi_cpu(cp, grhs, gsolid, gh, dir_xmax, omega, sweeps); rep.report_parity("ch_jacobi", back(dp, gNP), cp); }
		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *dp = dev(gp); ch_subtract_gradient_gpu(du, dv, dw, dp, gds, gd, cbc, rho, dt, dir_xmax);
			std::vector<double> cu = gu, cv = gv, cw = gw; ch_subtract_gradient_cpu(cu, cv, cw, gp, gsolid, gh, cbc, rho, dt, dir_xmax); rep.report_parity("ch_subtract_gradient_u", back(du, gNU), cu); }
		{ double* du = dev(gu), *dv = dev(gv), *dw = dev(gw), *ds = dev_zero(gNP); double a = ch_max_div_gpu(du, dv, dw, ds, gds, gd);
			double c = ch_max_div_cpu(gu, gv, gw, gsolid, gh); rep.report_parity("ch_max_div", {a}, {c}); }
	}

	std::printf("-- MMS: FV Laplacian order-of-accuracy under grading (task 5.2) --\n");
	{
		// Manufactured p = sin(kx)sin(ky)sin(kz); the FV operator Ap ≈ (kx²+ky²+kz²)·p. Measure the
		// L2 error over the interior (away from the box-Neumann boundary) at two h_fine and estimate
		// the observed order. Graceful degradation on stretched cells ⇒ expect ≳1st order.
		const double PI = 3.14159265358979323846, k = 2.0 * PI, k2 = 3.0 * k * k;
		double err[2], hf[2] = {0.04, 0.02};
		for (int r = 0; r < 2; ++r)
		{
			FineCoreSpec spec; spec.Lx = spec.Ly = spec.Lz = 1.0;
			spec.x0 = spec.y0 = spec.z0 = 0.3; spec.x1 = spec.y1 = spec.z1 = 0.7;
			spec.h_fine = hf[r]; spec.growth = 1.12;
			GridMetrics gm = GridMetrics::generate(spec);
			MacGrid gh = gm.host_view();
			int NPg = gh.p_count();
			std::vector<double> pp(NPg), Ap;
			for (int kk = 0; kk < gh.nz; ++kk) for (int jj = 0; jj < gh.ny; ++jj) for (int ii = 0; ii < gh.nx; ++ii)
				pp[gh.pidx(ii, jj, kk)] = std::sin(k * gh.xc(ii)) * std::sin(k * gh.yc(jj)) * std::sin(k * gh.zc(kk));
			poisson_apply_cpu(pp, Ap, gh);
			double num = 0.0, den = 0.0;
			for (int kk = 0; kk < gh.nz; ++kk) for (int jj = 0; jj < gh.ny; ++jj) for (int ii = 0; ii < gh.nx; ++ii)
			{
				double cx = gh.xc(ii), cy = gh.yc(jj), cz = gh.zc(kk);
				if (cx < 0.15 || cx > 0.85 || cy < 0.15 || cy > 0.85 || cz < 0.15 || cz > 0.85) continue; // skip boundary layers
				double exact = k2 * pp[gh.pidx(ii, jj, kk)];
				double e = Ap[gh.pidx(ii, jj, kk)] - exact;
				num += e * e; den += exact * exact;
			}
			err[r] = std::sqrt(num / den);
		}
		double order = std::log(err[0] / err[1]) / std::log(hf[0] / hf[1]);
		char d[96]; std::snprintf(d, sizeof d, "err %.2e->%.2e  order=%.2f", err[0], err[1], order);
		rep.check("mms_laplacian_order", order >= 0.9 && err[1] < err[0], d);
	}

	std::printf("-- graded MGPCG solve convergence w/ decimated coarse metrics (task 2.5) --\n");
	{
		// Even per-axis cell counts ⇒ a real multi-level hierarchy ⇒ the per-level decimation is
		// exercised. Solve A x = b (Dirichlet outlet ⇒ non-singular) on an empty graded channel;
		// convergence proves the graded coarse operators + transfers form a working preconditioner.
		FineCoreSpec spec;
		spec.Lx = spec.Ly = spec.Lz = 2.0;
		spec.x0 = spec.y0 = spec.z0 = 0.75; spec.x1 = spec.y1 = spec.z1 = 1.25;
		spec.h_fine = 0.0625; spec.growth = 1.15;
		GridMetrics gm = GridMetrics::generate(spec);
		MacGrid gd = gm.device_view(), gh = gm.host_view();
		int NPg = gh.p_count();
		std::vector<unsigned char> gsolid(NPg, 0);
		unsigned char* gds = dev_u8(gsolid);
		ChannelMgpcg solver(gd);
		solver.build_masks(gds);
		std::vector<double> b = filled(NPg, 301, -1.0, 1.0);
		double* db = dev(b); double* dxs = dev_zero(NPg);
		SolveResult res = solver.solve(dxs, db, 1e-5, 500, false);
		int levels = (int)solver.levels().size();
		char d[96]; std::snprintf(d, sizeof d, "%dx%dx%d %d lvls %d it relres %.1e", gh.nx, gh.ny, gh.nz, levels, res.iters, res.relres);
		rep.check("graded_mgpcg_solve", res.converged, d);
	}

	std::printf("-- windloads bbox-in-uniform-core guard (task 4.3) --\n");
	{
		FineCoreSpec spec;
		spec.Lx = spec.Ly = spec.Lz = 1.0;
		spec.x0 = spec.y0 = spec.z0 = 0.4; spec.x1 = spec.y1 = spec.z1 = 0.6;
		spec.h_fine = 0.03; spec.growth = 1.15;
		GridMetrics gm = GridMetrics::generate(spec);
		MacGrid gh = gm.host_view();
		int NPg = gh.p_count();
		std::vector<double> pf = filled(NPg, 401, -1.0, 1.0);
		WindLoadParams prm; prm.rho = 1.225; prm.u_ref = 5.0;
		auto make_solid = [&](double x0, double x1) {
			std::vector<unsigned char> s(NPg, 0);
			for (int k = 0; k < gh.nz; ++k) for (int j = 0; j < gh.ny; ++j) for (int i = 0; i < gh.nx; ++i)
			{
				double cx = gh.xc(i), cy = gh.yc(j), cz = gh.zc(k);
				if (cx > x0 && cx < x1 && cy > 0.45 && cy < 0.55 && cz > 0.45 && cz < 0.55) s[gh.pidx(i, j, k)] = 1;
			}
			return s;
		};
		// (a) building entirely inside the uniform fine core → integrates without error.
		bool ok_inside = false;
		try { std::vector<unsigned char> s = make_solid(0.45, 0.55); WindLoads L = compute_wind_loads(pf.data(), s.data(), gh, prm, nullptr); ok_inside = L.solid_cells > 0; }
		catch (...) { ok_inside = false; }
		rep.check("windloads_in_core_ok", ok_inside, "building in fine core integrates");
		// (b) building reaching into the graded transition → guard throws (flagged, not silently wrong).
		bool threw = false;
		try { std::vector<unsigned char> s = make_solid(0.45, 0.80); (void)compute_wind_loads(pf.data(), s.data(), gh, prm, nullptr); }
		catch (const std::exception&) { threw = true; }
		rep.check("windloads_transition_throws", threw, "building in transition rejected");
	}

	std::printf("-- voxelize_building world→index on a graded grid (task 4.1) --\n");
	{
		// A square wall centerline around world (1,1) must voxelize to solid cells AT world (1,1) — the
		// graded world→index mapping. The old floor(x/h_fine) would mis-place it (the core doesn't start
		// at x=0). Also assert every solid cell is in the uniform core (dx=dy=h_fine).
		FineCoreSpec spec; spec.Lx = spec.Ly = 2.0; spec.Lz = 1.0;
		spec.x0 = spec.y0 = 0.8; spec.x1 = spec.y1 = 1.2; spec.z0 = 0.0; spec.z1 = 0.6;
		spec.h_fine = 0.04; spec.growth = 1.2;
		GridMetrics gm = GridMetrics::generate(spec);
		MacGrid gh = gm.host_view();
		Footprint fp;
		Loop2D loop; loop.push_back({0.9, 0.9}); loop.push_back({1.1, 0.9}); loop.push_back({1.1, 1.1}); loop.push_back({0.9, 1.1});
		fp.loops.push_back(loop); fp.bbox_min = {0.9, 0.9}; fp.bbox_max = {1.1, 1.1};
		BuildingParams bp; bp.wall_thickness = 0.06; bp.wall_height = 0.4; bp.corner_radius = 0.0; bp.roof_thickness = 0.0; bp.roof_overhang = 0.0; bp.base_z = 0.0;
		int sc = 0;
		std::vector<unsigned char> mask = voxelize_building(fp, bp, gh, &sc);
		double xmn = 1e9, xmx = -1e9, ymn = 1e9, ymx = -1e9;
		bool allcore = true;
		for (int k = 0; k < gh.nz; ++k) for (int j = 0; j < gh.ny; ++j) for (int i = 0; i < gh.nx; ++i)
			if (mask[gh.pidx(i, j, k)])
			{
				xmn = std::min(xmn, gh.xc(i)); xmx = std::max(xmx, gh.xc(i));
				ymn = std::min(ymn, gh.yc(j)); ymx = std::max(ymx, gh.yc(j));
				if (std::fabs(gh.dx(i) - spec.h_fine) > 1e-9 || std::fabs(gh.dy(j) - spec.h_fine) > 1e-9) allcore = false;
			}
		// Solid bbox must be CENTRED on the footprint (1,1) and span ≈ the footprint width (0.2 m). The old
		// floor(x/h_fine) mapping would place it at a different world location (the core doesn't start at 0).
		bool located = sc > 0 && std::fabs(0.5 * (xmn + xmx) - 1.0) < 0.03 && std::fabs(0.5 * (ymn + ymx) - 1.0) < 0.03
			&& (xmx - xmn) > 0.15 && (xmx - xmn) < 0.25 && (ymx - ymn) > 0.15 && (ymx - ymn) < 0.25;
		char d[96]; std::snprintf(d, sizeof d, "%d cells x[%.2f,%.2f] y[%.2f,%.2f]", sc, xmn, xmx, ymn, ymx);
		rep.check("voxelize_building_located", located && allcore, d);
	}

	std::printf("-- seal_enclosed_voids: fill a sealed interior, leave an open structure fluid --\n");
	{
		// A 4x4 footprint [2..5] with a perimeter wall ring over k=[0..4] (sits on the ground k=0), optionally
		// capped by a full roof slab at k=5. Interior = (3..4)x(3..4)x(0..4) = 20 cells. Sealed (walls+roof+
		// ground) ⇒ the 20 interior cells are enclosed ⇒ filled. Open-top (no roof) ⇒ interior reaches the
		// sky ⇒ nothing filled. This also exercises the ground-exclusion (interior floor touches k==0).
		MacGrid g; g.nx = 8; g.ny = 8; g.nz = 8; g.h = 1.0;
		const int N = g.p_count();
		auto ring_walls = [&](std::vector<unsigned char>& m, bool roof)
		{
			for (int k = 0; k <= 4; ++k)
				for (int j = 2; j <= 5; ++j)
					for (int i = 2; i <= 5; ++i)
						if (i == 2 || i == 5 || j == 2 || j == 5) m[g.pidx(i, j, k)] = 1;
			if (roof)
				for (int j = 2; j <= 5; ++j)
					for (int i = 2; i <= 5; ++i)
						m[g.pidx(i, j, 5)] = 1;
		};
		std::vector<unsigned char> sealed((std::size_t)N, 0); ring_walls(sealed, true);
		const int filled_sealed = seal_enclosed_voids(sealed, g);
		std::vector<unsigned char> open((std::size_t)N, 0); ring_walls(open, false);
		const int filled_open = seal_enclosed_voids(open, g);
		const bool ok = (filled_sealed == 20) && (filled_open == 0);
		char d[96]; std::snprintf(d, sizeof d, "sealed +%d (want 20), open +%d (want 0)", filled_sealed, filled_open);
		rep.check("seal_enclosed_voids", ok, d);
	}

	std::printf("-- graded scene round-trip: regenerate from the persisted recipe (tasks 4.5/5.4) --\n");
	{
		// A .scn persists the fine-core SPEC and regenerates the metrics on load (recipe-as-source-of-
		// truth). The round-trip invariant is that GridMetrics::generate is deterministic: the same spec
		// yields an identical grid — same widths, centres, and hence world↔index mapping.
		FineCoreSpec spec;
		spec.Lx = 4.0; spec.Ly = 3.0; spec.Lz = 2.0;
		spec.x0 = 1.5; spec.x1 = 2.5; spec.y0 = 1.0; spec.y1 = 2.0; spec.z0 = 0.5; spec.z1 = 1.5;
		spec.h_fine = 0.05; spec.growth = 1.18;
		GridMetrics a = GridMetrics::generate(spec);
		GridMetrics b = GridMetrics::generate(spec); // "reload": regenerate from the same persisted spec
		MacGrid ga = a.host_view(), gb = b.host_view();
		bool same = ga.nx == gb.nx && ga.ny == gb.ny && ga.nz == gb.nz && std::fabs(ga.h - gb.h) < 1e-15 && std::fabs(a.h_min() - b.h_min()) < 1e-15;
		for (int i = 0; i < ga.nx && same; ++i)
			if (std::fabs(ga.dx(i) - gb.dx(i)) > 1e-15 || std::fabs(ga.xc(i) - gb.xc(i)) > 1e-15 || std::fabs(ga.xf(i) - gb.xf(i)) > 1e-15) same = false;
		for (int k = 0; k < ga.nz && same; ++k)
			if (std::fabs(ga.dz(k) - gb.dz(k)) > 1e-15 || std::fabs(ga.zf(k) - gb.zf(k)) > 1e-15) same = false;
		char d[64]; std::snprintf(d, sizeof d, "%dx%dx%d regenerated identical", ga.nx, ga.ny, ga.nz);
		rep.check("scene_regen_deterministic", same, d);
	}

	std::printf("-- h_fine sizing rule + resolution floor (task 6.1) --\n");
	{
		bool r10 = std::fabs(recommended_h_fine(0.3) - 0.03) < 1e-12;  // r/10
		bool floor_hi = below_resolution_floor(0.08, 0.3);             // 0.08 ≥ r/4=0.075 ⇒ under-resolved
		bool floor_ok = !below_resolution_floor(0.03, 0.3);            // 0.03 < 0.075 ⇒ comparison-grade
		bool sharp = !below_resolution_floor(0.08, 0.0);               // sharp (r=0) never flagged
		rep.check("hfine_sizing_rule", r10 && floor_hi && floor_ok && sharp, "h_fine≈r/10, floor at r/4");
	}

	free_all();
	std::printf("\nparity_probe: %d/%d checks passed (%d failed)\n", rep.total - rep.failed, rep.total, rep.failed);
	return rep.failed == 0 ? 0 : 1;
}
