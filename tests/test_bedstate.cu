// test_bedstate.cu — M5 bed state + morphodynamics: GPU-vs-CPU kernel parity (rel. max-norm 1e-5,
// CLAUDE.md) for f_pack reconstruction, box-filter, interface geometry, per-column geometry,
// bedload flux + upwind divergence, the Exner column update, and the avalanche sweep; plus
// sed_physics.h M5 value checks against RESEARCH §5/§6/§7 (slope θ_cr, Wong-Parker, Engelund-
// Fredsøe, Winterwerp) and avalanche mass conservation / repose relaxation.
#include "core/sediment/bedstate.h"
#include "core/sediment/avalanche.h"
#include "core/sediment/sed_physics.h"
#include "gpu_testutil.h"

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

using namespace scour::core;

namespace
{
	const double PI = 3.14159265358979323846;
	const double PHI32 = 32.0 * PI / 180.0;
	MacGrid grid() { MacGrid g; g.nx = 6; g.ny = 7; g.nz = 12; g.h = 0.05; return g; }
	std::vector<double> uni(int n, unsigned seed, double lo, double hi)
	{
		std::mt19937 rng(seed); std::uniform_real_distribution<double> d(lo, hi);
		std::vector<double> v(n); for (auto& x : v) x = d(rng); return v;
	}
	// physically-monotone f_pack from a random per-column grain thickness G.
	std::vector<double> fpack_from_randomG(MacGrid g, unsigned seed, double cpack)
	{
		auto G = uni(g.nx * g.ny, seed, 0.0, g.nz * g.h * cpack * 0.8);
		std::vector<double> f; fpack_from_G_cpu(G, f, g, cpack); return f;
	}
}

// ===================== sed_physics.h M5 value checks =====================================
TEST(SedPhysicsM5, SlopeThetaCrSanity)
{
	// RESEARCH §5 / research/01 §6: φ=32°, β=20° ⇒ up 1.49, down 0.39, transverse 0.76 (×θ_cr,flat).
	double b = 20.0 * PI / 180.0;
	EXPECT_NEAR(theta_cr_slope(1.0, b, 0.0, PHI32), 1.487, 0.01);       // upslope
	EXPECT_NEAR(theta_cr_slope(1.0, b, PI, PHI32), 0.392, 0.01);        // downslope
	EXPECT_NEAR(theta_cr_slope(1.0, b, PI / 2, PHI32), 0.764, 0.01);    // transverse
	// flat bed ⇒ factor 1 exactly.
	EXPECT_NEAR(theta_cr_slope(0.0492, 0.0, 0.0, PHI32), 0.0492, 1e-12);
	// overhang clamp ≥ 0.1·flat (β past φ, downslope).
	EXPECT_GE(theta_cr_slope(1.0, 50.0 * PI / 180.0, PI, PHI32), 0.1 - 1e-12);
}

TEST(SedPhysicsM5, WongParkerBedload)
{
	// research/03 §5 cross-check: θ=1.2050 ⇒ Φ=3.97·(1.1555)^1.5=4.93.
	EXPECT_NEAR(bedload_phi_wp(1.2050), 4.93, 0.02);
	EXPECT_EQ(bedload_phi_wp(0.0495), 0.0);       // at reference ⇒ 0
	EXPECT_EQ(bedload_phi_wp(0.03), 0.0);         // below ⇒ 0
	EXPECT_GT(bedload_phi_wp(0.2), 0.0);
}

TEST(SedPhysicsM5, EngelundFredsoeBedload)
{
	// gated below θ_cr; positive and monotone above; near-threshold matches the 18.74 closed form
	// to <1% (research/03 §Verify): at θ=0.10, θ_cr=0.05 ⇒ 18.74·(θ−θcr)(√θ−0.7√θcr).
	EXPECT_EQ(bedload_phi_ef(0.05, 0.05), 0.0);
	EXPECT_EQ(bedload_phi_ef(0.04, 0.05), 0.0);
	double phi = bedload_phi_ef(0.10, 0.05);
	double closed = 18.74 * (0.10 - 0.05) * (std::sqrt(0.10) - 0.7 * std::sqrt(0.05));
	EXPECT_NEAR(phi, closed, 0.01 * closed);
	EXPECT_GT(bedload_phi_ef(0.3, 0.05), bedload_phi_ef(0.15, 0.05));
}

TEST(SedPhysicsM5, WinterwerpErosionGate)
{
	double rho = 1025, rho_s = 2650, nu = 1.05e-6, d = 0.2e-3;
	double tcr = theta_cr_sw(grain_Dstar(d, rho, rho_s, nu));
	EXPECT_EQ(winterwerp_ulift(0.9 * tcr, tcr, d, rho, rho_s, nu), 0.0); // θ<θcr ⇒ no lift
	EXPECT_EQ(winterwerp_ulift(tcr, tcr, d, rho, rho_s, nu), 0.0);
	double u1 = winterwerp_ulift(2.0 * tcr, tcr, d, rho, rho_s, nu);
	double u2 = winterwerp_ulift(4.0 * tcr, tcr, d, rho, rho_s, nu);
	EXPECT_GT(u1, 0.0); EXPECT_GT(u2, u1);
}

TEST(SedPhysicsM5, NormalizedFillStates)
{
	EXPECT_NEAR(normalized_fill(0.32, 0.0), 0.5, 1e-12);
	EXPECT_NEAR(normalized_fill(SED_CPACK, 0.0), 1.0, 1e-12); // F≥1 ⇒ PACKED
	EXPECT_NEAR(normalized_fill(0.0, 0.0), 0.0, 1e-12);       // FLUID
}

// ===================== kernel GPU-vs-CPU parity =========================================
TEST(BedState, FpackFromGParity)
{
	MacGrid g = grid();
	auto hG = uni(g.nx * g.ny, 1, 0.0, g.nz * g.h * SED_CPACK);
	std::vector<double> ref; fpack_from_G_cpu(hG, ref, g, SED_CPACK);
	DevVec<double> G(hG), f((size_t)g.p_count());
	fpack_from_G_gpu(G.p, f.p, g, SED_CPACK); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(f.download(), ref), 1e-5);
}

TEST(BedState, BoxFilterParity)
{
	MacGrid g = grid();
	auto f = fpack_from_randomG(g, 2, SED_CPACK);
	std::vector<double> ref; bed_boxfilter_cpu(f, nullptr, ref, g, SED_CPACK);
	DevVec<double> fp(f), Ft((size_t)g.p_count());
	bed_boxfilter_gpu(fp.p, nullptr, Ft.p, g, SED_CPACK); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(Ft.download(), ref), 1e-5);
}

TEST(BedState, InterfaceGeomParity)
{
	MacGrid g = grid();
	auto f = fpack_from_randomG(g, 3, SED_CPACK);
	std::vector<double> Ft; bed_boxfilter_cpu(f, nullptr, Ft, g, SED_CPACK);
	std::vector<double> nbx, nby, nbz, Ab; bed_interface_geom_cpu(Ft, nbx, nby, nbz, Ab, g);
	DevVec<double> dFt(Ft), dnbx((size_t)g.p_count()), dnby((size_t)g.p_count()), dnbz((size_t)g.p_count()), dAb((size_t)g.p_count());
	bed_interface_geom_gpu(dFt.p, dnbx.p, dnby.p, dnbz.p, dAb.p, g); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(dnbx.download(), nbx), 1e-5);
	EXPECT_LT(rel_maxnorm(dnby.download(), nby), 1e-5);
	EXPECT_LT(rel_maxnorm(dnbz.download(), nbz), 1e-5);
	EXPECT_LT(rel_maxnorm(dAb.download(), Ab), 1e-5);
}

TEST(BedState, InterfaceGeomFlatAndRamp)
{
	// Flat bed ⇒ interface normal is +z (β=0). Linear ramp z_b=s·x ⇒ β≈atan(s) at the interface.
	MacGrid g; g.nx = 24; g.ny = 6; g.nz = 24; g.h = 0.05;
	double cpack = SED_CPACK;
	// flat
	std::vector<double> Gflat(g.nx * g.ny, 8.0 * g.h * cpack), fflat, Ftf, nbx, nby, nbz, Ab;
	fpack_from_G_cpu(Gflat, fflat, g, cpack);
	bed_boxfilter_cpu(fflat, nullptr, Ftf, g, cpack);
	bed_interface_geom_cpu(Ftf, nbx, nby, nbz, Ab, g);
	// interface cell in a flat column is where 0<F<1; here z_b=8h exactly ⇒ cell 8 boundary; use ramp
	// for the β value test and just assert the flat field has |n_bx|,|n_by| ~ 0 on filled interior.
	int idx = g.pidx(g.nx / 2, g.ny / 2, 7);
	EXPECT_NEAR(nbx[idx], 0.0, 1e-9); EXPECT_NEAR(nby[idx], 0.0, 1e-9);
	// ramp
	double s = std::tan(20.0 * PI / 180.0);
	std::vector<double> Gr(g.nx * g.ny), fr, Ftr, rnbx, rnby, rnbz, rAb, beta, upx, upy;
	for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
	{
		double zb = 6.0 * g.h + s * (i + 0.5) * g.h; // gentle ramp centred in z
		Gr[j * g.nx + i] = zb * cpack;
	}
	fpack_from_G_cpu(Gr, fr, g, cpack);
	bed_boxfilter_cpu(fr, nullptr, Ftr, g, cpack);
	bed_interface_geom_cpu(Ftr, rnbx, rnby, rnbz, rAb, g);
	bed_column_geom_cpu(fr, rnbx, rnby, rnbz, beta, upx, upy, g, cpack);
	// interior column β within a few degrees of atan(s)=20°; up-slope points in +x (upx≈+1).
	int col = (g.ny / 2) * g.nx + g.nx / 2;
	EXPECT_NEAR(beta[col] / PI * 180.0, 20.0, 5.0);
	EXPECT_GT(upx[col], 0.7);
}

TEST(BedState, ColumnGeomParity)
{
	MacGrid g = grid();
	auto f = fpack_from_randomG(g, 4, SED_CPACK);
	std::vector<double> Ft, nbx, nby, nbz, Ab;
	bed_boxfilter_cpu(f, nullptr, Ft, g, SED_CPACK);
	bed_interface_geom_cpu(Ft, nbx, nby, nbz, Ab, g);
	std::vector<double> beta, upx, upy; bed_column_geom_cpu(f, nbx, nby, nbz, beta, upx, upy, g, SED_CPACK);
	DevVec<double> df(f), dnbx(nbx), dnby(nby), dnbz(nbz);
	DevVec<double> dbeta((size_t)(g.nx * g.ny)), dupx((size_t)(g.nx * g.ny)), dupy((size_t)(g.nx * g.ny));
	bed_column_geom_gpu(df.p, dnbx.p, dnby.p, dnbz.p, dbeta.p, dupx.p, dupy.p, g, SED_CPACK); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(dbeta.download(), beta), 1e-5);
	EXPECT_LT(rel_maxnorm(dupx.download(), upx), 1e-5);
	EXPECT_LT(rel_maxnorm(dupy.download(), upy), 1e-5);
}

TEST(BedState, BedloadFluxParity)
{
	MacGrid g = grid(); int ncol = g.nx * g.ny;
	auto tx = uni(ncol, 5, -3.0, 3.0), ty = uni(ncol, 6, -3.0, 3.0);
	auto beta = uni(ncol, 7, 0.0, 0.3), upx = uni(ncol, 8, -1.0, 1.0), upy = uni(ncol, 9, -1.0, 1.0);
	MorphoParams p; p.d50 = 0.35e-3; p.rho = 1027; p.rho_s = 2650; p.nu = 1.36e-6;
	std::vector<double> qx, qy; bedload_flux_cpu(tx, ty, beta, upx, upy, qx, qy, g, p);
	DevVec<double> dtx(tx), dty(ty), dbeta(beta), dupx(upx), dupy(upy), dqx((size_t)ncol), dqy((size_t)ncol);
	bedload_flux_gpu(dtx.p, dty.p, dbeta.p, dupx.p, dupy.p, dqx.p, dqy.p, g, p); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(dqx.download(), qx), 1e-5);
	EXPECT_LT(rel_maxnorm(dqy.download(), qy), 1e-5);
}

TEST(BedState, BedloadDivParityAndValue)
{
	MacGrid g = grid(); int ncol = g.nx * g.ny;
	auto qx = uni(ncol, 10, -1e-4, 1e-4), qy = uni(ncol, 11, -1e-4, 1e-4);
	auto dx = uni(ncol, 12, -1.0, 1.0), dy = uni(ncol, 13, -1.0, 1.0); // random transport dirs
	for (int per : {0, 1})
	{
		std::vector<double> ref; bedload_div_cpu(qx, qy, dx, dy, ref, g, per);
		DevVec<double> dqx(qx), dqy(qy), ddx(dx), ddy(dy), ddiv((size_t)ncol);
		bedload_div_gpu(dqx.p, dqy.p, ddx.p, ddy.p, ddiv.p, g, per); cudaDeviceSynchronize();
		EXPECT_LT(rel_maxnorm(ddiv.download(), ref), 1e-5) << "periodic=" << per;
	}
	// value: qx = a·x (all positive), transport +x ⇒ upwind backward diff = a exactly; qy=0.
	double a = 1e-4; std::vector<double> lx(ncol), ly(ncol, 0.0), dirx(ncol, 1.0), diry(ncol, 0.0);
	for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) lx[j * g.nx + i] = a * (i + 0.5) * g.h;
	std::vector<double> dv; bedload_div_cpu(lx, ly, dirx, diry, dv, g, 0);
	for (int j = 0; j < g.ny; ++j) for (int i = 1; i < g.nx - 1; ++i) // interior (skip both wall cells)
		EXPECT_NEAR(dv[j * g.nx + i], a, 1e-9);
	// periodic ⇒ backward diff holds in every cell (no walls).
	std::vector<double> dvp; bedload_div_cpu(lx, ly, dirx, diry, dvp, g, 1);
	for (int j = 0; j < g.ny; ++j) for (int i = 1; i < g.nx; ++i)
		EXPECT_NEAR(dvp[j * g.nx + i], a, 1e-9);
	// open streamwise faces (xopen=1, non-periodic): the inlet cell gets a zero-gradient ghost ⇒ div=0
	// (no spurious erosion), and every downstream cell — including the outlet (free export) — ⇒ a.
	std::vector<double> dvo; bedload_div_cpu(lx, ly, dirx, diry, dvo, g, 0, /*xopen*/ 1);
	for (int j = 0; j < g.ny; ++j)
	{
		EXPECT_NEAR(dvo[j * g.nx + 0], 0.0, 1e-9);
		for (int i = 1; i < g.nx; ++i) EXPECT_NEAR(dvo[j * g.nx + i], a, 1e-9);
	}
	// GPU-vs-CPU parity for the open path (random fields).
	std::vector<double> refo; bedload_div_cpu(qx, qy, dx, dy, refo, g, 0, 1);
	DevVec<double> oqx(qx), oqy(qy), odx(dx), ody(dy), odiv((size_t)ncol);
	bedload_div_gpu(oqx.p, oqy.p, odx.p, ody.p, odiv.p, g, 0, 1); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(odiv.download(), refo), 1e-5) << "xopen parity";
}

TEST(BedState, MorphoExnerParity)
{
	MacGrid g = grid(); int ncol = g.nx * g.ny, np = g.p_count();
	auto hG = uni(ncol, 20, 2.0 * g.h * SED_CPACK, 8.0 * g.h * SED_CPACK);
	auto hc = uni(np, 21, 0.0, 0.02);
	auto tx = uni(ncol, 22, 0.0, 1.5), ty = uni(ncol, 23, -0.5, 0.5);
	auto beta = uni(ncol, 24, 0.0, 0.2), upx = uni(ncol, 25, -1.0, 1.0), upy = uni(ncol, 26, -1.0, 1.0);
	auto div = uni(ncol, 27, -1e-6, 1e-6);
	MorphoParams p; p.d50 = 0.2e-3; p.rho = 1027; p.rho_s = 2650; p.nu = 1.36e-6;
	p.ws0 = settling_ws(p.d50, p.rho, p.rho_s, p.nu); p.morfac = 1.0; p.hindered = 1;
	p.deposition_on = 1; p.erosion_on = 1; p.bedload_on = 1;
	double dt = 0.05;
	std::vector<double> Gc = hG, cc = hc; std::vector<double> clip(ncol, 0.0);
	std::vector<double> depc(ncol, 0.0), eroc(ncol, 0.0); // applied deposition/pickup grain (GUI diag)
	morpho_exner_cpu(Gc, &cc, tx, ty, beta, upx, upy, div, &clip, g, p, dt, depc.data(), eroc.data());
	DevVec<double> dG(hG), dc(hc), dtx(tx), dty(ty), dbeta(beta), dupx(upx), dupy(upy), ddiv(div), dclip((size_t)ncol);
	DevVec<double> ddep((size_t)ncol), dero((size_t)ncol);
	morpho_exner_gpu(dG.p, dc.p, dtx.p, dty.p, dbeta.p, dupx.p, dupy.p, ddiv.p, dclip.p, g, p, dt, ddep.p, dero.p); cudaDeviceSynchronize();
	EXPECT_LT(rel_maxnorm(dG.download(), Gc), 1e-5);
	EXPECT_LT(rel_maxnorm(dc.download(), cc), 1e-5);
	EXPECT_LT(rel_maxnorm(ddep.download(), depc), 1e-5) << "dep_out parity";
	EXPECT_LT(rel_maxnorm(dero.download(), eroc), 1e-5) << "ero_out parity";
	double sdep = 0.0, sero = 0.0; for (int n = 0; n < ncol; ++n) { sdep += depc[n]; sero += eroc[n]; }
	EXPECT_GT(sdep, 0.0) << "deposition should occur"; // genuine signal, not all-zero
	EXPECT_GT(sero, 0.0) << "erosion should occur";
}

TEST(BedState, MorphoExnerConservesGrain)
{
	// deposition+erosion move grain between G and c by exactly ∓dG_de ⇒ V_bed+V_susp invariant.
	MacGrid g; g.nx = 5; g.ny = 5; g.nz = 16; g.h = 0.05; int ncol = g.nx * g.ny, np = g.p_count();
	std::vector<double> hG(ncol, 4.0 * g.h * SED_CPACK), hc(np, 0.005);
	auto tx = uni(ncol, 30, 0.0, 1.0); std::vector<double> ty(ncol, 0.0);
	std::vector<double> beta(ncol, 0.0), upx(ncol, 1.0), upy(ncol, 0.0), div(ncol, 0.0);
	MorphoParams p; p.d50 = 0.2e-3; p.rho = 1027; p.rho_s = 2650; p.nu = 1.36e-6;
	p.ws0 = settling_ws(p.d50, p.rho, p.rho_s, p.nu); p.hindered = 1; p.bedload_on = 0;
	double dt = 0.02, Vc = g.h * g.h * g.h, Va = g.h * g.h;
	auto total = [&](const std::vector<double>& G, const std::vector<double>& c) {
		double s = 0; for (double v : G) s += v * Va; for (double v : c) s += v * Vc; return s;
	};
	double m0 = total(hG, hc);
	for (int s = 0; s < 50; ++s) { std::vector<double> clip; morpho_exner_cpu(hG, &hc, tx, ty, beta, upx, upy, div, nullptr, g, p, dt); }
	double m1 = total(hG, hc);
	EXPECT_LT(std::fabs(m1 - m0) / m0, 1e-12);
}

// ===================== avalanche =========================================================
TEST(Avalanche, SweepParity)
{
	MacGrid g; g.nx = 11; g.ny = 9; g.nz = 1; g.h = 0.05; int ncol = g.nx * g.ny;
	auto z = uni(ncol, 40, 0.0, 1.0);
	for (int per : {0, 1})
	{
		AvalancheParams ap; ap.periodic = per;
		std::vector<double> ref; avalanche_sweep_cpu(z, ref, g, ap);
		DevVec<double> dz(z), dout((size_t)ncol);
		avalanche_sweep_gpu(dz.p, dout.p, g, ap); cudaDeviceSynchronize();
		EXPECT_LT(rel_maxnorm(dout.download(), ref), 1e-5) << "periodic=" << per;
	}
}

TEST(Avalanche, MassConservationAndRepose)
{
	// A rough field relaxes to slopes ≤ tan(30°) with grain mass exactly conserved.
	MacGrid g; g.nx = 31; g.ny = 31; g.nz = 1; g.h = 0.05; int ncol = g.nx * g.ny;
	auto z = uni(ncol, 41, 0.0, 2.0);
	AvalancheParams ap; ap.tan_repose = std::tan(30.0 * PI / 180.0); ap.relax = 0.1; ap.periodic = 0;
	double m0 = 0; for (double v : z) m0 += v;
	std::vector<double> cur = z, nxt;
	for (int s = 0; s < 20000; ++s) { avalanche_sweep_cpu(cur, nxt, g, ap); cur.swap(nxt); }
	double m1 = 0; for (double v : cur) m1 += v;
	EXPECT_LT(std::fabs(m1 - m0) / m0, 1e-12);
	EXPECT_LE(avalanche_max_tanslope(cur, g, 0), ap.tan_repose * 1.03);
}
