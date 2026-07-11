// bedshear.h — log-law wall function producing the bed shear stress τ_b for the Stam
// fluid core (RESEARCH §4, research/04; cavity rules RESEARCH §4 / research/13 §8).
//
// The sediment module consumes τ_b through {u,v,w,ν_t,τ_b} only (PLAN §3): this file is
// the τ_b provider. Everything physical is a pure __host__ __device__ inline so the CUDA
// kernels and their CPU reference twins share identical arithmetic (GPU-vs-CPU parity
// tests at rel. max-norm 1e-5, CLAUDE.md hard rule).
//
// Core wall function (research/04 §5):  u* = κ·U_p / ln(z_p/z0),  τ_b = ρ·u*²,
// with U_p the tangential velocity probed at a FIXED normal distance z_p above the bed
// (never a stair-step first cell), z0 = d50/12 (grain roughness ks = 2.5·d50). Rough
// regime is explicit; transitional (ks+ = u*ks/ν in 5..70) uses the Christoffersen–
// Jonsson z0 fixed-point (research/04 §4,§6). Anti-stair-step mitigations (research/04):
// interpolate u* (not u), 3×3 tangential smoothing of τ_b, rate-limit wall updates to
// [0.5×,2×]/step, and EMA time-filter τ_b over 1..3 s before erosion/transport.
// Units: SI (m, s, m/s, Pa). κ, ks, z0 from RESEARCH.md verbatim.
#pragma once

#include "core/fluid/mac_grid.h"

#include <cmath>
#include <vector>

namespace scour::core
{
	// Roughness/regime selection for the u* extraction.
	enum WallRegime : int
	{
		WALL_AUTO = 0,         // pick rough/transitional by ks+ (research/04 §2) — physical default
		WALL_ROUGH = 1,        // force fully-rough explicit branch (z0 = ks/30 = d50/12)
		WALL_TRANSITIONAL = 2  // force Christoffersen–Jonsson fixed-point (all regimes)
	};

	// Wall-function parameters. ks/z0 are DERIVED from d50 at use sites (never hard-coded).
	struct WallParams
	{
		double kappa = 0.40;   // von Kármán (RESEARCH §4 verbatim)
		double d50 = 0.2e-3;   // m — grain size
		double nu = 1.36e-6;   // m²/s — molecular viscosity (config, RESEARCH §1)
		double rho = 1027.0;   // kg/m³ (RESEARCH §1)
		int regime = WALL_AUTO;
		int cj_iters = 5;      // Christoffersen–Jonsson fixed-point sweeps (research/04 §6: 2–5)
		double t_avg = 2.0;    // s — EMA window for τ_b (research/04 LES note, RESEARCH §4: 1–3 s)
		int smooth = 1;        // 3×3 tangential smoothing passes of τ_b (0 = off)
		double rate_lo = 0.5;  // per-step wall-update clamp (OpenFOAM practice, research/04 §7)
		double rate_hi = 2.0;
	};

	SCOUR_HD inline double wall_ks(double d50) { return 2.5 * d50; }        // research/04 §3
	SCOUR_HD inline double wall_z0_grain(double d50) { return d50 / 12.0; } // = ks/30

	// Guarded rough inversion: u* = κ·U_p/ln(z_p/z0), clamp ln arg ≥ e so u* ≤ κ·U_p
	// (research/04 §6 guard). z0 must be > 0.
	SCOUR_HD inline double wall_ustar_rough(double Up, double zp, double z0, double kappa)
	{
		double arg = zp / z0;
		if (arg < 2.718281828459045) arg = 2.718281828459045;
		return kappa * Up / log(arg);
	}

	// Christoffersen–Jonsson transitional roughness length given the current u*
	// (research/04 §4): z0 = (ks/30)(1 − e^(−u*·ks/27ν)) + ν/(9u*).
	SCOUR_HD inline double wall_z0_cj(double ustar, double ks, double nu)
	{
		double u = ustar > 1e-12 ? ustar : 1e-12;
		return (ks / 30.0) * (1.0 - exp(-u * ks / (27.0 * nu))) + nu / (9.0 * u);
	}

	// Full u* extraction from a probed tangential speed U_p at height z_p (research/04
	// §5,§6). AUTO selects rough vs transitional by ks+ from a rough first guess.
	SCOUR_HD inline double wall_ustar(double Up, double zp, double d50, double nu, double kappa, int regime, int cj_iters)
	{
		double ks = wall_ks(d50);
		double us = wall_ustar_rough(Up, zp, ks / 30.0, kappa); // rough guess u*(0)
		if (regime == WALL_ROUGH) return us;
		if (regime == WALL_AUTO)
		{
			double ksplus = us * ks / nu;            // grain Reynolds number
			if (ksplus > 70.0) return us;            // fully rough → explicit (research/04 §2)
		}
		// Transitional (or AUTO with ks+ ≤ 70): Christoffersen–Jonsson fixed point.
		for (int it = 0; it < cj_iters; ++it)
		{
			double z0 = wall_z0_cj(us, ks, nu);
			us = wall_ustar_rough(Up, zp, z0, kappa);
		}
		return us;
	}

	// Probe height z_p = max(1.5·h, 2·ks) above the bed (research/04 §practical).
	SCOUR_HD inline double wall_zp(double h, double ks) { return fmax(1.5 * h, 2.0 * ks); }

	// ---- Cavity clearance four-branch rule (RESEARCH §4 / research/13 §8) ----------
	// Given clearance H_c (raycast bed-normal distance to the first solid), the probe
	// height z_p0 = max(1.5h,2ks), the mid-gap speed U_gap, and z0, return u* and set
	// `freeze` when erosion must be frozen (H_c < 3h → sub-grid screen model regime).
	enum CavityBranch : int
	{
		CAV_STANDARD = 0,   // H_c ≥ 4·z_p0
		CAV_QUARTER = 1,    // 6h ≤ H_c < 4·z_p0 : probe at H_c/4
		CAV_GAP_TURB = 2,   // 3h ≤ H_c < 6h, Re_gap > 2800 : mid-gap log law
		CAV_GAP_LAMINAR = 3,// 3h ≤ H_c < 6h, Re_gap ≤ 2800 : plane-Poiseuille τ = 6μU/H_c
		CAV_FROZEN = 4      // H_c < 3h : freeze erosion, screen model
	};

	// Classify the cavity branch from clearance alone (probe placement is done by caller).
	SCOUR_HD inline int cavity_branch(double Hc, double h, double zp0, double Ugap, double nu)
	{
		if (Hc >= 4.0 * zp0) return CAV_STANDARD;
		if (Hc >= 6.0 * h) return CAV_QUARTER;
		if (Hc >= 3.0 * h)
		{
			double Re_gap = Ugap * Hc / nu;
			return Re_gap > 2800.0 ? CAV_GAP_TURB : CAV_GAP_LAMINAR;
		}
		return CAV_FROZEN;
	}

	// τ_b [Pa] for the turbulent/laminar mid-gap branches (RESEARCH §4). For the log-law
	// branches the caller uses wall_ustar()+τ=ρu*²; this covers the gap closures only.
	SCOUR_HD inline double cavity_tau_gap_turb(double Ugap, double Hc, double z0, double rho, double kappa)
	{
		double arg = Hc / (2.0 * z0);
		if (arg < 2.718281828459045) arg = 2.718281828459045;
		double cd = kappa / log(arg);
		return rho * cd * cd * Ugap * Ugap;
	}
	SCOUR_HD inline double cavity_tau_gap_laminar(double Ugap, double Hc, double rho, double nu)
	{
		return 6.0 * (rho * nu) * Ugap / Hc; // μ = ρ·ν ; plane-Poiseuille (research/13 §8)
	}

	// ---- Vertical sampling of the tangential (u,v) speed at height zp above a FLAT bed.
	// Bilinear-in-xy at the column center is unnecessary for the flat M3 bed (homogeneous
	// in xy for the periodic gate); we average the two x-faces / y-faces of cell (i,j) and
	// linearly interpolate in z. Returns u_c,v_c at zp. (Voxel-bed normal probing arrives
	// with the M5 interface geometry.)
	SCOUR_HD inline void bed_sample_uv(const double* u, const double* v, MacGrid g, int i, int j, double zp, double& uc, double& vc)
	{
		double gz = zp / g.h - 0.5;
		int k0 = (int)floor(gz);
		if (k0 < 0) k0 = 0;
		if (k0 > g.nz - 2) k0 = g.nz - 2;
		double fz = gz - k0;
		if (fz < 0.0) fz = 0.0; else if (fz > 1.0) fz = 1.0;
		double u0 = 0.5 * (u[g.uidx(i, j, k0)] + u[g.uidx(i + 1, j, k0)]);
		double u1 = 0.5 * (u[g.uidx(i, j, k0 + 1)] + u[g.uidx(i + 1, j, k0 + 1)]);
		double v0 = 0.5 * (v[g.vidx(i, j, k0)] + v[g.vidx(i, j + 1, k0)]);
		double v1 = 0.5 * (v[g.vidx(i, j, k0 + 1)] + v[g.vidx(i, j + 1, k0 + 1)]);
		uc = u0 * (1.0 - fz) + u1 * fz;
		vc = v0 * (1.0 - fz) + v1 * fz;
	}
}
