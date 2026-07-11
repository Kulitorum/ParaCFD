// channel_bc.h — open-channel boundary descriptor and BC-aware sampling for M2.
// Extends the M1 closed-box MAC grid to an OPEN channel with an obstacle mask:
//   * xmin = INLET  (Dirichlet velocity: uniform U or rough log-law profile; cross-
//                    flow v=w=0). Pressure Neumann.  RESEARCH §8.1.
//   * xmax = OUTLET (Orlanski-type convective BC + global flux rescale; pressure
//                    Dirichlet p=0). RESEARCH §3, §8, research/08 §8.
//   * y,z faces = free-slip lateral / rigid free-slip lid & floor (RESEARCH §8.3-8.6).
//   * interior voxel solids: normal velocity u·n=0 (both surface BC modes); tangential
//     no-slip (resolved cylinder validation) or free-slip (production, RESEARCH §3).
//
// Sampling helpers mirror mac_grid.h but resolve the x-normal ghosts as inlet/outlet.
// Solid faces are held at 0 in the stored field (see channel_ops apply_solid_bc), so
// the trilinear samplers need no mask; the mask enters advection (near-solid reversion
// + backtrace clipping) and diffusion (tangential no-slip reflection) explicitly.
// Units: SI (m, s, m/s, Pa). Fields are double (M2 grids are ≤ a few M cells).
#pragma once

#include "core/fluid/mac_grid.h"

#include <cmath>

namespace windcfd::core
{
	enum InletMode : int
	{
		INLET_UNIFORM = 0, // u(z) = U_inlet (top-hat) — cylinder benchmark
		INLET_LOGLAW = 1   // u(z) = (u*/κ) ln(z/z0)  — open-channel / scour runs
	};

	// Interior-solid tangential condition.
	enum SolidBC : int
	{
		SOLID_FREESLIP = 0, // ∂u_t/∂n = 0 (RESEARCH §3 production default)
		SOLID_NOSLIP = 1    // u_t = 0 at the surface (resolved-BL cylinder validation)
	};

	struct ChannelBC
	{
		// Lateral / vertical faces (x is always inlet/outlet). Free-slip everywhere for
		// the M2 gate; WALL_NOSLIP is available for completeness (reflection ghost).
		int ymin = WALL_FREESLIP, ymax = WALL_FREESLIP;
		int zmin = WALL_FREESLIP, zmax = WALL_FREESLIP;

		// Inlet profile (Dirichlet). uniform -> U_inlet; loglaw -> (ustar/kappa) ln((z-bed_datum)/z0).
		int inlet_mode = INLET_UNIFORM;
		double U_inlet = 1.0; // m/s (uniform mode; also the reference bulk speed)
		double ustar = 0.0;   // m/s (loglaw mode)
		double z0 = 1.0;      // m   (loglaw mode roughness length = d50/12)
		double kappa = 0.40;  // von Kármán
		// Bed-top elevation at the inlet plane [m]: the log-law is referenced to it so u→0 AT the sand
		// surface (not at z=0). 0 ⇒ bed at the domain floor (M3 flat-bed channel; keeps loglaw identical).
		double bed_datum = 0.0;

		// Orlanski convective outflow speed (bulk). Clamped so Uc*dt/h ≤ 1 by the driver.
		double Uc = 1.0; // m/s

		int solid_mode = SOLID_NOSLIP; // interior-obstacle tangential BC

		// Streamwise flow direction for TIDAL REVERSAL (M9, RESEARCH §8). +1 (default): inlet on the
		// xmin face (i=0), Orlanski outlet on xmax (i=nx) — the M2/M3 open channel, unchanged. −1: the
		// faces SWAP — the inlet Dirichlet profile is imposed on xmax (u<0, flowing −x) and the outlet
		// is on xmin. U_inlet/ustar stay MAGNITUDES; this sign carries the direction. The tidal driver
		// flips it at slack (where |U|→0), so the swap is smooth. All flow_sign≥0 branches reduce to the
		// original code ⇒ M2/M3 byte-identical (they never set it; the gates depend on that).
		int flow_sign = 1;
	};

	// Inlet u on the x-face plane i=0, at vertical node z=(k+0.5)h. For loglaw the height is taken
	// ABOVE the inlet bed top (bed_datum), so the profile rises from ~0 at the sand surface — this is
	// what prevents the top-hat's full-U-at-the-bed spurious leading-edge scour (a developed BL inflow).
	WINDCFD_HD inline double channel_inlet_u(ChannelBC bc, MacGrid g, int k)
	{
		double up;
		if (bc.inlet_mode == INLET_UNIFORM) up = bc.U_inlet;
		else
		{
			double z = (k + 0.5) * g.h - bc.bed_datum;         // height above the inlet bed top
			if (z < bc.z0 * 1.0000001) z = bc.z0 * 1.0000001;  // guard log domain; ⇒ ~0 at/below the bed
			up = (bc.ustar / bc.kappa) * log(z / bc.z0);
		}
		return bc.flow_sign < 0 ? -up : up; // reversed tide: inlet drives −x (U_inlet stays a magnitude)
	}

	// Pick u* so the log-law u(z)=(u*/κ)ln(z/z0) is FLUX-matched to a target depth-averaged current U
	// over a fluid depth H above the bed: U = (u*/κ)(ln(H/z0) − 1 + z0/H). Guards a thin/invalid depth.
	WINDCFD_HD inline double loglaw_ustar_for_U(double U, double H, double z0, double kappa)
	{
		if (z0 <= 0.0 || H <= z0 * 2.718281828459045) return kappa * U; // degenerate: safe fallback
		double denom = log(H / z0) - 1.0 + z0 / H;
		if (denom < 0.1) denom = 0.1;
		return kappa * U / denom;
	}

	// ---- BC-aware ghost fetch of the three velocity components -------------------
	// Normal x-ghosts: i<0 -> inlet (u=profile; v,w=0), i>range -> outlet (zero-gradient
	// copy). Tangential y/z-ghosts: free-slip mirror (or no-slip reflection).

	WINDCFD_HD inline double ch_fetch_u(const double* u, MacGrid g, ChannelBC bc, int i, int j, int k)
	{
		double factor = 1.0;
		if (j < 0) { if (bc.ymin == WALL_NOSLIP) factor = -factor; j = 0; }
		else if (j > g.ny - 1) { if (bc.ymax == WALL_NOSLIP) factor = -factor; j = g.ny - 1; }
		if (k < 0) { if (bc.zmin == WALL_NOSLIP) factor = -factor; k = 0; }
		else if (k > g.nz - 1) { if (bc.zmax == WALL_NOSLIP) factor = -factor; k = g.nz - 1; }
		if (bc.flow_sign < 0)
		{
			if (i > g.nx) return factor * channel_inlet_u(bc, g, k); // reversed: Dirichlet inlet on xmax
			if (i < 0) i = 0;                                        // reversed: outlet zero-gradient copy at xmin
		}
		else
		{
			if (i < 0) return factor * channel_inlet_u(bc, g, k); // Dirichlet inlet profile at xmin
			if (i > g.nx) i = g.nx;                               // outlet: zero-gradient copy at xmax
		}
		return factor * u[g.uidx(i, j, k)];
	}

	WINDCFD_HD inline double ch_fetch_v(const double* v, MacGrid g, ChannelBC bc, int i, int j, int k)
	{
		double factor = 1.0;
		bool inlet_ghost = false;
		if (bc.flow_sign < 0)
		{
			if (i > g.nx - 1) { inlet_ghost = true; i = g.nx - 1; } // reversed: inlet on xmax
			else if (i < 0) { i = 0; }                             // reversed: outlet copy at xmin
		}
		else
		{
			if (i < 0) { inlet_ghost = true; i = 0; }
			else if (i > g.nx - 1) { i = g.nx - 1; } // outlet zero-gradient copy
		}
		if (k < 0) { if (bc.zmin == WALL_NOSLIP) factor = -factor; k = 0; }
		else if (k > g.nz - 1) { if (bc.zmax == WALL_NOSLIP) factor = -factor; k = g.nz - 1; }
		if (j < 0) j = 0; else if (j > g.ny) j = g.ny; // v normal
		if (inlet_ghost) return 0.0;                   // cross-flow = 0 at inlet
		return factor * v[g.vidx(i, j, k)];
	}

	WINDCFD_HD inline double ch_fetch_w(const double* w, MacGrid g, ChannelBC bc, int i, int j, int k)
	{
		double factor = 1.0;
		bool inlet_ghost = false;
		if (bc.flow_sign < 0)
		{
			if (i > g.nx - 1) { inlet_ghost = true; i = g.nx - 1; } // reversed: inlet on xmax
			else if (i < 0) { i = 0; }                             // reversed: outlet copy at xmin
		}
		else
		{
			if (i < 0) { inlet_ghost = true; i = 0; }
			else if (i > g.nx - 1) { i = g.nx - 1; } // outlet zero-gradient copy
		}
		if (j < 0) { if (bc.ymin == WALL_NOSLIP) factor = -factor; j = 0; }
		else if (j > g.ny - 1) { if (bc.ymax == WALL_NOSLIP) factor = -factor; j = g.ny - 1; }
		if (k < 0) k = 0; else if (k > g.nz) k = g.nz; // w normal
		if (inlet_ghost) return 0.0;                   // cross-flow = 0 at inlet
		return factor * w[g.widx(i, j, k)];
	}

	// ---- Trilinear samplers (component node spaces identical to mac_grid.h) ------
	WINDCFD_HD inline double ch_trilerp_u(const double* u, MacGrid g, ChannelBC bc, double x, double y, double z, double* mn, double* mx)
	{
		double gx = x / g.h, gy = y / g.h - 0.5, gz = z / g.h - 0.5;
		int i0 = (int)floor(gx), j0 = (int)floor(gy), k0 = (int)floor(gz);
		double fx = gx - i0, fy = gy - j0, fz = gz - k0;
		double c[2][2][2];
		for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int cc = 0; cc < 2; ++cc)
			c[a][b][cc] = ch_fetch_u(u, g, bc, i0 + a, j0 + b, k0 + cc);
		if (mn && mx)
		{
			double lo = c[0][0][0], hi = c[0][0][0];
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int cc = 0; cc < 2; ++cc)
			{ lo = c[a][b][cc] < lo ? c[a][b][cc] : lo; hi = c[a][b][cc] > hi ? c[a][b][cc] : hi; }
			*mn = lo; *mx = hi;
		}
		double c00 = c[0][0][0] * (1 - fx) + c[1][0][0] * fx, c10 = c[0][1][0] * (1 - fx) + c[1][1][0] * fx;
		double c01 = c[0][0][1] * (1 - fx) + c[1][0][1] * fx, c11 = c[0][1][1] * (1 - fx) + c[1][1][1] * fx;
		double c0 = c00 * (1 - fy) + c10 * fy, c1 = c01 * (1 - fy) + c11 * fy;
		return c0 * (1 - fz) + c1 * fz;
	}
	WINDCFD_HD inline double ch_trilerp_v(const double* v, MacGrid g, ChannelBC bc, double x, double y, double z, double* mn, double* mx)
	{
		double gx = x / g.h - 0.5, gy = y / g.h, gz = z / g.h - 0.5;
		int i0 = (int)floor(gx), j0 = (int)floor(gy), k0 = (int)floor(gz);
		double fx = gx - i0, fy = gy - j0, fz = gz - k0;
		double c[2][2][2];
		for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int cc = 0; cc < 2; ++cc)
			c[a][b][cc] = ch_fetch_v(v, g, bc, i0 + a, j0 + b, k0 + cc);
		if (mn && mx)
		{
			double lo = c[0][0][0], hi = c[0][0][0];
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int cc = 0; cc < 2; ++cc)
			{ lo = c[a][b][cc] < lo ? c[a][b][cc] : lo; hi = c[a][b][cc] > hi ? c[a][b][cc] : hi; }
			*mn = lo; *mx = hi;
		}
		double c00 = c[0][0][0] * (1 - fx) + c[1][0][0] * fx, c10 = c[0][1][0] * (1 - fx) + c[1][1][0] * fx;
		double c01 = c[0][0][1] * (1 - fx) + c[1][0][1] * fx, c11 = c[0][1][1] * (1 - fx) + c[1][1][1] * fx;
		double c0 = c00 * (1 - fy) + c10 * fy, c1 = c01 * (1 - fy) + c11 * fy;
		return c0 * (1 - fz) + c1 * fz;
	}
	WINDCFD_HD inline double ch_trilerp_w(const double* w, MacGrid g, ChannelBC bc, double x, double y, double z, double* mn, double* mx)
	{
		double gx = x / g.h - 0.5, gy = y / g.h - 0.5, gz = z / g.h;
		int i0 = (int)floor(gx), j0 = (int)floor(gy), k0 = (int)floor(gz);
		double fx = gx - i0, fy = gy - j0, fz = gz - k0;
		double c[2][2][2];
		for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int cc = 0; cc < 2; ++cc)
			c[a][b][cc] = ch_fetch_w(w, g, bc, i0 + a, j0 + b, k0 + cc);
		if (mn && mx)
		{
			double lo = c[0][0][0], hi = c[0][0][0];
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int cc = 0; cc < 2; ++cc)
			{ lo = c[a][b][cc] < lo ? c[a][b][cc] : lo; hi = c[a][b][cc] > hi ? c[a][b][cc] : hi; }
			*mn = lo; *mx = hi;
		}
		double c00 = c[0][0][0] * (1 - fx) + c[1][0][0] * fx, c10 = c[0][1][0] * (1 - fx) + c[1][1][0] * fx;
		double c01 = c[0][0][1] * (1 - fx) + c[1][0][1] * fx, c11 = c[0][1][1] * (1 - fx) + c[1][1][1] * fx;
		double c0 = c00 * (1 - fy) + c10 * fy, c1 = c01 * (1 - fy) + c11 * fy;
		return c0 * (1 - fz) + c1 * fz;
	}

	// ---- Solid-mask helpers (mask is a cell field: 1=solid, 0=fluid) -------------
	WINDCFD_HD inline bool ch_is_solid(const unsigned char* solid, MacGrid g, int i, int j, int k)
	{
		if (i < 0 || i > g.nx - 1 || j < 0 || j > g.ny - 1 || k < 0 || k > g.nz - 1) return false;
		return solid[g.pidx(i, j, k)] != 0;
	}
	// A u-face is "solid" (interior to the obstacle, x-normal) when both x-cells solid.
	WINDCFD_HD inline bool ch_usolid(const unsigned char* s, MacGrid g, int i, int j, int k)
	{ return ch_is_solid(s, g, i - 1, j, k) && ch_is_solid(s, g, i, j, k); }
	WINDCFD_HD inline bool ch_vsolid(const unsigned char* s, MacGrid g, int i, int j, int k)
	{ return ch_is_solid(s, g, i, j - 1, k) && ch_is_solid(s, g, i, j, k); }
	WINDCFD_HD inline bool ch_wsolid(const unsigned char* s, MacGrid g, int i, int j, int k)
	{ return ch_is_solid(s, g, i, j, k - 1) && ch_is_solid(s, g, i, j, k); }
}
