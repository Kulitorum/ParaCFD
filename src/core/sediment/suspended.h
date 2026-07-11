// suspended.h — suspended-sediment concentration transport (RESEARCH §6.1–6.2, research/02).
// One volumetric concentration field c per active grain size on the MAC cell centers:
//
//   ∂c/∂t + ∇·(u c) − ∂(w_s c)/∂z = ∇·((ν_t/σ_s)∇c)                         (RESEARCH §6.2)
//
// implemented operator-split as (1) MacCormack semi-Lagrangian advection with the settling
// velocity FOLDED INTO the advection velocity (vz = w_fluid − w_s; RESEARCH §6.2), (2) explicit
// variable-coefficient diffusion with D_c = ν_t/σ_s, and (3) the flux-form bed exchange
// E/ρs (pickup, θ_cr-gated) − w_s·c_b (deposition, NEVER thresholded — RESEARCH §5 asymmetry).
//
// Sediment consumes only {u,v,w,ν_t,τ_b} from the FluidCore (PLAN §3): u/v/w are the fluid face
// velocities, τ′_b the grain-related bed shear (bedshear.*), ν_t the eddy viscosity. Every kernel
// has a CPU reference twin exercised by a GPU-vs-CPU gtest at rel. max-norm 1e-5 (CLAUDE.md).
#pragma once

#include "core/fluid/mac_grid.h"
#include "core/sediment/sed_physics.h"

#include <vector>

namespace scour::core
{
	// Per-face scalar boundary condition for the concentration field c.
	enum ScalarFace : int
	{
		SC_NEUMANN = 0,  // zero-gradient (ghost = nearest interior). Walls / homogeneous directions.
		SC_DIRICHLET = 1 // fixed far-field value (clear-water c=0 at the surface; pinned c_a for Rouse)
	};

	// Cell-centered scalar BC on the six domain faces, with the Dirichlet values.
	struct ScalarBC
	{
		int xmin = SC_NEUMANN, xmax = SC_NEUMANN;
		int ymin = SC_NEUMANN, ymax = SC_NEUMANN;
		int zmin = SC_NEUMANN, zmax = SC_NEUMANN;
		double v_xmin = 0.0, v_xmax = 0.0, v_ymin = 0.0, v_ymax = 0.0, v_zmin = 0.0, v_zmax = 0.0;
	};

	// Bed-exchange parameters (RESEARCH §6.2). d50/ρ/ρs/ν give D*, τ_cr, θ_cr at runtime.
	struct SedBedParams
	{
		double d50 = 0.2e-3;   // m
		double rho = 1027.0;   // kg/m³
		double rho_s = 2650.0; // kg/m³
		double nu = 1.36e-6;   // m²/s
		double alpha = 0.00033;// van Rijn pickup coefficient (±30%), RESEARCH §6.2
		double ws0 = 0.0;      // clear-water settling velocity w_s(d50) [m/s] (precomputed from config)
		int hindered = 1;      // apply hindered w_s in the deposition flux (RESEARCH §6.1)
	};

	// ---- (0) Effective settling field w_s,eff = w_s0·(1−c)^4.7 (hindered) at cell centers -------
	void suspended_effective_ws_gpu(const double* c, double* ws, MacGrid g, double ws0, int hindered);
	void suspended_effective_ws_cpu(const std::vector<double>& c, std::vector<double>& ws, MacGrid g, double ws0, int hindered);

	// ---- (1a) MacCormack semi-Lagrangian advection of c, settling folded into vz (RESEARCH §6.2) -
	// u/v/w = fluid face velocities (vbc = fluid velocity BC for sampling); ws = cell-centered
	// effective settling [m/s] (subtracted from w_fluid). scratch ≥ p_count doubles (phiHat).
	// band = # cells near a wall to revert to 1st order (0 for obstacle-free columns).
	// NOTE: MacCormack SL is 2nd-order accurate but NOT strictly mass-conservative (~0.1–0.3% drift
	// floor, independent of resolution). For the sand mass-budget path use suspended_advect_cons.
	void suspended_advect_gpu(const double* cIn, const double* u, const double* v, const double* w,
		const double* ws, double* cOut, double* scratch, MacGrid g, BC vbc, ScalarBC sbc, double dt, int band);
	void suspended_advect_cpu(const std::vector<double>& cIn, const std::vector<double>& u, const std::vector<double>& v,
		const std::vector<double>& w, const std::vector<double>& ws, std::vector<double>& cOut,
		MacGrid g, BC vbc, ScalarBC sbc, double dt, int band);

	// ---- (1b) CONSERVATIVE flux-form advection of c, settling folded into the z-flux velocity ----
	// Dimensionally-split (x→y→z) TVD van-Leer upwind: F_face = v_face·c_face(limited); c updated as
	// c −= (dt/h)·ΔF. z-face velocity = w_fluid − w_s,face (settling, RESEARCH §6.2). ALL six domain
	// faces carry ZERO advective flux (settling cannot cross the bed/surface; the bed exchange sink/
	// source is applied separately). This conserves suspended sand mass to machine precision — the
	// property the sand-trapping KPIs require — and, being flux-form, correctly depletes the surface
	// so the Rouse profile develops from any initial state (unlike SL, whose uniform state is a
	// spurious fixed point). scratch ≥ p_count doubles. Every sweep has a CPU twin (parity 1e-5).
	void suspended_advect_cons_gpu(const double* cIn, const double* u, const double* v, const double* w,
		const double* ws, double* cOut, double* scratch, MacGrid g, double dt);
	void suspended_advect_cons_cpu(const std::vector<double>& cIn, const std::vector<double>& u, const std::vector<double>& v,
		const std::vector<double>& w, const std::vector<double>& ws, std::vector<double>& cOut, MacGrid g, double dt);

	// ---- (2) Explicit variable-coefficient diffusion  c += dt·∇·(D_c ∇c) (RESEARCH §6.2) ---------
	// Dc = cell-centered diffusivity [m²/s] (= ν_t/σ_s + molecular). Face coeff = ½(Dc[c]+Dc[nb]);
	// zero flux across SC_NEUMANN faces, Dirichlet ghost value across SC_DIRICHLET faces.
	void suspended_diffuse_gpu(const double* cIn, double* cOut, const double* Dc, MacGrid g, ScalarBC sbc, double dt);
	void suspended_diffuse_cpu(const std::vector<double>& cIn, std::vector<double>& cOut, const std::vector<double>& Dc,
		MacGrid g, ScalarBC sbc, double dt);

	// ---- (3) Flux-form bed exchange at the bed-adjacent cells (k=0) (RESEARCH §6.2) --------------
	// Applies dc = (dt/h)·(E/ρs − w_s·c_b): pickup E (van Rijn 2019, θ_cr-gated, τ_prime = grain
	// skin stress plane) and deposition w_s·c_b (hindered, NEVER thresholded). Writes the applied
	// NET volumetric bed flux per column (E/ρs − w_s·c_b)·dt [m] to `bedflux` (size nx·ny) so the
	// driver can close the sand-mass budget (M_susp + M_bed = const). c is clamped ≥ 0.
	void suspended_bed_exchange_gpu(double* c, const double* tau_prime, double* bedflux, MacGrid g, SedBedParams p, double dt);
	void suspended_bed_exchange_cpu(std::vector<double>& c, const std::vector<double>& tau_prime, std::vector<double>& bedflux,
		MacGrid g, SedBedParams p, double dt);

	// ---- (4) OPEN-X boundary advective flux (RESEARCH §6; open-sea inflow / outflow) -------------
	// suspended_advect_cons leaves ALL six faces closed (zero-flux). This re-opens the two x-faces
	// as a mass-EXACT operator-split source applied AFTER the conservative advection, so the gate
	// kernel stays byte-identical and the closed-domain (M4/M5, seabed_flat) budget is untouched
	// unless this is called. First-order upwind per x-face: on INFLOW the face carries the prescribed
	// far-field ghost concentration; on OUTFLOW it carries the interior cell — so the equilibrium load
	// enters at the inlet and the developed load leaves at the outlet, and TIDAL REVERSAL is handled
	// automatically (each face upwinds on its own velocity sign). The other four faces stay closed
	// (settling cannot cross the bed/surface; y walls are homogeneous). ghost_xmin / ghost_xmax are
	// the far-field concentration on the yz-plane, indexed [k·ny + j] (size ny·nz): the driver fills
	// them with the Rouse equilibrium profile (OPEN mode) or the recycled opposite plane (RECYCLE
	// mode). If `net_flux` (size 2·ny·nz, xmin block [0,ny·nz) then xmax block) is non-null it receives
	// the signed sand volume [m³] added to the domain at each boundary cell this step (Σ = net domain
	// sand gain) so the driver can close the OPEN budget ΔM_bed + ΔM_susp = ∫boundary. c updated in
	// place; each (j,k) touches only its own two boundary cells (cell 0 and cell nx−1), race-free.
	void suspended_boundary_flux_gpu(double* c, const double* u, const double* ghost_xmin, const double* ghost_xmax,
		double* net_flux, MacGrid g, double dt);
	void suspended_boundary_flux_cpu(std::vector<double>& c, const std::vector<double>& u,
		const std::vector<double>& ghost_xmin, const std::vector<double>& ghost_xmax,
		std::vector<double>* net_flux, MacGrid g, double dt);
}
