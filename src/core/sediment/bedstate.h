// bedstate.h — M5 bed state + morphodynamics (RESEARCH §5, §7; research/13, /03, /01, /16).
//
// The bed is a single 3D packed-grain-fraction field `f_pack` [-] on the MAC cell centers
// (grains only, pore water excluded; RESEARCH §7). Cell state from the normalized fill
// F = f_pack/(c_pack·(1−phi_obs)): PRINTED / PACKED (F≥1) / INTERFACE (0<F<1) / FLUID (F=0).
//
// For an OPEN BED (no printed overhangs — the M5 gate regime, and the Exner-equivalence limit
// of RESEARCH §7 eq. 6) the vertical fill is monotone: packed cells below, one interface cell,
// fluid above. The complete morphodynamic state per (x,y) column is then the grain thickness
//   G(x,y) = Σ_k f_pack·h   [m of solid]   ⇔   bed elevation  z_b = G/(1−p) = G/c_pack,
// and f_pack is reconstructed from G by `fpack_from_G`. A column sum of every source below
// reduces exactly to the Exner equation (1−p)·∂z_b/∂t = −∇·q_b − u_lift·c_pack + D
// (RESEARCH §7; research/13 §6). Overhang/layered columns (needed from M6/M8) generalise this;
// the box-filtered interface geometry here is already the general (n_b from ∇F̃) form.
//
// Every kernel has a CPU reference twin exercised by a GPU-vs-CPU gtest at rel. max-norm 1e-5
// (CLAUDE.md). Physics constants come from sed_physics.h (RESEARCH §5/§7 verbatim, config-derived).
#pragma once

#include "core/fluid/mac_grid.h"
#include "core/sediment/sed_physics.h"

#include <vector>

namespace scour::core
{
	// Cell state (diagnostic; RESEARCH §7 / research/13 §1).
	enum BedState : int
	{
		BS_FLUID = 0,
		BS_INTERFACE = 1,
		BS_PACKED = 2,
		BS_PRINTED = 3
	};

	// Morphodynamic parameters (all physics config-derived per RESEARCH §5/§7).
	struct MorphoParams
	{
		double d50 = 0.2e-3;   // m
		double rho = 1027.0;   // kg/m³
		double rho_s = 2650.0; // kg/m³
		double nu = 1.36e-6;   // m²/s
		double ws0 = 0.0;      // clear-water settling velocity w_s(d50) [m/s] (precomputed)
		double phi_repose = 0.55850536; // repose angle [rad] = 32° (RESEARCH §1; tan32=0.6249)
		double cpack = SED_CPACK;       // 0.64
		double morfac = 1.0;            // MORFAC M (RESEARCH §7; M≤10 steady, ≤5 reversing)
		double dz_limit_frac = 0.05;    // per-update |Δz_b| ≤ frac·h (RESEARCH §7 / Delft3D DzMax)
		int hindered = 1;      // hindered w_s in deposition (RESEARCH §6.1)
		double erosion_coeff = 0.018; // Winterwerp u_lift coefficient (RESEARCH §7 = 0.018); the M6
		                              // supervised erosion-rate calibration knob (default byte-identical)
		int bedload_formula = 0; // 0 = Wong–Parker (default), 1 = Engelund–Fredsøe (RESEARCH §6.3)
		int erosion_on = 1, deposition_on = 1, bedload_on = 1; // per-process toggles (gate drivers)
	};

	// ---- (1) Reconstruct the 3D f_pack field from the per-column grain thickness G ------------
	// Open-bed monotone fill: f_pack[k] = c_pack·clamp(z_b/h − k, 0, 1), z_b = G/c_pack.
	// G is size nx·ny [m]; fpack is size p_count [-]. (Exner-equivalent; RESEARCH §7 eq. 6.)
	void fpack_from_G_gpu(const double* G, double* fpack, MacGrid g, double cpack);
	void fpack_from_G_cpu(const std::vector<double>& G, std::vector<double>& fpack, MacGrid g, double cpack);

	// ---- (2a) Box-filtered normalized fill  F̃ = 3×3×3 mean of F = f_pack/(c_pack(1−phi_obs)) --
	// Ghost fill: below the floor (k<0) F=1 (packed floor), above the top (k≥nz) F=0 (fluid),
	// x/y edges replicate. phi_obs may be null (⇒ 0). RESEARCH §7 / research/13 §2.
	void bed_boxfilter_gpu(const double* fpack, const double* phi_obs, double* Ftilde, MacGrid g, double cpack);
	void bed_boxfilter_cpu(const std::vector<double>& fpack, const std::vector<double>* phi_obs, std::vector<double>& Ftilde, MacGrid g, double cpack);

	// ---- (2b) Interface normal & area from the box-filtered fill gradient ---------------------
	// n_b = −∇F̃/|∇F̃| (points into the fluid), A_b = |∇F̃|·V_cell (RESEARCH §7). Central diffs
	// with the same ghost rule as (2a). Degenerate |∇F̃|≈0 ⇒ n_b=(0,0,1), A_b=0. Outputs nbx,
	// nby,nbz (unit) and Ab [m²]; per-cell (size p_count).
	void bed_interface_geom_gpu(const double* Ftilde, double* nbx, double* nby, double* nbz, double* Ab, MacGrid g);
	void bed_interface_geom_cpu(const std::vector<double>& Ftilde, std::vector<double>& nbx, std::vector<double>& nby,
		std::vector<double>& nbz, std::vector<double>& Ab, MacGrid g);

	// ---- (2c) Reduce per-cell interface geometry to per-column β and up-slope direction -------
	// For each (x,y) column, pick the INTERFACE cell (0<F<1 nearest the bed top) and read its
	// normal: β = acos(clamp(nbz)) [rad]; up-slope horizontal unit = normalize(−nbx,−nby)
	// (fallback (1,0) if flat). Outputs beta_col, upx_col, upy_col (size nx·ny).
	void bed_column_geom_gpu(const double* fpack, const double* nbx, const double* nby, const double* nbz,
		double* beta_col, double* upx_col, double* upy_col, MacGrid g, double cpack);
	void bed_column_geom_cpu(const std::vector<double>& fpack, const std::vector<double>& nbx, const std::vector<double>& nby,
		const std::vector<double>& nbz, std::vector<double>& beta_col, std::vector<double>& upx_col, std::vector<double>& upy_col,
		MacGrid g, double cpack);

	// ---- (3a) Bedload volumetric flux vector q_b(x,y) on the columns -------------------------
	// From the per-column grain-skin bed shear τ_b=(taubx,tauby): θ=|τ|/((ρs−ρ)g d); slope-
	// corrected θ_cr (Soulsby, β/ψ from up-slope vs τ direction); Φ (Wong–Parker or Engelund–
	// Fredsøe), GATED by θ_cr,slope; q_b = Φ·√((s−1)g d³) aligned with τ (RESEARCH §6.3, §5).
	// Outputs qx,qy [m²/s] (size nx·ny). τ fields are per-column (2D).
	void bedload_flux_gpu(const double* taubx, const double* tauby, const double* beta_col,
		const double* upx_col, const double* upy_col, double* qx, double* qy, MacGrid g, MorphoParams p);
	void bedload_flux_cpu(const std::vector<double>& taubx, const std::vector<double>& tauby, const std::vector<double>& beta_col,
		const std::vector<double>& upx_col, const std::vector<double>& upy_col, std::vector<double>& qx, std::vector<double>& qy,
		MacGrid g, MorphoParams p);

	// ---- (3b) Upwind divergence of the bedload flux, along the transport direction -----------
	// ∇·q = (F_{i+1/2}−F_{i−1/2})/h + …, face flux F upwinded ALONG THE TRANSPORT DIRECTION
	// (dirx,diry — the flow/τ tangent; only its sign is used): F_{i+1/2} = q_i if transport is +x
	// at that face else q_{i+1} (RESEARCH §7 / research/03: "upwind ∇·q_b along the transport
	// direction"; central is unstable — sedExnerFoam). This is the physically-correct upwind for a
	// morphodynamic flux whose direction is set by the flow, and it is single-valued per face ⇒
	// conservative (telescoping). `periodic` wraps x/y (else zero flux across the domain walls).
	// Output div [m/s] (size nx·ny). For an aligned load q = q_b·t̂ (q_b≥0) the transport sign is
	// sign(q); pass q as dir to recover the Godunov flux-split, or a fixed current for a 1D Exner test.
	// `xopen` (default 0 = the legacy closed walls, so all existing callers stay byte-identical): when
	// set AND non-periodic, the streamwise (x) domain faces use a zero-gradient (open) bedload flux —
	// Fxl = q at the inlet, Fxr = q at the outlet — so the far-field bed feeds/drains its equilibrium
	// bedload instead of the wall's spurious inlet erosion / outlet pile. y-faces stay lateral walls.
	void bedload_div_gpu(const double* qx, const double* qy, const double* dirx, const double* diry, double* div, MacGrid g, int periodic, int xopen = 0);
	void bedload_div_cpu(const std::vector<double>& qx, const std::vector<double>& qy, const std::vector<double>& dirx,
		const std::vector<double>& diry, std::vector<double>& div, MacGrid g, int periodic, int xopen = 0);

	// ---- (4) Exner column update: deposition + Winterwerp erosion + bedload divergence -------
	// Per column: D = w_s,hindered·c_b (NOT thresholded), E_v = u_lift·c_pack (Winterwerp, GATED),
	// bedload = −div_qb. dG = M·dt·(D − E_v − div_qb), limited so |Δz_b|≤dz_limit_frac·h; G += dG.
	// The deposition/erosion grain is exchanged with the suspended field c at the bed-adjacent
	// cell k_bed=clamp(⌊z_b/h⌋): dc[k_bed] = −(applied grain change from D,E)/h ⇒ V_bed+V_susp
	// conserved to machine precision (bedload is bed-internal). `clip_count` (size nx·ny, may be
	// null) is set to 1 where the limiter clipped. c may be null (clear-water: no exchange).
	// RESEARCH §5 (asymmetry: deposition never thresholded), §7 (Exner units: D volumetric, E via
	// u_lift·c_pack; MORFAC; per-step limiter).
	// `dep_out`/`ero_out` (size nx·ny, may be null — GUI diagnostic only) receive the ACTUAL applied
	// deposition and erosion grain thickness this step [m of solid], post-limiter/post-MORFAC. They are
	// pure WRITE-ONLY outputs — never read by the physics — so passing them is byte-identical to null
	// (existing gates unaffected). Net bed-exchange rate = (dep_out − ero_out)/(c_pack·M·dt).
	void morpho_exner_gpu(double* G, double* c, const double* taubx, const double* tauby, const double* beta_col,
		const double* upx_col, const double* upy_col, const double* div_qb, double* clip_count, MacGrid g, MorphoParams p, double dt,
		double* dep_out = nullptr, double* ero_out = nullptr);
	void morpho_exner_cpu(std::vector<double>& G, std::vector<double>* c, const std::vector<double>& taubx, const std::vector<double>& tauby,
		const std::vector<double>& beta_col, const std::vector<double>& upx_col, const std::vector<double>& upy_col,
		const std::vector<double>& div_qb, std::vector<double>* clip_count, MacGrid g, MorphoParams p, double dt,
		double* dep_out = nullptr, double* ero_out = nullptr);
}
