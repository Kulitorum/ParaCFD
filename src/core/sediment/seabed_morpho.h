// seabed_morpho.h — the three thin device helpers the live erodible-seabed scenario needs on
// TOP of the gate-verified sediment kernels (suspended.*, bedstate.*, avalanche.*, bedshear.*),
// which are reused UNCHANGED. Everything physical is shared with the CPU twin through the
// __host__ __device__ inlines in bedshear.h / sed_physics.h (GPU-vs-CPU parity, CLAUDE.md rule).
//
// The scenario is a fixed 3-D grid with a DEEP sand reservoir: the bottom `sand_depth` metres are
// packed grain (f_pack from the per-column grain thickness G, z_b = G/c_pack), fluid above, and a
// rigid printed structure seated on the bed. The bed is SOLID to the flow; as z_b changes the flow
// mask is rebuilt (owner: the engine). These helpers bridge the elevated/moving bed to the kernels
// that assume the bed sits at k=0:
//   (1) seabed_bedshear — τ_b sampled at z_p ABOVE the current bed top z_b (not the domain floor).
//   (2) seabed_confine_c — keep suspended sand in the fluid column (no leak into the reservoir).
//   (3) seabed_bed_post — clamp G ∈ [0,Gmax] and pin the rigid structure-footprint columns.
// Units: SI (m, s, Pa). τ plane fields are nx·ny (pl(i,j)=j·nx+i); c/f_pack are p_count.
#pragma once

#include "core/fluid/bedshear.h" // WallParams + wall-function inlines
#include "core/fluid/mac_grid.h"

#include <vector>

namespace scour::core
{
	// ---- (1) Per-column grain-skin bed shear τ_b from the LIVE flow -----------------------------
	// For each (x,y) column: z_b = G/c_pack; probe the tangential (u,v) speed at height z_b + z_p
	// (z_p = wall_zp(h,ks)); invert the log law (bedshear.h) for u* referenced to z_p above the bed;
	// τx = ρ u*² û, τy = ρ u*² v̂. Outputs plane fields taux,tauy,ustar (size nx·ny). Reuses the
	// bedshear.h inlines VERBATIM (only the sample height is offset by z_b). RESEARCH §4, research/04.
	void seabed_bedshear_gpu(const double* u, const double* v, const double* G,
		double* taux, double* tauy, double* ustar, MacGrid g, WallParams wp, double cpack);
	void seabed_bedshear_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& G,
		std::vector<double>& taux, std::vector<double>& tauy, std::vector<double>& ustar, MacGrid g, WallParams wp, double cpack);

	// ---- (2) Confine the suspended field c to the fluid region above the bed --------------------
	// The gate suspended kernels assume the bed is at k=0; here it is at k_bed = ⌊z_b/h⌋ with a deep
	// reservoir below. Settling would otherwise flux eroded/settled sand DOWN into the reservoir.
	// Per column: find k_fluid = first cell ≥ k_bed that is NOT structure (`structure` mask, size
	// p_count; sand below k_bed is identified from G, not the mask); sum c over [0,k_fluid) and move
	// it into c[k_fluid] (mass-conserving), zeroing the below/blocked cells. `structure` may be null
	// (⇒ k_fluid = k_bed everywhere). One thread per column.
	void seabed_confine_c_gpu(double* c, const double* G, const unsigned char* structure, MacGrid g, double cpack);
	void seabed_confine_c_cpu(std::vector<double>& c, const std::vector<double>& G, const std::vector<unsigned char>* structure, MacGrid g, double cpack);

	// ---- (3) Bed post-process: clamp + pin the rigid structure footprint -----------------------
	// Clamp every column's G to [0, Gmax] (no scour through bedrock; bounded deposition), then pin
	// G = G_fixed on frozen columns (the structure footprint — rigid, must not erode or accrete).
	// `frozen`/`G_fixed` size nx·ny; `frozen` may be null (⇒ clamp only). One thread per column.
	void seabed_bed_post_gpu(double* G, const double* G_fixed, const unsigned char* frozen, MacGrid g, double Gmax);
	void seabed_bed_post_cpu(std::vector<double>& G, const std::vector<double>& G_fixed, const std::vector<unsigned char>* frozen, MacGrid g, double Gmax);

	// ---- (4) HSV shear multiplier (M6 supervised calibration crutch, RESEARCH §11 #1) -----------
	// Scale the per-column bed shear by a spatial multiplier field `mult` (size nx·ny): τ → mult·τ
	// (taux,tauy) and u* → √mult·u* (τ = ρ u*²). The diffusive SL solver under-resolves the horseshoe
	// vortex and underpredicts its bed shear (~30% with body-fitted RANS; far more here — see gate_M6);
	// this is the sanctioned OPTIONAL knob to restore the near-pile amplification. The VALUE and REGION
	// are the orchestrator's + MH's calibration choice, never the model's default (mult ≡ 1 unless a
	// field is set). Applied to the raw τ before smoothing/EMA.
	void seabed_apply_shear_mult_gpu(double* taux, double* tauy, double* ustar, const double* mult, MacGrid g);
	void seabed_apply_shear_mult_cpu(std::vector<double>& taux, std::vector<double>& tauy, std::vector<double>& ustar, const std::vector<double>& mult, MacGrid g);
}
