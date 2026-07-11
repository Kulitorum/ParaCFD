// avalanche.h — M5 layered-column avalanching / sand-slide (RESEARCH §7; research/03 §8, /13 §7).
//
// Submerged angle of repose φ_s ≈ 32° (sand). Roulund et al. (2005) hysteresis: a slide ACTIVATES
// where the local bed slope exceeds 32° and RELAXES the slope down to 30° (32/30° hysteresis). On a
// voxel grid this is the standard mass-conserving 8-neighbour relaxation on the bed-elevation field
// z_b (= G/c_pack): for each ordered column pair (a→b) with (z_a−z_b) > L_ab·tan(θ_repose),
// transfer m = relax·((z_a−z_b) − L_ab·tan(θ_repose)) of elevation from a to b (antisymmetric ⇒
// grain mass exactly conserved). L_ab = h (orthogonal) or h·√2 (diagonal). Iterate Jacobi sweeps
// until every slope ≤ tan(θ_repose). Open-bed z_b here is the single-segment reduction of the
// general layered-column (multi-segment per (x,y)) representation of RESEARCH §7 / research/13 §7.
//
// GPU kernel + CPU reference twin, GPU-vs-CPU parity at rel. max-norm 1e-5 (CLAUDE.md).
#pragma once

#include "core/fluid/mac_grid.h"

#include <vector>

namespace scour::core
{
	struct AvalancheParams
	{
		double tan_trigger = 0.62487; // tan(32°) — activation threshold (hysteresis, driver-level)
		double tan_repose = 0.57735;  // tan(30°) — relax/stop angle (the sweep's target)
		double relax = 0.1;           // Jacobi under-relaxation r (8·r<1 for stability)
		int periodic = 0;             // wrap x/y neighbours (else domain-boundary walls block transfer)
	};

	// One Jacobi avalanche sweep on the bed-elevation field z (size nx·ny): z_out = z + Δz, Δz the
	// net antisymmetric pairwise transfer over the 8 neighbours (mass-conserving to machine
	// precision). z_in and z_out must not alias.
	void avalanche_sweep_gpu(const double* z_in, double* z_out, MacGrid g, AvalancheParams ap);
	void avalanche_sweep_cpu(const std::vector<double>& z_in, std::vector<double>& z_out, MacGrid g, AvalancheParams ap);

	// Maximum bed slope tan(β) over all 8-neighbour column pairs (for the hysteresis
	// activation test and convergence check). Host-side reduction over the field.
	double avalanche_max_tanslope(const std::vector<double>& z, MacGrid g, int periodic);
}
