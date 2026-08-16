// bedshear_ops.h — host-callable launchers for the bed-shear wall-function kernels and
// their serial CPU reference twins (GPU-vs-CPU gtest, rel. max-norm 1e-5). RESEARCH §4,
// research/04. The bed τ_b lives on the flat-bed plane (nx·ny, cell-centred, xy-index
// pl(i,j)=j·nx+i); the momentum-sink writes back onto the k=0 MAC faces.
#pragma once

#include "core/fluid/bedshear.h"

#include <vector>

namespace paracfd::core
{
	// Raw τ_b from the velocity field over the flat bed: probe (u,v) at z_p, invert the
	// log law for u* (regime per WallParams), τx=ρu*²·û, τy=ρu*²·v̂. Outputs plane fields
	// taux,tauy,ustar of length nx·ny. (research/04 §5, anti-stair-step: interpolated u*.)
	void bedshear_compute_gpu(const double* u, const double* v, double* taux, double* tauy, double* ustar,
		MacGrid g, WallParams wp);
	void bedshear_compute_cpu(const std::vector<double>& u, const std::vector<double>& v,
		std::vector<double>& taux, std::vector<double>& tauy, std::vector<double>& ustar,
		MacGrid g, WallParams wp);

	// 3×3 tangential box smoothing of a bed-plane field (Neumann edges). `passes` sweeps
	// using `scratch` (length nx·ny). Anti-stair-step mitigation (research/04 §mitigations).
	void bedshear_smooth_gpu(double* fx, double* fy, double* scratch, MacGrid g, int passes);
	void bedshear_smooth_cpu(std::vector<double>& fx, std::vector<double>& fy, MacGrid g, int passes);

	// Per-step rate-limit ([0.5×,2×] of previous magnitude) then EMA time-filter into the
	// persistent τ_ema plane: τ_ema ← (1−a)·τ_ema + a·τ_lim, a = dt/t_avg (research/04).
	void bedshear_filter_gpu(const double* taux, const double* tauy, double* emax, double* emay,
		MacGrid g, WallParams wp, double dt, int first);
	void bedshear_filter_cpu(const std::vector<double>& taux, const std::vector<double>& tauy,
		std::vector<double>& emax, std::vector<double>& emay, MacGrid g, WallParams wp, double dt, int first);

	// Apply the wall shear as a momentum sink on the bottom (k=0) MAC faces:
	// u -= (τx_face/(ρ·h))·dt, v -= (τy_face/(ρ·h))·dt, τ_face = xy-average of the plane
	// τ across the face's two cells (periodic in x/y). This is the equilibrium wall model
	// that balances the body force in the periodic-channel driver (research/04 LES note).
	void bedshear_apply_sink_gpu(double* u, double* v, const double* emax, const double* emay,
		MacGrid g, double rho, double dt);
	void bedshear_apply_sink_cpu(std::vector<double>& u, std::vector<double>& v,
		const std::vector<double>& emax, const std::vector<double>& emay, MacGrid g, double rho, double dt);
}
