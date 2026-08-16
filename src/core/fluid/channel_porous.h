// channel_porous.h — sub-grid thin-screen / porous momentum sink for perforated structures whose
// holes are below the ≥ 8-cell resolution floor (RESEARCH §8, research/13 §8, research/07). A
// perforated face imposes a pressure jump Δp = ½·ρ·k·u_n², k = 1/β² − 1 (β = open-area ratio,
// Steiros–Hultmark 2018 / Taylor–Davies). Over a screen of one cell thickness h that is a
// volumetric deceleration a = ½·k·u_n²/h on the through-screen (face-normal) velocity component.
//
// Applied as a linearized-implicit Forchheimer drag on the provisional (post-advection/diffusion)
// velocity, BEFORE the projection restores incompressibility:
//   u_face ← u_face / (1 + dt·½·k_face·|u_face|/h)
// k_face is taken as the max porous-k of the two cells sharing the face, so a face fully inside a
// porous wall gets its wall's resistance and a face straddling porous/fluid still resists. β = 1
// (fully open) ⇒ k = 0 ⇒ no-op; a null porous field ⇒ the launcher is never called ⇒ M2/M3/seabed
// byte-identical. Suspended sediment passes through with its settling intact (the sink is on the
// fluid only) — exactly the sub-grid perforation behaviour research/13 §8 prescribes.
//
// Pure __host__ __device__ inline shared by the kernel and its CPU twin (GPU-vs-CPU parity, hard rule).
#pragma once

#include "core/fluid/mac_grid.h"

#include <cmath>
#include <vector>

namespace paracfd::core
{
	// Linearized-implicit quadratic (Forchheimer) drag on one face-normal velocity component.
	// kface = screen resistance k = 1/β²−1 (≥0); h = cell size; returns the decelerated component.
	PARACFD_HD inline double porous_face_drag(double uc, double kface, double dt, double h)
	{
		if (kface <= 0.0 || h <= 0.0) return uc;
		double denom = 1.0 + dt * 0.5 * kface * fabs(uc) / h;
		return uc / denom;
	}

	// Apply the porous drag to all three face velocities in place. porous_k is a per-cell field
	// (size p_count, 0 = no screen). Device launcher; CPU twin `ch_porous_drag_cpu` mirrors it.
	void ch_porous_drag_gpu(double* u, double* v, double* w, const double* porous_k, MacGrid g, double dt);
	void ch_porous_drag_cpu(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w,
		const std::vector<double>& porous_k, MacGrid g, double dt);
}
