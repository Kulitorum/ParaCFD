// periodic_ops.h — testable kernels for the M3 periodic channel: the open-channel
// Prandtl mixing-length eddy viscosity and the wall-normal (z) explicit diffusion with
// free-slip top/bottom (the bed shear is injected separately as a wall-model sink). Each
// has a serial CPU twin (GPU-vs-CPU gtest, rel. max-norm 1e-5). RESEARCH §4, research/04.
#pragma once

#include "core/fluid/mac_grid.h"

#include <vector>

namespace paracfd::core
{
	// Mixing-length eddy viscosity at cell centres: ν_t = l²·|dU/dz|,
	// l = κ·z·√(1−z/Lz), z = (k+0.5)·h, |dU/dz| = |d/dz √(u_c²+v_c²)| via the horizontal
	// components (free-slip mirror in z). This algebraic closure makes the body-forced
	// steady profile an exact log law (research/04); Smagorinsky Δ=h cannot fill 5 m.
	void periodic_mixlen_gpu(const double* u, const double* v, double* nut, MacGrid g, double kappa, double Lz);
	void periodic_mixlen_cpu(const std::vector<double>& u, const std::vector<double>& v,
		std::vector<double>& nut, MacGrid g, double kappa, double Lz);

	// Explicit wall-normal diffusion of u,v: f += dt/h²·[F_up − F_dn], F = (ν+ν_t_face)Δf,
	// with ZERO flux through the bed (k=0) and lid (k=nz) faces (free-slip; the bed shear
	// is the wall-model sink). ν_t_face = ½(ν_t[k]+ν_t[k+1]).
	void periodic_zdiffuse_gpu(const double* uIn, const double* vIn, double* uOut, double* vOut,
		const double* nut, MacGrid g, double dt, double nu);
	void periodic_zdiffuse_cpu(const std::vector<double>& uIn, const std::vector<double>& vIn,
		std::vector<double>& uOut, std::vector<double>& vOut, const std::vector<double>& nut,
		MacGrid g, double dt, double nu);
}
