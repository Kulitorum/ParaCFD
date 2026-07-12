// channel_ops.h — host-callable launchers for the M2 open-channel MAC kernels and
// their serial CPU reference twins (GPU-vs-CPU gtest at rel. max-norm 1e-5, per the
// CLAUDE.md hard rule). Everything here is mask + open-BC aware:
//   * MacCormack advection with near-solid/near-boundary 1st-order reversion and
//     backtrace clipping out of solids (RESEARCH §3.2).
//   * Smagorinsky ν_t and explicit diffusion with solid-tangential no-slip/free-slip.
//   * Inlet Dirichlet, Orlanski-type outlet, global flux rescale, solid u·n=0.
//   * Masked + Dirichlet-outlet Poisson operator pieces (RESEARCH §3.3, §8).
// "cell field" = nx*ny*nz; face fields = u/v/w counts (mac_grid.h). solid is a cell
// field of unsigned char (1=solid,0=fluid); nearsolid is its `band`-dilation.
#pragma once

#include "core/fluid/channel_bc.h"

#include <vector>

namespace windcfd::core
{
	// ---- Advection (masked MacCormack, RK2 backtrace) ---------------------------
	void ch_advect_gpu(const double* uIn, const double* vIn, const double* wIn,
		double* uOut, double* vOut, double* wOut, double* scratch,
		const unsigned char* solid, const unsigned char* nearsolid,
		MacGrid g, ChannelBC bc, double dt);
	void ch_advect_cpu(const std::vector<double>& uIn, const std::vector<double>& vIn, const std::vector<double>& wIn,
		std::vector<double>& uOut, std::vector<double>& vOut, std::vector<double>& wOut,
		const std::vector<unsigned char>& solid, const std::vector<unsigned char>& nearsolid,
		MacGrid g, ChannelBC bc, double dt);

	// ---- Smagorinsky ν_t at cell centres (solids read as 0 velocity) ------------
	void ch_smagorinsky_gpu(const double* u, const double* v, const double* w,
		double* nut, const unsigned char* solid, MacGrid g, ChannelBC bc, double Cs);
	void ch_smagorinsky_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		std::vector<double>& nut, const std::vector<unsigned char>& solid, MacGrid g, ChannelBC bc, double Cs);

	// ---- Explicit diffusion  f_out = f_in + dt*(ν+ν_t)*Lap(f_in) ----------------
	// Solid tangential faces use no-slip reflection or free-slip mirror per bc.solid_mode.
	void ch_diffuse_gpu(const double* uIn, const double* vIn, const double* wIn,
		double* uOut, double* vOut, double* wOut, const double* nut, const unsigned char* solid,
		MacGrid g, ChannelBC bc, double dt, double nu);
	void ch_diffuse_cpu(const std::vector<double>& uIn, const std::vector<double>& vIn, const std::vector<double>& wIn,
		std::vector<double>& uOut, std::vector<double>& vOut, std::vector<double>& wOut,
		const std::vector<double>* nut, const std::vector<unsigned char>& solid,
		MacGrid g, ChannelBC bc, double dt, double nu);

	// ---- Boundary application ---------------------------------------------------
	// Inlet: set u-face plane i=0 to the Dirichlet profile; v,w at i=0 plane left as is
	// (cross-flow handled by fetch). Also zeroes v,w on the inlet cell column for cleanliness.
	void ch_apply_inlet_gpu(double* u, MacGrid g, ChannelBC bc);
	// Orlanski convective outlet on the i=nx u-plane: u_b += -(Uc*dt/dx_out)*(u[nx]-u[nx-1]), where
	// dx_out is the LOCAL outlet cell width g.dx(nx-1) (graded-correct; == h on a uniform grid). The
	// coefficient is formed inside the kernel so the device metric arrays are dereferenced on-device.
	void ch_orlanski_gpu(double* u, MacGrid g, ChannelBC bc, double dt);
	// Zero the interior solid faces (u·n=0 on obstacle surfaces & inside).
	void ch_apply_solid_bc_gpu(double* u, double* v, double* w, const unsigned char* solid, MacGrid g);
	// Volumetric flux [m^3/s] through the x-face plane i_plane: Σ u·dy(j)·dz(k) (true face areas, so it is
	// correct on a graded grid; == Σu·h² on a uniform grid).
	double ch_uplane_flux_gpu(const double* u, double* planeScratch, MacGrid g, int i_plane);
	// Rescale the outlet u-plane (i=nx) by factor s.
	void ch_scale_uplane_gpu(double* u, MacGrid g, int i_plane, double s);

	// ---- Masked + Dirichlet-outlet Poisson pieces (dir_xmax=1 => p=0 at i=nx) ----
	// rhs[c] = -(rho/dt) div(u)[c] for fluid cells, 0 for solid cells.
	void ch_poisson_rhs_gpu(const double* u, const double* v, const double* w,
		double* rhs, const unsigned char* solid, MacGrid g, double rho, double dt);
	void ch_poisson_rhs_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		std::vector<double>& rhs, const std::vector<unsigned char>& solid, MacGrid g, double rho, double dt);
	// Ap = (diag*p - Σ fluid-nb p)/h^2 ; diag = #fluid-nb (+1 if dir_xmax and i==nx-1).
	void ch_poisson_apply_gpu(const double* p, double* Ap, const unsigned char* solid, MacGrid g, int dir_xmax);
	void ch_poisson_apply_cpu(const std::vector<double>& p, std::vector<double>& Ap, const std::vector<unsigned char>& solid, MacGrid g, int dir_xmax);
	// residual r = rhs - A p.
	void ch_poisson_residual_gpu(const double* p, const double* rhs, double* r, const unsigned char* solid, MacGrid g, int dir_xmax);
	// damped-Jacobi sweeps (solids skipped): p += ω (rhs - A p)/diag.
	void ch_jacobi_gpu(double* p, const double* rhs, double* scratch, const unsigned char* solid, MacGrid g, int dir_xmax, double omega, int sweeps);
	void ch_jacobi_cpu(std::vector<double>& p, const std::vector<double>& rhs, const std::vector<unsigned char>& solid, MacGrid g, int dir_xmax, double omega, int sweeps);
	// red-black Gauss-Seidel over the boundary band (solids skipped).
	void ch_gs_band_gpu(double* p, const double* rhs, const unsigned char* solid, MacGrid g, int dir_xmax, int band, int sweeps, bool forward);
	// u -= (dt/rho) grad p on interior fluid faces; inlet face skipped, outlet uses p_ghost=0.
	void ch_subtract_gradient_gpu(double* u, double* v, double* w, const double* p, const unsigned char* solid, MacGrid g, ChannelBC bc, double rho, double dt, int dir_xmax);
	void ch_subtract_gradient_cpu(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, const std::vector<double>& p, const std::vector<unsigned char>& solid, MacGrid g, ChannelBC bc, double rho, double dt, int dir_xmax);
	// max |div u| over fluid cells [1/s].
	double ch_max_div_gpu(const double* u, const double* v, const double* w, double* divscratch, const unsigned char* solid, MacGrid g);
	double ch_max_div_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<unsigned char>& solid, MacGrid g);
}
