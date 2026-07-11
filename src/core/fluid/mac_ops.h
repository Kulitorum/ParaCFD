// mac_ops.h — host-callable launchers for the MAC-solver CUDA kernels and their
// serial CPU reference twins. Each kernel/reference pair is exercised by a
// GPU-vs-CPU gtest at rel. tol 1e-5 (CLAUDE.md hard rule). Fields are double
// precision throughout M1 (2M cells at 128^3 is tiny; double makes the Ghia <5%
// match and the 1e-6 divergence gate robust — RESEARCH §3, PLAN M1).
//
// Device functions take raw device pointers; CPU references take std::vector<double>
// of matching length. "cell field" = nx*ny*nz; face fields = u/v/w counts (mac_grid.h).
#pragma once

#include "core/fluid/mac_grid.h"

#include <vector>

namespace windcfd::core
{
	// ---- Advection (MacCormack semi-Lagrangian, RK2 backtrace) ------------------
	// Advect (uIn,vIn,wIn) by its own velocity into (uOut,vOut,wOut). scratch must
	// hold at least max(u_count,v_count,w_count) doubles. Reverts to 1st order within
	// `band` cells of any wall and clamps to the 8-corner min/max (RESEARCH §3.2).
	void mac_advect_gpu(const double* uIn, const double* vIn, const double* wIn,
		double* uOut, double* vOut, double* wOut, double* scratch,
		MacGrid g, BC bc, double dt, int band = 1);

	void mac_advect_cpu(const std::vector<double>& uIn, const std::vector<double>& vIn, const std::vector<double>& wIn,
		std::vector<double>& uOut, std::vector<double>& vOut, std::vector<double>& wOut,
		MacGrid g, BC bc, double dt, int band = 1);

	// ---- Smagorinsky eddy viscosity nu_t at cell centers -----------------------
	// nu_t = (Cs*h)^2 * |S|, |S| = sqrt(2 S_ij S_ij) (RESEARCH §3.4). Cs=0 -> zero.
	void smagorinsky_nut_gpu(const double* u, const double* v, const double* w,
		double* nut, MacGrid g, BC bc, double Cs);
	void smagorinsky_nut_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		std::vector<double>& nut, MacGrid g, BC bc, double Cs);

	// ---- Explicit diffusion  f_out = f_in + dt*(nu + nu_t)*Lap(f_in) ------------
	// nut is a cell field (may be null -> molecular only). RESEARCH §3.4.
	void mac_diffuse_gpu(const double* uIn, const double* vIn, const double* wIn,
		double* uOut, double* vOut, double* wOut, const double* nut,
		MacGrid g, BC bc, double dt, double nu);
	void mac_diffuse_cpu(const std::vector<double>& uIn, const std::vector<double>& vIn, const std::vector<double>& wIn,
		std::vector<double>& uOut, std::vector<double>& vOut, std::vector<double>& wOut, const std::vector<double>* nut,
		MacGrid g, BC bc, double dt, double nu);

	// ---- Projection pieces ------------------------------------------------------
	// rhs[c] = -(rho/dt) * div(u)[c]  (cell field). RESEARCH §3.3.
	void poisson_rhs_gpu(const double* u, const double* v, const double* w,
		double* rhs, MacGrid g, double rho, double dt);
	void poisson_rhs_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		std::vector<double>& rhs, MacGrid g, double rho, double dt);

	// A*p where A = -Laplacian (positive semi-definite, all-Neumann box): drop
	// out-of-box neighbours (dp/dn=0). Ap[c] = (cnt*p_c - sum p_nb)/h^2.
	void poisson_apply_gpu(const double* p, double* Ap, MacGrid g);
	void poisson_apply_cpu(const std::vector<double>& p, std::vector<double>& Ap, MacGrid g);

	// Damped-Jacobi smoothing sweep(s): p += omega * (rhs - A p)/diag, diag=cnt/h^2.
	void jacobi_smooth_gpu(double* p, const double* rhs, double* scratch, MacGrid g, double omega, int sweeps);
	void jacobi_smooth_cpu(std::vector<double>& p, const std::vector<double>& rhs, MacGrid g, double omega, int sweeps);

	// Red-black Gauss-Seidel over the boundary band (cells within `band` of the box
	// surface). `forward` selects colour order (for a symmetric V-cycle). RESEARCH §3.3.
	void gs_band_gpu(double* p, const double* rhs, MacGrid g, int band, int sweeps, bool forward);
	void gs_band_cpu(std::vector<double>& p, const std::vector<double>& rhs, MacGrid g, int band, int sweeps, bool forward);

	// residual r = rhs - A p (cell field).
	void poisson_residual_gpu(const double* p, const double* rhs, double* r, MacGrid g);
	void poisson_residual_cpu(const std::vector<double>& p, const std::vector<double>& rhs, std::vector<double>& r, MacGrid g);

	// Full-weighting restriction (fine -> coarse, factor-2). Coarse dims = fine/2.
	void restrict_gpu(const double* fine, double* coarse, MacGrid gf, MacGrid gc);
	void restrict_cpu(const std::vector<double>& fine, std::vector<double>& coarse, MacGrid gf, MacGrid gc);

	// Trilinear prolongation with accumulation: fine += P * coarse.
	void prolong_add_gpu(const double* coarse, double* fine, MacGrid gc, MacGrid gf);
	void prolong_add_cpu(const std::vector<double>& coarse, std::vector<double>& fine, MacGrid gc, MacGrid gf);

	// u -= (dt/rho) grad(p) on interior faces (wall-normal faces untouched). RESEARCH §3.3.
	void subtract_gradient_gpu(double* u, double* v, double* w, const double* p, MacGrid g, double rho, double dt);
	void subtract_gradient_cpu(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, const std::vector<double>& p, MacGrid g, double rho, double dt);

	// max |div u| over all cells [1/s] (divergence-check kernel, RESEARCH §3.3 gate).
	double max_abs_divergence_gpu(const double* u, const double* v, const double* w, double* divscratch, MacGrid g);
	double max_abs_divergence_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, MacGrid g);

	// ---- small device reductions (thrust-backed) --------------------------------
	double reduce_sum_gpu(const double* x, int n);
	double reduce_max_abs_gpu(const double* x, int n);
	double dot_gpu(const double* a, const double* b, int n);
	void axpy_gpu(double alpha, const double* x, double* y, int n); // y += alpha*x
	void scale_add_gpu(double* y, double alpha, const double* x, double beta, int n); // y = alpha*x + beta*y
	void add_const_gpu(double* x, double c, int n);
}
