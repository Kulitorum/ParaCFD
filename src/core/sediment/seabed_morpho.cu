// seabed_morpho.cu — see seabed_morpho.h. The three seabed-scenario helper kernels + CPU twins.
// Each kernel shares its arithmetic with the CPU twin through the bedshear.h / mac_grid.h inlines
// so the GPU-vs-CPU parity test (tests/test_seabed.cu) certifies only launch/index/memory.
#include "core/sediment/seabed_morpho.h"

#include <cuda_runtime.h>

#include <cmath>
#include <vector>

namespace scour::core
{
	namespace
	{
		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }
		SCOUR_HD inline int pl(MacGrid g, int i, int j) { return j * g.nx + i; }

		// ---- (1) τ_b from the flow at z_p above the CURRENT bed top z_b -----------------------
		// Identical to fluid/bedshear.cu's bedshear_cell but the sample height is offset by z_b.
		SCOUR_HD inline void seabed_bedshear_cell(const double* u, const double* v, const double* G, MacGrid g,
			WallParams wp, double cpack, int i, int j, double& taux, double& tauy, double& us)
		{
			double ks = wall_ks(wp.d50);
			double zp = wall_zp(g.h, ks);                 // probe height ABOVE the bed
			double zb = (cpack > 1e-12) ? G[pl(g, i, j)] / cpack : 0.0;
			double uc, vc;
			bed_sample_uv(u, v, g, i, j, zb + zp, uc, vc); // sample z_p above the bed top
			double Up = sqrt(uc * uc + vc * vc);
			us = (Up > 1e-12) ? wall_ustar(Up, zp, wp.d50, wp.nu, wp.kappa, wp.regime, wp.cj_iters) : 0.0;
			double tau = wp.rho * us * us;
			if (Up > 1e-12) { taux = tau * uc / Up; tauy = tau * vc / Up; }
			else { taux = 0.0; tauy = 0.0; }
		}

		__global__ void k_bedshear(const double* u, const double* v, const double* G,
			double* taux, double* tauy, double* ustar, MacGrid g, WallParams wp, double cpack, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = t / g.nx;
			double tx, ty, us; seabed_bedshear_cell(u, v, G, g, wp, cpack, i, j, tx, ty, us);
			taux[t] = tx; tauy[t] = ty; ustar[t] = us;
		}

		// ---- (4) HSV shear multiplier: τ → m·τ, u* → √m·u* (per column) ----------------------
		SCOUR_HD inline void shear_mult_col(double* taux, double* tauy, double* ustar, const double* mult, int t)
		{
			double m = mult[t];
			taux[t] *= m; tauy[t] *= m;
			ustar[t] *= sqrt(m > 0.0 ? m : 0.0);
		}
		__global__ void k_shear_mult(double* taux, double* tauy, double* ustar, const double* mult, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			shear_mult_col(taux, tauy, ustar, mult, t);
		}

		// ---- (2) confine c to the fluid column ----------------------------------------------
		SCOUR_HD inline void confine_column(double* c, const double* G, const unsigned char* structure,
			MacGrid g, double cpack, int i, int j)
		{
			double zb = (cpack > 1e-12) ? G[pl(g, i, j)] / cpack : 0.0;
			int kbed = (int)floor(zb / g.h);
			if (kbed < 0) kbed = 0; if (kbed > g.nz - 1) kbed = g.nz - 1;
			// first cell ≥ kbed that is not structure (sand below kbed is handled by kbed itself).
			int kf = kbed;
			while (kf < g.nz && structure && structure[g.pidx(i, j, kf)]) ++kf;
			if (kf >= g.nz) kf = g.nz - 1; // fully blocked column (should not occur): keep last cell
			double sum = 0.0;
			for (int k = 0; k < kf; ++k) { int idx = g.pidx(i, j, k); sum += c[idx]; c[idx] = 0.0; }
			c[g.pidx(i, j, kf)] += sum;
		}
		__global__ void k_confine(double* c, const double* G, const unsigned char* structure, MacGrid g, double cpack, int ncol)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= ncol) return;
			int i = t % g.nx, j = t / g.nx;
			confine_column(c, G, structure, g, cpack, i, j);
		}

		// ---- (3) clamp + pin ----------------------------------------------------------------
		__global__ void k_bed_post(double* G, const double* G_fixed, const unsigned char* frozen, double Gmax, int ncol)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= ncol) return;
			double gg = G[t];
			if (gg < 0.0) gg = 0.0; else if (gg > Gmax) gg = Gmax;
			if (frozen && frozen[t]) gg = G_fixed[t];
			G[t] = gg;
		}
	}

	// ================= GPU launchers =====================================
	void seabed_bedshear_gpu(const double* u, const double* v, const double* G,
		double* taux, double* tauy, double* ustar, MacGrid g, WallParams wp, double cpack)
	{
		int n = g.nx * g.ny;
		k_bedshear<<<gsz(n), 256>>>(u, v, G, taux, tauy, ustar, g, wp, cpack, n);
	}

	void seabed_apply_shear_mult_gpu(double* taux, double* tauy, double* ustar, const double* mult, MacGrid g)
	{
		int n = g.nx * g.ny;
		k_shear_mult<<<gsz(n), 256>>>(taux, tauy, ustar, mult, n);
	}

	void seabed_confine_c_gpu(double* c, const double* G, const unsigned char* structure, MacGrid g, double cpack)
	{
		int ncol = g.nx * g.ny;
		k_confine<<<gsz(ncol), 256>>>(c, G, structure, g, cpack, ncol);
	}

	void seabed_bed_post_gpu(double* G, const double* G_fixed, const unsigned char* frozen, MacGrid g, double Gmax)
	{
		int ncol = g.nx * g.ny;
		k_bed_post<<<gsz(ncol), 256>>>(G, G_fixed, frozen, Gmax, ncol);
	}

	// ================= CPU reference twins ===============================
	void seabed_bedshear_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& G,
		std::vector<double>& taux, std::vector<double>& tauy, std::vector<double>& ustar, MacGrid g, WallParams wp, double cpack)
	{
		int n = g.nx * g.ny; taux.assign(n, 0.0); tauy.assign(n, 0.0); ustar.assign(n, 0.0);
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
		{
			double tx, ty, us; seabed_bedshear_cell(u.data(), v.data(), G.data(), g, wp, cpack, i, j, tx, ty, us);
			int p = pl(g, i, j); taux[p] = tx; tauy[p] = ty; ustar[p] = us;
		}
	}

	void seabed_apply_shear_mult_cpu(std::vector<double>& taux, std::vector<double>& tauy, std::vector<double>& ustar, const std::vector<double>& mult, MacGrid g)
	{
		int n = g.nx * g.ny;
		for (int t = 0; t < n; ++t) shear_mult_col(taux.data(), tauy.data(), ustar.data(), mult.data(), t);
	}

	void seabed_confine_c_cpu(std::vector<double>& c, const std::vector<double>& G, const std::vector<unsigned char>* structure, MacGrid g, double cpack)
	{
		const unsigned char* s = structure ? structure->data() : nullptr;
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			confine_column(c.data(), G.data(), s, g, cpack, i, j);
	}

	void seabed_bed_post_cpu(std::vector<double>& G, const std::vector<double>& G_fixed, const std::vector<unsigned char>* frozen, MacGrid g, double Gmax)
	{
		int ncol = g.nx * g.ny;
		for (int t = 0; t < ncol; ++t)
		{
			double gg = G[t];
			if (gg < 0.0) gg = 0.0; else if (gg > Gmax) gg = Gmax;
			if (frozen && (*frozen)[t]) gg = G_fixed[t];
			G[t] = gg;
		}
	}
}
