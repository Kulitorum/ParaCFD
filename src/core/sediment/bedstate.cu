// bedstate.cu — CUDA kernels + CPU reference twins for the M5 bed state + morphodynamics
// (RESEARCH §5, §7; research/13, /03, /01, /16). The physically-shared arithmetic lives in the
// anonymous-namespace SCOUR_HD inlines so the GPU kernel and its CPU twin are bit-for-bit the same
// code (GPU-vs-CPU parity at rel. max-norm 1e-5, CLAUDE.md). See bedstate.h for the API contract.
#include "core/sediment/bedstate.h"
#include "core/fluid/mac_grid.h"
#include "core/sediment/sed_physics.h"

#include <cuda_runtime.h>

#include <vector>

namespace scour::core
{
	namespace
	{
		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }

		SCOUR_HD inline int clampi(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }

		// ---- normalized fill F with the interface-geometry ghost rule ------------------------
		// k<0 → below the floor is PACKED (F=1); k≥nz → above the top is FLUID (F=0); x/y replicate.
		SCOUR_HD inline double fetch_F(const double* fpack, const double* phi_obs, MacGrid g, double cpack, int i, int j, int k)
		{
			if (k < 0) return 1.0;
			if (k > g.nz - 1) return 0.0;
			i = clampi(i, 0, g.nx - 1);
			j = clampi(j, 0, g.ny - 1);
			int idx = g.pidx(i, j, k);
			double po = phi_obs ? phi_obs[idx] : 0.0;
			return normalized_fill(fpack[idx], po);
		}

		// generic ghost fetch on a cell field with prescribed below/above fill (for ∇F̃).
		SCOUR_HD inline double fetch_ghost(const double* f, MacGrid g, int i, int j, int k, double below, double above)
		{
			if (k < 0) return below;
			if (k > g.nz - 1) return above;
			i = clampi(i, 0, g.nx - 1);
			j = clampi(j, 0, g.ny - 1);
			return f[g.pidx(i, j, k)];
		}

		SCOUR_HD inline void fpack_column_node(const double* G, double* fpack, MacGrid g, double cpack, int i, int j)
		{
			double zb = G[j * g.nx + i] / cpack; // bed elevation [m]
			for (int k = 0; k < g.nz; ++k)
			{
				double f = zb / g.h - (double)k; // fractional fill of cell k
				if (f < 0.0) f = 0.0; else if (f > 1.0) f = 1.0;
				fpack[g.pidx(i, j, k)] = cpack * f;
			}
		}

		SCOUR_HD inline double boxfilter_node(const double* fpack, const double* phi_obs, MacGrid g, double cpack, int i, int j, int k)
		{
			double s = 0.0;
			for (int dk = -1; dk <= 1; ++dk)
				for (int dj = -1; dj <= 1; ++dj)
					for (int di = -1; di <= 1; ++di)
						s += fetch_F(fpack, phi_obs, g, cpack, i + di, j + dj, k + dk);
			return s * (1.0 / 27.0);
		}

		SCOUR_HD inline void geom_node(const double* Ft, MacGrid g, int i, int j, int k,
			double& nbx, double& nby, double& nbz, double& Ab)
		{
			double inv2h = 1.0 / (2.0 * g.h);
			double gx = (fetch_ghost(Ft, g, i + 1, j, k, 1.0, 0.0) - fetch_ghost(Ft, g, i - 1, j, k, 1.0, 0.0)) * inv2h;
			double gy = (fetch_ghost(Ft, g, i, j + 1, k, 1.0, 0.0) - fetch_ghost(Ft, g, i, j - 1, k, 1.0, 0.0)) * inv2h;
			double gz = (fetch_ghost(Ft, g, i, j, k + 1, 1.0, 0.0) - fetch_ghost(Ft, g, i, j, k - 1, 1.0, 0.0)) * inv2h;
			double mag = sqrt(gx * gx + gy * gy + gz * gz);
			double Vcell = g.h * g.h * g.h;
			if (mag < 1e-12)
			{
				nbx = 0.0; nby = 0.0; nbz = 1.0; Ab = 0.0; return; // flat/degenerate → up-normal, no area
			}
			nbx = -gx / mag; nby = -gy / mag; nbz = -gz / mag; // n_b = −∇F̃/|∇F̃| (into the fluid)
			Ab = mag * Vcell;
		}

		// per-column: pick the interface cell (lowest cell with 0<F<1) and derive β, up-slope dir.
		SCOUR_HD inline void column_geom_node(const double* fpack, const double* nbx, const double* nby, const double* nbz,
			double* beta_col, double* upx_col, double* upy_col, MacGrid g, double cpack, int i, int j)
		{
			int kif = -1;
			for (int k = 0; k < g.nz; ++k)
			{
				double f = fpack[g.pidx(i, j, k)] / cpack; // = F (phi_obs=0 open bed)
				if (f > 1e-9 && f < 1.0 - 1e-9) { kif = k; break; }
			}
			int col = j * g.nx + i;
			if (kif < 0) { beta_col[col] = 0.0; upx_col[col] = 1.0; upy_col[col] = 0.0; return; }
			int idx = g.pidx(i, j, kif);
			double nz = nbz[idx]; if (nz > 1.0) nz = 1.0; else if (nz < -1.0) nz = -1.0;
			beta_col[col] = acos(nz);
			double hx = -nbx[idx], hy = -nby[idx]; // up-slope = horizontal proj of ∇z_b = −(n_bx,n_by)
			double hm = sqrt(hx * hx + hy * hy);
			if (hm < 1e-12) { upx_col[col] = 1.0; upy_col[col] = 0.0; }
			else { upx_col[col] = hx / hm; upy_col[col] = hy / hm; }
		}

		// per-column bedload flux vector from the grain-skin bed shear (RESEARCH §6.3, §5).
		SCOUR_HD inline void bedload_flux_node(const double* taubx, const double* tauby, const double* beta_col,
			const double* upx_col, const double* upy_col, double* qx, double* qy, MacGrid g, MorphoParams p, int i, int j)
		{
			int col = j * g.nx + i;
			double tx = taubx[col], ty = tauby[col];
			double tm = sqrt(tx * tx + ty * ty);
			if (tm < 1e-30) { qx[col] = 0.0; qy[col] = 0.0; return; }
			double gam = sed_submerged_gamma(p.rho, p.rho_s);
			double theta = tm / (gam * p.d50);
			double Ds = grain_Dstar(p.d50, p.rho, p.rho_s, p.nu);
			double tcr_flat = theta_cr_sw(Ds);
			// ψ = angle between flow (τ̂) and up-slope
			double thx = tx / tm, thy = ty / tm;
			double ux = upx_col[col], uy = upy_col[col];
			double cpsi = thx * ux + thy * uy;
			double spsi = thx * uy - thy * ux;
			double psi = atan2(spsi, cpsi);
			double tcr = theta_cr_slope(tcr_flat, beta_col[col], psi, p.phi_repose);
			if (theta <= tcr) { qx[col] = 0.0; qy[col] = 0.0; return; } // θ_cr gates bedload
			double phi = p.bedload_formula == 1 ? bedload_phi_ef(theta, tcr) : bedload_phi_wp(theta);
			double qb = bedload_qb(phi, p.d50, p.rho, p.rho_s);
			qx[col] = qb * thx; qy[col] = qb * thy;
		}

		SCOUR_HD inline double bedload_div_node(const double* qx, const double* qy, const double* dirx, const double* diry,
			MacGrid g, int periodic, int xopen, int i, int j)
		{
			// Upwind ALONG THE TRANSPORT DIRECTION (research/03): face flux = q of the upwind cell,
			// chosen by the sign of the transport (dir) at the face. Single-valued per face ⇒
			// conservative. Non-periodic x-faces: xopen ⇒ zero-gradient (open-sea bedload feed/drain,
			// ghost q = qC), else zero bedload flux (wall). y-faces are always lateral walls.
			int col = j * g.nx + i;
			double qC = qx[col], dC = dirx[col];
			// x-faces
			double Fxr, Fxl;
			if (i == g.nx - 1 && !periodic) Fxr = xopen ? qC : 0.0; // outlet: free export at the interior rate / wall
			else { int ip = (i + 1 > g.nx - 1) ? 0 : i + 1; double qR = qx[j * g.nx + ip], dR = dirx[j * g.nx + ip]; Fxr = (dC + dR >= 0.0) ? qC : qR; }
			if (i == 0 && !periodic) Fxl = xopen ? qC : 0.0; // inlet: equilibrium feed (zero-gradient) / wall
			else { int im = (i - 1 < 0) ? g.nx - 1 : i - 1; double qL = qx[j * g.nx + im], dL = dirx[j * g.nx + im]; Fxl = (dL + dC >= 0.0) ? qL : qC; }
			double div = (Fxr - Fxl) / g.h;
			if (g.ny > 1)
			{
				double qCy = qy[col], dCy = diry[col];
				double Fyr, Fyl;
				if (j == g.ny - 1 && !periodic) Fyr = 0.0;
				else { int jp = (j + 1 > g.ny - 1) ? 0 : j + 1; double qR = qy[jp * g.nx + i], dR = diry[jp * g.nx + i]; Fyr = (dCy + dR >= 0.0) ? qCy : qR; }
				if (j == 0 && !periodic) Fyl = 0.0;
				else { int jm = (j - 1 < 0) ? g.ny - 1 : j - 1; double qL = qy[jm * g.nx + i], dL = diry[jm * g.nx + i]; Fyl = (dL + dCy >= 0.0) ? qL : qCy; }
				div += (Fyr - Fyl) / g.h;
			}
			return div;
		}

		// per-column Exner update: deposition (not thresholded) + Winterwerp erosion (gated) +
		// bedload divergence, with MORFAC and the per-step |Δz_b| limiter. Suspended exchange is
		// the exact negative of the applied deposition/erosion grain change (mass conservative).
		SCOUR_HD inline void morpho_exner_node(double* G, double* c, const double* taubx, const double* tauby,
			const double* beta_col, const double* upx_col, const double* upy_col, const double* div_qb,
			double* clip_count, MacGrid g, MorphoParams p, double dt, int i, int j,
			double* dep_out = nullptr, double* ero_out = nullptr)
		{
			int col = j * g.nx + i;
			double Gc = G[col];
			double zb = Gc / p.cpack;
			int kbed = clampi((int)floor(zb / g.h), 0, g.nz - 1);
			// Burial: suspended sand the growing bed has fully covered (cells entirely below the bed
			// surface z_b) is reclassified as bed grain — conservative, and keeps the still-water
			// deposition rate = w_s·c0 as the bed rises past a cell (RESEARCH §7 Exner equivalence).
			if (c)
			{
				for (int kk = 0; kk < kbed; ++kk)
				{
					if ((double)(kk + 1) * g.h <= zb)
					{
						int bci = g.pidx(i, j, kk);
						Gc += c[bci] * g.h; c[bci] = 0.0;
					}
				}
			}
			int cidx = g.pidx(i, j, kbed);
			double cb = c ? c[cidx] : 0.0;
			double ws = p.ws0 * (p.hindered ? hindered_factor(cb) : 1.0);
			double M = p.morfac;

			// deposition (NOT thresholded, RESEARCH §5 asymmetry): grain thickness to bed.
			double dep_grain = 0.0;
			if (p.deposition_on && c)
			{
				dep_grain = M * dt * ws * cb;          // D = w_s·c_b [m/s] × M·dt
				double avail = cb * g.h;               // grain available in the bed cell
				if (dep_grain > avail) dep_grain = avail;
			}

			// Winterwerp erosion (θ_cr-gated): grain thickness from bed to suspension.
			double ero_grain = 0.0;
			if (p.erosion_on)
			{
				double tx = taubx ? taubx[col] : 0.0, ty = tauby ? tauby[col] : 0.0;
				double tm = sqrt(tx * tx + ty * ty);
				if (tm > 1e-30)
				{
					double gam = sed_submerged_gamma(p.rho, p.rho_s);
					double theta = tm / (gam * p.d50);
					double Ds = grain_Dstar(p.d50, p.rho, p.rho_s, p.nu);
					double tcr_flat = theta_cr_sw(Ds);
					double thx = tx / tm, thy = ty / tm;
					double ux = upx_col[col], uy = upy_col[col];
					double psi = atan2(thx * uy - thy * ux, thx * ux + thy * uy);
					double tcr = theta_cr_slope(tcr_flat, beta_col[col], psi, p.phi_repose);
					double ulift = winterwerp_ulift(theta, tcr, p.d50, p.rho, p.rho_s, p.nu, p.erosion_coeff);
					ero_grain = M * dt * ulift * p.cpack; // E_v = u_lift·c_pack [m/s] × M·dt
				}
			}
			double bedavail = Gc + dep_grain;
			if (ero_grain > bedavail) ero_grain = bedavail; // cannot erode more than present

			// bedload (bed-internal): net grain into the column = −∇·q_b.
			double bl_grain = 0.0;
			if (p.bedload_on && div_qb) bl_grain = -M * dt * div_qb[col];

			double dG_de = dep_grain - ero_grain; // suspended-coupled part
			double dG = dG_de + bl_grain;

			// per-step limiter: |Δz_b| ≤ dz_limit_frac·h (RESEARCH §7 / Delft3D DzMax).
			double lim = p.dz_limit_frac * g.h;
			double dz = dG / p.cpack;
			if (fabs(dz) > lim && fabs(dG) > 1e-300)
			{
				double s = lim * p.cpack / fabs(dG);
				dG_de *= s; bl_grain *= s; dep_grain *= s; ero_grain *= s; dG = dG_de + bl_grain;
				if (clip_count) clip_count[col] = 1.0;
			}
			double Gn = Gc + dG;
			if (Gn < 0.0) Gn = 0.0;
			G[col] = Gn;
			if (c) c[cidx] = cb - dG_de / g.h; // exact negative of the bed grain change
			// GUI diagnostic outputs: the ACTUAL applied suspended-exchange grain this step (post-limiter,
			// post-MORFAC). Delivery = dep_grain, pickup = ero_grain [m of solid]. Write-only ⇒ no physics.
			if (dep_out) dep_out[col] = dep_grain;
			if (ero_out) ero_out[col] = ero_grain;
		}

		// ---- kernels -------------------------------------------------------------------------
		__global__ void k_fpack_from_G(const double* G, double* fpack, MacGrid g, double cpack, int ncol)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= ncol) return;
			int i = t % g.nx, j = t / g.nx;
			fpack_column_node(G, fpack, g, cpack, i, j);
		}
		__global__ void k_boxfilter(const double* fpack, const double* phi_obs, double* Ft, MacGrid g, double cpack, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			Ft[t] = boxfilter_node(fpack, phi_obs, g, cpack, i, j, k);
		}
		__global__ void k_geom(const double* Ft, double* nbx, double* nby, double* nbz, double* Ab, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			geom_node(Ft, g, i, j, k, nbx[t], nby[t], nbz[t], Ab[t]);
		}
		__global__ void k_column_geom(const double* fpack, const double* nbx, const double* nby, const double* nbz,
			double* beta_col, double* upx_col, double* upy_col, MacGrid g, double cpack, int ncol)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= ncol) return;
			int i = t % g.nx, j = t / g.nx;
			column_geom_node(fpack, nbx, nby, nbz, beta_col, upx_col, upy_col, g, cpack, i, j);
		}
		__global__ void k_bedload_flux(const double* taubx, const double* tauby, const double* beta_col,
			const double* upx_col, const double* upy_col, double* qx, double* qy, MacGrid g, MorphoParams p, int ncol)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= ncol) return;
			int i = t % g.nx, j = t / g.nx;
			bedload_flux_node(taubx, tauby, beta_col, upx_col, upy_col, qx, qy, g, p, i, j);
		}
		__global__ void k_bedload_div(const double* qx, const double* qy, const double* dirx, const double* diry, double* div, MacGrid g, int periodic, int xopen, int ncol)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= ncol) return;
			int i = t % g.nx, j = t / g.nx;
			div[t] = bedload_div_node(qx, qy, dirx, diry, g, periodic, xopen, i, j);
		}
		__global__ void k_morpho_exner(double* G, double* c, const double* taubx, const double* tauby, const double* beta_col,
			const double* upx_col, const double* upy_col, const double* div_qb, double* clip_count, MacGrid g, MorphoParams p, double dt, int ncol,
			double* dep_out, double* ero_out)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= ncol) return;
			int i = t % g.nx, j = t / g.nx;
			morpho_exner_node(G, c, taubx, tauby, beta_col, upx_col, upy_col, div_qb, clip_count, g, p, dt, i, j, dep_out, ero_out);
		}
	}

	// ================= launchers =================================================
	void fpack_from_G_gpu(const double* G, double* fpack, MacGrid g, double cpack)
	{
		int ncol = g.nx * g.ny; k_fpack_from_G<<<gsz(ncol), 256>>>(G, fpack, g, cpack, ncol);
	}
	void bed_boxfilter_gpu(const double* fpack, const double* phi_obs, double* Ftilde, MacGrid g, double cpack)
	{
		int n = g.p_count(); k_boxfilter<<<gsz(n), 256>>>(fpack, phi_obs, Ftilde, g, cpack, n);
	}
	void bed_interface_geom_gpu(const double* Ftilde, double* nbx, double* nby, double* nbz, double* Ab, MacGrid g)
	{
		int n = g.p_count(); k_geom<<<gsz(n), 256>>>(Ftilde, nbx, nby, nbz, Ab, g, n);
	}
	void bed_column_geom_gpu(const double* fpack, const double* nbx, const double* nby, const double* nbz,
		double* beta_col, double* upx_col, double* upy_col, MacGrid g, double cpack)
	{
		int ncol = g.nx * g.ny; k_column_geom<<<gsz(ncol), 256>>>(fpack, nbx, nby, nbz, beta_col, upx_col, upy_col, g, cpack, ncol);
	}
	void bedload_flux_gpu(const double* taubx, const double* tauby, const double* beta_col,
		const double* upx_col, const double* upy_col, double* qx, double* qy, MacGrid g, MorphoParams p)
	{
		int ncol = g.nx * g.ny; k_bedload_flux<<<gsz(ncol), 256>>>(taubx, tauby, beta_col, upx_col, upy_col, qx, qy, g, p, ncol);
	}
	void bedload_div_gpu(const double* qx, const double* qy, const double* dirx, const double* diry, double* div, MacGrid g, int periodic, int xopen)
	{
		int ncol = g.nx * g.ny; k_bedload_div<<<gsz(ncol), 256>>>(qx, qy, dirx, diry, div, g, periodic, xopen, ncol);
	}
	void morpho_exner_gpu(double* G, double* c, const double* taubx, const double* tauby, const double* beta_col,
		const double* upx_col, const double* upy_col, const double* div_qb, double* clip_count, MacGrid g, MorphoParams p, double dt,
		double* dep_out, double* ero_out)
	{
		int ncol = g.nx * g.ny; k_morpho_exner<<<gsz(ncol), 256>>>(G, c, taubx, tauby, beta_col, upx_col, upy_col, div_qb, clip_count, g, p, dt, ncol, dep_out, ero_out);
	}

	// ================= CPU reference twins =======================================
	void fpack_from_G_cpu(const std::vector<double>& G, std::vector<double>& fpack, MacGrid g, double cpack)
	{
		fpack.assign(g.p_count(), 0.0);
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			fpack_column_node(G.data(), fpack.data(), g, cpack, i, j);
	}
	void bed_boxfilter_cpu(const std::vector<double>& fpack, const std::vector<double>* phi_obs, std::vector<double>& Ftilde, MacGrid g, double cpack)
	{
		Ftilde.assign(g.p_count(), 0.0);
		const double* po = phi_obs ? phi_obs->data() : nullptr;
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			Ftilde[g.pidx(i, j, k)] = boxfilter_node(fpack.data(), po, g, cpack, i, j, k);
	}
	void bed_interface_geom_cpu(const std::vector<double>& Ftilde, std::vector<double>& nbx, std::vector<double>& nby,
		std::vector<double>& nbz, std::vector<double>& Ab, MacGrid g)
	{
		int n = g.p_count(); nbx.assign(n, 0.0); nby.assign(n, 0.0); nbz.assign(n, 0.0); Ab.assign(n, 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
		{
			int idx = g.pidx(i, j, k);
			geom_node(Ftilde.data(), g, i, j, k, nbx[idx], nby[idx], nbz[idx], Ab[idx]);
		}
	}
	void bed_column_geom_cpu(const std::vector<double>& fpack, const std::vector<double>& nbx, const std::vector<double>& nby,
		const std::vector<double>& nbz, std::vector<double>& beta_col, std::vector<double>& upx_col, std::vector<double>& upy_col,
		MacGrid g, double cpack)
	{
		int ncol = g.nx * g.ny; beta_col.assign(ncol, 0.0); upx_col.assign(ncol, 0.0); upy_col.assign(ncol, 0.0);
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			column_geom_node(fpack.data(), nbx.data(), nby.data(), nbz.data(), beta_col.data(), upx_col.data(), upy_col.data(), g, cpack, i, j);
	}
	void bedload_flux_cpu(const std::vector<double>& taubx, const std::vector<double>& tauby, const std::vector<double>& beta_col,
		const std::vector<double>& upx_col, const std::vector<double>& upy_col, std::vector<double>& qx, std::vector<double>& qy,
		MacGrid g, MorphoParams p)
	{
		int ncol = g.nx * g.ny; qx.assign(ncol, 0.0); qy.assign(ncol, 0.0);
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			bedload_flux_node(taubx.data(), tauby.data(), beta_col.data(), upx_col.data(), upy_col.data(), qx.data(), qy.data(), g, p, i, j);
	}
	void bedload_div_cpu(const std::vector<double>& qx, const std::vector<double>& qy, const std::vector<double>& dirx,
		const std::vector<double>& diry, std::vector<double>& div, MacGrid g, int periodic, int xopen)
	{
		int ncol = g.nx * g.ny; div.assign(ncol, 0.0);
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			div[j * g.nx + i] = bedload_div_node(qx.data(), qy.data(), dirx.data(), diry.data(), g, periodic, xopen, i, j);
	}
	void morpho_exner_cpu(std::vector<double>& G, std::vector<double>* c, const std::vector<double>& taubx, const std::vector<double>& tauby,
		const std::vector<double>& beta_col, const std::vector<double>& upx_col, const std::vector<double>& upy_col,
		const std::vector<double>& div_qb, std::vector<double>* clip_count, MacGrid g, MorphoParams p, double dt,
		double* dep_out, double* ero_out)
	{
		double* cp = c ? c->data() : nullptr;
		double* cc = clip_count ? clip_count->data() : nullptr;
		const double* dq = div_qb.empty() ? nullptr : div_qb.data();
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			morpho_exner_node(G.data(), cp, taubx.data(), tauby.data(), beta_col.data(), upx_col.data(), upy_col.data(), dq, cc, g, p, dt, i, j, dep_out, ero_out);
	}
}
