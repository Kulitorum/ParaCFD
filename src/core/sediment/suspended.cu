// suspended.cu — CUDA kernels + CPU reference twins for suspended-sediment transport
// (RESEARCH §6.1–6.2, research/02). The physically-shared arithmetic lives in the anonymous-
// namespace SCOUR_HD inlines so the GPU kernel and the CPU twin are bit-for-bit the same code
// (GPU-vs-CPU parity at rel. max-norm 1e-5, CLAUDE.md hard rule). See suspended.h for the API.
#include "core/sediment/suspended.h"
#include "core/fluid/mac_grid.h"
#include "core/sediment/sed_physics.h"

#include <cuda_runtime.h>

#include <vector>

namespace scour::core
{
	namespace
	{
		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }

		// ---- cell-centered scalar sampling with ScalarBC ghosts ------------------------------
		SCOUR_HD inline double fetch_scalar(const double* c, MacGrid g, ScalarBC bc, int i, int j, int k)
		{
			bool dir = false; double dval = 0.0;
			if (i < 0) { if (bc.xmin == SC_DIRICHLET) { dir = true; dval = bc.v_xmin; } i = 0; }
			else if (i > g.nx - 1) { if (bc.xmax == SC_DIRICHLET) { dir = true; dval = bc.v_xmax; } i = g.nx - 1; }
			if (j < 0) { if (bc.ymin == SC_DIRICHLET) { dir = true; dval = bc.v_ymin; } j = 0; }
			else if (j > g.ny - 1) { if (bc.ymax == SC_DIRICHLET) { dir = true; dval = bc.v_ymax; } j = g.ny - 1; }
			if (k < 0) { if (bc.zmin == SC_DIRICHLET) { dir = true; dval = bc.v_zmin; } k = 0; }
			else if (k > g.nz - 1) { if (bc.zmax == SC_DIRICHLET) { dir = true; dval = bc.v_zmax; } k = g.nz - 1; }
			if (dir) return dval;
			return c[g.pidx(i, j, k)];
		}

		SCOUR_HD inline double trilerp_scalar(const double* c, MacGrid g, ScalarBC bc, double x, double y, double z, double* mn, double* mx)
		{
			double gx = x / g.h - 0.5, gy = y / g.h - 0.5, gz = z / g.h - 0.5;
			int i0 = (int)floor(gx), j0 = (int)floor(gy), k0 = (int)floor(gz);
			double fx = gx - i0, fy = gy - j0, fz = gz - k0;
			double cc[2][2][2];
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int d = 0; d < 2; ++d)
				cc[a][b][d] = fetch_scalar(c, g, bc, i0 + a, j0 + b, k0 + d);
			if (mn && mx)
			{
				double lo = cc[0][0][0], hi = cc[0][0][0];
				for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int d = 0; d < 2; ++d)
				{ lo = cc[a][b][d] < lo ? cc[a][b][d] : lo; hi = cc[a][b][d] > hi ? cc[a][b][d] : hi; }
				*mn = lo; *mx = hi;
			}
			double c00 = cc[0][0][0] * (1 - fx) + cc[1][0][0] * fx;
			double c10 = cc[0][1][0] * (1 - fx) + cc[1][1][0] * fx;
			double c01 = cc[0][0][1] * (1 - fx) + cc[1][0][1] * fx;
			double c11 = cc[0][1][1] * (1 - fx) + cc[1][1][1] * fx;
			double c0 = c00 * (1 - fy) + c10 * fy;
			double c1 = c01 * (1 - fy) + c11 * fy;
			return c0 * (1 - fz) + c1 * fz;
		}

		SCOUR_HD inline ScalarBC neumann_bc()
		{
			ScalarBC b; // all SC_NEUMANN by default
			return b;
		}

		// ---- advection velocity for c: (u,v,w_fluid − w_s) with settling folded in -----------
		SCOUR_HD inline void vel_at_c(const double* u, const double* v, const double* w, const double* ws,
			MacGrid g, BC vbc, double x, double y, double z, double& vx, double& vy, double& vz)
		{
			vx = trilerp_u(u, g, vbc, x, y, z, nullptr, nullptr);
			vy = trilerp_v(v, g, vbc, x, y, z, nullptr, nullptr);
			double wf = trilerp_w(w, g, vbc, x, y, z, nullptr, nullptr);
			ScalarBC wsbc = neumann_bc();
			double wsv = trilerp_scalar(ws, g, wsbc, x, y, z, nullptr, nullptr);
			vz = wf - wsv; // settling ⇒ extra downward advection (RESEARCH §6.2)
		}

		SCOUR_HD inline void rk2_trace_c(const double* u, const double* v, const double* w, const double* ws,
			MacGrid g, BC vbc, double dt, double x, double y, double z, double& xo, double& yo, double& zo)
		{
			double vx, vy, vz;
			vel_at_c(u, v, w, ws, g, vbc, x, y, z, vx, vy, vz);
			double xm = x - 0.5 * dt * vx, ym = y - 0.5 * dt * vy, zm = z - 0.5 * dt * vz;
			clamp_to_domain(g, xm, ym, zm);
			double vmx, vmy, vmz;
			vel_at_c(u, v, w, ws, g, vbc, xm, ym, zm, vmx, vmy, vmz);
			xo = x - dt * vmx; yo = y - dt * vmy; zo = z - dt * vmz;
			clamp_to_domain(g, xo, yo, zo);
		}

		SCOUR_HD inline void cell_pos(MacGrid g, int i, int j, int k, double& x, double& y, double& z)
		{ x = (i + 0.5) * g.h; y = (j + 0.5) * g.h; z = (k + 0.5) * g.h; }

		SCOUR_HD inline bool near_wall_z(MacGrid g, int i, int j, int k, int band)
		{
			// Column tests are x/y-homogeneous; reversion is meaningful only in z. Guard by the
			// wall-normal (z) proximity so a thin x/y strip does not force 1st order everywhere.
			return k < band || k > g.nz - 1 - band;
		}

		SCOUR_HD inline void forward_c_node(const double* cIn, const double* u, const double* v, const double* w,
			const double* ws, double* phiHat, MacGrid g, BC vbc, ScalarBC sbc, double dt, int i, int j, int k)
		{
			int idx = g.pidx(i, j, k);
			double x, y, z; cell_pos(g, i, j, k, x, y, z);
			double xb, yb, zb; rk2_trace_c(u, v, w, ws, g, vbc, dt, x, y, z, xb, yb, zb);
			phiHat[idx] = trilerp_scalar(cIn, g, sbc, xb, yb, zb, nullptr, nullptr);
		}

		SCOUR_HD inline void correct_c_node(const double* cIn, const double* phiHat, const double* u, const double* v,
			const double* w, const double* ws, double* out, MacGrid g, BC vbc, ScalarBC sbc, double dt, int band, int i, int j, int k)
		{
			int idx = g.pidx(i, j, k);
			double phf = phiHat[idx];
			if (band > 0 && near_wall_z(g, i, j, k, band)) { out[idx] = phf; return; } // 1st-order reversion
			double x, y, z; cell_pos(g, i, j, k, x, y, z);
			double xb, yb, zb; rk2_trace_c(u, v, w, ws, g, vbc, dt, x, y, z, xb, yb, zb);
			double mn, mx;
			trilerp_scalar(cIn, g, sbc, xb, yb, zb, &mn, &mx); // 8-corner bounds
			double xf, yf, zf; rk2_trace_c(u, v, w, ws, g, vbc, -dt, x, y, z, xf, yf, zf);
			double phn = trilerp_scalar(phiHat, g, sbc, xf, yf, zf, nullptr, nullptr);
			double corrected = phf + 0.5 * (cIn[idx] - phn);
			if (corrected < mn) corrected = mn; else if (corrected > mx) corrected = mx; // clamp
			out[idx] = corrected;
		}

		// ---- explicit variable-coefficient z/x/y diffusion of c -------------------------------
		SCOUR_HD inline double face_D(const double* Dc, MacGrid g, int i, int j, int k, int ni, int nj, int nk)
		{
			// arithmetic mean of the two straddling cell diffusivities; the diffusivity field
			// itself always uses a zero-gradient (Neumann) ghost at the domain faces.
			double d0 = Dc[g.pidx(i, j, k)];
			double d1 = fetch_scalar(Dc, g, neumann_bc(), ni, nj, nk);
			return 0.5 * (d0 + d1);
		}

		SCOUR_HD inline double diffuse_c_node(const double* cIn, const double* Dc, MacGrid g, ScalarBC bc, double dt, int i, int j, int k)
		{
			double c0 = cIn[g.pidx(i, j, k)];
			double inv_h2 = 1.0 / (g.h * g.h);
			double flux = 0.0;
			// six faces; flux = D_face·(c_nb − c0). Ghost per ScalarBC for c; Dc mirrors (Neumann).
			int di[6] = { +1,-1, 0, 0, 0, 0 }, dj[6] = { 0, 0,+1,-1, 0, 0 }, dk[6] = { 0, 0, 0, 0,+1,-1 };
			for (int f = 0; f < 6; ++f)
			{
				int in = i + di[f], jn = j + dj[f], kn = k + dk[f];
				double cn = fetch_scalar(cIn, g, bc, in, jn, kn);
				double Df = face_D(Dc, g, i, j, k, in, jn, kn);
				flux += Df * (cn - c0);
			}
			return c0 + dt * inv_h2 * flux;
		}

		// ---- flux-form bed exchange at k=0 (RESEARCH §6.2) ------------------------------------
		SCOUR_HD inline void bed_exchange_node(double* c, const double* tau_prime, double* bedflux, MacGrid g, SedBedParams p, double dt, int i, int j)
		{
			int idx = g.pidx(i, j, 0);
			double cb = c[idx];
			double ws = p.ws0 * (p.hindered ? hindered_factor(cb) : 1.0);
			double tau = tau_prime ? tau_prime[g.pidx(i, j, 0)] : 0.0;
			double E = vanrijn_pickup(tau, p.d50, p.rho, p.rho_s, p.nu, p.alpha); // kg/m²/s, θ_cr-gated
			double Fbed = E / p.rho_s - ws * cb;         // net upward volumetric flux [m/s]
			double dc = (dt / g.h) * Fbed;               // source in the bed cell (flux per cell height)
			double cnew = cb + dc;
			if (cnew < 0.0) cnew = 0.0;                  // never negative concentration
			c[idx] = cnew;
			if (bedflux) bedflux[j * g.nx + i] = (cnew - cb) * g.h; // ACTUAL applied Δ(c·h) per column [m]
		}

		// ---- open-x boundary advective flux at cells i=0 (xmin face) and i=nx-1 (xmax face) ----
		// One (j,k) column-face pair per call. idh = dt/h; the two omitted faces of advx_node are:
		//   xmin face of cell 0  : u-face uidx(0,j,k)     ⇒ cell 0    gains  +idh·F_xmin
		//   xmax face of cell nx-1: u-face uidx(nx,j,k)   ⇒ cell nx-1 gains  −idh·F_xmax
		// First-order upwind face value: inflow ⇒ ghost, outflow ⇒ interior. Signed net domain sand
		// [m³] per boundary cell (post-clamp, so it matches the ACTUAL applied Δc) → net_flux[t] (xmin)
		// and net_flux[ny·nz + t] (xmax), t = k·ny + j. c clamped ≥ 0 (CFL keeps it non-binding).
		SCOUR_HD inline void boundary_flux_node(double* c, const double* u, const double* gmin, const double* gmax,
			double* net_flux, MacGrid g, double dt, int j, int k)
		{
			double idh = dt / g.h;
			double h3 = g.h * g.h * g.h;
			int t = k * g.ny + j; // yz-plane index (matches ghost / net_flux layout)
			int nyz = g.ny * g.nz;
			// xmin face (left face of cell i=0)
			{
				double uf = u[g.uidx(0, j, k)];
				int idx0 = g.pidx(0, j, k);
				double cface = uf >= 0.0 ? gmin[t] : c[idx0]; // inflow: ghost; outflow: interior (upwind)
				double cold = c[idx0];
				double cnew = cold + idh * (uf * cface);      // += idh·F_xmin
				if (cnew < 0.0) cnew = 0.0;
				c[idx0] = cnew;
				if (net_flux) net_flux[t] = (cnew - cold) * h3;
			}
			// xmax face (right face of cell i=nx-1)
			{
				double uf = u[g.uidx(g.nx, j, k)];
				int idxN = g.pidx(g.nx - 1, j, k);
				double cface = uf >= 0.0 ? c[idxN] : gmax[t]; // outflow: interior; inflow (reversal): ghost
				double cold = c[idxN];
				double cnew = cold - idh * (uf * cface);      // −= idh·F_xmax
				if (cnew < 0.0) cnew = 0.0;
				c[idxN] = cnew;
				if (net_flux) net_flux[nyz + t] = (cnew - cold) * h3;
			}
		}

		SCOUR_HD inline double eff_ws_node(double c, double ws0, int hindered)
		{
			return ws0 * (hindered ? hindered_factor(c) : 1.0);
		}

		// ---- conservative TVD van-Leer flux-form advection (mass-conserving) -----------------
		SCOUR_HD inline double vl_phi(double a, double b) // van Leer limiter φ(r), r=a/b
		{
			if (fabs(b) < 1e-300) return 0.0;
			double r = a / b;
			return (r + fabs(r)) / (1.0 + fabs(r));
		}
		// flux v·c_face with a Sweby/van-Leer limited face value from the 4-cell stencil. The
		// (1−|ν|) factor (ν = v·dt/h = face Courant number) is the flux-limited Lax–Wendroff blend
		// that keeps the single-stage forward-Euler update TVD/stable up to Courant 1 (a plain MUSCL
		// reconstruction is stable only to ~0.5). Smooth regions ⇒ φ→1 ⇒ 2nd-order Lax–Wendroff.
		SCOUR_HD inline double vanleer_flux(double cLL, double cL, double cR, double cRR, double vf, double nu_face)
		{
			double corr = 1.0 - fabs(nu_face); if (corr < 0.0) corr = 0.0;
			if (vf >= 0.0) { double b = cR - cL; double phi = vl_phi(cL - cLL, b); return vf * (cL + 0.5 * phi * corr * b); }
			double b = cL - cR; double phi = vl_phi(cR - cRR, b); return vf * (cR + 0.5 * phi * corr * b);
		}
		SCOUR_HD inline double cget(const double* c, MacGrid g, int i, int j, int k) // zero-gradient clamp
		{
			if (i < 0) i = 0; else if (i > g.nx - 1) i = g.nx - 1;
			if (j < 0) j = 0; else if (j > g.ny - 1) j = g.ny - 1;
			if (k < 0) k = 0; else if (k > g.nz - 1) k = g.nz - 1;
			return c[g.pidx(i, j, k)];
		}
		SCOUR_HD inline double advx_node(const double* c, const double* u, MacGrid g, double dt, int i, int j, int k)
		{
			double idh = dt / g.h;
			double uf, Fl = 0.0, Fr = 0.0;
			if (i != 0) { uf = u[g.uidx(i, j, k)]; Fl = vanleer_flux(cget(c, g, i - 2, j, k), cget(c, g, i - 1, j, k), cget(c, g, i, j, k), cget(c, g, i + 1, j, k), uf, uf * idh); }
			if (i != g.nx - 1) { uf = u[g.uidx(i + 1, j, k)]; Fr = vanleer_flux(cget(c, g, i - 1, j, k), cget(c, g, i, j, k), cget(c, g, i + 1, j, k), cget(c, g, i + 2, j, k), uf, uf * idh); }
			return c[g.pidx(i, j, k)] - idh * (Fr - Fl);
		}
		SCOUR_HD inline double advy_node(const double* c, const double* v, MacGrid g, double dt, int i, int j, int k)
		{
			double idh = dt / g.h;
			double vf, Fl = 0.0, Fr = 0.0;
			if (j != 0) { vf = v[g.vidx(i, j, k)]; Fl = vanleer_flux(cget(c, g, i, j - 2, k), cget(c, g, i, j - 1, k), cget(c, g, i, j, k), cget(c, g, i, j + 1, k), vf, vf * idh); }
			if (j != g.ny - 1) { vf = v[g.vidx(i, j + 1, k)]; Fr = vanleer_flux(cget(c, g, i, j - 1, k), cget(c, g, i, j, k), cget(c, g, i, j + 1, k), cget(c, g, i, j + 2, k), vf, vf * idh); }
			return c[g.pidx(i, j, k)] - idh * (Fr - Fl);
		}
		SCOUR_HD inline double advz_node(const double* c, const double* w, const double* ws, MacGrid g, double dt, int i, int j, int k)
		{
			double idh = dt / g.h;
			double Fb = 0.0, Ft = 0.0;
			if (k != 0) { double wsf = 0.5 * (ws[g.pidx(i, j, k - 1)] + ws[g.pidx(i, j, k)]); double vf = w[g.widx(i, j, k)] - wsf;
				Fb = vanleer_flux(cget(c, g, i, j, k - 2), cget(c, g, i, j, k - 1), cget(c, g, i, j, k), cget(c, g, i, j, k + 1), vf, vf * idh); }
			if (k != g.nz - 1) { double wsf = 0.5 * (ws[g.pidx(i, j, k)] + ws[g.pidx(i, j, k + 1)]); double vf = w[g.widx(i, j, k + 1)] - wsf;
				Ft = vanleer_flux(cget(c, g, i, j, k - 1), cget(c, g, i, j, k), cget(c, g, i, j, k + 1), cget(c, g, i, j, k + 2), vf, vf * idh); }
			return c[g.pidx(i, j, k)] - idh * (Ft - Fb);
		}

		// ---- kernels -------------------------------------------------------------------------
		__global__ void k_eff_ws(const double* c, double* ws, MacGrid g, double ws0, int hindered, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			ws[t] = eff_ws_node(c[t], ws0, hindered);
		}
		__global__ void k_fwd(const double* cIn, const double* u, const double* v, const double* w, const double* ws,
			double* phiHat, MacGrid g, BC vbc, ScalarBC sbc, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			forward_c_node(cIn, u, v, w, ws, phiHat, g, vbc, sbc, dt, i, j, k);
		}
		__global__ void k_corr(const double* cIn, const double* phiHat, const double* u, const double* v, const double* w,
			const double* ws, double* out, MacGrid g, BC vbc, ScalarBC sbc, double dt, int band, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			correct_c_node(cIn, phiHat, u, v, w, ws, out, g, vbc, sbc, dt, band, i, j, k);
		}
		__global__ void k_advx(const double* c, const double* u, double* out, MacGrid g, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			out[t] = advx_node(c, u, g, dt, i, j, k);
		}
		__global__ void k_advy(const double* c, const double* v, double* out, MacGrid g, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			out[t] = advy_node(c, v, g, dt, i, j, k);
		}
		__global__ void k_advz(const double* c, const double* w, const double* ws, double* out, MacGrid g, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			out[t] = advz_node(c, w, ws, g, dt, i, j, k);
		}
		__global__ void k_diffuse(const double* cIn, double* cOut, const double* Dc, MacGrid g, ScalarBC sbc, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			cOut[t] = diffuse_c_node(cIn, Dc, g, sbc, dt, i, j, k);
		}
		__global__ void k_bed(double* c, const double* tau_prime, double* bedflux, MacGrid g, SedBedParams p, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return; // n = nx*ny
			int i = t % g.nx, j = t / g.nx;
			bed_exchange_node(c, tau_prime, bedflux, g, p, dt, i, j);
		}
		__global__ void k_bflux(double* c, const double* u, const double* gmin, const double* gmax, double* net_flux, MacGrid g, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return; // n = ny*nz
			int j = t % g.ny, k = t / g.ny;
			boundary_flux_node(c, u, gmin, gmax, net_flux, g, dt, j, k);
		}
	}

	// ================= launchers =================================================
	void suspended_effective_ws_gpu(const double* c, double* ws, MacGrid g, double ws0, int hindered)
	{
		int n = g.p_count(); k_eff_ws<<<gsz(n), 256>>>(c, ws, g, ws0, hindered, n);
	}

	void suspended_advect_gpu(const double* cIn, const double* u, const double* v, const double* w,
		const double* ws, double* cOut, double* scratch, MacGrid g, BC vbc, ScalarBC sbc, double dt, int band)
	{
		int n = g.p_count();
		k_fwd<<<gsz(n), 256>>>(cIn, u, v, w, ws, scratch, g, vbc, sbc, dt, n);
		k_corr<<<gsz(n), 256>>>(cIn, scratch, u, v, w, ws, cOut, g, vbc, sbc, dt, band, n);
	}

	void suspended_advect_cons_gpu(const double* cIn, const double* u, const double* v, const double* w,
		const double* ws, double* cOut, double* scratch, MacGrid g, double dt)
	{
		int n = g.p_count();
		// dimensional split cIn →(x) cOut →(y) scratch →(z) cOut.
		k_advx<<<gsz(n), 256>>>(cIn, u, cOut, g, dt, n);
		k_advy<<<gsz(n), 256>>>(cOut, v, scratch, g, dt, n);
		k_advz<<<gsz(n), 256>>>(scratch, w, ws, cOut, g, dt, n);
	}

	void suspended_diffuse_gpu(const double* cIn, double* cOut, const double* Dc, MacGrid g, ScalarBC sbc, double dt)
	{
		int n = g.p_count(); k_diffuse<<<gsz(n), 256>>>(cIn, cOut, Dc, g, sbc, dt, n);
	}

	void suspended_bed_exchange_gpu(double* c, const double* tau_prime, double* bedflux, MacGrid g, SedBedParams p, double dt)
	{
		int n = g.nx * g.ny; k_bed<<<gsz(n), 256>>>(c, tau_prime, bedflux, g, p, dt, n);
	}

	void suspended_boundary_flux_gpu(double* c, const double* u, const double* ghost_xmin, const double* ghost_xmax,
		double* net_flux, MacGrid g, double dt)
	{
		int n = g.ny * g.nz; k_bflux<<<gsz(n), 256>>>(c, u, ghost_xmin, ghost_xmax, net_flux, g, dt, n);
	}

	// ================= CPU reference twins =======================================
	void suspended_effective_ws_cpu(const std::vector<double>& c, std::vector<double>& ws, MacGrid g, double ws0, int hindered)
	{
		ws.assign(g.p_count(), 0.0);
		for (int t = 0; t < g.p_count(); ++t) ws[t] = eff_ws_node(c[t], ws0, hindered);
	}

	void suspended_advect_cpu(const std::vector<double>& cIn, const std::vector<double>& u, const std::vector<double>& v,
		const std::vector<double>& w, const std::vector<double>& ws, std::vector<double>& cOut,
		MacGrid g, BC vbc, ScalarBC sbc, double dt, int band)
	{
		int n = g.p_count();
		std::vector<double> phiHat(n, 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			forward_c_node(cIn.data(), u.data(), v.data(), w.data(), ws.data(), phiHat.data(), g, vbc, sbc, dt, i, j, k);
		cOut.assign(n, 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			correct_c_node(cIn.data(), phiHat.data(), u.data(), v.data(), w.data(), ws.data(), cOut.data(), g, vbc, sbc, dt, band, i, j, k);
	}

	void suspended_advect_cons_cpu(const std::vector<double>& cIn, const std::vector<double>& u, const std::vector<double>& v,
		const std::vector<double>& w, const std::vector<double>& ws, std::vector<double>& cOut, MacGrid g, double dt)
	{
		int n = g.p_count();
		std::vector<double> a(n, 0.0), b(n, 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			a[g.pidx(i, j, k)] = advx_node(cIn.data(), u.data(), g, dt, i, j, k);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			b[g.pidx(i, j, k)] = advy_node(a.data(), v.data(), g, dt, i, j, k);
		cOut.assign(n, 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			cOut[g.pidx(i, j, k)] = advz_node(b.data(), w.data(), ws.data(), g, dt, i, j, k);
	}

	void suspended_diffuse_cpu(const std::vector<double>& cIn, std::vector<double>& cOut, const std::vector<double>& Dc,
		MacGrid g, ScalarBC sbc, double dt)
	{
		cOut.assign(g.p_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			cOut[g.pidx(i, j, k)] = diffuse_c_node(cIn.data(), Dc.data(), g, sbc, dt, i, j, k);
	}

	void suspended_bed_exchange_cpu(std::vector<double>& c, const std::vector<double>& tau_prime, std::vector<double>& bedflux,
		MacGrid g, SedBedParams p, double dt)
	{
		bedflux.assign((size_t)g.nx * g.ny, 0.0);
		const double* tp = tau_prime.empty() ? nullptr : tau_prime.data();
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			bed_exchange_node(c.data(), tp, bedflux.data(), g, p, dt, i, j);
	}

	void suspended_boundary_flux_cpu(std::vector<double>& c, const std::vector<double>& u,
		const std::vector<double>& ghost_xmin, const std::vector<double>& ghost_xmax,
		std::vector<double>* net_flux, MacGrid g, double dt)
	{
		int nyz = g.ny * g.nz;
		if (net_flux) net_flux->assign((size_t)2 * nyz, 0.0);
		double* nf = net_flux ? net_flux->data() : nullptr;
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j)
			boundary_flux_node(c.data(), u.data(), ghost_xmin.data(), ghost_xmax.data(), nf, g, dt, j, k);
	}
}
