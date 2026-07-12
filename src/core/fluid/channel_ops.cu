// channel_ops.cu — M2 open-channel MAC operators (masked + open BCs). See channel_ops.h.
// Every kernel has a serial CPU twin exercised by a GPU-vs-CPU gtest (rel. max-norm
// 1e-5). RESEARCH §3 (solver), §8 (BCs), research/08, research/14.
#include "core/fluid/channel_ops.h"
#include "core/fluid/mac_ops.h" // reductions (reduce_sum_gpu, reduce_max_abs_gpu)

#include <cuda_runtime.h>

#include <cmath>
#include <vector>

namespace windcfd::core
{
	namespace
	{
		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }

		// component geometry ----------------------------------------------------
		WINDCFD_HD inline void comp_extent(int comp, MacGrid g, int& ni, int& nj, int& nk)
		{
			if (comp == 0) { ni = g.nx + 1; nj = g.ny; nk = g.nz; }
			else if (comp == 1) { ni = g.nx; nj = g.ny + 1; nk = g.nz; }
			else { ni = g.nx; nj = g.ny; nk = g.nz + 1; }
		}
		WINDCFD_HD inline int comp_idx(int comp, MacGrid g, int i, int j, int k)
		{ return comp == 0 ? g.uidx(i, j, k) : comp == 1 ? g.vidx(i, j, k) : g.widx(i, j, k); }
		WINDCFD_HD inline void node_pos(int comp, MacGrid g, int i, int j, int k, double& x, double& y, double& z)
		{
			if (comp == 0) { x = g.xf(i); y = g.yc(j); z = g.zc(k); }
			else if (comp == 1) { x = g.xc(i); y = g.yf(j); z = g.zc(k); }
			else { x = g.xc(i); y = g.yc(j); z = g.zf(k); }
		}
		// advection updates the interior faces; inlet(i=0)/outlet(i=nx) u-faces and the
		// free-slip normal v/w faces are set by the BC steps, so copy them through.
		WINDCFD_HD inline bool comp_interior(int comp, MacGrid g, int i, int j, int k)
		{
			if (comp == 0) return i >= 1 && i <= g.nx - 1 && j >= 0 && j <= g.ny - 1 && k >= 0 && k <= g.nz - 1;
			if (comp == 1) return i >= 0 && i <= g.nx - 1 && j >= 1 && j <= g.ny - 1 && k >= 0 && k <= g.nz - 1;
			return i >= 0 && i <= g.nx - 1 && j >= 0 && j <= g.ny - 1 && k >= 1 && k <= g.nz - 1;
		}
		WINDCFD_HD inline bool comp_solidface(int comp, const unsigned char* s, MacGrid g, int i, int j, int k)
		{ return comp == 0 ? ch_usolid(s, g, i, j, k) : comp == 1 ? ch_vsolid(s, g, i, j, k) : ch_wsolid(s, g, i, j, k); }

		WINDCFD_HD inline double chf(int comp, const double* f, MacGrid g, ChannelBC bc, int i, int j, int k)
		{ return comp == 0 ? ch_fetch_u(f, g, bc, i, j, k) : comp == 1 ? ch_fetch_v(f, g, bc, i, j, k) : ch_fetch_w(f, g, bc, i, j, k); }
		WINDCFD_HD inline double chsample(int comp, const double* f, MacGrid g, ChannelBC bc, double x, double y, double z, double* mn, double* mx)
		{ return comp == 0 ? ch_trilerp_u(f, g, bc, x, y, z, mn, mx) : comp == 1 ? ch_trilerp_v(f, g, bc, x, y, z, mn, mx) : ch_trilerp_w(f, g, bc, x, y, z, mn, mx); }

		// ---- advection helpers ------------------------------------------------
		WINDCFD_HD inline bool cell_of_point_solid(MacGrid g, const unsigned char* solid, double x, double y, double z)
		{
			int i = (int)floor(grid_fx(g, x)), j = (int)floor(grid_fy(g, y)), k = (int)floor(grid_fz(g, z));
			return ch_is_solid(solid, g, i, j, k);
		}
		// Pull a backtrace/forward-trace endpoint out of any solid it landed in, by
		// bisecting toward the (fluid) origin node (RESEARCH §3.2 ray clipping).
		WINDCFD_HD inline void clip_solid(MacGrid g, const unsigned char* solid, double x0, double y0, double z0, double& x, double& y, double& z)
		{
			if (!cell_of_point_solid(g, solid, x, y, z)) return;
			double lo = 0.0, hi = 1.0; // fraction from origin(0) to endpoint(1)
			for (int it = 0; it < 12; ++it)
			{
				double m = 0.5 * (lo + hi);
				double xm = x0 + m * (x - x0), ym = y0 + m * (y - y0), zm = z0 + m * (z - z0);
				if (cell_of_point_solid(g, solid, xm, ym, zm)) hi = m; else lo = m;
			}
			x = x0 + lo * (x - x0); y = y0 + lo * (y - y0); z = z0 + lo * (z - z0);
		}
		WINDCFD_HD inline void vel_at(const double* u, const double* v, const double* w, MacGrid g, ChannelBC bc,
			double x, double y, double z, double& vx, double& vy, double& vz)
		{
			vx = ch_trilerp_u(u, g, bc, x, y, z, nullptr, nullptr);
			vy = ch_trilerp_v(v, g, bc, x, y, z, nullptr, nullptr);
			vz = ch_trilerp_w(w, g, bc, x, y, z, nullptr, nullptr);
		}
		WINDCFD_HD inline void rk2_trace(const double* u, const double* v, const double* w, MacGrid g, ChannelBC bc,
			const unsigned char* solid, double dt, double x, double y, double z, double& xo, double& yo, double& zo)
		{
			double vx, vy, vz; vel_at(u, v, w, g, bc, x, y, z, vx, vy, vz);
			double xm = x - 0.5 * dt * vx, ym = y - 0.5 * dt * vy, zm = z - 0.5 * dt * vz;
			clamp_to_domain(g, xm, ym, zm); clip_solid(g, solid, x, y, z, xm, ym, zm);
			double vmx, vmy, vmz; vel_at(u, v, w, g, bc, xm, ym, zm, vmx, vmy, vmz);
			xo = x - dt * vmx; yo = y - dt * vmy; zo = z - dt * vmz;
			clamp_to_domain(g, xo, yo, zo); clip_solid(g, solid, x, y, z, xo, yo, zo);
		}
		// 1st-order reversion band: near a domain edge OR near a solid (via nearsolid).
		WINDCFD_HD inline bool near_revert(int comp, MacGrid g, const unsigned char* nearsolid, int i, int j, int k)
		{
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk);
			if (i < 1 || i > ni - 2 || j < 1 || j > nj - 2 || k < 1 || k > nk - 2) return true;
			// nearsolid is a cell field; a face is near-solid if either straddling cell is.
			int ca_i = i, cb_i = i, ca_j = j, cb_j = j, ca_k = k, cb_k = k;
			if (comp == 0) { ca_i = i - 1; cb_i = i; }
			else if (comp == 1) { ca_j = j - 1; cb_j = j; }
			else { ca_k = k - 1; cb_k = k; }
			auto ns = [&](int ii, int jj, int kk) -> bool {
				if (ii < 0 || ii > g.nx - 1 || jj < 0 || jj > g.ny - 1 || kk < 0 || kk > g.nz - 1) return false;
				return nearsolid[g.pidx(ii, jj, kk)] != 0; };
			return ns(ca_i, ca_j, ca_k) || ns(cb_i, cb_j, cb_k);
		}

		WINDCFD_HD inline void forward_node(int comp, const double* field, const double* u, const double* v, const double* w,
			double* phiHat, const unsigned char* solid, MacGrid g, ChannelBC bc, double dt, int i, int j, int k)
		{
			int idx = comp_idx(comp, g, i, j, k);
			if (comp_solidface(comp, solid, g, i, j, k)) { phiHat[idx] = 0.0; return; }
			if (!comp_interior(comp, g, i, j, k)) { phiHat[idx] = field[idx]; return; }
			double x, y, z; node_pos(comp, g, i, j, k, x, y, z);
			double xb, yb, zb; rk2_trace(u, v, w, g, bc, solid, dt, x, y, z, xb, yb, zb);
			phiHat[idx] = chsample(comp, field, g, bc, xb, yb, zb, nullptr, nullptr);
		}
		WINDCFD_HD inline void correct_node(int comp, const double* field, const double* phiHat,
			const double* u, const double* v, const double* w, double* out,
			const unsigned char* solid, const unsigned char* nearsolid, MacGrid g, ChannelBC bc, double dt, int i, int j, int k)
		{
			int idx = comp_idx(comp, g, i, j, k);
			if (comp_solidface(comp, solid, g, i, j, k)) { out[idx] = 0.0; return; }
			if (!comp_interior(comp, g, i, j, k)) { out[idx] = field[idx]; return; }
			double phf = phiHat[idx];
			if (near_revert(comp, g, nearsolid, i, j, k)) { out[idx] = phf; return; } // 1st order
			double x, y, z; node_pos(comp, g, i, j, k, x, y, z);
			double xb, yb, zb; rk2_trace(u, v, w, g, bc, solid, dt, x, y, z, xb, yb, zb);
			double mn, mx; chsample(comp, field, g, bc, xb, yb, zb, &mn, &mx);
			double xf, yf, zf; rk2_trace(u, v, w, g, bc, solid, -dt, x, y, z, xf, yf, zf);
			double phn = chsample(comp, phiHat, g, bc, xf, yf, zf, nullptr, nullptr);
			double corrected = phf + 0.5 * (field[idx] - phn);
			if (corrected < mn) corrected = mn; else if (corrected > mx) corrected = mx;
			out[idx] = corrected;
		}

		__global__ void k_forward(int comp, const double* field, const double* u, const double* v, const double* w,
			double* phiHat, const unsigned char* solid, MacGrid g, ChannelBC bc, double dt, int total)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= total) return;
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk);
			int i = t % ni, j = (t / ni) % nj, k = t / (ni * nj);
			forward_node(comp, field, u, v, w, phiHat, solid, g, bc, dt, i, j, k);
		}
		__global__ void k_correct(int comp, const double* field, const double* phiHat,
			const double* u, const double* v, const double* w, double* out,
			const unsigned char* solid, const unsigned char* nearsolid, MacGrid g, ChannelBC bc, double dt, int total)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= total) return;
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk);
			int i = t % ni, j = (t / ni) % nj, k = t / (ni * nj);
			correct_node(comp, field, phiHat, u, v, w, out, solid, nearsolid, g, bc, dt, i, j, k);
		}
		void advect_one_gpu(int comp, const double* field, const double* u, const double* v, const double* w,
			double* out, double* phiHat, const unsigned char* solid, const unsigned char* nearsolid, MacGrid g, ChannelBC bc, double dt)
		{
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk);
			int total = ni * nj * nk;
			k_forward<<<gsz(total), 256>>>(comp, field, u, v, w, phiHat, solid, g, bc, dt, total);
			k_correct<<<gsz(total), 256>>>(comp, field, phiHat, u, v, w, out, solid, nearsolid, g, bc, dt, total);
		}

		// ---- Smagorinsky ------------------------------------------------------
		WINDCFD_HD inline double uc_(const double* u, MacGrid g, ChannelBC bc, int i, int j, int k)
		{ return 0.5 * (ch_fetch_u(u, g, bc, i, j, k) + ch_fetch_u(u, g, bc, i + 1, j, k)); }
		WINDCFD_HD inline double vc_(const double* v, MacGrid g, ChannelBC bc, int i, int j, int k)
		{ return 0.5 * (ch_fetch_v(v, g, bc, i, j, k) + ch_fetch_v(v, g, bc, i, j + 1, k)); }
		WINDCFD_HD inline double wc_(const double* w, MacGrid g, ChannelBC bc, int i, int j, int k)
		{ return 0.5 * (ch_fetch_w(w, g, bc, i, j, k) + ch_fetch_w(w, g, bc, i, j, k + 1)); }
		WINDCFD_HD inline double smag_nut(const double* u, const double* v, const double* w, const unsigned char* solid,
			MacGrid g, ChannelBC bc, double Cs, int i, int j, int k)
		{
			if (ch_is_solid(solid, g, i, j, k)) return 0.0;
			double invx = 1.0 / g.dx(i), invy = 1.0 / g.dy(j), invz = 1.0 / g.dz(k), sx = 1.0 / (g.gapx(i - 1) + g.gapx(i)), sy = 1.0 / (g.gapy(j - 1) + g.gapy(j)), sz = 1.0 / (g.gapz(k - 1) + g.gapz(k));
			double dudx = (ch_fetch_u(u, g, bc, i + 1, j, k) - ch_fetch_u(u, g, bc, i, j, k)) * invx;
			double dvdy = (ch_fetch_v(v, g, bc, i, j + 1, k) - ch_fetch_v(v, g, bc, i, j, k)) * invy;
			double dwdz = (ch_fetch_w(w, g, bc, i, j, k + 1) - ch_fetch_w(w, g, bc, i, j, k)) * invz;
			double dudy = (uc_(u, g, bc, i, j + 1, k) - uc_(u, g, bc, i, j - 1, k)) * sy;
			double dudz = (uc_(u, g, bc, i, j, k + 1) - uc_(u, g, bc, i, j, k - 1)) * sz;
			double dvdx = (vc_(v, g, bc, i + 1, j, k) - vc_(v, g, bc, i - 1, j, k)) * sx;
			double dvdz = (vc_(v, g, bc, i, j, k + 1) - vc_(v, g, bc, i, j, k - 1)) * sz;
			double dwdx = (wc_(w, g, bc, i + 1, j, k) - wc_(w, g, bc, i - 1, j, k)) * sx;
			double dwdy = (wc_(w, g, bc, i, j + 1, k) - wc_(w, g, bc, i, j - 1, k)) * sy;
			double Sxx = dudx, Syy = dvdy, Szz = dwdz;
			double Sxy = 0.5 * (dudy + dvdx), Sxz = 0.5 * (dudz + dwdx), Syz = 0.5 * (dvdz + dwdy);
			double Smag = sqrt(2.0 * (Sxx * Sxx + Syy * Syy + Szz * Szz + 2.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz)));
			double cd = Cs * cbrt(g.dx(i) * g.dy(j) * g.dz(k)); return cd * cd * Smag;
		}
		__global__ void k_nut(const double* u, const double* v, const double* w, double* nut, const unsigned char* solid,
			MacGrid g, ChannelBC bc, double Cs, int total)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= total) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			nut[g.pidx(i, j, k)] = smag_nut(u, v, w, solid, g, bc, Cs, i, j, k);
		}

		// ---- diffusion --------------------------------------------------------
		WINDCFD_HD inline double wall_ghost(double c, ChannelBC bc) { return bc.solid_mode == SOLID_NOSLIP ? -c : c; }
		WINDCFD_HD inline int normal_axis(int comp) { return comp; }
		WINDCFD_HD inline void diffuse_node(int comp, const double* field, const double* nut, const unsigned char* solid,
			double* out, MacGrid g, ChannelBC bc, double dt, double nu, int i, int j, int k)
		{
			int idx = comp_idx(comp, g, i, j, k);
			if (comp_solidface(comp, solid, g, i, j, k)) { out[idx] = 0.0; return; }
			if (!comp_interior(comp, g, i, j, k)) { out[idx] = field[idx]; return; }
			double c = field[idx]; double hL, hR;
			int na = normal_axis(comp);
			double lap = 0.0;
			// x pair
			{
				double np, nm;
				if (na == 0) { np = chf(comp, field, g, bc, i + 1, j, k); nm = chf(comp, field, g, bc, i - 1, j, k); }
				else { np = comp_solidface(comp, solid, g, i + 1, j, k) ? wall_ghost(c, bc) : chf(comp, field, g, bc, i + 1, j, k);
					nm = comp_solidface(comp, solid, g, i - 1, j, k) ? wall_ghost(c, bc) : chf(comp, field, g, bc, i - 1, j, k); }
				hL = (na == 0) ? g.dx(i - 1) : g.gapx(i - 1); hR = (na == 0) ? g.dx(i) : g.gapx(i);
				lap += d2_axis(nm, c, np, hL, hR);
			}
			// y pair
			{
				double np, nm;
				if (na == 1) { np = chf(comp, field, g, bc, i, j + 1, k); nm = chf(comp, field, g, bc, i, j - 1, k); }
				else { np = comp_solidface(comp, solid, g, i, j + 1, k) ? wall_ghost(c, bc) : chf(comp, field, g, bc, i, j + 1, k);
					nm = comp_solidface(comp, solid, g, i, j - 1, k) ? wall_ghost(c, bc) : chf(comp, field, g, bc, i, j - 1, k); }
				hL = (na == 1) ? g.dy(j - 1) : g.gapy(j - 1); hR = (na == 1) ? g.dy(j) : g.gapy(j);
				lap += d2_axis(nm, c, np, hL, hR);
			}
			// z pair
			{
				double np, nm;
				if (na == 2) { np = chf(comp, field, g, bc, i, j, k + 1); nm = chf(comp, field, g, bc, i, j, k - 1); }
				else { np = comp_solidface(comp, solid, g, i, j, k + 1) ? wall_ghost(c, bc) : chf(comp, field, g, bc, i, j, k + 1);
					nm = comp_solidface(comp, solid, g, i, j, k - 1) ? wall_ghost(c, bc) : chf(comp, field, g, bc, i, j, k - 1); }
				hL = (na == 2) ? g.dz(k - 1) : g.gapz(k - 1); hR = (na == 2) ? g.dz(k) : g.gapz(k);
				lap += d2_axis(nm, c, np, hL, hR);
			}
			double nut_f = 0.0;
			if (nut)
			{
				int a0i = i, a0j = j, a0k = k, a1i = i, a1j = j, a1k = k;
				if (comp == 0) { a0i = i - 1; }
				else if (comp == 1) { a0j = j - 1; }
				else { a0k = k - 1; }
				double n0 = (a0i >= 0 && a0j >= 0 && a0k >= 0 && a0i < g.nx && a0j < g.ny && a0k < g.nz) ? nut[g.pidx(a0i, a0j, a0k)] : 0.0;
				double n1 = (a1i < g.nx && a1j < g.ny && a1k < g.nz) ? nut[g.pidx(a1i, a1j, a1k)] : 0.0;
				nut_f = 0.5 * (n0 + n1);
			}
			out[idx] = c + dt * (nu + nut_f) * lap;
		}
		__global__ void k_diffuse(int comp, const double* field, const double* nut, const unsigned char* solid, double* out,
			MacGrid g, ChannelBC bc, double dt, double nu, int ni, int nj, int nk)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; int total = ni * nj * nk; if (t >= total) return;
			int i = t % ni, j = (t / ni) % nj, k = t / (ni * nj);
			diffuse_node(comp, field, nut, solid, out, g, bc, dt, nu, i, j, k);
		}
		void diffuse_one_gpu(int comp, const double* field, const double* nut, const unsigned char* solid, double* out, MacGrid g, ChannelBC bc, double dt, double nu)
		{
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk);
			int total = ni * nj * nk;
			k_diffuse<<<gsz(total), 256>>>(comp, field, nut, solid, out, g, bc, dt, nu, ni, nj, nk);
		}

		// ---- boundary kernels -------------------------------------------------
		__global__ void k_inlet(double* u, MacGrid g, ChannelBC bc, int njk)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= njk) return;
			int j = t % g.ny, k = t / g.ny;
			int iplane = bc.flow_sign < 0 ? g.nx : 0; // reversed tide: Dirichlet inlet on the xmax face
			u[g.uidx(iplane, j, k)] = channel_inlet_u(bc, g, k);
		}
		__global__ void k_orlanski(double* u, MacGrid g, ChannelBC bc, double Uc_dt, int njk)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= njk) return;
			int j = t % g.ny, k = t / g.ny;
			// coef = Uc·dt/Δx with Δx the LOCAL outlet cell width — NOT the scalar h. On a graded grid the
			// outlet sits in the coarse far field (dx(nx-1) ≫ h_fine), so dividing by h would over-convect the
			// outlet by dx/h (≈6× here) and pump an instability upstream from the exit. Uniform ⇒ dx==h (same).
			if (bc.flow_sign < 0)
			{
				// Reversed tide: outlet on the xmin face; outflow is −x, so clamp to ≤ 0 (an incoming +x
				// jet at this exit is what the flux-rescale would otherwise amplify).
				double coef = Uc_dt / g.dx(0); if (coef > 1.0) coef = 1.0; else if (coef < 0.0) coef = 0.0;
				double ub = u[g.uidx(0, j, k)], uin = u[g.uidx(1, j, k)];
				double uo = ub - coef * (ub - uin); u[g.uidx(0, j, k)] = uo < 0.0 ? uo : 0.0;
			}
			else
			{
				double coef = Uc_dt / g.dx(g.nx - 1); if (coef > 1.0) coef = 1.0; else if (coef < 0.0) coef = 0.0;
				double ub = u[g.uidx(g.nx, j, k)], uin = u[g.uidx(g.nx - 1, j, k)];
				double uo = ub - coef * (ub - uin); u[g.uidx(g.nx, j, k)] = uo > 0.0 ? uo : 0.0; // outflow-only: clamp reversed inflow at the +X exit (else the flux-rescale amplifies a reversed jet — the seabed exit runaway)
			}
		}
		__global__ void k_solidbc(double* u, double* v, double* w, const unsigned char* solid, MacGrid g)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x;
			int nu = g.u_count(), nv = g.v_count(), nw = g.w_count();
			if (t < nu)
			{
				int i = t % (g.nx + 1), j = (t / (g.nx + 1)) % g.ny, k = t / ((g.nx + 1) * g.ny);
				if (i >= 1 && i <= g.nx - 1 && (ch_is_solid(solid, g, i - 1, j, k) || ch_is_solid(solid, g, i, j, k))) u[t] = 0.0;
			}
			else if (t < nu + nv)
			{
				int s = t - nu; int i = s % g.nx, j = (s / g.nx) % (g.ny + 1), k = s / (g.nx * (g.ny + 1));
				if (j >= 1 && j <= g.ny - 1 && (ch_is_solid(solid, g, i, j - 1, k) || ch_is_solid(solid, g, i, j, k))) v[s] = 0.0;
			}
			else if (t < nu + nv + nw)
			{
				int s = t - nu - nv; int i = s % g.nx, j = (s / g.nx) % g.ny, k = s / (g.nx * g.ny);
				if (k >= 1 && k <= g.nz - 1 && (ch_is_solid(solid, g, i, j, k - 1) || ch_is_solid(solid, g, i, j, k))) w[s] = 0.0;
			}
		}
		__global__ void k_gather_uplane(const double* u, double* out, MacGrid g, int iplane, int njk)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= njk) return;
			int j = t % g.ny, k = t / g.ny;
			// Weight by the TRUE face area dy(j)·dz(k) — NOT a uniform h². On a graded grid the inlet/outlet
			// plane has non-uniform cells; a plain Σu·h² over-counts the fine cells, so q_in/q_out ≠ 1 even
			// when mass IS conserved, and the outlet mass-rescale then injects a spurious kick every step
			// (the corner blow-up). dy·dz == h² on a uniform grid, so this is byte-identical there.
			out[t] = u[g.uidx(iplane, j, k)] * g.dy(j) * g.dz(k);
		}
		__global__ void k_scale_uplane(double* u, MacGrid g, int iplane, double s, int njk)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= njk) return;
			int j = t % g.ny, k = t / g.ny;
			u[g.uidx(iplane, j, k)] *= s;
		}

		// ---- masked Poisson pieces -------------------------------------------
		// dir_xmax encodes the Dirichlet-p=0 OUTLET face: >0 ⇒ xmax (i=nx−1 cell, the M2/M3 default),
		// <0 ⇒ xmin (i=0 cell, the reversed-tide outlet), 0 ⇒ none. The velocity-Dirichlet inlet face
		// stays Neumann (no p-neighbour). All dir_xmax>0 paths are identical to the original code.
		WINDCFD_HD inline void nb_stencil(const double* p, const unsigned char* solid, MacGrid g, int dir_xmax, int i, int j, int k, int& count, double& nbsum)
		{
			count = 0; nbsum = 0.0;
			if (i > 0 && !ch_is_solid(solid, g, i - 1, j, k)) { ++count; nbsum += p[g.pidx(i - 1, j, k)]; }
			else if (i == 0 && dir_xmax < 0) { ++count; } // reversed outlet Dirichlet p=0 at xmin
			if (i < g.nx - 1) { if (!ch_is_solid(solid, g, i + 1, j, k)) { ++count; nbsum += p[g.pidx(i + 1, j, k)]; } }
			else if (dir_xmax > 0) { ++count; } // outlet Dirichlet p=0 at xmax
			if (j > 0 && !ch_is_solid(solid, g, i, j - 1, k)) { ++count; nbsum += p[g.pidx(i, j - 1, k)]; }
			if (j < g.ny - 1 && !ch_is_solid(solid, g, i, j + 1, k)) { ++count; nbsum += p[g.pidx(i, j + 1, k)]; }
			if (k > 0 && !ch_is_solid(solid, g, i, j, k - 1)) { ++count; nbsum += p[g.pidx(i, j, k - 1)]; }
			if (k < g.nz - 1 && !ch_is_solid(solid, g, i, j, k + 1)) { ++count; nbsum += p[g.pidx(i, j, k + 1)]; }
		}
		// Variable-coefficient FV Poisson stencil (masked + Dirichlet-outlet), the metric-aware form of
		// nb_stencil. Per real fluid face f: w_f = 1/(d_centre-to-centre · cell_width_normal). The
		// Dirichlet-outlet ghost (p=0) on the domain x-face uses d_f = the edge cell width (the existing
		// mirror-cell convention, weight 1/dx² uniform). diag = Σ w_f, wnb = Σ w_f·p_nb ⇒ Ap = diag·p_c − wnb.
		// Uniform spacings ⇒ (count·p_c − nbsum)/h², identical to nb_stencil. Stays SPD.
		WINDCFD_HD inline void fv_stencil_ch(const double* p, const unsigned char* solid, MacGrid g, int dir_xmax, int i, int j, int k, double& diag, double& wnb)
		{
			diag = 0.0; wnb = 0.0;
			if (i > 0 && !ch_is_solid(solid, g, i - 1, j, k)) { double w = 1.0 / (g.dxc(i) * g.dx(i)); diag += w; wnb += w * p[g.pidx(i - 1, j, k)]; }
			else if (i == 0 && dir_xmax < 0) diag += 1.0 / (g.dx(0) * g.dx(0)); // reversed outlet Dirichlet p=0 at xmin
			if (i < g.nx - 1) { if (!ch_is_solid(solid, g, i + 1, j, k)) { double w = 1.0 / (g.dxc(i + 1) * g.dx(i)); diag += w; wnb += w * p[g.pidx(i + 1, j, k)]; } }
			else if (dir_xmax > 0) diag += 1.0 / (g.dx(g.nx - 1) * g.dx(g.nx - 1)); // outlet Dirichlet p=0 at xmax
			if (j > 0 && !ch_is_solid(solid, g, i, j - 1, k)) { double w = 1.0 / (g.dyc(j) * g.dy(j)); diag += w; wnb += w * p[g.pidx(i, j - 1, k)]; }
			if (j < g.ny - 1 && !ch_is_solid(solid, g, i, j + 1, k)) { double w = 1.0 / (g.dyc(j + 1) * g.dy(j)); diag += w; wnb += w * p[g.pidx(i, j + 1, k)]; }
			if (k > 0 && !ch_is_solid(solid, g, i, j, k - 1)) { double w = 1.0 / (g.dzc(k) * g.dz(k)); diag += w; wnb += w * p[g.pidx(i, j, k - 1)]; }
			if (k < g.nz - 1 && !ch_is_solid(solid, g, i, j, k + 1)) { double w = 1.0 / (g.dzc(k + 1) * g.dz(k)); diag += w; wnb += w * p[g.pidx(i, j, k + 1)]; }
		}
		WINDCFD_HD inline double div_cell(const double* u, const double* v, const double* w, MacGrid g, int i, int j, int k)
		{
			return (u[g.uidx(i + 1, j, k)] - u[g.uidx(i, j, k)]) / g.dx(i)
				+ (v[g.vidx(i, j + 1, k)] - v[g.vidx(i, j, k)]) / g.dy(j)
				+ (w[g.widx(i, j, k + 1)] - w[g.widx(i, j, k)]) / g.dz(k);
		}
		__global__ void k_rhs(const double* u, const double* v, const double* w, double* rhs, const unsigned char* solid, MacGrid g, double rho, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			rhs[t] = ch_is_solid(solid, g, i, j, k) ? 0.0 : -(rho / dt) * div_cell(u, v, w, g, i, j, k);
		}
		__global__ void k_apply(const double* p, double* Ap, const unsigned char* solid, MacGrid g, int dir_xmax, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			if (ch_is_solid(solid, g, i, j, k)) { Ap[t] = 0.0; return; }
			double diag, wnb; fv_stencil_ch(p, solid, g, dir_xmax, i, j, k, diag, wnb);
			Ap[t] = diag * p[t] - wnb;
		}
		__global__ void k_residual(const double* p, const double* rhs, double* r, const unsigned char* solid, MacGrid g, int dir_xmax, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			if (ch_is_solid(solid, g, i, j, k)) { r[t] = 0.0; return; }
			double diag, wnb; fv_stencil_ch(p, solid, g, dir_xmax, i, j, k, diag, wnb);
			r[t] = rhs[t] - (diag * p[t] - wnb);
		}
		__global__ void k_jacobi(const double* p, const double* rhs, double* pout, const unsigned char* solid, MacGrid g, int dir_xmax, double omega, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			if (ch_is_solid(solid, g, i, j, k)) { pout[t] = 0.0; return; }
			double diag, wnb; fv_stencil_ch(p, solid, g, dir_xmax, i, j, k, diag, wnb);
			if (diag == 0.0) { pout[t] = p[t]; return; }
			/* diag already set by fv_stencil_ch */
			double Ap = diag * p[t] - wnb;
			pout[t] = p[t] + omega * (rhs[t] - Ap) / diag;
		}
		__global__ void k_gs(double* p, const double* rhs, const unsigned char* solid, MacGrid g, int dir_xmax, int band, int color, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			if (((i + j + k) & 1) != color) return;
			bool inband = i < band || i >= g.nx - band || j < band || j >= g.ny - band || k < band || k >= g.nz - band;
			if (!inband) return;
			if (ch_is_solid(solid, g, i, j, k)) return;
			double diag, wnb; fv_stencil_ch(p, solid, g, dir_xmax, i, j, k, diag, wnb);
			if (diag == 0.0) return;
			p[t] = (rhs[t] + wnb) / diag;
		}
		// ---- zebra line smoother (graded anisotropy) --------------------------
		// Tridiagonal row of the FV operator at (i,j,k) restricted to a line along `axis`:
		// a/c = −(coupling to the prev/next line cell), b = the FULL FV diagonal, d = rhs +
		// Σ transverse w·p_nb. The face weights REPRODUCE fv_stencil_ch exactly (same
		// expressions, fv_stencil_ch itself untouched to keep the uniform golden bitwise).
		// Solid cells and isolated fluid cells (diag==0) become identity rows, matching
		// k_jacobi/k_zero_solids.
		WINDCFD_HD inline void fv_line_row(const double* p, const double* rhs, const unsigned char* solid, MacGrid g, int dir_xmax,
			int i, int j, int k, int axis, double& a, double& b, double& c, double& d)
		{
			int idx = g.pidx(i, j, k);
			if (ch_is_solid(solid, g, i, j, k)) { a = 0.0; b = 1.0; c = 0.0; d = 0.0; return; }
			double diag = 0.0, tr = 0.0, wlo = 0.0, wup = 0.0;
			if (i > 0 && !ch_is_solid(solid, g, i - 1, j, k)) { double w = 1.0 / (g.dxc(i) * g.dx(i)); diag += w; if (axis == 0) wlo = w; else tr += w * p[g.pidx(i - 1, j, k)]; }
			else if (i == 0 && dir_xmax < 0) diag += 1.0 / (g.dx(0) * g.dx(0)); // reversed outlet Dirichlet p=0 at xmin
			if (i < g.nx - 1) { if (!ch_is_solid(solid, g, i + 1, j, k)) { double w = 1.0 / (g.dxc(i + 1) * g.dx(i)); diag += w; if (axis == 0) wup = w; else tr += w * p[g.pidx(i + 1, j, k)]; } }
			else if (dir_xmax > 0) diag += 1.0 / (g.dx(g.nx - 1) * g.dx(g.nx - 1)); // outlet Dirichlet p=0 at xmax
			if (j > 0 && !ch_is_solid(solid, g, i, j - 1, k)) { double w = 1.0 / (g.dyc(j) * g.dy(j)); diag += w; if (axis == 1) wlo = w; else tr += w * p[g.pidx(i, j - 1, k)]; }
			if (j < g.ny - 1 && !ch_is_solid(solid, g, i, j + 1, k)) { double w = 1.0 / (g.dyc(j + 1) * g.dy(j)); diag += w; if (axis == 1) wup = w; else tr += w * p[g.pidx(i, j + 1, k)]; }
			if (k > 0 && !ch_is_solid(solid, g, i, j, k - 1)) { double w = 1.0 / (g.dzc(k) * g.dz(k)); diag += w; if (axis == 2) wlo = w; else tr += w * p[g.pidx(i, j, k - 1)]; }
			if (k < g.nz - 1 && !ch_is_solid(solid, g, i, j, k + 1)) { double w = 1.0 / (g.dzc(k + 1) * g.dz(k)); diag += w; if (axis == 2) wup = w; else tr += w * p[g.pidx(i, j, k + 1)]; }
			if (diag == 0.0) { a = 0.0; b = 1.0; c = 0.0; d = p[idx]; return; }
			a = -wlo; b = diag; c = -wup; d = rhs[idx] + tr;
		}
		// Solve one whole line exactly (Thomas), in place in p. In-place is race-free under
		// zebra colouring: a line's transverse reads all land on OPPOSITE-colour lines, and a
		// row's own old p is consumed before its slot is reused for the swept RHS d'. Solid
		// faces zero the along-line coupling (a=0), so the chain segments itself across
		// obstacles. A pure-Neumann segment (no transverse fluid coupling anywhere, no
		// Dirichlet end) is singular — its last pivot vanishes; that row degrades to identity.
		// The trigger depends only on matrix coefficients (never on rhs/p), so the smoother
		// stays one fixed linear operator.
		WINDCFD_HD inline void line_solve_ch(double* p, const double* rhs, double* cprime, const unsigned char* solid,
			MacGrid g, int dir_xmax, int axis, int line)
		{
			int n, stride, i = 0, j = 0, k = 0;
			if (axis == 0) { n = g.nx; stride = 1; j = line % g.ny; k = line / g.ny; }
			else if (axis == 1) { n = g.ny; stride = g.nx; i = line % g.nx; k = line / g.nx; }
			else { n = g.nz; stride = g.nx * g.ny; i = line % g.nx; j = line / g.nx; }
			int base = g.pidx(i, j, k); // axis coordinate is 0 here
			double cp = 0.0, dp = 0.0;
			for (int m = 0; m < n; ++m)
			{
				int ii = i, jj = j, kk = k;
				if (axis == 0) ii = m; else if (axis == 1) jj = m; else kk = m;
				int idx = base + m * stride;
				double a, b, c, d;
				fv_line_row(p, rhs, solid, g, dir_xmax, ii, jj, kk, axis, a, b, c, d);
				double denom = b - a * cp;
				if (fabs(denom) > 1e-12 * (fabs(b) + fabs(a * cp))) { cp = c / denom; dp = (d - a * dp) / denom; }
				else { cp = 0.0; dp = p[idx]; } // singular pivot: keep this cell (identity row)
				cprime[idx] = cp;
				p[idx] = dp;
			}
			double x = p[base + (n - 1) * stride];
			for (int m = n - 2; m >= 0; --m)
			{
				int idx = base + m * stride;
				x = p[idx] - cprime[idx] * x;
				p[idx] = x;
			}
		}
		__global__ void k_line_smooth(double* p, const double* rhs, double* cprime, const unsigned char* solid, MacGrid g, int dir_xmax, int axis, int color, int nlines)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= nlines) return;
			int t1 = (axis == 0) ? t % g.ny : t % g.nx;
			int t2 = (axis == 0) ? t / g.ny : t / g.nx;
			if (((t1 + t2) & 1) != color) return;
			line_solve_ch(p, rhs, cprime, solid, g, dir_xmax, axis, t);
		}
		__global__ void k_subgrad(double* u, double* v, double* w, const double* p, const unsigned char* solid, MacGrid g, double coef, int dir_xmax)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x;
			int nu = g.u_count(), nv = g.v_count(), nw = g.w_count();
			if (t < nu)
			{
				int i = t % (g.nx + 1), j = (t / (g.nx + 1)) % g.ny, k = t / ((g.nx + 1) * g.ny);
				if (i >= 1 && i <= g.nx - 1)
				{
					if (!ch_is_solid(solid, g, i - 1, j, k) && !ch_is_solid(solid, g, i, j, k))
						u[t] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i - 1, j, k)]) / g.dxc(i);
				}
				else if (i == g.nx && dir_xmax > 0 && !ch_is_solid(solid, g, g.nx - 1, j, k))
					u[t] -= coef * (0.0 - p[g.pidx(g.nx - 1, j, k)]) / g.dx(g.nx - 1); // outlet ghost p=0 at xmax
				else if (i == 0 && dir_xmax < 0 && !ch_is_solid(solid, g, 0, j, k))
					u[t] -= coef * (p[g.pidx(0, j, k)] - 0.0) / g.dx(0); // reversed outlet ghost p=0 at xmin
			}
			else if (t < nu + nv)
			{
				int s = t - nu; int i = s % g.nx, j = (s / g.nx) % (g.ny + 1), k = s / (g.nx * (g.ny + 1));
				if (j >= 1 && j <= g.ny - 1 && !ch_is_solid(solid, g, i, j - 1, k) && !ch_is_solid(solid, g, i, j, k))
					v[s] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i, j - 1, k)]) / g.dyc(j);
			}
			else if (t < nu + nv + nw)
			{
				int s = t - nu - nv; int i = s % g.nx, j = (s / g.nx) % g.ny, k = s / (g.nx * g.ny);
				if (k >= 1 && k <= g.nz - 1 && !ch_is_solid(solid, g, i, j, k - 1) && !ch_is_solid(solid, g, i, j, k))
					w[s] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i, j, k - 1)]) / g.dzc(k);
			}
		}
		__global__ void k_divabs(const double* u, const double* v, const double* w, double* out, const unsigned char* solid, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			out[t] = ch_is_solid(solid, g, i, j, k) ? 0.0 : fabs(div_cell(u, v, w, g, i, j, k));
		}
	} // anonymous namespace

	// ================= public launchers ==================================
	void ch_advect_gpu(const double* uIn, const double* vIn, const double* wIn,
		double* uOut, double* vOut, double* wOut, double* scratch,
		const unsigned char* solid, const unsigned char* nearsolid, MacGrid g, ChannelBC bc, double dt)
	{
		advect_one_gpu(0, uIn, uIn, vIn, wIn, uOut, scratch, solid, nearsolid, g, bc, dt);
		advect_one_gpu(1, vIn, uIn, vIn, wIn, vOut, scratch, solid, nearsolid, g, bc, dt);
		advect_one_gpu(2, wIn, uIn, vIn, wIn, wOut, scratch, solid, nearsolid, g, bc, dt);
	}
	void ch_smagorinsky_gpu(const double* u, const double* v, const double* w, double* nut, const unsigned char* solid, MacGrid g, ChannelBC bc, double Cs)
	{ int total = g.p_count(); k_nut<<<gsz(total), 256>>>(u, v, w, nut, solid, g, bc, Cs, total); }
	void ch_diffuse_gpu(const double* uIn, const double* vIn, const double* wIn,
		double* uOut, double* vOut, double* wOut, const double* nut, const unsigned char* solid, MacGrid g, ChannelBC bc, double dt, double nu)
	{
		diffuse_one_gpu(0, uIn, nut, solid, uOut, g, bc, dt, nu);
		diffuse_one_gpu(1, vIn, nut, solid, vOut, g, bc, dt, nu);
		diffuse_one_gpu(2, wIn, nut, solid, wOut, g, bc, dt, nu);
	}
	void ch_apply_inlet_gpu(double* u, MacGrid g, ChannelBC bc)
	{ int njk = g.ny * g.nz; k_inlet<<<gsz(njk), 256>>>(u, g, bc, njk); }
	void ch_orlanski_gpu(double* u, MacGrid g, ChannelBC bc, double dt)
	{ int njk = g.ny * g.nz; k_orlanski<<<gsz(njk), 256>>>(u, g, bc, bc.Uc * dt, njk); } // coef=Uc·dt/dx(outlet) computed in-kernel (device metric deref; graded-correct)
	void ch_apply_solid_bc_gpu(double* u, double* v, double* w, const unsigned char* solid, MacGrid g)
	{ int n = g.u_count() + g.v_count() + g.w_count(); k_solidbc<<<gsz(n), 256>>>(u, v, w, solid, g); }
	double ch_uplane_flux_gpu(const double* u, double* planeScratch, MacGrid g, int i_plane)
	{
		int njk = g.ny * g.nz;
		k_gather_uplane<<<gsz(njk), 256>>>(u, planeScratch, g, i_plane, njk);
		return reduce_sum_gpu(planeScratch, njk); // face-area (dy·dz) folded into the gather; uniform ⇒ Σu·h² as before
	}
	void ch_scale_uplane_gpu(double* u, MacGrid g, int i_plane, double s)
	{ int njk = g.ny * g.nz; k_scale_uplane<<<gsz(njk), 256>>>(u, g, i_plane, s, njk); }

	void ch_poisson_rhs_gpu(const double* u, const double* v, const double* w, double* rhs, const unsigned char* solid, MacGrid g, double rho, double dt)
	{ int n = g.p_count(); k_rhs<<<gsz(n), 256>>>(u, v, w, rhs, solid, g, rho, dt, n); }
	void ch_poisson_apply_gpu(const double* p, double* Ap, const unsigned char* solid, MacGrid g, int dir_xmax)
	{ int n = g.p_count(); k_apply<<<gsz(n), 256>>>(p, Ap, solid, g, dir_xmax, n); }
	void ch_poisson_residual_gpu(const double* p, const double* rhs, double* r, const unsigned char* solid, MacGrid g, int dir_xmax)
	{ int n = g.p_count(); k_residual<<<gsz(n), 256>>>(p, rhs, r, solid, g, dir_xmax, n); }
	void ch_jacobi_gpu(double* p, const double* rhs, double* scratch, const unsigned char* solid, MacGrid g, int dir_xmax, double omega, int sweeps)
	{
		int n = g.p_count(); double* src = p; double* dst = scratch;
		for (int s = 0; s < sweeps; ++s) { k_jacobi<<<gsz(n), 256>>>(src, rhs, dst, solid, g, dir_xmax, omega, n); double* tmp = src; src = dst; dst = tmp; }
		if (src != p) cudaMemcpy(p, src, sizeof(double) * n, cudaMemcpyDeviceToDevice);
	}
	void ch_gs_band_gpu(double* p, const double* rhs, const unsigned char* solid, MacGrid g, int dir_xmax, int band, int sweeps, bool forward)
	{
		int n = g.p_count();
		for (int s = 0; s < sweeps; ++s)
		{
			int c0 = forward ? 0 : 1, c1 = forward ? 1 : 0;
			k_gs<<<gsz(n), 256>>>(p, rhs, solid, g, dir_xmax, band, c0, n);
			k_gs<<<gsz(n), 256>>>(p, rhs, solid, g, dir_xmax, band, c1, n);
		}
	}
	void ch_line_smooth_gpu(double* p, const double* rhs, double* cprime, const unsigned char* solid, MacGrid g, int dir_xmax, int axis, int color)
	{
		int nlines = axis == 0 ? g.ny * g.nz : axis == 1 ? g.nx * g.nz : g.nx * g.ny;
		k_line_smooth<<<gsz(nlines), 256>>>(p, rhs, cprime, solid, g, dir_xmax, axis, color, nlines);
	}
	void ch_subtract_gradient_gpu(double* u, double* v, double* w, const double* p, const unsigned char* solid, MacGrid g, ChannelBC bc, double rho, double dt, int dir_xmax)
	{ int n = g.u_count() + g.v_count() + g.w_count(); k_subgrad<<<gsz(n), 256>>>(u, v, w, p, solid, g, dt / rho, dir_xmax); }
	double ch_max_div_gpu(const double* u, const double* v, const double* w, double* divscratch, const unsigned char* solid, MacGrid g)
	{ int n = g.p_count(); k_divabs<<<gsz(n), 256>>>(u, v, w, divscratch, solid, g, n); return reduce_max_abs_gpu(divscratch, n); }

	// ================= CPU references ====================================
	void ch_advect_cpu(const std::vector<double>& uIn, const std::vector<double>& vIn, const std::vector<double>& wIn,
		std::vector<double>& uOut, std::vector<double>& vOut, std::vector<double>& wOut,
		const std::vector<unsigned char>& solid, const std::vector<unsigned char>& nearsolid, MacGrid g, ChannelBC bc, double dt)
	{
		auto one = [&](int comp, const std::vector<double>& field, std::vector<double>& out) {
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk); int total = ni * nj * nk;
			std::vector<double> phiHat(total);
			for (int k = 0; k < nk; ++k) for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i)
				forward_node(comp, field.data(), uIn.data(), vIn.data(), wIn.data(), phiHat.data(), solid.data(), g, bc, dt, i, j, k);
			out.assign(total, 0.0);
			for (int k = 0; k < nk; ++k) for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i)
				correct_node(comp, field.data(), phiHat.data(), uIn.data(), vIn.data(), wIn.data(), out.data(), solid.data(), nearsolid.data(), g, bc, dt, i, j, k);
		};
		one(0, uIn, uOut); one(1, vIn, vOut); one(2, wIn, wOut);
	}
	void ch_smagorinsky_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		std::vector<double>& nut, const std::vector<unsigned char>& solid, MacGrid g, ChannelBC bc, double Cs)
	{
		nut.assign(g.p_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			nut[g.pidx(i, j, k)] = smag_nut(u.data(), v.data(), w.data(), solid.data(), g, bc, Cs, i, j, k);
	}
	void ch_diffuse_cpu(const std::vector<double>& uIn, const std::vector<double>& vIn, const std::vector<double>& wIn,
		std::vector<double>& uOut, std::vector<double>& vOut, std::vector<double>& wOut,
		const std::vector<double>* nut, const std::vector<unsigned char>& solid, MacGrid g, ChannelBC bc, double dt, double nu)
	{
		const double* pn = nut ? nut->data() : nullptr;
		uOut.assign(g.u_count(), 0.0); vOut.assign(g.v_count(), 0.0); wOut.assign(g.w_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i)
			diffuse_node(0, uIn.data(), pn, solid.data(), uOut.data(), g, bc, dt, nu, i, j, k);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j <= g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			diffuse_node(1, vIn.data(), pn, solid.data(), vOut.data(), g, bc, dt, nu, i, j, k);
		for (int k = 0; k <= g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			diffuse_node(2, wIn.data(), pn, solid.data(), wOut.data(), g, bc, dt, nu, i, j, k);
	}
	void ch_poisson_rhs_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		std::vector<double>& rhs, const std::vector<unsigned char>& solid, MacGrid g, double rho, double dt)
	{
		rhs.assign(g.p_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			rhs[g.pidx(i, j, k)] = ch_is_solid(solid.data(), g, i, j, k) ? 0.0 : -(rho / dt) * div_cell(u.data(), v.data(), w.data(), g, i, j, k);
	}
	void ch_poisson_apply_cpu(const std::vector<double>& p, std::vector<double>& Ap, const std::vector<unsigned char>& solid, MacGrid g, int dir_xmax)
	{
		Ap.assign(g.p_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
		{
			if (ch_is_solid(solid.data(), g, i, j, k)) { Ap[g.pidx(i, j, k)] = 0.0; continue; }
			double diag, wnb; fv_stencil_ch(p.data(), solid.data(), g, dir_xmax, i, j, k, diag, wnb);
			Ap[g.pidx(i, j, k)] = diag * p[g.pidx(i, j, k)] - wnb;
		}
	}
	void ch_jacobi_cpu(std::vector<double>& p, const std::vector<double>& rhs, const std::vector<unsigned char>& solid, MacGrid g, int dir_xmax, double omega, int sweeps)
	{
		std::vector<double> pout(p.size());
		for (int s = 0; s < sweeps; ++s)
		{
			for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			{
				int c = g.pidx(i, j, k);
				if (ch_is_solid(solid.data(), g, i, j, k)) { pout[c] = 0.0; continue; }
				double diag, wnb; fv_stencil_ch(p.data(), solid.data(), g, dir_xmax, i, j, k, diag, wnb);
				if (diag == 0.0) { pout[c] = p[c]; continue; }
				/* diag already set by fv_stencil_ch */
				double Ap = diag * p[c] - wnb;
				pout[c] = p[c] + omega * (rhs[c] - Ap) / diag;
			}
			p = pout;
		}
	}
	void ch_line_smooth_cpu(std::vector<double>& p, const std::vector<double>& rhs, const std::vector<unsigned char>& solid, MacGrid g, int dir_xmax, int axis, int color)
	{
		std::vector<double> cprime(p.size(), 0.0);
		int nlines = axis == 0 ? g.ny * g.nz : axis == 1 ? g.nx * g.nz : g.nx * g.ny;
		for (int t = 0; t < nlines; ++t)
		{
			int t1 = (axis == 0) ? t % g.ny : t % g.nx;
			int t2 = (axis == 0) ? t / g.ny : t / g.nx;
			if (((t1 + t2) & 1) != color) continue;
			line_solve_ch(p.data(), rhs.data(), cprime.data(), solid.data(), g, dir_xmax, axis, t);
		}
	}
	void ch_subtract_gradient_cpu(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, const std::vector<double>& p, const std::vector<unsigned char>& solid, MacGrid g, ChannelBC bc, double rho, double dt, int dir_xmax)
	{
		double coef = dt / rho;
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i)
		{
			if (i >= 1 && i <= g.nx - 1)
			{ if (!ch_is_solid(solid.data(), g, i - 1, j, k) && !ch_is_solid(solid.data(), g, i, j, k)) u[g.uidx(i, j, k)] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i - 1, j, k)]) / g.dxc(i); }
			else if (i == g.nx && dir_xmax > 0 && !ch_is_solid(solid.data(), g, g.nx - 1, j, k))
				u[g.uidx(i, j, k)] -= coef * (0.0 - p[g.pidx(g.nx - 1, j, k)]) / g.dx(g.nx - 1);
			else if (i == 0 && dir_xmax < 0 && !ch_is_solid(solid.data(), g, 0, j, k))
				u[g.uidx(i, j, k)] -= coef * (p[g.pidx(0, j, k)] - 0.0) / g.dx(0);
		}
		for (int k = 0; k < g.nz; ++k) for (int j = 1; j <= g.ny - 1; ++j) for (int i = 0; i < g.nx; ++i)
			if (!ch_is_solid(solid.data(), g, i, j - 1, k) && !ch_is_solid(solid.data(), g, i, j, k)) v[g.vidx(i, j, k)] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i, j - 1, k)]) / g.dyc(j);
		for (int k = 1; k <= g.nz - 1; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			if (!ch_is_solid(solid.data(), g, i, j, k - 1) && !ch_is_solid(solid.data(), g, i, j, k)) w[g.widx(i, j, k)] -= coef * (p[g.pidx(i, j, k)] - p[g.pidx(i, j, k - 1)]) / g.dzc(k);
	}
	double ch_max_div_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<unsigned char>& solid, MacGrid g)
	{
		double m = 0.0;
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			if (!ch_is_solid(solid.data(), g, i, j, k)) m = std::max(m, std::fabs(div_cell(u.data(), v.data(), w.data(), g, i, j, k)));
		return m;
	}
}
