// advect.cu — MacCormack semi-Lagrangian advection of the MAC velocity field.
// Two SL sweeps: phi^{n+1} = phi_hat^{n+1} + 0.5*(phi^n - phi_hat^n), clamped to the
// 8-corner min/max of the forward interpolation, reverting to 1st order within `band`
// cells of a wall. RK2 backtrace, rays clipped to the domain. RESEARCH §3.2 (Selle
// et al. 2008), verified in research/08.
#include "core/fluid/mac_grid.h"
#include "core/fluid/mac_ops.h"

#include <cuda_runtime.h>
#include <vector>

namespace windcfd::core
{
	namespace
	{
		WINDCFD_HD inline void node_pos(int comp, MacGrid g, int i, int j, int k, double& x, double& y, double& z)
		{
			if (comp == 0) { x = g.xf(i); y = g.yc(j); z = g.zc(k); }
			else if (comp == 1) { x = g.xc(i); y = g.yf(j); z = g.zc(k); }
			else { x = g.xc(i); y = g.yc(j); z = g.zf(k); }
		}

		WINDCFD_HD inline double sample_comp(int comp, const double* f, MacGrid g, BC bc, double x, double y, double z, double* mn, double* mx)
		{
			if (comp == 0) return trilerp_u(f, g, bc, x, y, z, mn, mx);
			if (comp == 1) return trilerp_v(f, g, bc, x, y, z, mn, mx);
			return trilerp_w(f, g, bc, x, y, z, mn, mx);
		}

		WINDCFD_HD inline void vel_at(const double* u, const double* v, const double* w, MacGrid g, BC bc,
			double x, double y, double z, double& vx, double& vy, double& vz)
		{
			vx = trilerp_u(u, g, bc, x, y, z, nullptr, nullptr);
			vy = trilerp_v(v, g, bc, x, y, z, nullptr, nullptr);
			vz = trilerp_w(w, g, bc, x, y, z, nullptr, nullptr);
		}

		// RK2 backtrace by time-step `dt` (dt>0 traces upstream/backward in time).
		WINDCFD_HD inline void rk2_trace(const double* u, const double* v, const double* w, MacGrid g, BC bc, double dt,
			double x, double y, double z, double& xo, double& yo, double& zo)
		{
			double vx, vy, vz;
			vel_at(u, v, w, g, bc, x, y, z, vx, vy, vz);
			double xm = x - 0.5 * dt * vx, ym = y - 0.5 * dt * vy, zm = z - 0.5 * dt * vz;
			clamp_to_domain(g, xm, ym, zm);
			double vmx, vmy, vmz;
			vel_at(u, v, w, g, bc, xm, ym, zm, vmx, vmy, vmz);
			xo = x - dt * vmx; yo = y - dt * vmy; zo = z - dt * vmz;
			clamp_to_domain(g, xo, yo, zo);
		}

		WINDCFD_HD inline void comp_extent(int comp, MacGrid g, int& ni, int& nj, int& nk)
		{
			if (comp == 0) { ni = g.nx + 1; nj = g.ny; nk = g.nz; }
			else if (comp == 1) { ni = g.nx; nj = g.ny + 1; nk = g.nz; }
			else { ni = g.nx; nj = g.ny; nk = g.nz + 1; }
		}

		WINDCFD_HD inline bool comp_interior(int comp, MacGrid g, int i, int j, int k)
		{
			if (comp == 0) return i >= 1 && i <= g.nx - 1 && j >= 0 && j <= g.ny - 1 && k >= 0 && k <= g.nz - 1;
			if (comp == 1) return i >= 0 && i <= g.nx - 1 && j >= 1 && j <= g.ny - 1 && k >= 0 && k <= g.nz - 1;
			return i >= 0 && i <= g.nx - 1 && j >= 0 && j <= g.ny - 1 && k >= 1 && k <= g.nz - 1;
		}

		WINDCFD_HD inline bool near_wall(int comp, MacGrid g, int i, int j, int k, int band)
		{
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk);
			return i < band || i > ni - 1 - band || j < band || j > nj - 1 - band || k < band || k > nk - 1 - band;
		}

		WINDCFD_HD inline int comp_idx(int comp, MacGrid g, int i, int j, int k)
		{
			if (comp == 0) return g.uidx(i, j, k);
			if (comp == 1) return g.vidx(i, j, k);
			return g.widx(i, j, k);
		}

		// Forward SL sweep: phiHat = A(field). Interior nodes get the backtrace value,
		// boundary nodes are copied through (keeps wall-normal faces at their BC value).
		WINDCFD_HD inline void forward_node(int comp, const double* field, const double* u, const double* v, const double* w,
			double* phiHat, MacGrid g, BC bc, double dt, int i, int j, int k)
		{
			int idx = comp_idx(comp, g, i, j, k);
			if (!comp_interior(comp, g, i, j, k)) { phiHat[idx] = field[idx]; return; }
			double x, y, z; node_pos(comp, g, i, j, k, x, y, z);
			double xb, yb, zb; rk2_trace(u, v, w, g, bc, dt, x, y, z, xb, yb, zb);
			phiHat[idx] = sample_comp(comp, field, g, bc, xb, yb, zb, nullptr, nullptr);
		}

		// Correction sweep: MacCormack combine + clamp + near-wall reversion.
		WINDCFD_HD inline void correct_node(int comp, const double* field, const double* phiHat,
			const double* u, const double* v, const double* w, double* out,
			MacGrid g, BC bc, double dt, int band, int i, int j, int k)
		{
			int idx = comp_idx(comp, g, i, j, k);
			if (!comp_interior(comp, g, i, j, k)) { out[idx] = field[idx]; return; }
			double phf = phiHat[idx];
			if (near_wall(comp, g, i, j, k, band)) { out[idx] = phf; return; } // 1st-order reversion
			double x, y, z; node_pos(comp, g, i, j, k, x, y, z);
			double xb, yb, zb; rk2_trace(u, v, w, g, bc, dt, x, y, z, xb, yb, zb);
			double mn, mx;
			sample_comp(comp, field, g, bc, xb, yb, zb, &mn, &mx); // 8-corner bounds
			double xf, yf, zf; rk2_trace(u, v, w, g, bc, -dt, x, y, z, xf, yf, zf);
			double phn = sample_comp(comp, phiHat, g, bc, xf, yf, zf, nullptr, nullptr);
			double corrected = phf + 0.5 * (field[idx] - phn);
			if (corrected < mn) corrected = mn; else if (corrected > mx) corrected = mx; // clamp
			out[idx] = corrected;
		}

		__global__ void k_forward(int comp, const double* field, const double* u, const double* v, const double* w,
			double* phiHat, MacGrid g, BC bc, double dt, int total)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x;
			if (t >= total) return;
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk);
			int i = t % ni, j = (t / ni) % nj, k = t / (ni * nj);
			forward_node(comp, field, u, v, w, phiHat, g, bc, dt, i, j, k);
		}

		__global__ void k_correct(int comp, const double* field, const double* phiHat,
			const double* u, const double* v, const double* w, double* out,
			MacGrid g, BC bc, double dt, int band, int total)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x;
			if (t >= total) return;
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk);
			int i = t % ni, j = (t / ni) % nj, k = t / (ni * nj);
			correct_node(comp, field, phiHat, u, v, w, out, g, bc, dt, band, i, j, k);
		}

		void advect_one_gpu(int comp, const double* field, const double* u, const double* v, const double* w,
			double* out, double* phiHat, MacGrid g, BC bc, double dt, int band)
		{
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk);
			int total = ni * nj * nk;
			int block = 256, grid = (total + block - 1) / block;
			k_forward<<<grid, block>>>(comp, field, u, v, w, phiHat, g, bc, dt, total);
			k_correct<<<grid, block>>>(comp, field, phiHat, u, v, w, out, g, bc, dt, band, total);
		}
	}

	void mac_advect_gpu(const double* uIn, const double* vIn, const double* wIn,
		double* uOut, double* vOut, double* wOut, double* scratch,
		MacGrid g, BC bc, double dt, int band)
	{
		advect_one_gpu(0, uIn, uIn, vIn, wIn, uOut, scratch, g, bc, dt, band);
		advect_one_gpu(1, vIn, uIn, vIn, wIn, vOut, scratch, g, bc, dt, band);
		advect_one_gpu(2, wIn, uIn, vIn, wIn, wOut, scratch, g, bc, dt, band);
	}

	// ---- CPU reference ---------------------------------------------------------
	namespace
	{
		void advect_one_cpu(int comp, const std::vector<double>& field,
			const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
			std::vector<double>& out, MacGrid g, BC bc, double dt, int band)
		{
			int ni, nj, nk; comp_extent(comp, g, ni, nj, nk);
			int total = ni * nj * nk;
			std::vector<double> phiHat(total);
			for (int k = 0; k < nk; ++k) for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i)
				forward_node(comp, field.data(), u.data(), v.data(), w.data(), phiHat.data(), g, bc, dt, i, j, k);
			out.assign(total, 0.0);
			for (int k = 0; k < nk; ++k) for (int j = 0; j < nj; ++j) for (int i = 0; i < ni; ++i)
				correct_node(comp, field.data(), phiHat.data(), u.data(), v.data(), w.data(), out.data(), g, bc, dt, band, i, j, k);
		}
	}

	void mac_advect_cpu(const std::vector<double>& uIn, const std::vector<double>& vIn, const std::vector<double>& wIn,
		std::vector<double>& uOut, std::vector<double>& vOut, std::vector<double>& wOut,
		MacGrid g, BC bc, double dt, int band)
	{
		advect_one_cpu(0, uIn, uIn, vIn, wIn, uOut, g, bc, dt, band);
		advect_one_cpu(1, vIn, uIn, vIn, wIn, vOut, g, bc, dt, band);
		advect_one_cpu(2, wIn, uIn, vIn, wIn, wOut, g, bc, dt, band);
	}
}
