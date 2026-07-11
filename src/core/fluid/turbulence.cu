// turbulence.cu — explicit Smagorinsky eddy viscosity and explicit diffusion of the
// MAC velocity field. nu_t = (Cs*Δ)^2*|S|, |S| = sqrt(2 S_ij S_ij) at cell centres with
// per-cell filter width Δ = (dx·dy·dz)^(1/3) (graded-structured-grid change; ⇒ Cs·h on a
// uniform grid); diffusion is treated explicitly (h^2/6nu >> advective dt), Stam's implicit
// solve is dropped. RESEARCH §3.4 (verified research/08 item 6). M1 gate runs Cs=0. Strain
// gradients + the diffusion Laplacian use local metric spacings (d2_axis), collapsing exactly
// to the uniform stencil when every spacing equals h.
#include "core/fluid/mac_grid.h"
#include "core/fluid/mac_ops.h"

#include <cuda_runtime.h>
#include <vector>

namespace windcfd::core
{
	namespace
	{
		WINDCFD_HD inline double uc(const double* u, MacGrid g, BC bc, int i, int j, int k)
		{ return 0.5 * (fetch_u(u, g, bc, i, j, k) + fetch_u(u, g, bc, i + 1, j, k)); }
		WINDCFD_HD inline double vc(const double* v, MacGrid g, BC bc, int i, int j, int k)
		{ return 0.5 * (fetch_v(v, g, bc, i, j, k) + fetch_v(v, g, bc, i, j + 1, k)); }
		WINDCFD_HD inline double wc(const double* w, MacGrid g, BC bc, int i, int j, int k)
		{ return 0.5 * (fetch_w(w, g, bc, i, j, k) + fetch_w(w, g, bc, i, j, k + 1)); }

		WINDCFD_HD inline double smag_nut(const double* u, const double* v, const double* w,
			MacGrid g, BC bc, double Cs, int i, int j, int k)
		{
			// Diagonal strains over the local cell width; cross strains over the 2-cell centre span
			// (gap(i-1)+gap(i) = 2h uniform). OOB-safe via the gap accessors' mirror-ghost widths.
			double invx = 1.0 / g.dx(i), invy = 1.0 / g.dy(j), invz = 1.0 / g.dz(k);
			double sx = 1.0 / (g.gapx(i - 1) + g.gapx(i)), sy = 1.0 / (g.gapy(j - 1) + g.gapy(j)), sz = 1.0 / (g.gapz(k - 1) + g.gapz(k));
			double dudx = (fetch_u(u, g, bc, i + 1, j, k) - fetch_u(u, g, bc, i, j, k)) * invx;
			double dvdy = (fetch_v(v, g, bc, i, j + 1, k) - fetch_v(v, g, bc, i, j, k)) * invy;
			double dwdz = (fetch_w(w, g, bc, i, j, k + 1) - fetch_w(w, g, bc, i, j, k)) * invz;
			double dudy = (uc(u, g, bc, i, j + 1, k) - uc(u, g, bc, i, j - 1, k)) * sy;
			double dudz = (uc(u, g, bc, i, j, k + 1) - uc(u, g, bc, i, j, k - 1)) * sz;
			double dvdx = (vc(v, g, bc, i + 1, j, k) - vc(v, g, bc, i - 1, j, k)) * sx;
			double dvdz = (vc(v, g, bc, i, j, k + 1) - vc(v, g, bc, i, j, k - 1)) * sz;
			double dwdx = (wc(w, g, bc, i + 1, j, k) - wc(w, g, bc, i - 1, j, k)) * sx;
			double dwdy = (wc(w, g, bc, i, j + 1, k) - wc(w, g, bc, i, j - 1, k)) * sy;
			double Sxx = dudx, Syy = dvdy, Szz = dwdz;
			double Sxy = 0.5 * (dudy + dvdx), Sxz = 0.5 * (dudz + dwdx), Syz = 0.5 * (dvdz + dwdy);
			double Smag = sqrt(2.0 * (Sxx * Sxx + Syy * Syy + Szz * Szz + 2.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz)));
			double cd = Cs * cbrt(g.dx(i) * g.dy(j) * g.dz(k));
			return cd * cd * Smag;
		}

		__global__ void k_nut(const double* u, const double* v, const double* w, double* nut, MacGrid g, BC bc, double Cs, int total)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x;
			if (t >= total) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			nut[g.pidx(i, j, k)] = smag_nut(u, v, w, g, bc, Cs, i, j, k);
		}

		// diffuse one component (comp: 0=u,1=v,2=w) into out. The Laplacian is the sum of three
		// non-uniform second derivatives (d2_axis): the component's NORMAL axis is face-spaced (dx/dy/dz),
		// the two TANGENTIAL axes are centre-spaced (gap*). Uniform spacings ⇒ the (Σnb − 6c)/h² stencil.
		WINDCFD_HD inline void diffuse_node(int comp, const double* field, const double* nut,
			double* out, MacGrid g, BC bc, double dt, double nu, int i, int j, int k)
		{
			int idx;
			double val;
			if (comp == 0)
			{
				idx = g.uidx(i, j, k);
				bool interior = i >= 1 && i <= g.nx - 1 && j >= 0 && j <= g.ny - 1 && k >= 0 && k <= g.nz - 1;
				if (!interior) { out[idx] = field[idx]; return; }
				double c = field[idx];
				double lap = d2_axis(fetch_u(field, g, bc, i - 1, j, k), c, fetch_u(field, g, bc, i + 1, j, k), g.dx(i - 1), g.dx(i))
					+ d2_axis(fetch_u(field, g, bc, i, j - 1, k), c, fetch_u(field, g, bc, i, j + 1, k), g.gapy(j - 1), g.gapy(j))
					+ d2_axis(fetch_u(field, g, bc, i, j, k - 1), c, fetch_u(field, g, bc, i, j, k + 1), g.gapz(k - 1), g.gapz(k));
				double nut_f = nut ? 0.5 * (nut[g.pidx(i - 1, j, k)] + nut[g.pidx(i, j, k)]) : 0.0;
				val = c + dt * (nu + nut_f) * lap;
			}
			else if (comp == 1)
			{
				idx = g.vidx(i, j, k);
				bool interior = i >= 0 && i <= g.nx - 1 && j >= 1 && j <= g.ny - 1 && k >= 0 && k <= g.nz - 1;
				if (!interior) { out[idx] = field[idx]; return; }
				double c = field[idx];
				double lap = d2_axis(fetch_v(field, g, bc, i - 1, j, k), c, fetch_v(field, g, bc, i + 1, j, k), g.gapx(i - 1), g.gapx(i))
					+ d2_axis(fetch_v(field, g, bc, i, j - 1, k), c, fetch_v(field, g, bc, i, j + 1, k), g.dy(j - 1), g.dy(j))
					+ d2_axis(fetch_v(field, g, bc, i, j, k - 1), c, fetch_v(field, g, bc, i, j, k + 1), g.gapz(k - 1), g.gapz(k));
				double nut_f = nut ? 0.5 * (nut[g.pidx(i, j - 1, k)] + nut[g.pidx(i, j, k)]) : 0.0;
				val = c + dt * (nu + nut_f) * lap;
			}
			else
			{
				idx = g.widx(i, j, k);
				bool interior = i >= 0 && i <= g.nx - 1 && j >= 0 && j <= g.ny - 1 && k >= 1 && k <= g.nz - 1;
				if (!interior) { out[idx] = field[idx]; return; }
				double c = field[idx];
				double lap = d2_axis(fetch_w(field, g, bc, i - 1, j, k), c, fetch_w(field, g, bc, i + 1, j, k), g.gapx(i - 1), g.gapx(i))
					+ d2_axis(fetch_w(field, g, bc, i, j - 1, k), c, fetch_w(field, g, bc, i, j + 1, k), g.gapy(j - 1), g.gapy(j))
					+ d2_axis(fetch_w(field, g, bc, i, j, k - 1), c, fetch_w(field, g, bc, i, j, k + 1), g.dz(k - 1), g.dz(k));
				double nut_f = nut ? 0.5 * (nut[g.pidx(i, j, k - 1)] + nut[g.pidx(i, j, k)]) : 0.0;
				val = c + dt * (nu + nut_f) * lap;
			}
			out[idx] = val;
		}

		__global__ void k_diffuse(int comp, const double* field, const double* nut, double* out,
			MacGrid g, BC bc, double dt, double nu, int ni, int nj, int nk)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x;
			int total = ni * nj * nk;
			if (t >= total) return;
			int i = t % ni, j = (t / ni) % nj, k = t / (ni * nj);
			diffuse_node(comp, field, nut, out, g, bc, dt, nu, i, j, k);
		}

		void diffuse_one_gpu(int comp, const double* field, const double* nut, double* out, MacGrid g, BC bc, double dt, double nu)
		{
			int ni = (comp == 0) ? g.nx + 1 : g.nx;
			int nj = (comp == 1) ? g.ny + 1 : g.ny;
			int nk = (comp == 2) ? g.nz + 1 : g.nz;
			int total = ni * nj * nk, block = 256, grid = (total + block - 1) / block;
			k_diffuse<<<grid, block>>>(comp, field, nut, out, g, bc, dt, nu, ni, nj, nk);
		}
	}

	void smagorinsky_nut_gpu(const double* u, const double* v, const double* w, double* nut, MacGrid g, BC bc, double Cs)
	{
		int total = g.p_count(), block = 256, grid = (total + block - 1) / block;
		k_nut<<<grid, block>>>(u, v, w, nut, g, bc, Cs, total);
	}

	void smagorinsky_nut_cpu(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		std::vector<double>& nut, MacGrid g, BC bc, double Cs)
	{
		nut.assign(g.p_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			nut[g.pidx(i, j, k)] = smag_nut(u.data(), v.data(), w.data(), g, bc, Cs, i, j, k);
	}

	void mac_diffuse_gpu(const double* uIn, const double* vIn, const double* wIn,
		double* uOut, double* vOut, double* wOut, const double* nut,
		MacGrid g, BC bc, double dt, double nu)
	{
		diffuse_one_gpu(0, uIn, nut, uOut, g, bc, dt, nu);
		diffuse_one_gpu(1, vIn, nut, vOut, g, bc, dt, nu);
		diffuse_one_gpu(2, wIn, nut, wOut, g, bc, dt, nu);
	}

	void mac_diffuse_cpu(const std::vector<double>& uIn, const std::vector<double>& vIn, const std::vector<double>& wIn,
		std::vector<double>& uOut, std::vector<double>& vOut, std::vector<double>& wOut, const std::vector<double>* nut,
		MacGrid g, BC bc, double dt, double nu)
	{
		const double* pn = nut ? nut->data() : nullptr;
		uOut.assign(g.u_count(), 0.0); vOut.assign(g.v_count(), 0.0); wOut.assign(g.w_count(), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i)
			diffuse_node(0, uIn.data(), pn, uOut.data(), g, bc, dt, nu, i, j, k);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j <= g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			diffuse_node(1, vIn.data(), pn, vOut.data(), g, bc, dt, nu, i, j, k);
		for (int k = 0; k <= g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			diffuse_node(2, wIn.data(), pn, wOut.data(), g, bc, dt, nu, i, j, k);
	}
}
