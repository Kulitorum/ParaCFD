// channel_porous.cu — see channel_porous.h. Kernel + CPU twin for the sub-grid porous momentum sink.
#include "core/fluid/channel_porous.h"

#include <cuda_runtime.h>

namespace paracfd::core
{
	namespace
	{
		// k for a face = max porous-k of the two cells it separates (a face gets a screen if either
		// neighbour is porous). Out-of-domain neighbours contribute 0 (open).
		PARACFD_HD inline double kface_u(const double* pk, MacGrid g, int i, int j, int k)
		{
			double a = (i - 1 >= 0 && i - 1 < g.nx) ? pk[g.pidx(i - 1, j, k)] : 0.0;
			double b = (i >= 0 && i < g.nx) ? pk[g.pidx(i, j, k)] : 0.0;
			return a > b ? a : b;
		}
		PARACFD_HD inline double kface_v(const double* pk, MacGrid g, int i, int j, int k)
		{
			double a = (j - 1 >= 0 && j - 1 < g.ny) ? pk[g.pidx(i, j - 1, k)] : 0.0;
			double b = (j >= 0 && j < g.ny) ? pk[g.pidx(i, j, k)] : 0.0;
			return a > b ? a : b;
		}
		PARACFD_HD inline double kface_w(const double* pk, MacGrid g, int i, int j, int k)
		{
			double a = (k - 1 >= 0 && k - 1 < g.nz) ? pk[g.pidx(i, j, k - 1)] : 0.0;
			double b = (k >= 0 && k < g.nz) ? pk[g.pidx(i, j, k)] : 0.0;
			return a > b ? a : b;
		}

		__global__ void k_porous_u(double* u, const double* pk, MacGrid g, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int nxp = g.nx + 1;
			int i = t % nxp, j = (t / nxp) % g.ny, k = t / (nxp * g.ny);
			u[t] = porous_face_drag(u[t], kface_u(pk, g, i, j, k), dt, g.h);
		}
		__global__ void k_porous_v(double* v, const double* pk, MacGrid g, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int nyp = g.ny + 1;
			int i = t % g.nx, j = (t / g.nx) % nyp, k = t / (g.nx * nyp);
			v[t] = porous_face_drag(v[t], kface_v(pk, g, i, j, k), dt, g.h);
		}
		__global__ void k_porous_w(double* w, const double* pk, MacGrid g, double dt, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			w[t] = porous_face_drag(w[t], kface_w(pk, g, i, j, k), dt, g.h);
		}
	}

	void ch_porous_drag_gpu(double* u, double* v, double* w, const double* porous_k, MacGrid g, double dt)
	{
		if (!porous_k) return;
		int nu = g.u_count(), nv = g.v_count(), nw = g.w_count();
		k_porous_u<<<(nu + 255) / 256, 256>>>(u, porous_k, g, dt, nu);
		k_porous_v<<<(nv + 255) / 256, 256>>>(v, porous_k, g, dt, nv);
		k_porous_w<<<(nw + 255) / 256, 256>>>(w, porous_k, g, dt, nw);
	}

	void ch_porous_drag_cpu(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w,
		const std::vector<double>& porous_k, MacGrid g, double dt)
	{
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i)
			u[g.uidx(i, j, k)] = porous_face_drag(u[g.uidx(i, j, k)], kface_u(porous_k.data(), g, i, j, k), dt, g.h);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j <= g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			v[g.vidx(i, j, k)] = porous_face_drag(v[g.vidx(i, j, k)], kface_v(porous_k.data(), g, i, j, k), dt, g.h);
		for (int k = 0; k <= g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			w[g.widx(i, j, k)] = porous_face_drag(w[g.widx(i, j, k)], kface_w(porous_k.data(), g, i, j, k), dt, g.h);
	}
}
