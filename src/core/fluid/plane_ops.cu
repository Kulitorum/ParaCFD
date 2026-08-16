// plane_ops.cu — see plane_ops.h.
#include "core/fluid/plane_ops.h"

#include <cuda_runtime.h>

namespace paracfd::core
{
	namespace
	{
		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }
		__global__ void k_add_u(double* u, const double* up, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int j = t % g.ny, k = t / g.ny;
			u[g.uidx(0, j, k)] += up[t];
		}
		__global__ void k_set_v(double* v, const double* vp, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int j = t % (g.ny + 1), k = t / (g.ny + 1);
			v[g.vidx(0, j, k)] = vp[t];
		}
		__global__ void k_set_w(double* w, const double* wp, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int j = t % g.ny, k = t / g.ny;
			w[g.widx(0, j, k)] = wp[t];
		}
	}

	void plane_add_u_gpu(double* u, const double* up, MacGrid g)
	{
		int n = g.ny * g.nz; k_add_u<<<gsz(n), 256>>>(u, up, g, n);
	}
	void plane_set_vw_gpu(double* v, double* w, const double* vp, const double* wp, MacGrid g)
	{
		int nv = (g.ny + 1) * g.nz, nw = g.ny * (g.nz + 1);
		k_set_v<<<gsz(nv), 256>>>(v, vp, g, nv);
		k_set_w<<<gsz(nw), 256>>>(w, wp, g, nw);
	}
}
