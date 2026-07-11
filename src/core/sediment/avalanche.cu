// avalanche.cu — CUDA kernel + CPU reference twin for the 8-neighbour sand-slide sweep
// (RESEARCH §7; research/03 §8, /13 §7). Shared arithmetic in the anonymous-namespace SCOUR_HD
// node function (GPU-vs-CPU parity at rel. max-norm 1e-5, CLAUDE.md). See avalanche.h.
#include "core/sediment/avalanche.h"
#include "core/fluid/mac_grid.h"

#include <cuda_runtime.h>

#include <cmath>
#include <vector>

namespace scour::core
{
	namespace
	{
		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }

		SCOUR_HD inline bool nbr(MacGrid g, int periodic, int i, int j, int di, int dj, int& oi, int& oj)
		{
			oi = i + di; oj = j + dj;
			if (periodic)
			{
				if (oi < 0) oi += g.nx; else if (oi > g.nx - 1) oi -= g.nx;
				if (oj < 0) oj += g.ny; else if (oj > g.ny - 1) oj -= g.ny;
				return true;
			}
			return !(oi < 0 || oi > g.nx - 1 || oj < 0 || oj > g.ny - 1);
		}

		// Net antisymmetric transfer at column (i,j) over its 8 neighbours (Jacobi, from z_in).
		SCOUR_HD inline double avalanche_dz_node(const double* z, MacGrid g, AvalancheParams ap, int i, int j)
		{
			double zc = z[j * g.nx + i];
			double dz = 0.0;
			double sq2 = 1.4142135623730951;
			for (int dj = -1; dj <= 1; ++dj)
				for (int di = -1; di <= 1; ++di)
				{
					if (di == 0 && dj == 0) continue;
					int oi, oj;
					if (!nbr(g, ap.periodic, i, j, di, dj, oi, oj)) continue;
					double zn = z[oj * g.nx + oi];
					double L = (di != 0 && dj != 0) ? g.h * sq2 : g.h;
					double thr = L * ap.tan_repose;
					double diff = zc - zn;
					if (diff > thr) dz -= ap.relax * (diff - thr);      // self higher → lose to nb
					else if (-diff > thr) dz += ap.relax * (-diff - thr); // nb higher → gain from nb
				}
			return dz;
		}

		__global__ void k_sweep(const double* z_in, double* z_out, MacGrid g, AvalancheParams ap, int ncol)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= ncol) return;
			int i = t % g.nx, j = t / g.nx;
			z_out[t] = z_in[t] + avalanche_dz_node(z_in, g, ap, i, j);
		}
	}

	void avalanche_sweep_gpu(const double* z_in, double* z_out, MacGrid g, AvalancheParams ap)
	{
		int ncol = g.nx * g.ny; k_sweep<<<gsz(ncol), 256>>>(z_in, z_out, g, ap, ncol);
	}
	void avalanche_sweep_cpu(const std::vector<double>& z_in, std::vector<double>& z_out, MacGrid g, AvalancheParams ap)
	{
		int ncol = g.nx * g.ny; z_out.assign(ncol, 0.0);
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			z_out[j * g.nx + i] = z_in[j * g.nx + i] + avalanche_dz_node(z_in.data(), g, ap, i, j);
	}

	double avalanche_max_tanslope(const std::vector<double>& z, MacGrid g, int periodic)
	{
		double mx = 0.0, sq2 = 1.4142135623730951;
		for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
		{
			double zc = z[j * g.nx + i];
			for (int dj = -1; dj <= 1; ++dj) for (int di = -1; di <= 1; ++di)
			{
				if (di == 0 && dj == 0) continue;
				int oi, oj;
				if (!nbr(g, periodic, i, j, di, dj, oi, oj)) continue;
				double L = (di != 0 && dj != 0) ? g.h * sq2 : g.h;
				double s = std::fabs(zc - z[oj * g.nx + oi]) / L;
				if (s > mx) mx = s;
			}
		}
		return mx;
	}
}
