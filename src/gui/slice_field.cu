// slice_field.cu — see slice_field.h. Pure CUDA (no Qt/GL): one __host__ __device__
// sampler+colour map shared by the GPU kernel and the CPU reference, so the GPU-vs-CPU
// parity test only validates launch/index/memory correctness.
#include "gui/slice_field.h"

#include "gui/colormap.h"

#include <algorithm>
#include <cfloat>
#include <cmath>

namespace paracfd::gui
{
	using paracfd::core::MacGrid;

	namespace
	{
		PARACFD_HD inline float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }
		PARACFD_HD inline int clampi(int x, int lo, int hi) { return x < lo ? lo : (x > hi ? hi : x); }

		// Perceptual 5-stop gradient blue->cyan->green->yellow->red over t in [0,1]. Delegates to
		// the SHARED ramp (gui/colormap.h) so the slice, the flow arrows and the legend agree, and
		// so the GPU kernel and its CPU reference stay bit-for-bit identical (plain float lerps).
		PARACFD_HD inline float4 colormap(float t)
		{
			float4 c;
			scour_colormap(t, c.x, c.y, c.z);
			c.w = 1.0f;
			return c;
		}

		// Cell-centred scalar of the chosen field at MAC cell (i,j,k).
		PARACFD_HD inline float sample_scalar(const double* u, const double* v, const double* w, const double* p,
			MacGrid g, Field field, int i, int j, int k)
		{
			switch (field)
			{
			case Field::VelU:
				return (float)(0.5 * (u[g.uidx(i, j, k)] + u[g.uidx(i + 1, j, k)]));
			case Field::VelV:
				return (float)(0.5 * (v[g.vidx(i, j, k)] + v[g.vidx(i, j + 1, k)]));
			case Field::VelW:
				return (float)(0.5 * (w[g.widx(i, j, k)] + w[g.widx(i, j, k + 1)]));
			case Field::Pressure:
				return p ? (float)p[g.pidx(i, j, k)] : 0.0f;
			case Field::SpeedMag:
			default:
			{
				double uc = 0.5 * (u[g.uidx(i, j, k)] + u[g.uidx(i + 1, j, k)]);
				double vc = 0.5 * (v[g.vidx(i, j, k)] + v[g.vidx(i, j + 1, k)]);
				double wc = 0.5 * (w[g.widx(i, j, k)] + w[g.widx(i, j, k + 1)]);
				return (float)sqrt(uc * uc + vc * vc + wc * wc);
			}
			}
		}

		// Cell-centred speed magnitude |(u,v,w)| at MAC cell (i,j,k). Shared by the auto-range
		// reduction's GPU kernel and CPU reference so the arrow colour scale matches exactly.
		PARACFD_HD inline float cell_speed(const double* u, const double* v, const double* w, MacGrid g, int i, int j, int k)
		{
			double uc = 0.5 * (u[g.uidx(i, j, k)] + u[g.uidx(i + 1, j, k)]);
			double vc = 0.5 * (v[g.vidx(i, j, k)] + v[g.vidx(i, j + 1, k)]);
			double wc = 0.5 * (w[g.widx(i, j, k)] + w[g.widx(i, j, k + 1)]);
			return (float)sqrt(uc * uc + vc * vc + wc * wc);
		}

		// Full per-vertex evaluation: world position -> nearest cell -> scalar -> colour.
		PARACFD_HD inline float4 eval_vertex(const double* u, const double* v, const double* w, const double* p,
			const SliceParams& sp, int a, int b)
		{
			float x, y, z;
			slice_vertex_world(sp, a, b, x, y, z);
			MacGrid g = sp.grid;
			// World → cell index via the grid's world↔index map (graded-aware; floor(x/h) uniform). On a
			// graded grid sp.grid is the DEVICE view, so grid_fx reads device metric arrays on the device.
			int i = clampi((int)floor(paracfd::core::grid_fx(g, (double)x)), 0, g.nx - 1);
			int j = clampi((int)floor(paracfd::core::grid_fy(g, (double)y)), 0, g.ny - 1);
			int k = clampi((int)floor(paracfd::core::grid_fz(g, (double)z)), 0, g.nz - 1);
			float s = sample_scalar(u, v, w, p, g, sp.field, i, j, k);
			float denom = (sp.vmax > sp.vmin) ? (sp.vmax - sp.vmin) : 1.0f;
			return colormap((s - sp.vmin) / denom);
		}

		__global__ void k_slice(const double* u, const double* v, const double* w, const double* p,
			SliceParams sp, float4* out)
		{
			int idx = blockIdx.x * blockDim.x + threadIdx.x;
			int n = sp.nu * sp.nv;
			if (idx >= n) return;
			int a = idx % sp.nu, b = idx / sp.nu;
			out[idx] = eval_vertex(u, v, w, p, sp, a, b);
		}

		// --- Auto-range reduction ------------------------------------------------------------
		// Float atomic min/max via CAS (min/max are order-independent ⇒ the reduction is
		// deterministic and bit-exact against the CPU reference regardless of block scheduling).
		__device__ inline void atomic_min_f(float* addr, float val)
		{
			int* ia = reinterpret_cast<int*>(addr);
			int old = *ia, assumed;
			do
			{
				assumed = old;
				if (__int_as_float(assumed) <= val) break;
				old = atomicCAS(ia, assumed, __float_as_int(val));
			} while (assumed != old);
		}
		__device__ inline void atomic_max_f(float* addr, float val)
		{
			int* ia = reinterpret_cast<int*>(addr);
			int old = *ia, assumed;
			do
			{
				assumed = old;
				if (__int_as_float(assumed) >= val) break;
				old = atomicCAS(ia, assumed, __float_as_int(val));
			} while (assumed != old);
		}

		__global__ void k_range_init(float* out3)
		{
			out3[0] = FLT_MAX;  // field_min
			out3[1] = -FLT_MAX; // field_max
			out3[2] = 0.0f;     // speed_max (|u| >= 0)
		}

		// One block-reduced min/max/speed-max over a grid-stride slice of the pressure-cell index
		// space (cell == pidx(i,j,k) by construction of the MAC layout), then a single atomic per value.
		__global__ void k_range_reduce(const double* u, const double* v, const double* w, const double* p,
			const unsigned char* solid, MacGrid g, Field field, float* out3)
		{
			__shared__ float smin[256];
			__shared__ float smax[256];
			__shared__ float sspd[256];
			const int tid = threadIdx.x;
			const int ncell = g.nx * g.ny * g.nz;
			const int nxy = g.nx * g.ny;
			float lmin = FLT_MAX, lmax = -FLT_MAX, lspd = 0.0f;
			for (int cell = blockIdx.x * blockDim.x + tid; cell < ncell; cell += gridDim.x * blockDim.x)
			{
				if (solid && solid[cell]) continue; // fluid cells only (obstacle/structure excluded)
				const int i = cell % g.nx;
				const int j = (cell / g.nx) % g.ny;
				const int k = cell / nxy;
				float s = sample_scalar(u, v, w, p, g, field, i, j, k);
				float sp = cell_speed(u, v, w, g, i, j, k);
				if (isfinite(s)) { lmin = fminf(lmin, s); lmax = fmaxf(lmax, s); }
				if (isfinite(sp)) lspd = fmaxf(lspd, sp);
			}
			smin[tid] = lmin; smax[tid] = lmax; sspd[tid] = lspd;
			__syncthreads();
			for (int s = blockDim.x / 2; s > 0; s >>= 1)
			{
				if (tid < s)
				{
					smin[tid] = fminf(smin[tid], smin[tid + s]);
					smax[tid] = fmaxf(smax[tid], smax[tid + s]);
					sspd[tid] = fmaxf(sspd[tid], sspd[tid + s]);
				}
				__syncthreads();
			}
			if (tid == 0)
			{
				atomic_min_f(&out3[0], smin[0]);
				atomic_max_f(&out3[1], smax[0]);
				atomic_max_f(&out3[2], sspd[0]);
			}
		}
	} // namespace

	void slice_fill_gpu(const double* u, const double* v, const double* w, const double* p,
		const SliceParams& sp, float4* out, cudaStream_t stream)
	{
		int n = sp.nu * sp.nv;
		if (n <= 0) return;
		int block = 256, grid = (n + block - 1) / block;
		k_slice<<<grid, block, 0, stream>>>(u, v, w, p, sp, out);
	}

	void slice_fill_cpu(const double* u, const double* v, const double* w, const double* p,
		const SliceParams& sp, float4* out)
	{
		for (int b = 0; b < sp.nv; ++b)
			for (int a = 0; a < sp.nu; ++a)
				out[b * sp.nu + a] = eval_vertex(u, v, w, p, sp, a, b);
	}

	void slice_reduce_gpu(const double* u, const double* v, const double* w, const double* p,
		const unsigned char* solid, MacGrid g, Field field, float* out3_dev, cudaStream_t stream)
	{
		int ncell = g.p_count();
		if (ncell <= 0 || !out3_dev) return;
		k_range_init<<<1, 1, 0, stream>>>(out3_dev);
		const int block = 256;
		int grid = (ncell + block - 1) / block;
		if (grid > 1024) grid = 1024; // grid-stride covers the remainder; caps atomic contention
		k_range_reduce<<<grid, block, 0, stream>>>(u, v, w, p, solid, g, field, out3_dev);
	}

	void slice_reduce_cpu(const double* u, const double* v, const double* w, const double* p,
		const unsigned char* solid, MacGrid g, Field field, float out3[3])
	{
		float fmin = FLT_MAX, fmax = -FLT_MAX, smax = 0.0f;
		for (int k = 0; k < g.nz; ++k)
			for (int j = 0; j < g.ny; ++j)
				for (int i = 0; i < g.nx; ++i)
				{
					if (solid && solid[g.pidx(i, j, k)]) continue;
					float s = sample_scalar(u, v, w, p, g, field, i, j, k);
					float sp = cell_speed(u, v, w, g, i, j, k);
					if (std::isfinite(s)) { fmin = std::min(fmin, s); fmax = std::max(fmax, s); }
					if (std::isfinite(sp)) smax = std::max(smax, sp);
				}
		out3[0] = fmin; out3[1] = fmax; out3[2] = smax;
	}
}
