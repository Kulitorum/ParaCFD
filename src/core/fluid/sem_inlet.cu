// sem_inlet.cu — see sem_inlet.h. Jarrin (2006) SEM plane kernels + generator class.
// RESEARCH §8, research/14 §3. Physics inlines live in the header (parity test 1e-5).
#include "core/fluid/sem_inlet.h"
#include "core/fluid/plane_ops.h"

#include <cuda_runtime.h>

#include <cmath>
#include <vector>

namespace scour::core
{
	namespace
	{
		inline int gsz(int n, int b = 256) { return (n + b - 1) / b; }

		__global__ void k_uplane(const double* ex, const double* ey, const double* ez,
			const double* e1, const double* e2, const double* e3, double* up, MacGrid g, SemParams p, double pref, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int j = t % g.ny, k = t / g.ny;
			double y = (j + 0.5) * g.h, z = (k + 0.5) * g.h;
			double S1, S2, S3; sem_sums(y, z, ex, ey, ez, e1, e2, e3, p.N, p.sigma, pref, S1, S2, S3);
			double a11, a22, a31, a33; sem_cholesky(p.ustar, z, p.h_dom, p.gamma, a11, a22, a31, a33);
			up[t] = a11 * S1;
		}
		__global__ void k_vplane(const double* ex, const double* ey, const double* ez,
			const double* e1, const double* e2, const double* e3, double* vp, MacGrid g, SemParams p, double pref, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int j = t % (g.ny + 1), k = t / (g.ny + 1);
			double y = j * g.h, z = (k + 0.5) * g.h;
			double S1, S2, S3; sem_sums(y, z, ex, ey, ez, e1, e2, e3, p.N, p.sigma, pref, S1, S2, S3);
			double a11, a22, a31, a33; sem_cholesky(p.ustar, z, p.h_dom, p.gamma, a11, a22, a31, a33);
			vp[t] = a22 * S2;
		}
		__global__ void k_wplane(const double* ex, const double* ey, const double* ez,
			const double* e1, const double* e2, const double* e3, double* wp, MacGrid g, SemParams p, double pref, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int j = t % g.ny, k = t / g.ny;
			double y = (j + 0.5) * g.h, z = k * g.h;
			double S1, S2, S3; sem_sums(y, z, ex, ey, ez, e1, e2, e3, p.N, p.sigma, pref, S1, S2, S3);
			double a11, a22, a31, a33; sem_cholesky(p.ustar, z, p.h_dom, p.gamma, a11, a22, a31, a33);
			wp[t] = a31 * S1 + a33 * S3;
		}
	}

	void sem_planes_gpu(const double* ex, const double* ey, const double* ez,
		const double* e1, const double* e2, const double* e3,
		double* up, double* vp, double* wp, MacGrid g, SemParams p)
	{
		double pref = std::sqrt(sem_boxvol(p.sigma, p.Ly, p.Lz) / (p.sigma * p.sigma * p.sigma));
		int nu = g.ny * g.nz, nv = (g.ny + 1) * g.nz, nw = g.ny * (g.nz + 1);
		k_uplane<<<gsz(nu), 256>>>(ex, ey, ez, e1, e2, e3, up, g, p, pref, nu);
		k_vplane<<<gsz(nv), 256>>>(ex, ey, ez, e1, e2, e3, vp, g, p, pref, nv);
		k_wplane<<<gsz(nw), 256>>>(ex, ey, ez, e1, e2, e3, wp, g, p, pref, nw);
	}

	void sem_planes_cpu(const std::vector<double>& ex, const std::vector<double>& ey, const std::vector<double>& ez,
		const std::vector<double>& e1, const std::vector<double>& e2, const std::vector<double>& e3,
		std::vector<double>& up, std::vector<double>& vp, std::vector<double>& wp, MacGrid g, SemParams p)
	{
		double pref = std::sqrt(sem_boxvol(p.sigma, p.Ly, p.Lz) / (p.sigma * p.sigma * p.sigma));
		up.assign((size_t)g.ny * g.nz, 0.0);
		vp.assign((size_t)(g.ny + 1) * g.nz, 0.0);
		wp.assign((size_t)g.ny * (g.nz + 1), 0.0);
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j)
		{
			double y = (j + 0.5) * g.h, z = (k + 0.5) * g.h;
			double S1, S2, S3; sem_sums(y, z, ex.data(), ey.data(), ez.data(), e1.data(), e2.data(), e3.data(), p.N, p.sigma, pref, S1, S2, S3);
			double a11, a22, a31, a33; sem_cholesky(p.ustar, z, p.h_dom, p.gamma, a11, a22, a31, a33);
			up[k * g.ny + j] = a11 * S1;
		}
		for (int k = 0; k < g.nz; ++k) for (int j = 0; j <= g.ny; ++j)
		{
			double y = j * g.h, z = (k + 0.5) * g.h;
			double S1, S2, S3; sem_sums(y, z, ex.data(), ey.data(), ez.data(), e1.data(), e2.data(), e3.data(), p.N, p.sigma, pref, S1, S2, S3);
			double a11, a22, a31, a33; sem_cholesky(p.ustar, z, p.h_dom, p.gamma, a11, a22, a31, a33);
			vp[k * (g.ny + 1) + j] = a22 * S2;
		}
		for (int k = 0; k <= g.nz; ++k) for (int j = 0; j < g.ny; ++j)
		{
			double y = (j + 0.5) * g.h, z = k * g.h;
			double S1, S2, S3; sem_sums(y, z, ex.data(), ey.data(), ez.data(), e1.data(), e2.data(), e3.data(), p.N, p.sigma, pref, S1, S2, S3);
			double a11, a22, a31, a33; sem_cholesky(p.ustar, z, p.h_dom, p.gamma, a11, a22, a31, a33);
			wp[k * g.ny + j] = a31 * S1 + a33 * S3;
		}
	}

	// ---- generator class -------------------------------------------------
	namespace
	{
		double* dalloc(int n) { double* p = nullptr; cudaMalloc(&p, sizeof(double) * (size_t)n); return p; }
	}

	SemInlet::SemInlet(MacGrid g, SemParams p) : g_(g), p_(p), rng_(p.seed ? p.seed : 1u)
	{
		hex_.resize(p_.N); hey_.resize(p_.N); hez_.resize(p_.N);
		he1_.resize(p_.N); he2_.resize(p_.N); he3_.resize(p_.N);
		for (int k = 0; k < p_.N; ++k) regenerate(k, false); // uniform across the whole box
		ex_ = dalloc(p_.N); ey_ = dalloc(p_.N); ez_ = dalloc(p_.N);
		e1_ = dalloc(p_.N); e2_ = dalloc(p_.N); e3_ = dalloc(p_.N);
		up_ = dalloc(g_.ny * g_.nz); vp_ = dalloc((g_.ny + 1) * g_.nz); wp_ = dalloc(g_.ny * (g_.nz + 1));
		hup_.assign((size_t)g_.ny * g_.nz, 0.0);
	}

	SemInlet::~SemInlet()
	{
		for (double* p : {ex_, ey_, ez_, e1_, e2_, e3_, up_, vp_, wp_}) cudaFree(p);
	}

	// LCG uniform in [0,1) — deterministic eddy placement (quality non-critical, RESEARCH
	// determinism note: fixed seed).
	static inline double lcg01(unsigned& s) { s = s * 1664525u + 1013904223u; return ((s >> 8) & 0xFFFFFFu) * (1.0 / 16777216.0); }

	void SemInlet::regenerate(int k, bool upstream)
	{
		double sig = p_.sigma;
		hex_[k] = upstream ? -sig : (-sig + 2.0 * sig * lcg01(rng_));
		hey_[k] = -sig + (p_.Ly + 2.0 * sig) * lcg01(rng_);
		hez_[k] = -sig + (p_.Lz + 2.0 * sig) * lcg01(rng_);
		he1_[k] = lcg01(rng_) < 0.5 ? -1.0 : 1.0;
		he2_[k] = lcg01(rng_) < 0.5 ? -1.0 : 1.0;
		he3_[k] = lcg01(rng_) < 0.5 ? -1.0 : 1.0;
	}

	void SemInlet::advance(double dt)
	{
		double Uc = p_.U_d, sig = p_.sigma;
		for (int k = 0; k < p_.N; ++k)
		{
			hex_[k] += Uc * dt;
			if (hex_[k] > sig) regenerate(k, true); // re-inject at the upstream face
		}
		cudaMemcpy(ex_, hex_.data(), sizeof(double) * p_.N, cudaMemcpyHostToDevice);
		cudaMemcpy(ey_, hey_.data(), sizeof(double) * p_.N, cudaMemcpyHostToDevice);
		cudaMemcpy(ez_, hez_.data(), sizeof(double) * p_.N, cudaMemcpyHostToDevice);
		cudaMemcpy(e1_, he1_.data(), sizeof(double) * p_.N, cudaMemcpyHostToDevice);
		cudaMemcpy(e2_, he2_.data(), sizeof(double) * p_.N, cudaMemcpyHostToDevice);
		cudaMemcpy(e3_, he3_.data(), sizeof(double) * p_.N, cudaMemcpyHostToDevice);
		sem_planes_gpu(ex_, ey_, ez_, e1_, e2_, e3_, up_, vp_, wp_, g_, p_);
		cudaMemcpy(hup_.data(), up_, sizeof(double) * g_.ny * g_.nz, cudaMemcpyDeviceToHost);
	}

	void SemInlet::add_u(double* u, MacGrid g) { plane_add_u_gpu(u, up_, g); }
	void SemInlet::set_vw(double* v, double* w, MacGrid g) { plane_set_vw_gpu(v, w, vp_, wp_, g); }
}
