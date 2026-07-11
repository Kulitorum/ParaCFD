// cylinder.cu — see cylinder.h. Host driver + a cell-centre time-mean accumulation
// kernel + a host radix-2 FFT for the Strouhal number and a control-volume momentum
// balance for Cd (RESEARCH §10 V3, research/11 §4).
#include "core/fluid/cylinder.h"
#include "core/fluid/channel_core.h"
#include "core/fluid/channel_mask.h"

#include <cuda_runtime.h>

#include <cmath>
#include <cstdio>
#include <vector>

namespace scour::core
{
	namespace
	{
		constexpr double PI = 3.14159265358979323846;

		__global__ void k_accum(const double* u, const double* v, const double* p,
			double* Su, double* Sv, double* Sp, double* Suu, double* Suv, MacGrid g, int n)
		{
			int t = blockIdx.x * blockDim.x + threadIdx.x; if (t >= n) return;
			int i = t % g.nx, j = (t / g.nx) % g.ny, k = t / (g.nx * g.ny);
			double uc = 0.5 * (u[g.uidx(i, j, k)] + u[g.uidx(i + 1, j, k)]);
			double vc = 0.5 * (v[g.vidx(i, j, k)] + v[g.vidx(i, j + 1, k)]);
			double pc = p[g.pidx(i, j, k)];
			Su[t] += uc; Sv[t] += vc; Sp[t] += pc; Suu[t] += uc * uc; Suv[t] += uc * vc;
		}

		// In-place iterative radix-2 Cooley-Tukey FFT (n a power of two).
		void fft(std::vector<double>& re, std::vector<double>& im)
		{
			int n = (int)re.size();
			for (int i = 1, j = 0; i < n; ++i)
			{
				int bit = n >> 1;
				for (; j & bit; bit >>= 1) j ^= bit;
				j ^= bit;
				if (i < j) { std::swap(re[i], re[j]); std::swap(im[i], im[j]); }
			}
			for (int len = 2; len <= n; len <<= 1)
			{
				double ang = -2.0 * PI / len;
				double wr = std::cos(ang), wi = std::sin(ang);
				for (int i = 0; i < n; i += len)
				{
					double cwr = 1.0, cwi = 0.0;
					for (int k = 0; k < len / 2; ++k)
					{
						double ur = re[i + k], ui = im[i + k];
						double vr = re[i + k + len / 2] * cwr - im[i + k + len / 2] * cwi;
						double vi = re[i + k + len / 2] * cwi + im[i + k + len / 2] * cwr;
						re[i + k] = ur + vr; im[i + k] = ui + vi;
						re[i + k + len / 2] = ur - vr; im[i + k + len / 2] = ui - vi;
						double ncwr = cwr * wr - cwi * wi, ncwi = cwr * wi + cwi * wr;
						cwr = ncwr; cwi = ncwi;
					}
				}
			}
		}

		// Dominant frequency of a uniformly-sampled real series via Hann-windowed FFT
		// with parabolic peak interpolation. Returns f in cycles per (dt-unit).
		double dominant_freq(const std::vector<double>& sig, double dt)
		{
			int n = (int)sig.size();
			double mean = 0.0; for (double x : sig) mean += x; mean /= n;
			std::vector<double> re(n), im(n, 0.0);
			for (int i = 0; i < n; ++i)
			{
				double w = 0.5 * (1.0 - std::cos(2.0 * PI * i / (n - 1))); // Hann
				re[i] = (sig[i] - mean) * w;
			}
			fft(re, im);
			int kmax = 2; double best = -1.0;
			for (int k = 2; k < n / 2; ++k)
			{
				double mag = re[k] * re[k] + im[k] * im[k];
				if (mag > best) { best = mag; kmax = k; }
			}
			// parabolic interpolation on log-magnitude
			auto lm = [&](int k) { double m = re[k] * re[k] + im[k] * im[k]; return 0.5 * std::log(m + 1e-300); };
			double y0 = lm(kmax - 1), y1 = lm(kmax), y2 = lm(kmax + 1);
			double denom = (y0 - 2 * y1 + y2);
			double delta = (denom != 0.0) ? 0.5 * (y0 - y2) / denom : 0.0;
			if (delta > 0.5) delta = 0.5; else if (delta < -0.5) delta = -0.5;
			double df = 1.0 / (n * dt);
			return (kmax + delta) * df;
		}

		double zerocross_freq(const std::vector<double>& sig, double dt)
		{
			int n = (int)sig.size();
			double mean = 0.0; for (double x : sig) mean += x; mean /= n;
			int first = -1, last = -1, count = 0;
			for (int i = 1; i < n; ++i)
			{
				double a = sig[i - 1] - mean, b = sig[i] - mean;
				if (a <= 0.0 && b > 0.0) { if (first < 0) first = i; last = i; ++count; }
			}
			if (count < 2 || first < 0) return 0.0;
			double span = (last - first) * dt;
			return (count - 1) / span;
		}
	} // namespace

	CylinderResult run_cylinder(const CylinderConfig& cfg)
	{
		CylinderResult res;
		MacGrid g; g.h = cfg.h;
		g.nx = cfg.nx_D * cfg.D; g.ny = cfg.ny_D * cfg.D; g.nz = cfg.nz;
		res.nx = g.nx; res.ny = g.ny; res.nz = g.nz;

		double Dphys = cfg.D * cfg.h, R = 0.5 * Dphys;
		double nu = cfg.U * Dphys / cfg.Re;
		int icx = (int)std::llround(cfg.xc_D * cfg.D), jcy = g.ny / 2;
		double xc = icx * cfg.h, yc = jcy * cfg.h;

		std::vector<unsigned char> solid;
		build_cylinder_mask(g, xc, yc, R, solid);

		ChannelBC bc;
		bc.inlet_mode = INLET_UNIFORM; bc.U_inlet = cfg.U; bc.Uc = cfg.U;
		bc.solid_mode = SOLID_NOSLIP; // resolved BL => Kármán street & correct Cd
		bc.ymin = bc.ymax = WALL_FREESLIP; bc.zmin = bc.zmax = WALL_FREESLIP;

		// fixed dt (uniform sampling for the FFT); wake peak ~1.7U.
		double dt_adv = cfg.cfl * cfg.h / (1.7 * cfg.U);
		double dt_diff = cfg.h * cfg.h / (6.0 * nu);
		double dt = 0.9 * std::min(dt_adv, dt_diff);
		res.dt = dt; res.Re = cfg.Re; res.D = Dphys; res.U = cfg.U;

		ChannelParams pr;
		pr.rho = cfg.rho; pr.nu = nu; pr.Cs = cfg.Cs; pr.cfl = cfg.cfl; pr.safety = 0.9;
		pr.proj_tol = 1e-4; pr.proj_max_iter = 80; pr.fixed_dt = dt; pr.advect_band = 1;

		ChannelFluidCore core(g, bc, pr, solid);
		core.init_uniform(cfg.U, cfg.v_blip * cfg.U);

		double Lx = g.nx * cfg.h;
		int spinup = (int)std::llround(cfg.spinup_flowthroughs * Lx / (cfg.U * dt));
		res.spinup_steps = spinup;
		int rec = cfg.record_samples; res.record_steps = rec;

		std::fprintf(stderr, "[cyl Re=%.0f] grid %dx%dx%d dt=%.4f spinup=%d rec=%d nu=%.4f\n",
			cfg.Re, g.nx, g.ny, g.nz, dt, spinup, rec, nu);
		std::fflush(stderr);

		for (int s = 0; s < spinup; ++s)
		{
			core.step();
			res.run_max_iters = std::max(res.run_max_iters, core.last_solve_iters());
			res.run_max_relres = std::max(res.run_max_relres, core.last_solve_relres());
			if ((s + 1) % 2000 == 0) { std::fprintf(stderr, "  spinup %d/%d iters=%d relres=%.1e\n", s + 1, spinup, core.last_solve_iters(), core.last_solve_relres()); std::fflush(stderr); }
		}

		// recording
		int n = g.p_count();
		double *Su, *Sv, *Sp, *Suu, *Suv;
		cudaMalloc(&Su, sizeof(double) * n); cudaMemset(Su, 0, sizeof(double) * n);
		cudaMalloc(&Sv, sizeof(double) * n); cudaMemset(Sv, 0, sizeof(double) * n);
		cudaMalloc(&Sp, sizeof(double) * n); cudaMemset(Sp, 0, sizeof(double) * n);
		cudaMalloc(&Suu, sizeof(double) * n); cudaMemset(Suu, 0, sizeof(double) * n);
		cudaMalloc(&Suv, sizeof(double) * n); cudaMemset(Suv, 0, sizeof(double) * n);

		int iprobe = icx + 5 * cfg.D, kmid = g.nz / 2;
		std::vector<double> vsig(rec);
		int block = 256, grid = (n + block - 1) / block;
		for (int s = 0; s < rec; ++s)
		{
			core.step();
			res.run_max_iters = std::max(res.run_max_iters, core.last_solve_iters());
			res.run_max_relres = std::max(res.run_max_relres, core.last_solve_relres());
			double qi = core.inlet_flux(), qo = core.outlet_flux();
			if (std::fabs(qi) > 1e-30) res.mass_imbalance = std::max(res.mass_imbalance, std::fabs(qo - qi) / std::fabs(qi));
			vsig[s] = core.get_v(iprobe, jcy, kmid);
			k_accum<<<grid, block>>>(core.u_dev(), core.v_dev(), core.p_dev(), Su, Sv, Sp, Suu, Suv, g, n);
			if ((s + 1) % 2000 == 0) { std::fprintf(stderr, "  record %d/%d v=%.4f iters=%d\n", s + 1, rec, vsig[s], core.last_solve_iters()); std::fflush(stderr); }
		}
		cudaDeviceSynchronize();

		// --- Strouhal ---
		double vmean = 0.0; for (double x : vsig) vmean += x; vmean /= rec;
		double vrms = 0.0; for (double x : vsig) vrms += (x - vmean) * (x - vmean); vrms = std::sqrt(vrms / rec);
		res.v_rms = vrms / cfg.U;
		double f = dominant_freq(vsig, dt);
		res.f_peak = f; res.St_fft = f * Dphys / cfg.U;
		res.St_zerocross = zerocross_freq(vsig, dt) * Dphys / cfg.U;
		res.sheds = (res.v_rms > 0.02);

		// --- Cd via control-volume momentum balance from the time means ---
		std::vector<double> U(n), V(n), P(n), UU(n), UV(n);
		cudaMemcpy(U.data(), Su, sizeof(double) * n, cudaMemcpyDeviceToHost);
		cudaMemcpy(V.data(), Sv, sizeof(double) * n, cudaMemcpyDeviceToHost);
		cudaMemcpy(P.data(), Sp, sizeof(double) * n, cudaMemcpyDeviceToHost);
		cudaMemcpy(UU.data(), Suu, sizeof(double) * n, cudaMemcpyDeviceToHost);
		cudaMemcpy(UV.data(), Suv, sizeof(double) * n, cudaMemcpyDeviceToHost);
		double inv = 1.0 / rec;
		for (int t = 0; t < n; ++t) { U[t] *= inv; V[t] *= inv; P[t] *= inv; UU[t] *= inv; UV[t] *= inv; }
		auto U_ = [&](int i, int j, int k) { return U[g.pidx(i, j, k)]; };
		auto V_ = [&](int i, int j, int k) { return V[g.pidx(i, j, k)]; };
		auto P_ = [&](int i, int j, int k) { return P[g.pidx(i, j, k)]; };
		auto UU_ = [&](int i, int j, int k) { return UU[g.pidx(i, j, k)]; };
		auto UV_ = [&](int i, int j, int k) { return UV[g.pidx(i, j, k)]; };

		int i1 = icx - 3 * cfg.D, i2 = icx + 6 * cfg.D, j1 = jcy - 4 * cfg.D, j2 = jcy + 4 * cfg.D;
		i1 = std::max(2, i1); i2 = std::min(g.nx - 3, i2);
		j1 = std::max(2, j1); j2 = std::min(g.ny - 3, j2);
		double hh = cfg.h * cfg.h, rho = cfg.rho, Fx = 0.0;
		for (int k = 0; k < g.nz; ++k)
			for (int j = j1; j <= j2; ++j)
			{
				auto Xface = [&](int i) { // x-face between cells i-1 and i
					double UUf = 0.5 * (UU_(i - 1, j, k) + UU_(i, j, k));
					double Pf = 0.5 * (P_(i - 1, j, k) + P_(i, j, k));
					double dUdx = (U_(i, j, k) - U_(i - 1, j, k)) / cfg.h;
					return rho * UUf + Pf - 2.0 * rho * nu * dUdx;
				};
				Fx += Xface(i1) * hh;      // left  (+)
				Fx -= Xface(i2 + 1) * hh;  // right (-)
			}
		for (int k = 0; k < g.nz; ++k)
			for (int i = i1; i <= i2; ++i)
			{
				auto Yface = [&](int j) { // y-face between cells j-1 and j
					double UVf = 0.5 * (UV_(i, j - 1, k) + UV_(i, j, k));
					double dUdy = (U_(i, j, k) - U_(i, j - 1, k)) / cfg.h;
					double Vfip = 0.5 * (V_(i + 1, j - 1, k) + V_(i + 1, j, k));
					double Vfim = 0.5 * (V_(i - 1, j - 1, k) + V_(i - 1, j, k));
					double dVdx = (Vfip - Vfim) / (2.0 * cfg.h);
					return rho * UVf - rho * nu * (dUdy + dVdx);
				};
				Fx -= Yface(j2 + 1) * hh;  // top    (-)
				Fx += Yface(j1) * hh;      // bottom (+)
			}
		double Lz = g.nz * cfg.h;
		res.Cd = 2.0 * Fx / (rho * cfg.U * cfg.U * Dphys * Lz);

		cudaFree(Su); cudaFree(Sv); cudaFree(Sp); cudaFree(Suu); cudaFree(Suv);
		return res;
	}
}
