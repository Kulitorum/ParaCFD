// sem_inlet.h — Jarrin (2006) Synthetic Eddy Method inlet-turbulence generator
// (RESEARCH §8, research/14 §3). N compact tent eddies convect through a box straddling
// the inlet plane; their superposition, coloured by the Cholesky factor of the
// Nezu–Nakagawa (1993) open-channel Reynolds-stress tensor, gives divergence-carrying
// inflow fluctuations that the MGPCG projection cleans each step (research/14 §3). Ambient
// TI from the un-amplified profile is 5–10% over the lower half depth; γ (≤2.5) scales the
// whole R_ij up for near-pile 10–20% TI targets.
//
// The physical inlines (Nezu profile, Cholesky a_ij(z), tent shape) are __host__ __device__
// so the GPU plane kernel and its CPU reference twin share arithmetic (parity test 1e-5).
#pragma once

#include "core/fluid/inlet_fluct.h"
#include "core/fluid/mac_grid.h"

#include <cmath>
#include <vector>

namespace windcfd::core
{
	struct SemParams
	{
		double sigma = 1.0;   // eddy length scale [m] = h_dom/5 (research/14 §3)
		int N = 150;          // eddy count (≈ V_B/σ³)
		double ustar = 0.0345;// inlet friction velocity [m/s]
		double h_dom = 5.0;   // domain height [m] (Nezu z/h)
		double U_d = 1.0;     // depth-avg speed = convection speed U_c
		double gamma = 1.0;   // Reynolds-stress amplitude scaling (≤2.5)
		double Ly = 10.0, Lz = 5.0; // inlet-plane extent [m]
		unsigned seed = 20260707u;
	};

	// Nezu–Nakagawa (1993) rms profiles (research/14 §2), z from bed, h = h_dom.
	WINDCFD_HD inline double nezu_urms(double us, double z, double h, double g) { return g * 2.30 * us * exp(-z / h); }
	WINDCFD_HD inline double nezu_vrms(double us, double z, double h, double g) { return g * 1.63 * us * exp(-z / h); }
	WINDCFD_HD inline double nezu_wrms(double us, double z, double h, double g) { return g * 1.27 * us * exp(-z / h); }
	// shear stress −u'w' = u*²(1−z/h) ⇒ R13 = <u'w'> = −u*²(1−z/h), scaled by γ².
	WINDCFD_HD inline double nezu_uw(double us, double z, double h, double g) { return -(g * g) * us * us * (1.0 - z / h); }

	// Lower-triangular Cholesky of R_ij (research/14 §3): a21=a32=0.
	WINDCFD_HD inline void sem_cholesky(double us, double z, double h, double g,
		double& a11, double& a22, double& a31, double& a33)
	{
		a11 = nezu_urms(us, z, h, g);
		a22 = nezu_vrms(us, z, h, g);
		double R33 = nezu_wrms(us, z, h, g); R33 = R33 * R33;
		a31 = a11 > 1e-12 ? nezu_uw(us, z, h, g) / a11 : 0.0;
		double d = R33 - a31 * a31;
		a33 = d > 0.0 ? sqrt(d) : 0.0;
	}

	// Tent shape f(s)=√(3/2)(1−|s|), |s|<1 (research/14 §3); ∫f² over [-1,1] = 1.
	WINDCFD_HD inline double sem_tent(double s) { double a = fabs(s); return a < 1.0 ? 1.2247448713915890 * (1.0 - a) : 0.0; }

	// Convected SEM sums S_j(y,z) = (1/√N) Σ_k ε_j^k f_σ(0−x_k, y−y_k, z−z_k).
	WINDCFD_HD inline void sem_sums(double y, double z, const double* ex, const double* ey, const double* ez,
		const double* e1, const double* e2, const double* e3, int N, double sigma, double pref,
		double& S1, double& S2, double& S3)
	{
		double s1 = 0, s2 = 0, s3 = 0, invs = 1.0 / sigma;
		for (int k = 0; k < N; ++k)
		{
			double fx = sem_tent((0.0 - ex[k]) * invs);
			if (fx == 0.0) continue;
			double fy = sem_tent((y - ey[k]) * invs);
			if (fy == 0.0) continue;
			double fz = sem_tent((z - ez[k]) * invs);
			if (fz == 0.0) continue;
			double f = pref * fx * fy * fz;
			s1 += e1[k] * f; s2 += e2[k] * f; s3 += e3[k] * f;
		}
		double inv = 1.0 / sqrt((double)N);
		S1 = s1 * inv; S2 = s2 * inv; S3 = s3 * inv;
	}

	// V_B and the f_σ prefactor √(V_B/σ³).
	WINDCFD_HD inline double sem_boxvol(double sigma, double Ly, double Lz)
	{ return (2.0 * sigma) * (Ly + 2.0 * sigma) * (Lz + 2.0 * sigma); }

	// Free-function plane evaluators (shared by class + parity test). Fill inlet-face
	// fluctuation planes: up[ny·nz], vp[(ny+1)·nz], wp[ny·(nz+1)].
	void sem_planes_gpu(const double* ex, const double* ey, const double* ez,
		const double* e1, const double* e2, const double* e3,
		double* up, double* vp, double* wp, MacGrid g, SemParams p);
	void sem_planes_cpu(const std::vector<double>& ex, const std::vector<double>& ey, const std::vector<double>& ez,
		const std::vector<double>& e1, const std::vector<double>& e2, const std::vector<double>& e3,
		std::vector<double>& up, std::vector<double>& vp, std::vector<double>& wp, MacGrid g, SemParams p);

	// SEM generator as an InletFluct. Deterministic given the seed.
	class SemInlet final : public InletFluct
	{
	public:
		SemInlet(MacGrid g, SemParams p);
		~SemInlet() override;
		SemInlet(const SemInlet&) = delete;
		SemInlet& operator=(const SemInlet&) = delete;

		void advance(double dt) override;
		void add_u(double* u, MacGrid g) override;
		void set_vw(double* v, double* w, MacGrid g) override;

		// Host copy of the current fluctuation planes (for the TI diagnostic).
		const std::vector<double>& u_plane() const { return hup_; }

	private:
		void regenerate(int k, bool upstream);
		MacGrid g_;
		SemParams p_;
		std::vector<double> hex_, hey_, hez_, he1_, he2_, he3_; // host eddy state
		std::vector<double> hup_;                                // host u' plane (diagnostic)
		double *ex_, *ey_, *ez_, *e1_, *e2_, *e3_;               // device eddy state
		double *up_, *vp_, *wp_;                                 // device fluct planes
		unsigned rng_;
	};
}
