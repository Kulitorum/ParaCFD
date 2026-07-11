// cavity.cpp — see cavity.h. Pure host driver; all device work goes through
// StamFluidCore. Ghia, Ghia & Shin (1982) reference tables verbatim from research/11.
#include "core/fluid/cavity.h"

#include <cmath>
#include <cstdio>
#include <vector>

namespace windcfd::core
{
	namespace
	{
		// Ghia et al. (1982), 129x129 multigrid. u along the vertical centreline x=0.5.
		constexpr int NG = 17;
		const double GY[NG] = {1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000};
		const double GU100[NG] = {1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000};
		const double GU1000[NG] = {1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.06080, -0.10648, -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0.00000};
		// v along the horizontal centreline y=0.5.
		const double GX[NG] = {1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000};
		const double GV100[NG] = {0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000};
		const double GV1000[NG] = {0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, -0.42665, -0.31966, 0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.00000};

		// Linear interpolation on an ascending (pos,val) profile, clamped at the ends.
		double interp(const std::vector<double>& pos, const std::vector<double>& val, double t)
		{
			if (t <= pos.front()) return val.front();
			if (t >= pos.back()) return val.back();
			int lo = 0, hi = (int)pos.size() - 1;
			while (hi - lo > 1) { int m = (lo + hi) / 2; if (pos[m] <= t) lo = m; else hi = m; }
			double f = (t - pos[lo]) / (pos[hi] - pos[lo]);
			return val[lo] * (1 - f) + val[hi] * f;
		}
	}

	CavityMetrics run_cavity(int N, double Re, int max_steps, int check_interval, double steady_tol, bool probe)
	{
		CavityMetrics m;
		m.N = N; m.Re = Re;
		MacGrid g; g.nx = N; g.ny = N; g.nz = N; g.h = 1.0 / N;
		BC bc; // defaults = closed cavity: x no-slip, y free-slip, z-bottom no-slip, z-top lid
		bc.lid_u = 1.0;
		StamParams pr;
		pr.rho = 1.0;              // velocity is rho-independent for the cavity
		pr.nu = bc.lid_u * 1.0 / Re; // nu = U*L/Re, L=1 (RESEARCH §10 V1)
		pr.Cs = 0.0;               // LES off
		pr.cfl = 1.0; pr.safety = 0.9;
		pr.proj_tol = 1e-4; pr.proj_max_iter = 30; pr.use_mg = true;

		StamFluidCore core(g, bc, pr);

		int step = 0;
		double change = 1e30;
		core.max_velocity_change_and_remark(); // mark initial (returns sentinel)
		for (; step < max_steps; ++step)
		{
			core.step();
			m.run_max_iters = std::max(m.run_max_iters, core.last_solve_iters());
			m.run_max_relres = std::max(m.run_max_relres, core.last_solve_relres());
			if ((step + 1) % check_interval == 0)
			{
				change = core.max_velocity_change_and_remark();
				if (change < steady_tol) { ++step; break; }
			}
		}
		m.steps = step;
		m.final_change = change;
		m.dt_last = core.last_dt();

		// --- extract centreline profiles ---------------------------------------
		FluidSnapshot s = core.snapshot();
		double h = g.h, L = 1.0;
		int ic = N / 2, jc = N / 2, kc = N / 2;

		// u(z) at x=0.5, mid-span. Augment with wall/lid BC endpoints.
		std::vector<double> zs, uz;
		zs.push_back(0.0); uz.push_back(0.0); // no-slip bottom
		for (int k = 0; k < N; ++k) { zs.push_back((k + 0.5) * h); uz.push_back(s.u[g.uidx(ic, jc, k)]); }
		zs.push_back(L); uz.push_back(bc.lid_u); // lid

		// w(x) at z=0.5, mid-span (our vertical velocity == Ghia v).
		std::vector<double> xs, wx;
		xs.push_back(0.0); wx.push_back(0.0); // no-slip left wall (w tangential = 0)
		for (int i = 0; i < N; ++i) { xs.push_back((i + 0.5) * h); wx.push_back(s.w[g.widx(i, jc, kc)]); }
		xs.push_back(L); wx.push_back(0.0); // no-slip right wall

		const double* gu = (Re < 500.0) ? GU100 : GU1000;
		const double* gv = (Re < 500.0) ? GV100 : GV1000;
		double su = 0.0, sv = 0.0, umin = 0.0;
		for (int p = 0; p < NG; ++p)
		{
			double du = interp(zs, uz, GY[p]) - gu[p];
			double dv = interp(xs, wx, GX[p]) - gv[p];
			su += du * du; sv += dv * dv;
		}
		for (double val : uz) umin = std::min(umin, val);
		m.rms_u = std::sqrt(su / NG);
		m.rms_v = std::sqrt(sv / NG);
		m.umin_center = umin;

		// --- pressure / divergence sub-gate ------------------------------------
		if (probe)
		{
			m.have_probe = true;
			m.probe_prod = core.projection_probe(1e-4, 30, true);
			m.probe_valid = core.projection_probe(1e-9, 60, true); // tight enough that max|div|<1e-6 U/h
			m.probe_jacobi = core.projection_probe(1e-4, 2000, false);
		}
		return m;
	}
}
