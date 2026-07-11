// channel_empty.cpp — see channel_empty.h. Pure host driver over ChannelFluidCore.
#include "core/fluid/channel_empty.h"
#include "core/fluid/channel_core.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace scour::core
{
	ChannelEmptyResult run_channel_empty(const ChannelEmptyConfig& cfg)
	{
		ChannelEmptyResult res;
		MacGrid g; g.h = cfg.h;
		g.nx = (int)std::llround(cfg.Lx / cfg.h);
		g.ny = (int)std::llround(cfg.Ly / cfg.h);
		g.nz = (int)std::llround(cfg.Lz / cfg.h);
		res.nx = g.nx; res.ny = g.ny; res.nz = g.nz; res.h = cfg.h;

		double z0 = cfg.d50 / 12.0;
		double h_dom = g.nz * cfg.h;
		double ustar = 0.40 * cfg.U / (std::log(h_dom / z0) - 1.0); // u* = κ U_d/(ln(h/z0)-1)
		res.ustar = ustar; res.z0 = z0;

		ChannelBC bc;
		bc.inlet_mode = INLET_LOGLAW; bc.U_inlet = cfg.U; bc.Uc = cfg.U;
		bc.ustar = ustar; bc.z0 = z0; bc.kappa = 0.40;
		bc.solid_mode = SOLID_FREESLIP;
		bc.ymin = bc.ymax = WALL_FREESLIP; bc.zmin = bc.zmax = WALL_FREESLIP;

		ChannelParams pr;
		pr.rho = 1027.0; pr.nu = cfg.nu; pr.Cs = cfg.Cs; pr.cfl = 1.0; pr.safety = 0.9;
		pr.proj_tol = 1e-4; pr.proj_max_iter = 60; pr.fixed_dt = 0.0; pr.advect_band = 1;

		std::vector<unsigned char> nosolid((size_t)g.p_count(), 0);
		ChannelFluidCore core(g, bc, pr, nosolid);
		core.init_inlet_profile();

		auto mass_imb = [&]() {
			double qi = core.inlet_flux(), qo = core.outlet_flux();
			return std::fabs(qi) > 1e-30 ? std::fabs(qo - qi) / std::fabs(qi) : 0.0;
		};
		double Lx = g.nx * cfg.h;
		double T_target = cfg.flowthroughs * Lx / cfg.U;
		double t = 0.0; int steps = 0;
		while (t < T_target) { t += core.step(); ++steps; res.mass_imbalance = std::max(res.mass_imbalance, mass_imb()); }

		FluidSnapshot s1 = core.snapshot();
		double T_half = 0.5 * Lx / cfg.U, t2 = 0.0;
		while (t2 < T_half) { t2 += core.step(); ++steps; res.mass_imbalance = std::max(res.mass_imbalance, mass_imb()); }
		FluidSnapshot s2 = core.snapshot();
		res.sim_time = t + t2; res.steps = steps;

		// steadiness: max |u| change between s1 and s2 (all components), normalised by U
		double dmax = 0.0;
		for (size_t i = 0; i < s1.u.size(); ++i) dmax = std::max(dmax, std::fabs(s2.u[i] - s1.u[i]));
		for (size_t i = 0; i < s1.v.size(); ++i) dmax = std::max(dmax, std::fabs(s2.v[i] - s1.v[i]));
		for (size_t i = 0; i < s1.w.size(); ++i) dmax = std::max(dmax, std::fabs(s2.w[i] - s1.w[i]));
		res.steady_change = dmax / cfg.U;

		// profile preservation at mid-domain (x = Lx/2, mid-span j), per vertical node k
		int imid = g.nx / 2, jmid = g.ny / 2;
		double se = 0.0, mx = 0.0;
		for (int k = 0; k < g.nz; ++k)
		{
			double z = (k + 0.5) * cfg.h;
			double u_in = (ustar / 0.40) * std::log(z / z0);
			double u_mid = s2.u[g.uidx(imid, jmid, k)];
			double e = std::fabs(u_mid - u_in) / cfg.U;
			se += e * e; mx = std::max(mx, e);
		}
		res.profile_err_rms = std::sqrt(se / g.nz);
		res.profile_err_max = mx;
		return res;
	}
}
