// channel_sem.cpp — see channel_sem.h. Host driver over ChannelFluidCore + SemInlet.
#include "core/fluid/channel_sem.h"
#include "core/fluid/channel_core.h"
#include "core/fluid/channel_ops.h" // ch_max_div_gpu
#include "core/fluid/sem_inlet.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace windcfd::core
{
	ChannelSemResult run_channel_sem(const ChannelSemConfig& cfg)
	{
		ChannelSemResult res;
		MacGrid g; g.h = cfg.h;
		g.nx = (int)std::llround(cfg.Lx / cfg.h);
		g.ny = (int)std::llround(cfg.Ly / cfg.h);
		g.nz = (int)std::llround(cfg.Lz / cfg.h);
		res.nx = g.nx; res.ny = g.ny; res.nz = g.nz; res.h = cfg.h;

		double z0 = cfg.d50 / 12.0;
		double ustar = 0.40 * cfg.U / (std::log(cfg.Lz / z0) - 1.0);
		res.z0 = z0; res.ustar = ustar;

		ChannelBC bc;
		bc.inlet_mode = INLET_LOGLAW; bc.U_inlet = cfg.U; bc.Uc = cfg.U;
		bc.ustar = ustar; bc.z0 = z0; bc.kappa = 0.40;
		bc.solid_mode = SOLID_FREESLIP;
		bc.ymin = bc.ymax = WALL_FREESLIP; bc.zmin = bc.zmax = WALL_FREESLIP;

		ChannelParams pr;
		pr.rho = 1.225; pr.nu = cfg.nu; pr.Cs = cfg.Cs; pr.cfl = 1.0; pr.safety = 0.9;
		pr.proj_tol = 1e-4; pr.proj_max_iter = 60; pr.fixed_dt = 0.0; pr.advect_band = 1;

		std::vector<unsigned char> nosolid((size_t)g.p_count(), 0);
		ChannelFluidCore core(g, bc, pr, nosolid);
		core.init_inlet_profile();

		SemParams sp; sp.sigma = cfg.sigma; sp.N = cfg.N; sp.ustar = ustar; sp.h_dom = cfg.Lz;
		sp.U_d = cfg.U; sp.gamma = cfg.gamma; sp.Ly = cfg.Ly; sp.Lz = cfg.Lz; sp.seed = cfg.seed;
		SemInlet sem(g, sp);
		core.set_inlet_fluct(&sem);

		double* divscr = nullptr; cudaMalloc(&divscr, sizeof(double) * g.p_count());

		double Lx = g.nx * cfg.h;
		double t = 0.0; int steps = 0;
		double T_spin = cfg.spinup_flowthroughs * Lx / cfg.U;
		while (t < T_spin) { t += core.step(); ++steps; }

		// recording window: accumulate mid-x streamwise stats + divergence audit.
		int imid = g.nx / 2;
		std::vector<double> su((size_t)g.ny * g.nz, 0.0), su2((size_t)g.ny * g.nz, 0.0);
		int nsamp = 0;
		double T_rec = cfg.record_flowthroughs * Lx / cfg.U, t0 = t, thalf = t + 0.5 * T_rec;
		double dfirst = 0, dsecond = 0; int nfirst = 0, nsecond = 0;
		while (t - t0 < T_rec)
		{
			double dt = core.step(); t += dt; ++steps;
			if (steps % cfg.sample_every != 0) continue;
			FluidSnapshot s = core.snapshot();
			for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j)
			{
				double uc = 0.5 * (s.u[g.uidx(imid, j, k)] + s.u[g.uidx(imid + 1, j, k)]);
				int p = k * g.ny + j; su[p] += uc; su2[p] += uc * uc;
			}
			double dmax = ch_max_div_gpu(core.u_dev(), core.v_dev(), core.w_dev(), divscr, core.solid_dev(), g);
			res.div_max = std::max(res.div_max, dmax);
			if (t < thalf) { dfirst += dmax; ++nfirst; } else { dsecond += dmax; ++nsecond; }
			++nsamp;
		}
		res.steps = steps; res.sim_time = t; res.samples = nsamp;
		res.div_first = nfirst ? dfirst / nfirst : 0.0;
		res.div_second = nsecond ? dsecond / nsecond : 0.0;

		// TI(z) spanwise-averaged streamwise intensity at mid-x.
		res.zc.resize(g.nz); res.TI_profile.assign(g.nz, 0.0);
		for (int k = 0; k < g.nz; ++k)
		{
			res.zc[k] = (k + 0.5) * cfg.h;
			double acc = 0.0;
			for (int j = 0; j < g.ny; ++j)
			{
				int p = k * g.ny + j;
				double mean = su[p] / std::max(1, nsamp);
				double var = su2[p] / std::max(1, nsamp) - mean * mean;
				if (var < 0.0) var = 0.0;
				acc += std::sqrt(var) / cfg.U;
			}
			res.TI_profile[k] = acc / g.ny;
		}
		// TI at z_ref and mean over the lower half [0.1,0.5]·Lz.
		double z_ref = cfg.z_ref_frac * cfg.Lz;
		int kref = std::min(g.nz - 1, std::max(0, (int)std::llround(z_ref / cfg.h - 0.5)));
		res.TI_ref = res.TI_profile[kref];
		double bacc = 0.0; int bn = 0;
		for (int k = 0; k < g.nz; ++k)
		{
			double z = res.zc[k];
			if (z >= 0.1 * cfg.Lz && z <= 0.5 * cfg.Lz) { bacc += res.TI_profile[k]; ++bn; }
		}
		res.TI_band_mean = bn ? bacc / bn : 0.0;

		cudaFree(divscr);
		return res;
	}
}
