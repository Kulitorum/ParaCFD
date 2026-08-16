// channel_sem.h — M3 SEM open-channel validation (gate V2, second half). A 10×10×5 m
// open channel (log-law inlet, Orlanski outlet, free-slip lid/laterals/bed, Smagorinsky
// LES) fed by the Jarrin SEM generator must (a) hold 5–10% ambient streamwise turbulence
// intensity at mid-domain and (b) show no divergence growth (the MGPCG projection cleans
// the SEM's non-solenoidal part each step). RESEARCH §8, research/14 §2–3.
#pragma once

#include <vector>

namespace paracfd::core
{
	struct ChannelSemConfig
	{
		double Lx = 10.0, Ly = 10.0, Lz = 5.0; // m
		double h = 0.1;                        // voxel [m]
		double U = 1.0;                        // depth-avg speed [m/s]
		double d50 = 0.2e-3;                   // grain size [m] → z0, inlet u*
		double nu = 1.5e-5;                   // m²/s
		double Cs = 0.11;                      // Smagorinsky (RESEARCH §3.4)
		double gamma = 1.0;                    // SEM Reynolds-stress amplitude (≤2.5)
		int N = 150;                           // SEM eddy count
		double sigma = 1.0;                    // SEM length scale [m] = h_dom/5
		double spinup_flowthroughs = 3.0;
		double record_flowthroughs = 3.0;
		int sample_every = 8;                  // TI/divergence sampling cadence [steps]
		double z_ref_frac = 0.1;               // TI report height z_ref = frac·Lz (Nezu ~7% at 0.1h)
		unsigned seed = 20260707u;
	};

	struct ChannelSemResult
	{
		int nx = 0, ny = 0, nz = 0, steps = 0, samples = 0;
		double h = 0, ustar = 0, z0 = 0, sim_time = 0;
		double TI_ref = 0;        // streamwise TI at z_ref, spanwise-averaged at mid-x
		double TI_band_mean = 0;  // mean over the lower half depth [0.1,0.5]·Lz
		double div_max = 0;       // max|∇·u| over the recording window (1/s)
		double div_first = 0, div_second = 0; // window-half means (growth check)
		std::vector<double> zc, TI_profile;   // TI(z) at mid-x
	};

	ChannelSemResult run_channel_sem(const ChannelSemConfig& cfg);
}
