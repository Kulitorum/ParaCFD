// channel_empty.h — empty open-channel validation for the M2 gate (V3, second half):
// a 10×10×5 m channel with a rough log-law inlet, free-slip lid/floor/laterals and an
// Orlanski outlet must run steady and preserve the inlet velocity profile at mid-domain
// within 5% (a sheared u(z) with v=w=0 is an exact steady Euler solution; this checks
// the outlet BC, flux balance and projection don't corrupt it). RESEARCH §8.1-8.3.
#pragma once

namespace windcfd::core
{
	struct ChannelEmptyConfig
	{
		double Lx = 10.0, Ly = 10.0, Lz = 5.0; // m
		double h = 0.1;                        // m (voxel)
		double U = 1.0;                        // m/s depth-averaged target
		double d50 = 0.35e-3;                  // m (=> z0 = d50/12)
		double nu = 1.5e-5;                   // m^2/s
		double flowthroughs = 4.0;             // spin-up length
		double Cs = 0.0;
	};

	struct ChannelEmptyResult
	{
		int nx = 0, ny = 0, nz = 0;
		double h = 0, ustar = 0, z0 = 0;
		double profile_err_max = 0; // max_k |u_mid(k)-u_in(k)| / U
		double profile_err_rms = 0;
		double steady_change = 0;   // max |u| change over the final half flow-through / U
		double mass_imbalance = 0;
		double sim_time = 0;
		int steps = 0;
	};

	ChannelEmptyResult run_channel_empty(const ChannelEmptyConfig& cfg);
}
