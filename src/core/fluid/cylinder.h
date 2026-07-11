// cylinder.h — flow-past-a-cylinder benchmark driver for the M2 gate (V3). Sets up a
// thin-spanwise open channel with a procedural circular cylinder (no-slip surface, so
// the resolved boundary layer at D=32 / Re=100 generates the Kármán street — free-slip
// suppresses vorticity generation and cannot reproduce Cd≈1.33 or shedding), spins up
// 10 flow-throughs, then records the cross-stream (v) probe at (5D,0) and the control-
// volume drag. Strouhal from an FFT of the probe; Cd from the time-mean CV momentum
// balance. RESEARCH §10 V3, research/11 §4.
#pragma once

#include "core/fluid/channel_bc.h"

namespace scour::core
{
	struct CylinderConfig
	{
		double Re = 100.0;
		int D = 32;               // cylinder diameter in voxels
		double U = 1.0;           // freestream speed [m/s]
		double rho = 1.0;         // density (nondimensional)
		double h = 1.0;           // voxel size
		int nx_D = 24;            // streamwise length / D   (≥ 24)
		int ny_D = 12;            // cross-stream width / D   (≥ 12)
		int nz = 4;               // spanwise thickness (voxels; free-slip => 2D)
		double xc_D = 8.0;        // cylinder centre from inlet / D
		int spinup_flowthroughs = 10;
		int record_samples = 16384; // power of two (uniform-dt series for the FFT); ≥20 shedding cycles
		double cfl = 1.0;
		double Cs = 0.0;          // LES off for the laminar benchmark
		double v_blip = 0.05;     // transverse seed to trigger shedding [·U]
	};

	struct CylinderResult
	{
		double Re = 0, D = 0, U = 0, dt = 0;
		int nx = 0, ny = 0, nz = 0, spinup_steps = 0, record_steps = 0;
		double St_fft = 0, St_zerocross = 0, f_peak = 0;
		double Cd = 0;
		double v_rms = 0;          // recording-window RMS(v_probe)/U
		double mass_imbalance = 0; // max |Qout-Qin|/|Qin| over recording
		int run_max_iters = 0;
		double run_max_relres = 0;
		bool sheds = false;
	};

	CylinderResult run_cylinder(const CylinderConfig& cfg);
}
