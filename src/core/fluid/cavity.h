// cavity.h — lid-driven cavity driver + Ghia (1982) comparison for the M1 gate (V1).
// RESEARCH §10 V1, tables in research/11. Closed unit cube, lid moving in +x on the
// top face; spanwise (y) free-slip keeps the flow 2D so the mid-span plane reproduces
// the 2-D Ghia centreline profiles. Mapping: Ghia u(y)->our u(z) at x=0.5; Ghia
// v(x)->our w(x) at z=0.5 (both vertical velocity), normalised by the lid speed.
#pragma once

#include "core/fluid/stam_fluid_core.h"

namespace windcfd::core
{
	struct CavityMetrics
	{
		int N = 0;
		double Re = 0.0;
		int steps = 0;
		double final_change = 0.0; // max|du| over the last check interval
		double dt_last = 0.0;
		double rms_u = 0.0;        // RMS error of u(z) vs Ghia u(y), normalised by U
		double rms_v = 0.0;        // RMS error of w(x) vs Ghia v(x), normalised by U
		double umin_center = 0.0;  // min of u along the vertical centreline
		int run_max_iters = 0;     // max MGPCG iterations over the (warm-started) run
		double run_max_relres = 0.0;

		// Pressure/divergence sub-gate (filled only when probe=true).
		bool have_probe = false;
		ProjectionProbe probe_prod;   // production tol 1e-4
		ProjectionProbe probe_valid;  // validation tol 1e-6-equivalent
		ProjectionProbe probe_jacobi; // Jacobi fallback (debug)
	};

	// Run a lid-driven cavity at Reynolds number Re on an N^3 grid (unit cube, U=1).
	// Stops early when the max velocity change over `check_interval` steps falls below
	// `steady_tol`, or after `max_steps`. `probe` toggles the pressure sub-gate.
	CavityMetrics run_cavity(int N, double Re, int max_steps, int check_interval, double steady_tol, bool probe);
}
