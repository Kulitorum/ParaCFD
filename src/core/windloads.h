// windloads.h — integrate the wind loads on a voxelized building from the CFD pressure field.
//
// The building is a solid-cell mask on the MAC grid; its surface is the set of faces between a
// solid cell and an adjacent fluid cell. Surface pressure at such a face is taken as the pressure
// of the first fluid cell (cell-centred, Pa — see mac_grid.h). Pressure force on the body is
//   F = -sum over exposed faces of (p_face - p_ref) * n_hat * h^2      (pressure pushes inward),
// n_hat = outward face normal (solid -> fluid). Reference pressure p_ref is the mean pressure over
// the upstream inlet slab (freestream). Coefficients use the dynamic pressure q = 1/2 rho U_ref^2:
//   Cp = (p - p_ref)/q ; Cd = Fx/(q*A_frontal) ; Cs = Fy/(q*A_frontal) ; Cl = Fz/(q*A_plan).
// Only the PRESSURE load is integrated (skin friction is small for a bluff body and not included).
//
// Host-only, OCC-free (lives in libwindcfd). Units: SI. Wind is assumed along +x (the inlet).
#pragma once

#include "core/fluid/mac_grid.h"

#include <vector>

namespace windcfd::core
{
	struct WindLoadParams
	{
		double rho = 1.225;  // kg/m^3 (air)
		double u_ref = 1.0;  // m/s free-stream wind speed
	};

	struct WindLoads
	{
		double Fx = 0, Fy = 0, Fz = 0;    // integrated pressure force on the building [N]
		double Cd = 0, Cs = 0, Cl = 0;    // drag (+x), side (+y), lift (+z) coefficients
		double Mx = 0, My = 0, Mz = 0;    // moment about the building centroid [N*m]
		double CMx = 0, CMy = 0, CMz = 0; // moment coefficients (/ q A_ref L_ref)
		double A_frontal = 0, A_plan = 0; // projected reference areas [m^2]
		double L_ref = 0;                 // reference length = building height [m]
		double p_ref = 0;                 // reference (free-stream) pressure used [Pa]
		double cp_min = 0, cp_max = 0;    // Cp range over the surface
		long long exposed_faces = 0;
		long long solid_cells = 0;
	};

	// Compute the wind loads from the cell-centred pressure field `p` (g.pidx order, Pa) and the
	// 1=solid/0=fluid mask. If `out_cell_cp` != null (size g.p_count()), it receives a per-solid-
	// surface-cell mean Cp for visualisation (non-surface cells set to NaN).
	WindLoads compute_wind_loads(const double* p, const unsigned char* solid, MacGrid g,
		const WindLoadParams& prm, std::vector<float>* out_cell_cp = nullptr);

	// Converged, time-averaged wind-load statistics over a user-controlled window (see LoadAverager).
	// A turbulent (LES) flow never settles instantaneously — only its statistics do — so an unaveraged
	// snapshot is one random draw; these are the trustworthy quantities. `*_mean` = time-mean coefficient;
	// `*_rms` = std-dev of the fluctuation (the gust/unsteady component); [`*_min`,`*_max`] = the observed
	// peaks (design-critical, especially the roof/corner suction). `cp_min`/`cp_max` are the range of the
	// TIME-AVERAGED per-cell surface Cp. `Cd_drift` is a convergence hint (fraction): how far the running
	// mean still moves between the first half of the window and the whole window (< 0 ⇒ too few samples).
	struct WindLoadStats
	{
		long long samples = 0;
		double duration = 0.0;          // sim-time span of the averaging window [s]
		double flow_through_time = 0.0; // one flow-through time L_x/U_ref [s] (filled in by the caller)
		double Cd_mean = 0, Cl_mean = 0, Cs_mean = 0;
		double Cd_rms = 0, Cl_rms = 0, Cs_rms = 0; // std-dev of the fluctuation about the mean
		double Cd_min = 0, Cd_max = 0;
		double Cl_min = 0, Cl_max = 0;
		double Cs_min = 0, Cs_max = 0;
		double cp_min = 0, cp_max = 0; // range of the time-averaged per-cell Cp over the surface
		double Cd_drift = -1.0;        // |mean_all - mean_firsthalf| / |mean_all| (< 0 ⇒ N/A)
	};

	// Accumulates instantaneous WindLoads samples (fed at the solver's load cadence) into converged
	// time-averaged statistics: a running mean + RMS(fluctuation) + peak (min/max) of Cd/Cl/Cs (Welford,
	// numerically stable), plus a running time-average of the per-cell surface Cp for the building
	// colouring. The averaging window is user-controlled (reset() begins it, e.g. once the flow looks
	// spun-up); the surface mask is assumed FIXED over a window (a rebuild starts a new one). Host-only.
	class LoadAverager
	{
	public:
		// Begin a fresh window: `p_count` = grid.p_count() (per-cell Cp length), `t_start` = the current
		// sim-time [s]. Zeroes every accumulator.
		void reset(int p_count, double t_start);
		// Fold one instantaneous sample: `L` from compute_wind_loads, `cell_cp` its per-cell Cp field
		// (g.p_count(), NaN off the surface), `sim_time` the current flow time [s].
		void add(const WindLoads& L, const std::vector<float>& cell_cp, double sim_time);
		long long samples() const { return n_; }
		// Snapshot the statistics. If `out_mean_cp` != null it is (re)sized to p_count and filled with the
		// time-averaged per-cell Cp (NaN off the surface) for visualisation.
		WindLoadStats result(std::vector<float>* out_mean_cp = nullptr) const;

	private:
		long long n_ = 0;
		double t_start_ = 0.0, t_last_ = 0.0;
		double m_cd_ = 0, m_cl_ = 0, m_cs_ = 0; // running means (Welford)
		double s_cd_ = 0, s_cl_ = 0, s_cs_ = 0; // running M2 (sum of squared deviations)
		double min_cd_ = 0, max_cd_ = 0, min_cl_ = 0, max_cl_ = 0, min_cs_ = 0, max_cs_ = 0;
		std::vector<double> cp_sum_; // per-cell running sum of Cp (surface cells)
		std::vector<int> cp_cnt_;    // per-cell sample count (0 ⇒ never on the surface)
		std::vector<float> cd_hist_; // running-mean Cd after each sample (for the convergence-drift hint)
	};
}
