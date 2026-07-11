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
}
