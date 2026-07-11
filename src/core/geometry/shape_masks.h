// shape_masks.h — procedural reference-shape voxel masks for the M8 shape-ranking harness
// (PLAN M8: "reference shapes generated procedurally … no CAD input needed"). OCC-free (in
// libscour): a shape is a cell solid-mask (1=solid, size p_count, same format as build_cylinder_mask
// / voxelize_mesh) so it drops straight into SeabedMorpho as the rigid `structure`. Each builder
// seats the shape ON the bed (z from bed_top upward), returns the solid-cell count, and centres it
// in x/y by argument. Shapes are sized to a common footprint so a ranking compares like with like.
//
// The set spans the trapping spectrum (RESEARCH §9): a flat PLATE (minimal enclosure, low trapping),
// a solid RING/annulus (encloses a plan area, reflective walls), and an open-top CUP/box (walls +
// floor, open top — the classic settling trap). Units: SI (m). A companion porous-k field (for
// perforated walls, M7) is built separately via screen_resistance_k + set_porous.
#pragma once

#include "core/fluid/mac_grid.h"

#include <cmath>
#include <vector>

namespace scour::core
{
	namespace shape_detail
	{
		inline int zrange(MacGrid g, double z, int& k) { k = (int)std::llround(z / g.h); return k; }
	}

	// Flat square plate: solid slab of half-extent `half` (m) and thickness `thick` (m) seated on the
	// bed. Minimal vertical obstruction — the low-trapping reference.
	inline int build_plate_mask(MacGrid g, double xc, double yc, double half, double bed_top, double thick, std::vector<unsigned char>& solid)
	{
		solid.assign((size_t)g.p_count(), 0);
		int k0 = (int)std::llround(bed_top / g.h), k1 = k0 + std::max(1, (int)std::llround(thick / g.h));
		int cnt = 0;
		for (int k = k0; k < k1 && k < g.nz; ++k)
			for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			{
				double x = (i + 0.5) * g.h, y = (j + 0.5) * g.h;
				if (std::fabs(x - xc) <= half && std::fabs(y - yc) <= half) { solid[g.pidx(i, j, k)] = 1; ++cnt; }
			}
		return cnt;
	}

	// Solid ring / annulus wall: solid where r_in ≤ r ≤ r_out, from the bed up by `height` (m).
	// Encloses a plan area with reflective walls (a mid-trapping, toe-scouring reference).
	inline int build_ring_mask(MacGrid g, double xc, double yc, double r_out, double r_in, double bed_top, double height, std::vector<unsigned char>& solid)
	{
		solid.assign((size_t)g.p_count(), 0);
		int k0 = (int)std::llround(bed_top / g.h), k1 = k0 + std::max(1, (int)std::llround(height / g.h));
		int cnt = 0;
		for (int k = k0; k < k1 && k < g.nz; ++k)
			for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			{
				double x = (i + 0.5) * g.h, y = (j + 0.5) * g.h, r = std::sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
				if (r >= r_in && r <= r_out) { solid[g.pidx(i, j, k)] = 1; ++cnt; }
			}
		return cnt;
	}

	// Open-bottom, open-top cup / box: four vertical walls of thickness `wall` (m) around a square of
	// half-extent `half`, from the bed up by `height`. NO floor — the interior bed is erodible so
	// sand entering over the rim is sheltered by the walls and ACCRETES on the bed inside (the classic
	// settling trap; matches Du's open-bottom units). A solid floor is deliberately omitted: this
	// single-interface f_pack bed cannot accrete sand ON a rigid surface (that needs the layered-column
	// overhang model, deferred), so an open bottom is required for the trap to register. Interior fluid.
	inline int build_cup_mask(MacGrid g, double xc, double yc, double half, double bed_top, double height, double wall, std::vector<unsigned char>& solid)
	{
		solid.assign((size_t)g.p_count(), 0);
		int k0 = (int)std::llround(bed_top / g.h), k1 = k0 + std::max(1, (int)std::llround(height / g.h));
		int cnt = 0;
		for (int k = k0; k < k1 && k < g.nz; ++k)
			for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i)
			{
				double x = (i + 0.5) * g.h, y = (j + 0.5) * g.h;
				double dx = std::fabs(x - xc), dy = std::fabs(y - yc);
				if (dx > half || dy > half) continue;                    // outside footprint
				bool wall_cell = (dx > half - wall || dy > half - wall); // perimeter walls only (open bottom)
				if (wall_cell) { solid[g.pidx(i, j, k)] = 1; ++cnt; }
			}
		return cnt;
	}
}
