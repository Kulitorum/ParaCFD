// channel_mask.h — procedural voxel obstacle masks for M2 (RESEARCH §3 obstacles,
// §7 packed-fraction bed arrives in M5). solid is a cell field (1=solid,0=fluid);
// nearsolid is its Chebyshev `band`-dilation used to trigger 1st-order advection
// reversion within CFL·h of a surface (RESEARCH §3.2). Host-side builders.
#pragma once

#include "core/fluid/mac_grid.h"

#include <vector>

namespace scour::core
{
	// Circular cylinder of radius R (m) centred at (xc,yc) in the x-y plane, extruded
	// through the full z extent (axis along z). Cell (i,j,k) solid iff its centre lies
	// inside the circle. Returns the enclosed solid-cell count.
	int build_cylinder_mask(MacGrid g, double xc, double yc, double R, std::vector<unsigned char>& solid);

	// Dilate `solid` by `band` cells (Chebyshev/∞-norm) into `nearsolid`.
	void dilate_mask(MacGrid g, const std::vector<unsigned char>& solid, int band, std::vector<unsigned char>& nearsolid);
}
