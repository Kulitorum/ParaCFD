// vti_synthetic.h — the single source of truth for the M0 synthetic test field.
//
// Both the C++ smoke test / CLI emitter and the Python validator
// (tools/validate_vti.py) compute this EXACT formula so the gate is a true
// producer/consumer round-trip. Keep the two implementations byte-for-byte in sync.
#pragma once

#include "core/vti_writer.h"

#include <cmath>

namespace paracfd::io
{
	// Deterministic smooth analytic field evaluated at point (i,j,k).
	// x=i*sx, y=j*sy, z=k*sz. Non-cubic dims are used on purpose to catch
	// index-ordering bugs. Value is finite for all inputs.
	inline double synthetic_field_value(int i, int j, int k, const VtiImage& g)
	{
		const double pi = 3.14159265358979323846;
		const double x = i * g.sx;
		const double y = j * g.sy;
		const double z = k * g.sz;
		const double Lx = g.nx * g.sx;
		const double Ly = g.ny * g.sy;
		const double Lz = g.nz * g.sz;
		const double fx = (Lx > 0.0) ? std::sin(2.0 * pi * x / Lx) : 0.0;
		const double fy = (Ly > 0.0) ? std::cos(2.0 * pi * y / Ly) : 0.0;
		const double fz = (Lz > 0.0) ? (z / Lz) : 0.0;
		return fx * fy + 0.5 * fz;
	}

	// Build the canonical M0 synthetic image: 17 x 13 x 9 points, h = 5 cm,
	// origin at the seabed corner, one scalar field named "synthetic".
	inline VtiImage make_synthetic_image()
	{
		VtiImage g;
		g.nx = 17;
		g.ny = 13;
		g.nz = 9;
		g.ox = 0.0;
		g.oy = 0.0;
		g.oz = 0.0;
		g.sx = 0.05;
		g.sy = 0.05;
		g.sz = 0.05;

		VtiField f;
		f.name = "synthetic";
		f.data.resize(g.num_points());
		for (int k = 0; k < g.nz; ++k)
			for (int j = 0; j < g.ny; ++j)
				for (int i = 0; i < g.nx; ++i)
				{
					const std::size_t idx = static_cast<std::size_t>(i) + static_cast<std::size_t>(g.nx) * (static_cast<std::size_t>(j) + static_cast<std::size_t>(g.ny) * static_cast<std::size_t>(k));
					f.data[idx] = static_cast<float>(synthetic_field_value(i, j, k, g));
				}
		g.point_fields.push_back(std::move(f));
		return g;
	}
}
