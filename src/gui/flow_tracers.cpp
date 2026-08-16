// flow_tracers.cpp — see flow_tracers.h. Inlet-seeded streamline integration + ribbon packing.
#include "gui/flow_tracers.h"

#include "gui/colormap.h"

#include <algorithm>
#include <cmath>

namespace paracfd::gui
{
	using paracfd::core::MacGrid;

	namespace
	{
		inline int clampi(int x, int lo, int hi) { return x < lo ? lo : (x > hi ? hi : x); }
	}

	// --- velocity sampling (mirrors FlowParticles: cell-centred trilinear MAC reconstruction) ------

	bool FlowTracers::is_solid(const FlowField& f, float x, float y, float z) const
	{
		const MacGrid& g = f.grid;
		if (!f.solid) return false;
		int i = clampi((int)std::floor(paracfd::core::grid_fx(g, x)), 0, g.nx - 1);
		int j = clampi((int)std::floor(paracfd::core::grid_fy(g, y)), 0, g.ny - 1);
		int k = clampi((int)std::floor(paracfd::core::grid_fz(g, z)), 0, g.nz - 1);
		return f.solid[(size_t)g.pidx(i, j, k)] != 0;
	}

	bool FlowTracers::crosses_fabric(const FlowField& f, float ax, float ay, float az, float bx, float by, float bz) const
	{
		return f.fabric && f.fabric->intersect_segment({ax, ay, az}, {bx, by, bz}, 1e-8).hit;
	}

	void FlowTracers::sample(const FlowField& f, float x, float y, float z, double& uu, double& vv, double& ww) const
	{
		const MacGrid& g = f.grid;
		// Continuous cell-centre index space (graded-aware world→centre-index map; x/h-0.5 uniform).
		float gx = (float)paracfd::core::grid_cx(g, x), gy = (float)paracfd::core::grid_cy(g, y), gz = (float)paracfd::core::grid_cz(g, z);
		int i0 = clampi((int)std::floor(gx), 0, g.nx - 1);
		int j0 = clampi((int)std::floor(gy), 0, g.ny - 1);
		int k0 = clampi((int)std::floor(gz), 0, g.nz - 1);
		int i1 = std::min(i0 + 1, g.nx - 1);
		int j1 = std::min(j0 + 1, g.ny - 1);
		int k1 = std::min(k0 + 1, g.nz - 1);
		float fx = std::clamp(gx - (float)i0, 0.0f, 1.0f);
		float fy = std::clamp(gy - (float)j0, 0.0f, 1.0f);
		float fz = std::clamp(gz - (float)k0, 0.0f, 1.0f);

		auto cc = [&](int i, int j, int k, double& cu, double& cv, double& cw)
		{
			if (f.solid && f.solid[(size_t)g.pidx(i, j, k)]) { cu = cv = cw = 0.0; return; }
			cu = 0.5 * (f.u[g.uidx(i, j, k)] + f.u[g.uidx(i + 1, j, k)]);
			cv = 0.5 * (f.v[g.vidx(i, j, k)] + f.v[g.vidx(i, j + 1, k)]);
			cw = 0.5 * (f.w[g.widx(i, j, k)] + f.w[g.widx(i, j, k + 1)]);
		};

		double U[2][2][2], V[2][2][2], W[2][2][2];
		const int ii[2] = { i0, i1 }, jj[2] = { j0, j1 }, kk[2] = { k0, k1 };
		for (int a = 0; a < 2; ++a)
			for (int b = 0; b < 2; ++b)
				for (int c = 0; c < 2; ++c)
					cc(ii[a], jj[b], kk[c], U[a][b][c], V[a][b][c], W[a][b][c]);

		auto trilerp = [&](double q[2][2][2]) -> double
		{
			double c00 = q[0][0][0] * (1 - fx) + q[1][0][0] * fx;
			double c10 = q[0][1][0] * (1 - fx) + q[1][1][0] * fx;
			double c01 = q[0][0][1] * (1 - fx) + q[1][0][1] * fx;
			double c11 = q[0][1][1] * (1 - fx) + q[1][1][1] * fx;
			double c0 = c00 * (1 - fy) + c10 * fy;
			double c1 = c01 * (1 - fy) + c11 * fy;
			return c0 * (1 - fz) + c1 * fz;
		};
		uu = trilerp(U); vv = trilerp(V); ww = trilerp(W);
	}

	// --- inlet seed lattice ------------------------------------------------------------------------

	void FlowTracers::build_seeds(const TracerView& view, const FlowField& f)
	{
		const MacGrid& g = f.grid;
		const float Lx = (float)g.Lx(), Ly = (float)g.Ly(), Lz = (float)g.Lz();
		const float h = (float)g.h;
		// Seed just inside the INLET face, which flips with the flow direction (a reversed tide drives
		// the current in −x from the x-max face). integrate_line() then steps along the local flow, so
		// the streamlines run downstream from whichever face is the inlet.
		const float x0 = (f.flow_sign < 0) ? (Lx - 0.5f * h) : (0.5f * h);
		const int density = std::max(1, view.density);
		auto lin = [](int i, int n, float L) { return (((float)i + 0.5f) / (float)n) * L; };

		seeds_.clear();
		if (view.three_d)
		{
			// Inlet plane (x=x0): a y×z grid at uniform spacing (density along the larger dimension).
			const float spacing = std::max(Ly, Lz) / (float)density;
			const int ny = clampi((int)std::lround(Ly / spacing), 1, 120);
			const int nz = clampi((int)std::lround(Lz / spacing), 1, 120);
			for (int k = 0; k < nz; ++k)
				for (int j = 0; j < ny; ++j)
					seeds_.push_back({ x0, lin(j, ny, Ly), lin(k, nz, Lz) });
		}
		else if (view.axis == 2) // Z-normal slice: seed the inlet edge (x0, *, plane), vary y
		{
			const int ny = clampi(density, 1, 200);
			for (int j = 0; j < ny; ++j)
				seeds_.push_back({ x0, lin(j, ny, Ly), view.plane_pos });
		}
		else if (view.axis == 1) // Y-normal slice: seed (x0, plane, *), vary z
		{
			const int nz = clampi(density, 1, 200);
			for (int k = 0; k < nz; ++k)
				seeds_.push_back({ x0, view.plane_pos, lin(k, nz, Lz) });
		}
		else // X-normal slice (the inlet-facing plane): a (y,z) grid at x=plane_pos, cross-flow streamlines
		{
			const float spacing = std::max(Ly, Lz) / (float)density;
			const int ny = clampi((int)std::lround(Ly / spacing), 1, 120);
			const int nz = clampi((int)std::lround(Lz / spacing), 1, 120);
			for (int k = 0; k < nz; ++k)
				for (int j = 0; j < ny; ++j)
					seeds_.push_back({ view.plane_pos, lin(j, ny, Ly), lin(k, nz, Lz) });
		}
		hold_.assign(seeds_.size(), 0.0f); // fresh holds for the new lattice
	}

	// --- one streamline ----------------------------------------------------------------------------

	int FlowTracers::integrate_line(const TracerView& view, const FlowField& f, float sx, float sy, float sz)
	{
		const MacGrid& g = f.grid;
		const float Lx = (float)g.Lx(), Ly = (float)g.Ly(), Lz = (float)g.Lz();
		const float ds = view.step_ds;

		line_.clear();
		float x = sx, y = sy, z = sz;
		for (int s = 0; s < view.max_points; ++s)
		{
			if (x < 0 || x > Lx || y < 0 || y > Ly || z < 0 || z > Lz) break; // left an edge
			if (is_solid(f, x, y, z)) break;                                    // met a solid
			double uu, vv, ww; sample(f, x, y, z, uu, vv, ww);
			double vx = uu, vy = vv, vz = ww;
			if (!view.three_d) // constrain to the slice plane: drop the out-of-plane component + pin
			{
				if (view.axis == 0) { vx = 0; x = view.plane_pos; }
				else if (view.axis == 1) { vy = 0; y = view.plane_pos; }
				else { vz = 0; z = view.plane_pos; }
			}
			double sp = std::sqrt(vx * vx + vy * vy + vz * vz);
			line_.push_back({ x, y, z, (float)sp });
			if (sp < 1e-5) break; // stalled (eddy core / dead water) — stop here
			// Step by a fixed arc length ALONG the flow so vertex spacing is uniform regardless of
			// speed (a clean line; speed is carried in the colour, not the spacing).
			const float inv = (float)(ds / sp);
			const float nx = x + (float)vx * inv, ny = y + (float)vy * inv, nz = z + (float)vz * inv;
			if (crosses_fabric(f, x, y, z, nx, ny, nz)) break;
			x = nx; y = ny; z = nz;
		}
		return (int)line_.size();
	}

	// --- ribbon geometry from the current line_ ----------------------------------------------------

	void FlowTracers::emit_ribbon(const TracerView& view, int n)
	{
		const float eps = 1e-6f;
		// Ribbon: two vertices (side ±1) per point, tangent from central differences. Coloured by the
		// LOCAL speed through the shared ramp (no animation); alpha fades only at the very tips.
		const float inv_scale = view.speed_scale > eps ? 1.0f / view.speed_scale : 1.0f;
		const float taper = std::max(3.0f, std::min(8.0f, 0.1f * (float)n));
		firsts_.push_back(vertexCount());
		counts_.push_back(2 * n);
		for (int i = 0; i < n; ++i)
		{
			const Pt& p = line_[i];
			const Pt& pa = line_[std::max(0, i - 1)];
			const Pt& pb = line_[std::min(n - 1, i + 1)];
			float tx = pb.x - pa.x, ty = pb.y - pa.y, tz = pb.z - pa.z;
			float tl = std::sqrt(tx * tx + ty * ty + tz * tz);
			if (tl > eps) { tx /= tl; ty /= tl; tz /= tl; }
			else { tx = 1.0f; ty = 0.0f; tz = 0.0f; }

			const float t = std::clamp(p.spd * inv_scale, 0.0f, 1.0f);
			float r, gc, b; scour_colormap(t, r, gc, b);
			const float fade_in = std::clamp((float)i / 3.0f, 0.0f, 1.0f);
			const float fade_out = std::clamp((float)(n - 1 - i) / taper, 0.0f, 1.0f);
			const float a = 0.95f * fade_in * fade_out;

			for (int side = -1; side <= 1; side += 2)
			{
				vbo_.push_back(p.x); vbo_.push_back(p.y); vbo_.push_back(p.z); vbo_.push_back((float)side);
				vbo_.push_back(tx); vbo_.push_back(ty); vbo_.push_back(tz);
				vbo_.push_back(r); vbo_.push_back(gc); vbo_.push_back(b); vbo_.push_back(a);
			}
		}
	}

	// --- per-frame: integrate every seed, apply the boring filter + temporal hold ------------------

	void FlowTracers::advance(const TracerView& view, const FlowField& f)
	{
		vbo_.clear(); firsts_.clear(); counts_.clear();
		if (!f.u || f.grid.nx <= 0) return;
		const MacGrid& g = f.grid;

		// Rebuild the seed lattice (and reset the holds) whenever the config changes. plane_pos only
		// matters in 2D, so it does not trigger a rebuild while orbiting a 3D view.
		const float plane_sig = view.three_d ? 0.0f : view.plane_pos;
		const bool changed = !cfg_valid_ || view.three_d != s_three_d_ || view.axis != s_axis_
			|| view.density != s_density_ || g.nx != s_nx_ || g.ny != s_ny_ || g.nz != s_nz_
			|| f.flow_sign != s_sign_ || std::fabs(plane_sig - s_plane_) > 1e-6f;
		if (changed)
		{
			build_seeds(view, f);
			cfg_valid_ = true;
			s_three_d_ = view.three_d; s_axis_ = view.axis; s_density_ = view.density;
			s_nx_ = g.nx; s_ny_ = g.ny; s_nz_ = g.nz; s_sign_ = f.flow_sign; s_plane_ = plane_sig;
		}

		const float dt = std::clamp(view.dt, 0.0f, 0.5f);
		for (size_t k = 0; k < seeds_.size(); ++k)
		{
			const int n = integrate_line(view, f, seeds_[k][0], seeds_[k][1], seeds_[k][2]);

			// "Interesting" = the streamline exists and its path length clears the boring threshold
			// (or the filter is off). Refresh the hold when interesting; otherwise let it decay — a
			// tracer keeps drawing until its hold runs out, so it can't flicker at the threshold.
			bool interesting = false;
			if (n >= 2)
			{
				if (view.min_length <= 0.0f) interesting = true;
				else
				{
					float arc = 0.0f;
					for (int i = 1; i < n; ++i)
					{
						const float dx = line_[i].x - line_[i - 1].x, dy = line_[i].y - line_[i - 1].y, dz = line_[i].z - line_[i - 1].z;
						arc += std::sqrt(dx * dx + dy * dy + dz * dz);
					}
					interesting = arc >= view.min_length;
				}
			}
			if (interesting) hold_[k] = view.hold_seconds;
			else hold_[k] = view.instant ? 0.0f : std::max(0.0f, hold_[k] - dt); // instant ⇒ no retention

			if (n >= 2 && (interesting || hold_[k] > 0.0f)) emit_ribbon(view, n);
		}
	}
}
