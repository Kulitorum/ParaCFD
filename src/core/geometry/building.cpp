// building.cpp — centerline surface -> 2D footprint -> thickened wall + flat roof solid mask.
// See building.h for the rationale (voxelize the centerline directly, no watertight solid).
#include "core/geometry/building.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace windcfd::core
{
	namespace
	{
		constexpr double kPi = 3.14159265358979323846;

		struct V2
		{
			double x = 0.0, y = 0.0;
		};
		inline V2 operator-(V2 a, V2 b) { return { a.x - b.x, a.y - b.y }; }
		inline V2 operator+(V2 a, V2 b) { return { a.x + b.x, a.y + b.y }; }
		inline V2 operator*(V2 a, double s) { return { a.x * s, a.y * s }; }
		inline double dot(V2 a, V2 b) { return a.x * b.x + a.y * b.y; }
		inline double cross(V2 a, V2 b) { return a.x * b.y - a.y * b.x; }
		inline double len(V2 a) { return std::sqrt(a.x * a.x + a.y * a.y); }
		inline V2 norm(V2 a)
		{
			double l = len(a);
			return l > 1e-15 ? V2{ a.x / l, a.y / l } : V2{ 0.0, 0.0 };
		}

		// Squared distance from point p to segment [a,b].
		double dist2_pt_seg(V2 p, V2 a, V2 b)
		{
			V2 ab = b - a, ap = p - a;
			double L2 = dot(ab, ab);
			double t = L2 > 1e-30 ? dot(ap, ab) / L2 : 0.0;
			t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
			V2 c = a + ab * t;
			V2 d = p - c;
			return dot(d, d);
		}

		// Min distance from (x,y) to any loop edge (loops treated as closed).
		double dist_to_loops(double x, double y, const std::vector<Loop2D>& loops)
		{
			V2 p{ x, y };
			double best = 1e300;
			for (const Loop2D& lp : loops)
			{
				const std::size_t m = lp.size();
				if (m < 2) continue;
				for (std::size_t i = 0, j = m - 1; i < m; j = i++)
				{
					V2 a{ lp[j][0], lp[j][1] }, b{ lp[i][0], lp[i][1] };
					double d2 = dist2_pt_seg(p, a, b);
					if (d2 < best) best = d2;
				}
			}
			return std::sqrt(best);
		}

		// Replace each sharp corner of a closed loop with a tangent circular arc of radius r
		// (clamped to the adjacent edge lengths). Nearly-straight vertices pass through. This
		// rounds the CENTERLINE; the +/- half-thickness thickening then yields an outer corner
		// radius ~ r + half. `arc_segs` points are emitted per arc.
		Loop2D round_loop(const Loop2D& in, double r, int arc_segs)
		{
			const std::size_t m = in.size();
			if (m < 3 || r <= 1e-9) return in;
			Loop2D out;
			out.reserve(m * (std::size_t)(arc_segs + 2));
			for (std::size_t i = 0; i < m; ++i)
			{
				V2 A{ in[(i + m - 1) % m][0], in[(i + m - 1) % m][1] };
				V2 B{ in[i][0], in[i][1] };
				V2 C{ in[(i + 1) % m][0], in[(i + 1) % m][1] };
				V2 e1 = B - A, e2 = C - B;
				double L1 = len(e1), L2 = len(e2);
				if (L1 < 1e-9 || L2 < 1e-9) { out.push_back({ B.x, B.y }); continue; }
				V2 a1 = norm(A - B); // from B back toward A
				V2 a2 = norm(C - B); // from B toward C
				double cosphi = std::max(-1.0, std::min(1.0, dot(a1, a2)));
				double phi = std::acos(cosphi); // interior angle at B
				if (phi > kPi - 0.05) { out.push_back({ B.x, B.y }); continue; } // ~straight
				double t = r / std::tan(phi * 0.5);
				t = std::min(t, std::min(0.49 * L1, 0.49 * L2));
				double reff = t * std::tan(phi * 0.5);
				V2 bis = norm(a1 + a2); // into the interior of the corner
				double sinh = std::sin(phi * 0.5);
				if (sinh < 1e-6) { out.push_back({ B.x, B.y }); continue; }
				V2 center = B + bis * (reff / sinh);
				V2 T1 = B + a1 * t, T2 = B + a2 * t;
				double ang1 = std::atan2(T1.y - center.y, T1.x - center.x);
				double ang2 = std::atan2(T2.y - center.y, T2.x - center.x);
				double dA = ang2 - ang1;
				while (dA > kPi) dA -= 2.0 * kPi;
				while (dA < -kPi) dA += 2.0 * kPi;
				for (int s = 0; s <= arc_segs; ++s)
				{
					double a = ang1 + dA * ((double)s / arc_segs);
					out.push_back({ center.x + reff * std::cos(a), center.y + reff * std::sin(a) });
				}
			}
			return out;
		}

		// Signed area of a closed polygon (shoelace); |area| distinguishes a real building
		// outline from an open/degenerate wall path.
		double poly_area(const Loop2D& p)
		{
			const std::size_t m = p.size();
			if (m < 3) return 0.0;
			double a = 0.0;
			for (std::size_t i = 0, j = m - 1; i < m; j = i++)
				a += p[j][0] * p[i][1] - p[i][0] * p[j][1];
			return 0.5 * a;
		}

		// Point-in-triangle (edge-sign test, winding-agnostic).
		bool point_in_tri(V2 p, V2 a, V2 b, V2 c)
		{
			double d1 = cross(p - a, b - a);
			double d2 = cross(p - b, c - b);
			double d3 = cross(p - c, a - c);
			bool neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
			bool pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);
			return !(neg && pos);
		}

		// 2D convex hull (Andrew's monotone chain), CCW, of every point across `loops`.
		Loop2D convex_hull(const std::vector<Loop2D>& loops)
		{
			std::vector<std::array<double, 2>> pts;
			for (const Loop2D& lp : loops)
				for (const auto& p : lp) pts.push_back(p);
			const int n = (int)pts.size();
			if (n < 3) return {};
			std::sort(pts.begin(), pts.end(), [](const std::array<double, 2>& a, const std::array<double, 2>& b)
				{ return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]); });
			auto crs = [](const std::array<double, 2>& o, const std::array<double, 2>& a, const std::array<double, 2>& b)
				{ return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]); };
			Loop2D h(2 * (std::size_t)n);
			int k = 0;
			for (int i = 0; i < n; ++i) { while (k >= 2 && crs(h[k - 2], h[k - 1], pts[i]) <= 0.0) --k; h[k++] = pts[i]; }
			const int lower = k + 1;
			for (int i = n - 2; i >= 0; --i) { while (k >= lower && crs(h[k - 2], h[k - 1], pts[i]) <= 0.0) --k; h[k++] = pts[i]; }
			h.resize(k > 0 ? (std::size_t)(k - 1) : 0);
			return h;
		}

		// Inside test for a CCW convex polygon (point is left of every edge).
		bool point_in_convex(double x, double y, const Loop2D& hpts)
		{
			const std::size_t m = hpts.size();
			if (m < 3) return false;
			for (std::size_t i = 0, j = m - 1; i < m; j = i++)
			{
				const double ex = hpts[i][0] - hpts[j][0], ey = hpts[i][1] - hpts[j][1];
				const double px = x - hpts[j][0], py = y - hpts[j][1];
				if (ex * py - ey * px < 0.0) return false; // right of an edge -> outside
			}
			return true;
		}
	} // namespace

	std::size_t Footprint::point_count() const
	{
		std::size_t n = 0;
		for (const Loop2D& l : loops) n += l.size();
		return n;
	}

	Footprint mesh_horizontal_section(const TriMesh& mesh, double z0)
	{
		Footprint fp;
		if (mesh.empty()) return fp;

		const double zmin = mesh.bbox_min[2], zmax = mesh.bbox_max[2];
		const double zrange = std::max(1e-9, zmax - zmin);
		if (std::isnan(z0)) z0 = zmin + 0.5 * zrange;
		z0 += 1.3e-4 * zrange; // tiny nudge off exact vertices/seams (cf. the slicer's slicing offset)

		const float* P = mesh.positions.data();
		const std::uint32_t* I = mesh.indices.data();
		const std::size_t ntri = mesh.triangle_count();

		// Raw section segments as (ax,ay,bx,by).
		std::vector<std::array<double, 4>> segs;
		segs.reserve(ntri / 4 + 4);
		for (std::size_t t = 0; t < ntri; ++t)
		{
			const std::uint32_t vi[3] = { I[3 * t], I[3 * t + 1], I[3 * t + 2] };
			double px[3], py[3], d[3];
			for (int c = 0; c < 3; ++c)
			{
				px[c] = P[3 * vi[c]];
				py[c] = P[3 * vi[c] + 1];
				d[c] = (double)P[3 * vi[c] + 2] - z0;
			}
			double cx[2], cy[2];
			int nc = 0;
			for (int e = 0; e < 3 && nc < 2; ++e)
			{
				int u = e, v = (e + 1) % 3;
				if (d[u] * d[v] < 0.0)
				{
					double tt = d[u] / (d[u] - d[v]);
					cx[nc] = px[u] + tt * (px[v] - px[u]);
					cy[nc] = py[u] + tt * (py[v] - py[u]);
					++nc;
				}
			}
			if (nc == 2) segs.push_back({ cx[0], cy[0], cx[1], cy[1] });
		}

		if (segs.empty()) return fp;

		// xy bbox + chaining tolerance.
		double xmn = 1e300, ymn = 1e300, xmx = -1e300, ymx = -1e300;
		for (const auto& s : segs)
		{
			xmn = std::min({ xmn, s[0], s[2] });
			ymn = std::min({ ymn, s[1], s[3] });
			xmx = std::max({ xmx, s[0], s[2] });
			ymx = std::max({ ymx, s[1], s[3] });
		}
		const double diag = std::hypot(xmx - xmn, ymx - ymn);
		const double tol = std::max(1e-9, 1e-5 * diag);
		const double tol2 = tol * tol;

		// Greedy endpoint chaining into loops.
		const std::size_t n = segs.size();
		std::vector<char> used(n, 0);
		auto near = [&](double ax, double ay, double bx, double by)
		{ double dx = ax - bx, dy = ay - by; return dx * dx + dy * dy <= tol2; };

		for (std::size_t s = 0; s < n; ++s)
		{
			if (used[s]) continue;
			used[s] = 1;
			Loop2D loop;
			loop.push_back({ segs[s][0], segs[s][1] });
			double sx = segs[s][0], sy = segs[s][1];
			double cxp = segs[s][2], cyp = segs[s][3];
			loop.push_back({ cxp, cyp });
			for (std::size_t guard = 0; guard < n; ++guard)
			{
				if (near(cxp, cyp, sx, sy) && loop.size() >= 3) break; // closed
				std::size_t nextk = n;
				double nx = 0.0, ny = 0.0;
				for (std::size_t k = 0; k < n; ++k)
				{
					if (used[k]) continue;
					if (near(segs[k][0], segs[k][1], cxp, cyp)) { nextk = k; nx = segs[k][2]; ny = segs[k][3]; break; }
					if (near(segs[k][2], segs[k][3], cxp, cyp)) { nextk = k; nx = segs[k][0]; ny = segs[k][1]; break; }
				}
				if (nextk == n) break; // open chain: stop
				used[nextk] = 1;
				cxp = nx;
				cyp = ny;
				if (near(cxp, cyp, sx, sy)) break; // closed (don't duplicate start)
				loop.push_back({ cxp, cyp });
			}
			if (loop.size() >= 3) fp.loops.push_back(std::move(loop));
		}

		// Footprint bbox over the kept loops.
		double bxmn = 1e300, bymn = 1e300, bxmx = -1e300, bymx = -1e300;
		for (const Loop2D& lp : fp.loops)
			for (const auto& p : lp)
			{
				bxmn = std::min(bxmn, p[0]);
				bymn = std::min(bymn, p[1]);
				bxmx = std::max(bxmx, p[0]);
				bymx = std::max(bymx, p[1]);
			}
		if (!fp.loops.empty())
		{
			fp.bbox_min = { bxmn, bymn };
			fp.bbox_max = { bxmx, bymx };
		}

		std::printf("[building] section z=%.4f m: %zu segments -> %zu loops, %zu pts, footprint %.3f x %.3f m\n",
			z0, segs.size(), fp.loops.size(), fp.point_count(),
			fp.bbox_max[0] - fp.bbox_min[0], fp.bbox_max[1] - fp.bbox_min[1]);
		return fp;
	}

	Footprint center_footprint(const Footprint& fp, double Lx, double Ly)
	{
		Footprint out = fp;
		if (fp.empty()) return out;
		const double cx = 0.5 * (fp.bbox_min[0] + fp.bbox_max[0]);
		const double cy = 0.5 * (fp.bbox_min[1] + fp.bbox_max[1]);
		const double dx = 0.5 * Lx - cx, dy = 0.5 * Ly - cy;
		for (Loop2D& lp : out.loops)
			for (auto& p : lp) { p[0] += dx; p[1] += dy; }
		out.bbox_min = { fp.bbox_min[0] + dx, fp.bbox_min[1] + dy };
		out.bbox_max = { fp.bbox_max[0] + dx, fp.bbox_max[1] + dy };
		return out;
	}

	std::vector<unsigned char> voxelize_building(const Footprint& fp, const BuildingParams& prm,
		MacGrid g, int* out_solid_count)
	{
		std::vector<unsigned char> mask((std::size_t)g.p_count(), 0);
		if (fp.empty() || g.p_count() == 0) { if (out_solid_count) *out_solid_count = 0; return mask; }

		const double half = 0.5 * prm.wall_thickness;
		const double top = prm.base_z + prm.wall_height;
		const double roof_top = top + prm.roof_thickness;

		// Rounded corners: pre-fillet the centerline so the outer corner radius ~= corner_radius.
		std::vector<Loop2D> wall_loops = fp.loops;
		if (prm.corner_radius > half + 1e-9)
			for (Loop2D& lp : wall_loops) lp = round_loop(lp, prm.corner_radius - half, 6);

		// True SHARP (square) corners at corner_radius 0: start from the width-t distance band
		// (rounded corners), then square off each CONVEX corner by adding an outward miter WEDGE
		// triangle. Local + robust on non-convex outlines (no global offset polygon to self-
		// intersect). Curved walls (small turn per facet) are left rounded via the angle threshold.
		const bool sharp = (prm.corner_radius <= 1e-6);
		struct Wedge
		{
			V2 a, b, c;
		};
		std::vector<Wedge> wedges;
		if (sharp)
		{
			for (const Loop2D& lp : fp.loops)
			{
				const std::size_t m = lp.size();
				if (m < 3) continue;
				const double s = poly_area(lp) >= 0.0 ? 1.0 : -1.0; // loop orientation (+1 = CCW)
				for (std::size_t i = 0; i < m; ++i)
				{
					V2 A{ lp[(i + m - 1) % m][0], lp[(i + m - 1) % m][1] };
					V2 B{ lp[i][0], lp[i][1] };
					V2 C{ lp[(i + 1) % m][0], lp[(i + 1) % m][1] };
					V2 din = norm(B - A), dout = norm(C - B);
					if (len(din) < 0.5 || len(dout) < 0.5) continue;
					if (cross(din, dout) * s <= 0.15) continue; // only genuine convex corners (skip curves/concave)
					const V2 uin = s > 0 ? V2{ din.y, -din.x } : V2{ -din.y, din.x };   // outward edge normals
					const V2 uout = s > 0 ? V2{ dout.y, -dout.x } : V2{ -dout.y, dout.x };
					const double denom = 1.0 + dot(uin, uout);
					V2 mm = denom < 1e-6 ? uout : (uin + uout) * (1.0 / denom);
					const double ml = len(mm);
					if (ml > 6.0) mm = mm * (6.0 / ml); // miter limit -> bevel very sharp spikes
					wedges.push_back({ B + uin * half, B + mm * half, B + uout * half });
				}
			}
		}
		auto wall_hit = [&](double x, double y) -> bool
		{
			if (!sharp) return dist_to_loops(x, y, wall_loops) <= half;
			if (dist_to_loops(x, y, fp.loops) <= half) return true; // the rounded band...
			V2 p{ x, y };
			for (const Wedge& w : wedges) // ...squared off at convex corners
				if (point_in_tri(p, w.a, w.b, w.c)) return true;
			return false;
		};

		// Roof coverage: a SOLID slab over the CONVEX HULL of the footprint, dilated by the
		// overhang. Using the hull guarantees a filled roof regardless of how the section split
		// the footprint into loops (a plain point-in-loops test can leave the roof hollow).
		const Loop2D roof_hull = convex_hull(fp.loops);
		const std::vector<Loop2D> roof_hull_v{ roof_hull };
		const double roof_reach = half + prm.roof_overhang;

		// xy region that can possibly be solid (footprint + wall/overhang reach), as cell indices.
		const double reach = std::max(half, half + prm.roof_overhang) + 1.5 * g.h;
		auto ci = [&](double v, int n) { int c = (int)std::floor(v / g.h); return c < 0 ? 0 : (c > n - 1 ? n - 1 : c); };
		const int i0 = ci(fp.bbox_min[0] - reach, g.nx), i1 = ci(fp.bbox_max[0] + reach, g.nx);
		const int j0 = ci(fp.bbox_min[1] - reach, g.ny), j1 = ci(fp.bbox_max[1] + reach, g.ny);

		int wall_cells = 0, roof_cells = 0;
		for (int k = 0; k < g.nz; ++k)
		{
			const double cz = (k + 0.5) * g.h;
			const bool inWall = (cz >= prm.base_z && cz <= top);
			const bool inRoof = (prm.roof_thickness > 0.0 && cz > top && cz <= roof_top);
			if (!inWall && !inRoof) continue;
			for (int j = j0; j <= j1; ++j)
			{
				const double cy = (j + 0.5) * g.h;
				for (int i = i0; i <= i1; ++i)
				{
					const double cx = (i + 0.5) * g.h;
					bool solid = false;
					if (inWall)
					{
						if (wall_hit(cx, cy)) { solid = true; ++wall_cells; }
					}
					else // inRoof: solid slab over the convex-hull footprint + overhang skirt
					{
						if (point_in_convex(cx, cy, roof_hull) || dist_to_loops(cx, cy, roof_hull_v) <= roof_reach)
						{ solid = true; ++roof_cells; }
					}
					if (solid) mask[(std::size_t)g.pidx(i, j, k)] = 1;
				}
			}
		}

		const int total = wall_cells + roof_cells;
		if (out_solid_count) *out_solid_count = total;
		std::printf("[building] voxelize: t=%.3f h=%.2f r=%.3f overhang=%.2f roof=%.2f -> %d solid cells (%d wall + %d roof, %.3f%% of domain)\n",
			prm.wall_thickness, prm.wall_height, prm.corner_radius, prm.roof_overhang, prm.roof_thickness,
			total, wall_cells, roof_cells, 100.0 * total / std::max(1, g.p_count()));
		return mask;
	}
}
