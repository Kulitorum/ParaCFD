#include "core/geometry/triangle_bvh.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace paracfd::core
{
	Vec3d normalized(Vec3d a)
	{
		const double l2 = length2(a);
		return l2 > 1e-60 ? a / std::sqrt(l2) : Vec3d{};
	}

	void Aabb3d::expand(Vec3d p)
	{
		lo.x = std::min(lo.x, p.x); lo.y = std::min(lo.y, p.y); lo.z = std::min(lo.z, p.z);
		hi.x = std::max(hi.x, p.x); hi.y = std::max(hi.y, p.y); hi.z = std::max(hi.z, p.z);
	}
	void Aabb3d::expand(const Aabb3d& b) { if (b.valid()) { expand(b.lo); expand(b.hi); } }
	bool Aabb3d::valid() const { return lo.x <= hi.x && lo.y <= hi.y && lo.z <= hi.z; }
	bool Aabb3d::overlaps(const Aabb3d& b) const
	{
		return lo.x <= b.hi.x && hi.x >= b.lo.x && lo.y <= b.hi.y && hi.y >= b.lo.y && lo.z <= b.hi.z && hi.z >= b.lo.z;
	}
	double Aabb3d::distance2(Vec3d p) const
	{
		double s = 0.0;
		for (int i = 0; i < 3; ++i)
		{
			const double d = p[i] < lo[i] ? lo[i] - p[i] : (p[i] > hi[i] ? p[i] - hi[i] : 0.0);
			s += d * d;
		}
		return s;
	}

	namespace
	{
		Aabb3d triangle_bounds(const BvhTriangle& t)
		{
			Aabb3d b; b.expand(t.a); b.expand(t.b); b.expand(t.c); return b;
		}

		bool separated_on_axis(Vec3d axis, Vec3d v0, Vec3d v1, Vec3d v2, Vec3d half)
		{
			const double a2 = length2(axis);
			if (a2 < 1e-60) return false;
			const double p0 = dot(v0, axis), p1 = dot(v1, axis), p2 = dot(v2, axis);
			const double mn = std::min({p0, p1, p2}), mx = std::max({p0, p1, p2});
			const double r = half.x * std::abs(axis.x) + half.y * std::abs(axis.y) + half.z * std::abs(axis.z);
			const double eps = 32.0 * std::numeric_limits<double>::epsilon() * (r + std::max(std::abs(mn), std::abs(mx)) + 1.0);
			return mn > r + eps || mx < -r - eps;
		}

		bool segment_triangle(Vec3d origin, Vec3d delta, const BvhTriangle& tri, double t_min, double t_max,
			double& t, double& u, double& v)
		{
			const Vec3d e1 = tri.b - tri.a, e2 = tri.c - tri.a;
			const Vec3d p = cross(delta, e2);
			const double det = dot(e1, p);
			const double scale = std::sqrt(std::max(0.0, length2(e1) * length2(e2) * length2(delta)));
			if (std::abs(det) <= std::max(1e-30, 64.0 * std::numeric_limits<double>::epsilon() * scale)) return false;
			const double inv = 1.0 / det;
			const Vec3d s = origin - tri.a;
			u = dot(s, p) * inv;
			const double beps = 64.0 * std::numeric_limits<double>::epsilon();
			if (u < -beps || u > 1.0 + beps) return false;
			const Vec3d q = cross(s, e1);
			v = dot(delta, q) * inv;
			if (v < -beps || u + v > 1.0 + beps) return false;
			t = dot(e2, q) * inv;
			return t >= t_min && t <= t_max;
		}

		bool segment_aabb(Vec3d a, Vec3d d, const Aabb3d& b, double t0, double t1)
		{
			for (int axis = 0; axis < 3; ++axis)
			{
				if (std::abs(d[axis]) < 1e-300)
				{
					if (a[axis] < b.lo[axis] || a[axis] > b.hi[axis]) return false;
					continue;
				}
				double q0 = (b.lo[axis] - a[axis]) / d[axis];
				double q1 = (b.hi[axis] - a[axis]) / d[axis];
				if (q0 > q1) std::swap(q0, q1);
				t0 = std::max(t0, q0); t1 = std::min(t1, q1);
				if (t0 > t1) return false;
			}
			return true;
		}

		Vec3d closest_point_triangle(Vec3d p, const BvhTriangle& t)
		{
			const Vec3d ab = t.b - t.a, ac = t.c - t.a, ap = p - t.a;
			const double d1 = dot(ab, ap), d2 = dot(ac, ap);
			if (d1 <= 0.0 && d2 <= 0.0) return t.a;
			const Vec3d bp = p - t.b;
			const double d3 = dot(ab, bp), d4 = dot(ac, bp);
			if (d3 >= 0.0 && d4 <= d3) return t.b;
			const double vc = d1 * d4 - d3 * d2;
			if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) return t.a + ab * (d1 / (d1 - d3));
			const Vec3d cp = p - t.c;
			const double d5 = dot(ab, cp), d6 = dot(ac, cp);
			if (d6 >= 0.0 && d5 <= d6) return t.c;
			const double vb = d5 * d2 - d1 * d6;
			if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) return t.a + ac * (d2 / (d2 - d6));
			const double va = d3 * d6 - d5 * d4;
			if (va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0)
				return t.b + (t.c - t.b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
			const double inv = 1.0 / (va + vb + vc);
			return t.a + ab * (vb * inv) + ac * (vc * inv);
		}
	}

	bool TriangleBvh::triangle_intersects_aabb(const BvhTriangle& tri, const Aabb3d& box)
	{
		if (!triangle_bounds(tri).overlaps(box)) return false;
		const Vec3d c = (box.lo + box.hi) * 0.5, half = (box.hi - box.lo) * 0.5;
		const Vec3d v0 = tri.a - c, v1 = tri.b - c, v2 = tri.c - c;
		if (separated_on_axis({1, 0, 0}, v0, v1, v2, half) || separated_on_axis({0, 1, 0}, v0, v1, v2, half) || separated_on_axis({0, 0, 1}, v0, v1, v2, half)) return false;
		const Vec3d e[3] = {v1 - v0, v2 - v1, v0 - v2};
		if (separated_on_axis(cross(e[0], e[1]), v0, v1, v2, half)) return false;
		for (Vec3d q : e)
		{
			if (separated_on_axis(cross(q, {1, 0, 0}), v0, v1, v2, half) ||
				separated_on_axis(cross(q, {0, 1, 0}), v0, v1, v2, half) ||
				separated_on_axis(cross(q, {0, 0, 1}), v0, v1, v2, half)) return false;
		}
		return true;
	}

	void TriangleBvh::build(const TriMesh& mesh, std::uint32_t leaf_size)
	{
		primitives_.clear(); nodes_.clear(); triangles_.clear(); triangles_by_id_.clear();
		leaf_size = std::max<std::uint32_t>(1, leaf_size);
		primitives_.reserve(mesh.triangle_count());
		triangles_by_id_.resize(mesh.triangle_count());
		for (std::uint32_t id = 0; id < mesh.triangle_count(); ++id)
		{
			const std::uint32_t i0 = mesh.indices[3 * id], i1 = mesh.indices[3 * id + 1], i2 = mesh.indices[3 * id + 2];
			if (3ull * std::max({i0, i1, i2}) + 2 >= mesh.positions.size()) continue;
			auto p = [&](std::uint32_t i) { return Vec3d{mesh.positions[3 * i], mesh.positions[3 * i + 1], mesh.positions[3 * i + 2]}; };
			BvhTriangle t{p(i0), p(i1), p(i2), id, mesh.has_face_provenance() ? mesh.source_face_ids[id] : id};
			if (length2(cross(t.b - t.a, t.c - t.a)) < 1e-60) continue;
			Primitive prim{t, triangle_bounds(t), (t.a + t.b + t.c) / 3.0};
			primitives_.push_back(prim);
			triangles_by_id_[id] = t;
		}
		if (primitives_.empty()) return;
		build_node(0, static_cast<std::uint32_t>(primitives_.size()), leaf_size);
		triangles_.reserve(primitives_.size());
		for (const Primitive& p : primitives_) triangles_.push_back(p.tri);
	}

	std::uint32_t TriangleBvh::build_node(std::uint32_t first, std::uint32_t count, std::uint32_t leaf_size)
	{
		Node node; Aabb3d cb;
		for (std::uint32_t i = first; i < first + count; ++i) { node.bounds.expand(primitives_[i].bounds); cb.expand(primitives_[i].centroid); }
		const std::uint32_t id = static_cast<std::uint32_t>(nodes_.size()); nodes_.push_back(node);
		if (count <= leaf_size)
		{
			nodes_[id].first = first; nodes_[id].count = count; return id;
		}
		const Vec3d extent = cb.hi - cb.lo;
		int axis = extent.y > extent.x ? 1 : 0; if (extent.z > extent[axis]) axis = 2;
		const std::uint32_t mid = first + count / 2;
		std::nth_element(primitives_.begin() + first, primitives_.begin() + mid, primitives_.begin() + first + count,
			[axis](const Primitive& a, const Primitive& b) { return a.centroid[axis] < b.centroid[axis]; });
		const std::uint32_t left = build_node(first, mid - first, leaf_size);
		const std::uint32_t right = build_node(mid, first + count - mid, leaf_size);
		nodes_[id].left = left; nodes_[id].right = right;
		return id;
	}

	void TriangleBvh::query_aabb(const Aabb3d& box, std::vector<std::uint32_t>& out) const
	{
		out.clear(); if (nodes_.empty() || !box.valid()) return;
		std::vector<std::uint32_t> stack{0};
		while (!stack.empty())
		{
			const Node& n = nodes_[stack.back()]; stack.pop_back();
			if (!n.bounds.overlaps(box)) continue;
			if (n.leaf())
			{
				for (std::uint32_t i = n.first; i < n.first + n.count; ++i)
					if (triangle_intersects_aabb(primitives_[i].tri, box)) out.push_back(primitives_[i].tri.triangle_id);
			}
			else { stack.push_back(n.left); stack.push_back(n.right); }
		}
		std::sort(out.begin(), out.end());
	}
	std::vector<std::uint32_t> TriangleBvh::query_aabb(const Aabb3d& box) const { std::vector<std::uint32_t> r; query_aabb(box, r); return r; }

	SegmentHit TriangleBvh::intersect_segment(Vec3d a, Vec3d b, double t_min, double t_max) const
	{
		SegmentHit best; if (nodes_.empty()) return best;
		const Vec3d d = b - a; double best_t = t_max;
		std::vector<std::uint32_t> stack{0};
		while (!stack.empty())
		{
			const Node& n = nodes_[stack.back()]; stack.pop_back();
			if (!segment_aabb(a, d, n.bounds, t_min, best_t)) continue;
			if (n.leaf())
			{
				for (std::uint32_t i = n.first; i < n.first + n.count; ++i)
				{
					double t, u, v;
					if (segment_triangle(a, d, primitives_[i].tri, t_min, best_t, t, u, v))
					{
						best_t = t; best.hit = true; best.t = t; best.u = u; best.v = v;
						best.position = a + d * t; best.geometric_normal = normalized(cross(primitives_[i].tri.b - primitives_[i].tri.a, primitives_[i].tri.c - primitives_[i].tri.a));
						best.triangle_id = primitives_[i].tri.triangle_id; best.source_face_id = primitives_[i].tri.source_face_id;
					}
				}
			}
			else { stack.push_back(n.left); stack.push_back(n.right); }
		}
		return best;
	}

	NearestSurfacePoint TriangleBvh::nearest(Vec3d p, double max_distance) const
	{
		NearestSurfacePoint best; if (nodes_.empty()) return best;
		double best2 = max_distance * max_distance;
		std::vector<std::uint32_t> stack{0};
		while (!stack.empty())
		{
			const std::uint32_t ni = stack.back(); stack.pop_back(); const Node& n = nodes_[ni];
			if (n.bounds.distance2(p) > best2) continue;
			if (n.leaf())
			{
				for (std::uint32_t i = n.first; i < n.first + n.count; ++i)
				{
					const Vec3d q = closest_point_triangle(p, primitives_[i].tri); const double d2 = length2(q - p);
					if (d2 < best2)
					{
						best2 = d2; best.found = true; best.point = q; best.distance = std::sqrt(d2);
						best.geometric_normal = normalized(cross(primitives_[i].tri.b - primitives_[i].tri.a, primitives_[i].tri.c - primitives_[i].tri.a));
						best.triangle_id = primitives_[i].tri.triangle_id; best.source_face_id = primitives_[i].tri.source_face_id;
					}
				}
			}
			else
			{
				const double dl = nodes_[n.left].bounds.distance2(p), dr = nodes_[n.right].bounds.distance2(p);
				if (dl < dr) { stack.push_back(n.right); stack.push_back(n.left); } else { stack.push_back(n.left); stack.push_back(n.right); }
			}
		}
		return best;
	}

	double TriangleBvh::distance(Vec3d p, double max_distance) const { return nearest(p, max_distance).distance; }
}
