#include "core/geometry/triangle_bvh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>

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
		void validate_topology_vertex_coordinates(const TriMesh& mesh)
		{
			if(mesh.has_malformed_topology_vertex_ids())
				throw std::invalid_argument("TriangleBvh topology vertex sidecar has "+
					std::to_string(mesh.topology_vertex_ids.size())+" IDs for "+
					std::to_string(mesh.vertex_count())+" vertices");
			if(!mesh.has_topology_vertex_ids())return;
			struct FirstUse
			{
				std::size_t vertex=0;
				std::array<double,3> position{};
			};
			std::map<std::uint32_t,FirstUse> first_use;
			for(std::size_t vertex=0;vertex<mesh.vertex_count();++vertex)
			{
				const std::uint32_t topology=mesh.topology_vertex_id(vertex);
				const std::array<double,3> position=mesh.vertex_position_double(vertex);
				const auto [found,inserted]=first_use.emplace(topology,FirstUse{vertex,position});
				if(inserted)continue;
				if(found->second.position==position)continue;
				throw std::invalid_argument("TriangleBvh topology vertex ID "+
					std::to_string(topology)+" has inconsistent coordinates at vertices "+
					std::to_string(found->second.vertex)+" and "+std::to_string(vertex));
			}
		}

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

		double segment_segment_distance2(Vec3d p0, Vec3d p1, Vec3d q0, Vec3d q1)
		{
			// Closest points of two closed segments, including parallel and degenerate cases.
			// Ericson's clamped formulation is deliberately kept in double precision because
			// this is static CAD preprocessing.
			const Vec3d d1=p1-p0,d2=q1-q0,r=p0-q0;
			const double a=dot(d1,d1),e=dot(d2,d2),f=dot(d2,r);
			double s=0.0,t=0.0;
			const double eps=64.0*std::numeric_limits<double>::epsilon();
			if(a<=eps&&e<=eps)return length2(p0-q0);
			if(a<=eps)t=std::clamp(f/e,0.0,1.0);
			else
			{
				const double c=dot(d1,r);
				if(e<=eps)s=std::clamp(-c/a,0.0,1.0);
				else
				{
					const double b=dot(d1,d2),denom=a*e-b*b;
					if(std::abs(denom)>eps*a*e)s=std::clamp((b*f-c*e)/denom,0.0,1.0);
					t=(b*s+f)/e;
					if(t<0.0){t=0.0;s=std::clamp(-c/a,0.0,1.0);}
					else if(t>1.0){t=1.0;s=std::clamp((b-c)/a,0.0,1.0);}
				}
			}
			return length2((p0+d1*s)-(q0+d2*t));
		}

		double segment_triangle_distance2(Vec3d a,Vec3d b,const BvhTriangle& tri)
		{
			double t=0,u=0,v=0;
			if(segment_triangle(a,b-a,tri,0.0,1.0,t,u,v))return 0.0;
			double best=std::min(length2(a-closest_point_triangle(a,tri)),length2(b-closest_point_triangle(b,tri)));
			best=std::min(best,segment_segment_distance2(a,b,tri.a,tri.b));
			best=std::min(best,segment_segment_distance2(a,b,tri.b,tri.c));
			return std::min(best,segment_segment_distance2(a,b,tri.c,tri.a));
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
		edge_certificate_by_triangle_.clear();edge_contact_tolerance_=0.0;edge_clearance_tolerance_=0.0;
		validate_topology_vertex_coordinates(mesh);
		if(mesh.has_malformed_cad_edge_atom_provenance())
			throw std::invalid_argument("TriangleBvh CAD edge atom sidecar has "+
				std::to_string(mesh.triangle_cad_edge_atom_ids.size())+" IDs for "+
				std::to_string(mesh.indices.size())+" triangle half-edges");
		leaf_size = std::max<std::uint32_t>(1, leaf_size);
		primitives_.reserve(mesh.triangle_count());
		triangles_by_id_.resize(mesh.triangle_count());
		edge_certificate_by_triangle_.resize(mesh.triangle_count());
		std::array<double,3> coordinate_lo{{std::numeric_limits<double>::infinity(),
			std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()}};
		std::array<double,3> coordinate_hi{{-std::numeric_limits<double>::infinity(),
			-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()}};
		auto scan_coordinate=[&](std::size_t coordinate,double value)
		{
			if(!std::isfinite(value))return;const int axis=static_cast<int>(coordinate%3);
			coordinate_lo[axis]=std::min(coordinate_lo[axis],value);
			coordinate_hi[axis]=std::max(coordinate_hi[axis],value);
		};
		if(mesh.has_fp64_positions())
		{
			for(std::size_t coordinate=0;coordinate<mesh.positions_fp64.size();++coordinate)
				scan_coordinate(coordinate,mesh.positions_fp64[coordinate]);
		}
		else
		{
			for(std::size_t coordinate=0;coordinate<mesh.positions.size();++coordinate)
				scan_coordinate(coordinate,static_cast<double>(mesh.positions[coordinate]));
		}
		double geometry_scale=0.0;for(int axis=0;axis<3;++axis)if(coordinate_lo[axis]<=coordinate_hi[axis])
			geometry_scale=std::max(geometry_scale,coordinate_hi[axis]-coordinate_lo[axis]);
		geometry_scale=std::max(geometry_scale,static_cast<double>(std::numeric_limits<float>::min()));
		// Float-only meshes cannot resolve geometry closer than a few coordinate ULPs,
		// so they retain a wide explicit ambiguity band.  Imported CAD carries an FP64
		// coordinate sidecar and per-edge BRep tolerances; its roundoff baseline can be
		// correspondingly narrow without pretending that every edge has one global CAD
		// tolerance.
		// Coincidence is a double-precision geometric predicate even when the source
		// coordinates were quantized to float.  FP32 uncertainty belongs in the wider
		// ambiguity/clearance band; using it as an attachment radius would silently close
		// a real, represented opening.  Extent rather than absolute world position keeps
		// the certificate invariant under a harmless translation.
		edge_contact_tolerance_=256.0*std::numeric_limits<double>::epsilon()*geometry_scale;
		edge_clearance_tolerance_=(mesh.has_fp64_positions()
			?4096.0*std::numeric_limits<double>::epsilon()
			:32.0*std::numeric_limits<float>::epsilon())*geometry_scale;
		for(auto& triangle:edge_certificate_by_triangle_)for(auto& edge:triangle)
		{
			edge.contact_tolerance=edge_contact_tolerance_;
			edge.clearance_tolerance=edge_clearance_tolerance_;
		}
		std::vector<std::uint8_t> valid_triangle(mesh.triangle_count(),0);
		for (std::uint32_t id = 0; id < mesh.triangle_count(); ++id)
		{
			const std::uint32_t i0 = mesh.indices[3 * id], i1 = mesh.indices[3 * id + 1], i2 = mesh.indices[3 * id + 2];
			if (3ull * std::max({i0, i1, i2}) + 2 >= mesh.positions.size()) continue;
			auto p = [&](std::uint32_t i)
			{
				const auto point=mesh.vertex_position_double(i);return Vec3d{point[0],point[1],point[2]};
			};
			BvhTriangle t{p(i0), p(i1), p(i2), id, mesh.has_face_provenance() ? mesh.source_face_ids[id] : id};
			if (length2(cross(t.b - t.a, t.c - t.a)) < 1e-60) continue;
			Primitive prim{t, triangle_bounds(t), (t.a + t.b + t.c) / 3.0};
			primitives_.push_back(prim);
			triangles_by_id_[id] = t;
			valid_triangle[id]=1;
			for(auto& edge:edge_certificate_by_triangle_[id])
				edge.unknown_reason=FabricEdgeUnknownReason::no_geometric_evidence;
		}
		if (primitives_.empty()) return;
		build_node(0, static_cast<std::uint32_t>(primitives_.size()), leaf_size);
		triangles_.reserve(primitives_.size());
		for (const Primitive& p : primitives_) triangles_.push_back(p.tri);

		struct IndexedHalfEdge
		{
			std::uint64_t vertex_key=0;
			std::uint64_t atom_id=TriMesh::kNoCadEdgeAtomId;
			std::uint32_t triangle=0;
			std::uint8_t local_edge=0;
			std::size_t flat_index=0;
		};
		std::vector<IndexedHalfEdge> half_edges;half_edges.reserve(3*primitives_.size());
		const std::size_t half_edge_count=3*mesh.triangle_count();
		std::vector<std::size_t> attachment_parent(half_edge_count);
		std::iota(attachment_parent.begin(),attachment_parent.end(),std::size_t{0});
		std::vector<std::vector<std::size_t>> claimed_peers(half_edge_count);
		std::vector<std::uint8_t> hard_unknown(half_edge_count,0);
		auto attachment_root=[&](std::size_t edge)
		{
			std::size_t root=edge;while(attachment_parent[root]!=root)root=attachment_parent[root];
			while(attachment_parent[edge]!=edge){const std::size_t next=attachment_parent[edge];attachment_parent[edge]=root;edge=next;}return root;
		};
		auto join_attachment=[&](std::size_t a,std::size_t b)
		{
			a=attachment_root(a);b=attachment_root(b);if(a!=b)attachment_parent[std::max(a,b)]=std::min(a,b);
		};
		auto certificate_for=[&](std::size_t flat)->FabricEdgeCertificate&
		{
			return edge_certificate_by_triangle_[flat/3][flat%3];
		};
		auto edge_points=[&](std::size_t flat)
		{
			const BvhTriangle& triangle=triangles_by_id_[flat/3];const Vec3d points[3]={triangle.a,triangle.b,triangle.c};
			return std::array<Vec3d,2>{points[flat%3],points[(flat%3+1)%3]};
		};
		std::map<std::uint32_t,double> cad_edge_max_tolerance;
		for(std::uint32_t triangle_id=0;triangle_id<mesh.triangle_count();++triangle_id)
		{
			if(!valid_triangle[triangle_id])continue;
			for(int edge=0;edge<3;++edge)
			{
				const std::uint32_t ia=mesh.indices[3*triangle_id+edge],ib=mesh.indices[3*triangle_id+(edge+1)%3];
				// Render vertices remain face-local so a rib and skin may retain independent
				// UVs/normals.  A conforming CAD import supplies the shared discrete identity
				// separately; legacy/programmatic meshes fall back to their render indices.
				const std::uint32_t topology_a=mesh.topology_vertex_id(ia);
				const std::uint32_t topology_b=mesh.topology_vertex_id(ib);
				const std::uint32_t lo=std::min(topology_a,topology_b),hi=std::max(topology_a,topology_b);
				const std::size_t flat=3*static_cast<std::size_t>(triangle_id)+edge;
				half_edges.push_back({(static_cast<std::uint64_t>(lo)<<32)|hi,
					mesh.cad_edge_atom_id(triangle_id,edge),triangle_id,
					static_cast<std::uint8_t>(edge),flat});
				hard_unknown[flat]=mesh.cad_edge_provenance_unknown(triangle_id,edge)?1u:0u;
				if(hard_unknown[flat])certificate_for(flat).unknown_reason=
					FabricEdgeUnknownReason::source_provenance;
				certificate_for(flat).cad_contact_id=mesh.cad_edge_contact_id(triangle_id,edge);
				const std::uint32_t cad_edge=mesh.cad_edge_id(triangle_id,edge);
				if(cad_edge!=TriMesh::kNoCadEdgeId)
					cad_edge_max_tolerance[cad_edge]=std::max(cad_edge_max_tolerance[cad_edge],
						std::max(0.0,mesh.cad_edge_tolerance(triangle_id,edge)));
			}
		}
		std::sort(half_edges.begin(),half_edges.end(),[](const IndexedHalfEdge& a,const IndexedHalfEdge& b)
		{
			if(a.vertex_key!=b.vertex_key)return a.vertex_key<b.vertex_key;
			if(a.atom_id!=b.atom_id)return a.atom_id<b.atom_id;
			if(a.triangle!=b.triangle)return a.triangle<b.triangle;
			return a.local_edge<b.local_edge;
		});

		for(std::size_t first=0;first<half_edges.size();)
		{
			std::size_t last=first+1;while(last<half_edges.size()
				&&half_edges[last].vertex_key==half_edges[first].vertex_key
				&&half_edges[last].atom_id==half_edges[first].atom_id)++last;
			if(last-first>1)
			{
				const std::size_t incidence=last-first;
				bool non_manifold=incidence!=2,provenance_unknown=false,declared_inconsistent=false;
				std::uint32_t declared_incidence=0;double group_contact=edge_contact_tolerance_;
				std::uint64_t declared_contact=TriMesh::kNoCadContactId;
				std::size_t contact_members=0;
				bool cad_boundary=false,cross_face=false;const std::uint32_t first_face=
					triangles_by_id_[half_edges[first].triangle].source_face_id;
				for(std::size_t q=first;q<last;++q)
				{
					const auto& edge=half_edges[q];const std::uint32_t cad_incidence=
						mesh.cad_edge_incident_face_count(edge.triangle,edge.local_edge);
					const std::uint64_t contact=mesh.cad_edge_contact_id(
						edge.triangle,edge.local_edge);
					const std::uint32_t contact_fan=contact==TriMesh::kNoCadContactId?0u:
						mesh.cad_edge_certified_fan_degree(edge.triangle,edge.local_edge);
					const std::uint32_t edge_declared_incidence=std::max({cad_incidence,
						mesh.cad_edge_is_periodic_seam(edge.triangle,edge.local_edge)?2u:0u,
						contact_fan});
					if(contact!=TriMesh::kNoCadContactId)
					{
						declared_inconsistent=declared_inconsistent||(declared_contact
							!=TriMesh::kNoCadContactId&&declared_contact!=contact)||contact_fan<2;
						declared_contact=contact;++contact_members;
					}
					provenance_unknown=provenance_unknown||hard_unknown[edge.flat_index]!=0;
					if(edge_declared_incidence>=2)
					{
						declared_inconsistent=declared_inconsistent||
							(declared_incidence!=0&&declared_incidence!=edge_declared_incidence);
						declared_incidence=std::max(declared_incidence,edge_declared_incidence);
					}
					group_contact=std::max(group_contact,
						std::max(0.0,mesh.cad_edge_tolerance(edge.triangle,edge.local_edge)));
					cad_boundary=cad_boundary||mesh.cad_edge_id(edge.triangle,edge.local_edge)!=TriMesh::kNoCadEdgeId;
					cross_face=cross_face||(mesh.has_face_provenance()&&
						triangles_by_id_[edge.triangle].source_face_id!=first_face);
				}
				const bool declaration_mismatch=declared_inconsistent||
					(declared_incidence>=2&&declared_incidence!=incidence)||
					(contact_members!=0&&contact_members!=incidence);
				const std::uint16_t fan=static_cast<std::uint16_t>(std::min<std::size_t>(
					std::numeric_limits<std::uint16_t>::max(),incidence));
				for(std::size_t q=first;q<last;++q)
				{
					const std::size_t flat=half_edges[q].flat_index;auto& output=certificate_for(flat);
					output.contact_tolerance=group_contact;
				if(provenance_unknown||declaration_mismatch)
					{
						hard_unknown[flat]=1;output.kind=FabricEdgeKind::unknown;
						output.role=FabricEdgeRole::unknown;output.incident_fan_degree=0;
						output.unknown_reason=provenance_unknown
							?FabricEdgeUnknownReason::source_provenance
							:FabricEdgeUnknownReason::indexed_fan_mismatch;
						continue;
					}
					output.kind=FabricEdgeKind::attached;
					output.unknown_reason=FabricEdgeUnknownReason::none;
					output.role=non_manifold?FabricEdgeRole::junction:
						((cad_boundary||cross_face)?FabricEdgeRole::manifold_seam:FabricEdgeRole::tessellation_interior);
					output.incident_fan_degree=non_manifold?std::max<std::uint16_t>(3,fan):2;
					for(std::size_t peer=first;peer<last;++peer)if(peer!=q)
						claimed_peers[flat].push_back(half_edges[peer].flat_index);
				}
				first=last;continue;
			}
			first=last;
		}

		struct BoundaryWitness
		{
			double lo=0.0,hi=0.0;
			std::size_t peer=0;
			double tolerance=0.0;
		};
		struct InteriorWitness{double lo=0.0,hi=0.0;std::uint32_t triangle=0;};
		auto edge_overlap_interval=[&](Vec3d a,Vec3d b,Vec3d c,Vec3d d,double tolerance,
			double& overlap_lo,double& overlap_hi)
		{
			const Vec3d ab=b-a,cd=d-c;const double ab2=length2(ab),cd2=length2(cd);
			if(!(ab2>0&&cd2>0))return false;const double ab_length=std::sqrt(ab2),cd_length=std::sqrt(cd2);
			const double angular=std::sqrt(length2(cross(ab/ab_length,cd/cd_length)));
			if(angular>std::max(512.0*std::numeric_limits<double>::epsilon(),
				4.0*tolerance/std::min(ab_length,cd_length)))return false;
			auto line_distance2=[](Vec3d point,Vec3d origin,Vec3d direction,double direction2)
			{
				return length2(cross(point-origin,direction))/direction2;
			};
			const double tolerance2=tolerance*tolerance;
			if(line_distance2(c,a,ab,ab2)>tolerance2||line_distance2(d,a,ab,ab2)>tolerance2||
				line_distance2(a,c,cd,cd2)>tolerance2||line_distance2(b,c,cd,cd2)>tolerance2)return false;
			const double tc=dot(c-a,ab)/ab2,td=dot(d-a,ab)/ab2;
			overlap_lo=std::max(0.0,std::min(tc,td));overlap_hi=std::min(1.0,std::max(tc,td));
			const double parameter_tolerance=std::max(512.0*std::numeric_limits<double>::epsilon(),
				tolerance/ab_length);
			if(overlap_hi-overlap_lo<=parameter_tolerance)return false;
			const Vec3d midpoint=a+ab*(0.5*(overlap_lo+overlap_hi));
			return segment_segment_distance2(midpoint,midpoint,c,d)<=tolerance2;
		};
		auto triangle_interior_interval=[&](Vec3d a,Vec3d b,const BvhTriangle& triangle,
			double tolerance,double& interval_lo,double& interval_hi)
		{
			const Vec3d ab=triangle.b-triangle.a,ac=triangle.c-triangle.a;
			const Vec3d raw_normal=cross(ab,ac);const double normal_length=std::sqrt(length2(raw_normal));
			if(!(normal_length>0))return false;const Vec3d normal=raw_normal/normal_length;
			const double distance_a=dot(a-triangle.a,normal),distance_b=dot(b-triangle.a,normal);
			if(std::max(std::abs(distance_a),std::abs(distance_b))>tolerance)return false;
			const Vec3d projected_a=a-normal*distance_a,projected_b=b-normal*distance_b;
			const double d00=dot(ab,ab),d01=dot(ab,ac),d11=dot(ac,ac),denominator=d00*d11-d01*d01;
			if(!(std::abs(denominator)>std::numeric_limits<double>::min()))return false;
			auto barycentric=[&](Vec3d point)
			{
				const Vec3d relative=point-triangle.a;const double d20=dot(relative,ab),d21=dot(relative,ac);
				const double u=(d11*d20-d01*d21)/denominator,v=(d00*d21-d01*d20)/denominator;
				return std::array<double,3>{1.0-u-v,u,v};
			};
			const auto weight_a=barycentric(projected_a),weight_b=barycentric(projected_b);
			const double triangle_scale=std::sqrt(std::max({length2(ab),length2(ac),length2(triangle.c-triangle.b)}));
			const double margin=std::max(512.0*std::numeric_limits<double>::epsilon(),tolerance/triangle_scale);
			interval_lo=0.0;interval_hi=1.0;
			for(int weight=0;weight<3;++weight)
			{
				const double first_weight=weight_a[weight],delta=weight_b[weight]-first_weight;
				if(std::abs(delta)<=std::numeric_limits<double>::epsilon())
				{
					if(first_weight<=margin)return false;continue;
				}
				// A transverse crossing of a triangle boundary is still interior to the
				// receiving sheet when the adjacent triangle continues it.  Clipping each
				// interval to a positive tolerance margin leaves a false uncovered strip at
				// every tessellation diagonal.  Use only a roundoff-sized negative allowance
				// for varying barycentric weights; the constant-weight branch above still
				// rejects a source edge that actually follows the candidate boundary.
				const double boundary_weight=-512.0*std::numeric_limits<double>::epsilon();
				const double crossing=(boundary_weight-first_weight)/delta;
				if(delta>0)interval_lo=std::max(interval_lo,crossing);
				else interval_hi=std::min(interval_hi,crossing);
			}
			return interval_hi-interval_lo>512.0*std::numeric_limits<double>::epsilon();
		};

		for(const IndexedHalfEdge& edge:half_edges)
		{
			auto& output=certificate_for(edge.flat_index);
			if(output.kind==FabricEdgeKind::attached||hard_unknown[edge.flat_index])continue;
			const std::uint32_t exact_cad_fan=mesh.cad_edge_certified_fan_degree(
				edge.triangle,edge.local_edge);
			if(exact_cad_fan>=2&&output.cad_contact_id!=TriMesh::kNoCadContactId)
			{
				output.kind=FabricEdgeKind::attached;
				output.role=exact_cad_fan>2?FabricEdgeRole::junction:FabricEdgeRole::manifold_seam;
				output.unknown_reason=FabricEdgeUnknownReason::none;
				output.incident_fan_degree=static_cast<std::uint16_t>(std::min<std::uint32_t>(
					exact_cad_fan,std::numeric_limits<std::uint16_t>::max()));
				output.contact_tolerance=std::max(edge_contact_tolerance_,
					std::max(0.0,mesh.cad_edge_tolerance(edge.triangle,edge.local_edge)));
				continue;
			}
			// Explicit topology IDs are the authoritative connectivity contract. Distinct
			// IDs may deliberately occupy the same coordinates (two coincident but
			// disconnected sheets), so the legacy geometric recovery below must never
			// weld them. A declared shared CAD edge that failed to acquire an indexed
			// peer is an error; an otherwise unclaimed edge is a confirmed opening.
			if(mesh.has_topology_vertex_ids())
			{
				const std::uint32_t cad_incidence=std::max<std::uint32_t>(
					mesh.cad_edge_incident_face_count(edge.triangle,edge.local_edge),
					mesh.cad_edge_is_periodic_seam(edge.triangle,edge.local_edge)?2u:0u);
				const bool invalid_contact=output.cad_contact_id!=TriMesh::kNoCadContactId;
				if(cad_incidence>=2||invalid_contact)
				{
					output.kind=FabricEdgeKind::unknown;output.role=FabricEdgeRole::unknown;
					output.incident_fan_degree=0;
					output.unknown_reason=FabricEdgeUnknownReason::indexed_fan_mismatch;
				}
				else
				{
					output.kind=FabricEdgeKind::confirmed_free;output.role=FabricEdgeRole::none;
					output.incident_fan_degree=1;
					output.unknown_reason=FabricEdgeUnknownReason::none;
				}
				continue;
			}
			const auto source_points=edge_points(edge.flat_index);const Vec3d a=source_points[0],b=source_points[1],delta=b-a;
			const std::uint32_t source_vertex_a=mesh.indices[3*edge.triangle+edge.local_edge];
			const std::uint32_t source_vertex_b=mesh.indices[3*edge.triangle+(edge.local_edge+1)%3];
			const std::uint32_t source_face=triangles_by_id_[edge.triangle].source_face_id;
			auto is_endpoint_neighbour=[&](std::uint32_t candidate,double tolerance)
			{
				const double tolerance2=tolerance*tolerance;
				for(unsigned corner=0;corner<3;++corner)
				{
					const std::uint32_t vertex=mesh.indices[3*candidate+corner];
					if(mesh.has_face_provenance()&&triangles_by_id_[candidate].source_face_id==source_face&&
						(vertex==source_vertex_a||vertex==source_vertex_b))return true;
					const auto position=mesh.vertex_position_double(vertex);
					const Vec3d point{position[0],position[1],position[2]};
					if(length2(point-a)<=tolerance2||length2(point-b)<=tolerance2)return true;
				}
				return false;
			};
			const double edge_length=std::sqrt(length2(delta));
			if(!(edge_length>edge_clearance_tolerance_))continue;
			const std::uint32_t source_cad_id=mesh.cad_edge_id(edge.triangle,edge.local_edge);
			const std::uint32_t reported_cad_incidence=mesh.cad_edge_incident_face_count(edge.triangle,edge.local_edge);
			const bool periodic_cad_seam=mesh.cad_edge_is_periodic_seam(edge.triangle,edge.local_edge);
			const std::uint32_t source_cad_incidence=std::max<std::uint32_t>(
				reported_cad_incidence,periodic_cad_seam?2u:0u);
			const double source_cad_tolerance=std::max(0.0,mesh.cad_edge_tolerance(edge.triangle,edge.local_edge));
			// A BRep edge tolerance certifies the representation of that edge; it does not
			// certify contact with a different edge or with the interior of another face.
			// Exact STEP preprocessing supplies explicit contact atoms for those cases.
			// The legacy/preview path may recover differently tessellated segments of the
			// *same* TopoDS_Edge, but must never turn proximity between distinct CAD
			// entities into an attachment and thereby close a real opening.
			const double source_contact_tolerance=std::max(edge_contact_tolerance_,source_cad_tolerance);
			double query_tolerance=edge_clearance_tolerance_;if(source_cad_id!=TriMesh::kNoCadEdgeId)
			{
				const auto found=cad_edge_max_tolerance.find(source_cad_id);
				if(found!=cad_edge_max_tolerance.end())query_tolerance=std::max(query_tolerance,found->second);
			}
			Aabb3d query;query.expand(a);query.expand(b);const Vec3d margin{query_tolerance,query_tolerance,query_tolerance};query.lo=query.lo-margin;query.hi=query.hi+margin;
			std::vector<std::uint32_t> candidate_triangles;query_aabb(query,candidate_triangles);
			std::sort(candidate_triangles.begin(),candidate_triangles.end());
			candidate_triangles.erase(std::unique(candidate_triangles.begin(),candidate_triangles.end()),candidate_triangles.end());
			std::vector<BoundaryWitness> boundary_witnesses;std::vector<InteriorWitness> interior_witnesses;
			std::vector<double> cuts{0.0,1.0};
			for(std::uint32_t candidate: candidate_triangles)
			{
				if(candidate==edge.triangle||candidate>=valid_triangle.size()||!valid_triangle[candidate])continue;
				const BvhTriangle& candidate_triangle=triangles_by_id_[candidate];
				double interior_lo=0,interior_hi=0;
				if(source_cad_id==TriMesh::kNoCadEdgeId&&
					triangle_interior_interval(a,b,candidate_triangle,source_contact_tolerance,
						interior_lo,interior_hi))
				{
					interior_witnesses.push_back({interior_lo,interior_hi,candidate});cuts.push_back(interior_lo);cuts.push_back(interior_hi);
				}
				for(int candidate_edge=0;candidate_edge<3;++candidate_edge)
				{
					const std::size_t peer=3*static_cast<std::size_t>(candidate)+candidate_edge;
					const auto peer_points=edge_points(peer);const std::uint32_t peer_cad_id=mesh.cad_edge_id(candidate,candidate_edge);
					if(source_cad_id!=TriMesh::kNoCadEdgeId&&peer_cad_id!=source_cad_id)
						continue;
					double pair_tolerance=source_contact_tolerance;
					if(peer_cad_id!=TriMesh::kNoCadEdgeId)
						pair_tolerance=std::max(pair_tolerance,
							std::max(0.0,mesh.cad_edge_tolerance(candidate,candidate_edge)));
					double overlap_lo=0,overlap_hi=0;
					if(!edge_overlap_interval(a,b,peer_points[0],peer_points[1],pair_tolerance,overlap_lo,overlap_hi))continue;
					boundary_witnesses.push_back({overlap_lo,overlap_hi,peer,pair_tolerance});cuts.push_back(overlap_lo);cuts.push_back(overlap_hi);
				}
			}
			std::sort(cuts.begin(),cuts.end());
			cuts.erase(std::unique(cuts.begin(),cuts.end(),[](double x,double y)
			{
				return std::abs(x-y)<=512.0*std::numeric_limits<double>::epsilon();
			}),cuts.end());
			bool saw_attached=false,saw_clear=false,saw_junction=false,saw_ambiguous=false;
			bool cad_coverage_ok=true;std::uint16_t maximum_fan=0;
			std::uint16_t minimum_observed_fan=std::numeric_limits<std::uint16_t>::max();
			std::uint16_t maximum_observed_fan=0;double used_contact=source_contact_tolerance;
			std::vector<std::size_t> used_peers;
			for(std::size_t cut=0;cut+1<cuts.size();++cut)
			{
				const double lo=cuts[cut],hi=cuts[cut+1],span=hi-lo;if(!(span>0))continue;
				// Predicate clipping can leave roundoff-sized intervals at a shared endpoint
				// or where a 1:N seam changes peer segment.  They have no representable
				// positive length and must not poison an otherwise complete reciprocal
				// attachment.  Real gaps in the wider representation-uncertainty band reach
				// the distance test below and remain explicitly unknown.
				if(span*edge_length<=64.0*edge_contact_tolerance_)continue;
				const double midpoint=0.5*(lo+hi);std::vector<std::size_t> peers;bool interior=false;
				for(const BoundaryWitness& witness:boundary_witnesses)if(midpoint>witness.lo&&midpoint<witness.hi)
				{
					peers.push_back(witness.peer);used_contact=std::max(used_contact,witness.tolerance);
				}
				std::sort(peers.begin(),peers.end());peers.erase(std::unique(peers.begin(),peers.end()),peers.end());
				for(const InteriorWitness& witness:interior_witnesses)
					interior=interior||(midpoint>witness.lo&&midpoint<witness.hi);
				const std::uint16_t observed_fan=static_cast<std::uint16_t>(std::min<std::size_t>(
					std::numeric_limits<std::uint16_t>::max(),1+peers.size()+(interior?2:0)));
				minimum_observed_fan=std::min(minimum_observed_fan,observed_fan);
				maximum_observed_fan=std::max(maximum_observed_fan,observed_fan);
				if(interior||peers.size()>1)
				{
					saw_junction=true;maximum_fan=std::max(maximum_fan,observed_fan);
					used_peers.insert(used_peers.end(),peers.begin(),peers.end());
				}
				else if(peers.size()==1)
				{
					saw_attached=true;maximum_fan=std::max<std::uint16_t>(maximum_fan,2);used_peers.push_back(peers.front());
				}
				else
				{
					// Adjacent triangles of the same open sheet normally meet this free edge
					// at a vertex.  Distance to the adjacent triangle can grow much more
					// slowly than distance travelled along the edge (for example at a
					// tessellation diagonal), so a two-clearance inset falsely labels an
					// ordinary perimeter as ambiguous.  Ignore only a small endpoint
					// neighbourhood; a genuinely near-parallel edge remains inside the
					// clearance band over the interior and therefore still fails closed.
					const double parameter_inset=std::min(0.24*span,std::max(
						512.0*std::numeric_limits<double>::epsilon(),16.0*edge_clearance_tolerance_/edge_length));
					const Vec3d segment_a=a+delta*(lo+parameter_inset),segment_b=a+delta*(hi-parameter_inset);
					double nearest2=std::numeric_limits<double>::infinity();
					for(std::uint32_t candidate:candidate_triangles)if(candidate!=edge.triangle&&candidate<valid_triangle.size()&&
						valid_triangle[candidate]&&!is_endpoint_neighbour(candidate,source_contact_tolerance))
						nearest2=std::min(nearest2,segment_triangle_distance2(segment_a,segment_b,triangles_by_id_[candidate]));
					if(nearest2>=edge_clearance_tolerance_*edge_clearance_tolerance_)saw_clear=true;
					else saw_ambiguous=true;
				}
				if(source_cad_id!=TriMesh::kNoCadEdgeId&&source_cad_incidence>=2)
				{
					bool has_same_cad_peer=false;for(std::size_t peer:peers)
						has_same_cad_peer=has_same_cad_peer||mesh.cad_edge_id(peer/3,static_cast<unsigned>(peer%3))==source_cad_id;
					cad_coverage_ok=cad_coverage_ok&&has_same_cad_peer;
				}
			}
			std::sort(used_peers.begin(),used_peers.end());used_peers.erase(std::unique(used_peers.begin(),used_peers.end()),used_peers.end());
			const bool mixed_open_attachment=saw_clear&&(saw_attached||saw_junction);
			const bool cad_mismatch=source_cad_id!=TriMesh::kNoCadEdgeId&&source_cad_incidence>=2&&
				(!cad_coverage_ok||saw_clear||minimum_observed_fan!=source_cad_incidence||
					maximum_observed_fan!=source_cad_incidence);
			if(saw_ambiguous||mixed_open_attachment||cad_mismatch||(!saw_clear&&!saw_attached&&!saw_junction))
			{
				output.kind=FabricEdgeKind::unknown;output.role=FabricEdgeRole::unknown;output.incident_fan_degree=0;
				if(saw_ambiguous)output.unknown_reason=FabricEdgeUnknownReason::proximity_ambiguity;
				else if(mixed_open_attachment)output.unknown_reason=FabricEdgeUnknownReason::mixed_open_attachment;
				else if(cad_mismatch)output.unknown_reason=FabricEdgeUnknownReason::cad_fan_mismatch;
				else output.unknown_reason=FabricEdgeUnknownReason::no_geometric_evidence;
			}
			else if(saw_junction||source_cad_incidence>2)
			{
				output.kind=FabricEdgeKind::attached;output.role=FabricEdgeRole::junction;
				output.unknown_reason=FabricEdgeUnknownReason::none;
				output.incident_fan_degree=std::max<std::uint16_t>(3,std::max<std::uint16_t>(maximum_fan,
					static_cast<std::uint16_t>(std::min<std::uint32_t>(source_cad_incidence,
						std::numeric_limits<std::uint16_t>::max()))));
			}
			else if(saw_attached)
			{
				output.kind=FabricEdgeKind::attached;output.role=FabricEdgeRole::manifold_seam;output.incident_fan_degree=2;
				output.unknown_reason=FabricEdgeUnknownReason::none;
			}
			else
			{
				output.kind=FabricEdgeKind::confirmed_free;output.role=FabricEdgeRole::none;output.incident_fan_degree=1;
				output.unknown_reason=FabricEdgeUnknownReason::none;
			}
			output.contact_tolerance=used_contact;output.clearance_tolerance=edge_clearance_tolerance_;
			if(output.kind==FabricEdgeKind::attached)
			{
				claimed_peers[edge.flat_index]=std::move(used_peers);
			}
		}

		// Promote a geometrically witnessed attachment onto an indexed reciprocal which was
		// certified before the general search.  This is required when a rib lies exactly on
		// an already-indexed skin diagonal: the rib sees both skin sectors, while the fast
		// indexed certificate initially sees only its same-index partner.
		for(const IndexedHalfEdge& edge:half_edges)
		{
			const auto& output=certificate_for(edge.flat_index);if(output.kind!=FabricEdgeKind::attached)continue;
			for(std::size_t peer:claimed_peers[edge.flat_index])
			{
				if(peer>=half_edge_count||certificate_for(peer).kind!=FabricEdgeKind::attached)continue;
				auto& reverse=claimed_peers[peer];if(std::find(reverse.begin(),reverse.end(),edge.flat_index)==reverse.end())
					reverse.push_back(edge.flat_index);
				auto& peer_output=certificate_for(peer);const std::uint16_t fan=std::max(
					output.incident_fan_degree,peer_output.incident_fan_degree);
				if(fan>2){peer_output.role=FabricEdgeRole::junction;peer_output.incident_fan_degree=fan;}
			}
		}

		// A boundary-edge attachment is a symmetric relation.  Per-edge interval coverage
		// can be asymmetric (for example a short segment covers only part of a long one),
		// and one side can carry contradictory CAD provenance.  Reconcile to a fixed point
		// before issuing any equality token: uncertainty on either side invalidates the
		// claimed pair on both sides.  Interior T-junction witnesses legitimately have no
		// reciprocal half-edge and therefore an empty claim list.
		bool changed=true;
		while(changed)
		{
			changed=false;
			for(const IndexedHalfEdge& edge:half_edges)
			{
				auto& output=certificate_for(edge.flat_index);
				if(output.kind!=FabricEdgeKind::attached)continue;
				bool compatible=true;
				for(std::size_t peer:claimed_peers[edge.flat_index])
				{
					if(peer>=half_edge_count||certificate_for(peer).kind!=FabricEdgeKind::attached||
						std::find(claimed_peers[peer].begin(),claimed_peers[peer].end(),edge.flat_index)
							==claimed_peers[peer].end())
					{
						compatible=false;break;
					}
				}
				if(!compatible)
				{
					output.kind=FabricEdgeKind::unknown;output.role=FabricEdgeRole::unknown;
					output.unknown_reason=FabricEdgeUnknownReason::nonreciprocal_attachment;
					output.incident_fan_degree=0;output.attachment_id=FabricEdgeCertificate::no_attachment;
					changed=true;
				}
			}
		}
		for(const IndexedHalfEdge& edge:half_edges)if(certificate_for(edge.flat_index).kind==FabricEdgeKind::attached)
			for(std::size_t peer:claimed_peers[edge.flat_index])
				if(peer<half_edge_count&&certificate_for(peer).kind==FabricEdgeKind::attached)
					join_attachment(edge.flat_index,peer);

		std::vector<std::uint64_t> minimum_attachment(half_edge_count,FabricEdgeCertificate::no_attachment);
		for(const IndexedHalfEdge& edge:half_edges)if(certificate_for(edge.flat_index).kind==FabricEdgeKind::attached)
		{
			const std::size_t root=attachment_root(edge.flat_index);minimum_attachment[root]=std::min(
				minimum_attachment[root],static_cast<std::uint64_t>(edge.flat_index)+1);
		}
		for(const IndexedHalfEdge& edge:half_edges)if(certificate_for(edge.flat_index).kind==FabricEdgeKind::attached)
			certificate_for(edge.flat_index).attachment_id=minimum_attachment[attachment_root(edge.flat_index)];
	}

	FabricEdgeCertificate TriangleBvh::edge_certificate(std::uint32_t original_triangle_id,
		int local_edge) const
	{
		if(original_triangle_id>=edge_certificate_by_triangle_.size()||local_edge<0||local_edge>=3)
		{
			FabricEdgeCertificate result;result.contact_tolerance=edge_contact_tolerance_;
			result.clearance_tolerance=edge_clearance_tolerance_;return result;
		}
		return edge_certificate_by_triangle_[original_triangle_id][local_edge];
	}

	FabricEdgeKind TriangleBvh::edge_kind(std::uint32_t original_triangle_id,int local_edge) const
	{
		return edge_certificate(original_triangle_id,local_edge).kind;
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
		std::array<std::uint32_t,128> stack{};std::size_t stack_size=1;stack[0]=0;
		while (stack_size)
		{
			const Node& n = nodes_[stack[--stack_size]];
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
			else { stack[stack_size++]=n.left;stack[stack_size++]=n.right; }
		}
		return best;
	}

	bool TriangleBvh::segment_touches_surface(Vec3d a,Vec3d b,double tolerance) const
	{
		if(nodes_.empty()||tolerance<0.0)return false;
		Aabb3d query;query.expand(a);query.expand(b);
		const Vec3d margin{tolerance,tolerance,tolerance};query.lo=query.lo-margin;query.hi=query.hi+margin;
		std::vector<std::uint32_t> candidates;query_aabb(query,candidates);
		const double limit2=tolerance*tolerance;
		for(std::uint32_t id:candidates)if(segment_triangle_distance2(a,b,triangle(id))<=limit2)return true;
		return false;
	}

	NearestSurfacePoint TriangleBvh::nearest(Vec3d p, double max_distance) const
	{
		NearestSurfacePoint best; if (nodes_.empty()) return best;
		double best2 = max_distance * max_distance;
		std::array<std::uint32_t,128> stack{};std::size_t stack_size=1;stack[0]=0;
		while (stack_size)
		{
			const std::uint32_t ni=stack[--stack_size];const Node& n=nodes_[ni];
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
				if(dl<dr){stack[stack_size++]=n.right;stack[stack_size++]=n.left;}else{stack[stack_size++]=n.left;stack[stack_size++]=n.right;}
			}
		}
		return best;
	}

	NearestSurfacePoint TriangleBvh::nearest_on_face(Vec3d p,std::uint32_t source_face_id,double max_distance) const
	{
		NearestSurfacePoint best;if(nodes_.empty())return best;
		double best2=max_distance*max_distance;std::vector<std::uint32_t> stack{0};
		while(!stack.empty())
		{
			const std::uint32_t ni=stack.back();stack.pop_back();const Node& n=nodes_[ni];
			if(n.bounds.distance2(p)>best2)continue;
			if(n.leaf())
			{
				for(std::uint32_t i=n.first;i<n.first+n.count;++i)
				{
					const BvhTriangle& triangle=primitives_[i].tri;if(triangle.source_face_id!=source_face_id)continue;
					const Vec3d q=closest_point_triangle(p,triangle);const double d2=length2(q-p);if(d2>=best2)continue;
					best2=d2;best.found=true;best.point=q;best.distance=std::sqrt(d2);
					best.geometric_normal=normalized(cross(triangle.b-triangle.a,triangle.c-triangle.a));
					best.triangle_id=triangle.triangle_id;best.source_face_id=triangle.source_face_id;
				}
			}
			else{stack.push_back(n.left);stack.push_back(n.right);}
		}
		return best;
	}

	double TriangleBvh::distance(Vec3d p, double max_distance) const { return nearest(p, max_distance).distance; }
}
