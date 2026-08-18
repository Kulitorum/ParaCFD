#include "core/geometry/local_surface_arrangement.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <utility>

namespace paracfd::core
{
	namespace
	{
		using Polygon = std::vector<Vec3d>;

		struct Plane
		{
			Vec3d normal{};
			double offset = 0.0;
			double distance(Vec3d point) const { return dot(normal, point) - offset; }
		};

		struct Polyhedron
		{
			std::vector<LocalArrangementFace> faces;
			std::vector<std::int8_t> side;
		};

		struct Measure
		{
			double value = 0.0;
			Vec3d centroid{};
		};

		struct PreparedTriangle
		{
			LocalSurfaceTriangle input;
			int input_index = -1;
			Polygon polygon;
			Vec3d normal{};
			double area = 0.0;
			Vec3d centroid{};
			int support_plane = -1;
			int boundary_face = -1;
		};

		struct AtomicPatch
		{
			int triangle = -1;
			int canonical_minus_atom = -1;
			int canonical_plus_atom = -1;
			Polygon polygon;
			double area = 0.0;
			Vec3d centroid{};
			double parent_facet_area = 0.0;
			double parent_fabric_area = 0.0;
			int support_plane = -1;
			// An arithmetic-scale carrier adjacent to a certified attached seam.  Its
			// finite area is retained, but its local atom pair can be a roundoff overlap
			// on the wrong sector of the seam; topology is inherited from a resolved
			// portion of the same source triangle below.
			bool attached_seam_carrier = false;
		};

		struct DisjointSet
		{
			explicit DisjointSet(int count) : parent(count), rank(count, 0) { std::iota(parent.begin(), parent.end(), 0); }
			int root(int value)
			{
				while(parent[value] != value) { parent[value] = parent[parent[value]]; value = parent[value]; }
				return value;
			}
			void join(int a, int b)
			{
				a = root(a); b = root(b); if(a == b) return;
				if(rank[a] < rank[b]) std::swap(a, b);
				parent[b] = a; if(rank[a] == rank[b]) ++rank[a];
			}
			std::vector<int> parent, rank;
		};

		constexpr std::uint64_t no_attachment=std::numeric_limits<std::uint64_t>::max();

		Vec3d triangle_vertex(const LocalSurfaceTriangle& triangle,int vertex)
		{
			if(vertex==0)return triangle.a;
			if(vertex==1)return triangle.b;
			return triangle.c;
		}

		Vec3d& triangle_vertex(LocalSurfaceTriangle& triangle,int vertex)
		{
			if(vertex==0)return triangle.a;
			if(vertex==1)return triangle.b;
			return triangle.c;
		}

		// CAD faces are tessellated independently, so two half-edges certified as the
		// same attachment can differ by a few representable coordinates.  Exact local
		// topology must use one geometric seam or that harmless CAD discrepancy becomes
		// a microscopic fluid passage.  Canonicalise only equality classes already
		// certified by the BVH: proximity by itself never closes an opening.
		bool canonicalize_certified_attachments(
			std::vector<LocalSurfaceTriangle>& triangles,double contact,double roundoff,
			std::string& error)
		{
			struct EdgeOccurrence{int triangle=-1,edge=-1;};
			struct Endpoint{int node=-1;double parameter=0.0;Vec3d projected{};};
			std::map<std::uint64_t,std::vector<EdgeOccurrence>> groups;
			for(int triangle=0;triangle<static_cast<int>(triangles.size());++triangle)
				for(int edge=0;edge<3;++edge)
					if(triangles[triangle].edge_kind[edge]==FabricEdgeKind::attached&&
						triangles[triangle].edge_attachment_id[edge]!=no_attachment)
						groups[triangles[triangle].edge_attachment_id[edge]].push_back({triangle,edge});

			DisjointSet equivalent_vertices(3*static_cast<int>(triangles.size()));
			std::vector<Vec3d> proposal_sum(3*triangles.size());
			std::vector<int> proposal_count(3*triangles.size(),0);
			const double tolerance=contact+16.0*roundoff;
			const double tolerance2=tolerance*tolerance;
			for(const auto& group:groups)
			{
				if(group.second.size()<2)continue;
				Vec3d origin{},direction{};double longest2=0.0;
				for(const EdgeOccurrence occurrence:group.second)
				{
					const auto& triangle=triangles[occurrence.triangle];
					const Vec3d a=triangle_vertex(triangle,occurrence.edge);
					const Vec3d b=triangle_vertex(triangle,(occurrence.edge+1)%3);
					if(length2(b-a)>longest2){longest2=length2(b-a);origin=a;direction=b-a;}
				}
				if(!(longest2>tolerance2))continue;
				direction=direction/std::sqrt(longest2);
				// Use the mean certified line rather than selecting either CAD face as
				// authoritative.  This keeps the result invariant to triangle order and,
				// for a constant-offset seam, moves both endpoints by the same amount.
				origin={};int origin_count=0;
				for(const EdgeOccurrence occurrence:group.second)
				{
					const auto& triangle=triangles[occurrence.triangle];
					origin=origin+triangle_vertex(triangle,occurrence.edge)+
						triangle_vertex(triangle,(occurrence.edge+1)%3);
					origin_count+=2;
				}
				origin=origin/static_cast<double>(origin_count);
				std::vector<Endpoint> endpoints;endpoints.reserve(2*group.second.size());
				double maximum_line_distance2=0.0;
				for(const EdgeOccurrence occurrence:group.second)
					for(int endpoint=0;endpoint<2;++endpoint)
					{
						const int vertex=(occurrence.edge+endpoint)%3;
						const Vec3d point=triangle_vertex(triangles[occurrence.triangle],vertex);
						const double parameter=dot(point-origin,direction);
						const Vec3d projected=origin+direction*parameter;
						const double line_distance2=length2(point-projected);
						maximum_line_distance2=std::max(maximum_line_distance2,line_distance2);
						if(line_distance2>tolerance2)
						{
							error="certified attachment token contains geometrically inconsistent half-edges";
							return false;
						}
						const int node=3*occurrence.triangle+vertex;
						endpoints.push_back({node,parameter,projected});
					}
				// Exact indexed edges and already-conforming 1:N subdivisions require no
				// movement.  Still union their endpoint occurrences so a correction from
				// another incident seam propagates across triangle-local vertex copies.
				const double movement_threshold=0.25*roundoff;
				if(maximum_line_distance2>movement_threshold*movement_threshold)
					for(const Endpoint& endpoint:endpoints)
					{
						proposal_sum[endpoint.node]=proposal_sum[endpoint.node]+endpoint.projected;
						++proposal_count[endpoint.node];
					}
				std::sort(endpoints.begin(),endpoints.end(),[](const Endpoint& a,const Endpoint& b)
				{
					if(a.parameter!=b.parameter)return a.parameter<b.parameter;
					return a.node<b.node;
				});
				for(std::size_t begin=0;begin<endpoints.size();)
				{
					std::size_t end=begin+1;
					while(end<endpoints.size()&&
						endpoints[end].parameter-endpoints[begin].parameter<=tolerance)++end;
					for(std::size_t index=begin+1;index<end;++index)
						equivalent_vertices.join(endpoints[begin].node,endpoints[index].node);
					begin=end;
				}
			}

			std::vector<Vec3d> root_sum(3*triangles.size());
			std::vector<int> root_count(3*triangles.size(),0);
			for(int node=0;node<static_cast<int>(proposal_count.size());++node)
				if(proposal_count[node]>0)
				{
					const int root=equivalent_vertices.root(node);
					root_sum[root]=root_sum[root]+proposal_sum[node];
					root_count[root]+=proposal_count[node];
				}
			for(int node=0;node<static_cast<int>(proposal_count.size());++node)
			{
				const int root=equivalent_vertices.root(node);
				if(root_count[root]>0)
				{
					triangle_vertex(triangles[node/3],node%3)=
						root_sum[root]/static_cast<double>(root_count[root]);
				}
			}
			return true;
		}

		double maximum_coordinate(const Aabb3d& box)
		{
			return std::max({1.0, std::abs(box.lo.x), std::abs(box.lo.y), std::abs(box.lo.z),
				std::abs(box.hi.x), std::abs(box.hi.y), std::abs(box.hi.z),
				box.hi.x-box.lo.x, box.hi.y-box.lo.y, box.hi.z-box.lo.z});
		}

		double polygon_coordinate_scale(const Polygon& polygon)
		{
			double scale=1.0;
			for(Vec3d point:polygon)
				scale=std::max({scale,std::abs(point.x),std::abs(point.y),std::abs(point.z)});
			return scale;
		}

		std::string diagnostic_number(double value)
		{
			std::ostringstream stream;
			stream<<std::scientific<<std::setprecision(17)<<value;
			return stream.str();
		}

		Plane canonical_plane(Vec3d normal, double offset)
		{
			const double magnitude = std::sqrt(length2(normal));
			if(!(magnitude > 0.0)) return {};
			normal = normal/magnitude;
			offset /= magnitude;
			// A largest-component pivot is discontinuous at 45 degrees: roundoff can
			// make n and -n select different components and therefore different signs.
			// A fixed lexicographic component order with a tiny zero band is tie-free.
			const double zero = 64.0*std::numeric_limits<double>::epsilon();
			int pivot = 2;
			if(std::abs(normal.x) > zero) pivot = 0;
			else if(std::abs(normal.y) > zero) pivot = 1;
			if(normal[pivot] < 0.0) { normal = normal * -1.0; offset = -offset; }
			return {normal, offset};
		}

		double snapped_distance(double distance,double snap_tolerance)
		{
			return std::abs(distance)<=snap_tolerance?0.0:distance;
		}

		Polygon clip_polygon(const Polygon& input,const Plane& plane,int keep_side,
			double boundary_tolerance,double snap_tolerance,double merge_tolerance)
		{
			Polygon output;
			if(input.empty()) return output;
			// If a tolerant half-space is requested, intersect its shifted boundary—not
			// d=0.  This keeps every transition parameter inside [0,1].  BSP geometry
			// passes zero here so its two half-clips share one exact boundary.
			const double boundary=keep_side<0?boundary_tolerance:-boundary_tolerance;
			auto inside = [&](double shifted_distance)
				{return keep_side<0?shifted_distance<=0.0:shifted_distance>=0.0;};
			auto append = [&](Vec3d point)
			{
				if(output.empty() || length2(output.back()-point) > merge_tolerance*merge_tolerance)
					output.push_back(point);
			};
			for(std::size_t index = 0; index < input.size(); ++index)
			{
				const Vec3d a = input[index], b = input[(index+1)%input.size()];
				const double raw_da=plane.distance(a)-boundary,raw_db=plane.distance(b)-boundary;
				const double da=snapped_distance(raw_da,snap_tolerance);
				const double db=snapped_distance(raw_db,snap_tolerance);
				const Vec3d snapped_a=da==0.0?a-plane.normal*raw_da:a;
				const Vec3d snapped_b=db==0.0?b-plane.normal*raw_db:b;
				const bool ia = inside(da), ib = inside(db);
				if(ia) append(snapped_a);
				if(ia != ib)
				{
					if(da == 0.0) append(snapped_a);
					else if(db == 0.0) append(snapped_b);
					else
					{
						const double denominator = da-db;
						if(std::abs(denominator) > std::numeric_limits<double>::min())
						{
							const double parameter = std::clamp(da/denominator,0.0,1.0);
							append(a+(b-a)*parameter);
						}
					}
				}
			}
			if(output.size() > 1 && length2(output.front()-output.back()) <= merge_tolerance*merge_tolerance)
				output.pop_back();
			return output;
		}

		void append_unique(std::vector<Vec3d>& points, Vec3d point, double tolerance)
		{
			const double limit = tolerance*tolerance;
			for(Vec3d existing : points) if(length2(existing-point) <= limit) return;
			points.push_back(point);
		}

		Polygon ordered_polygon(std::vector<Vec3d> points, Vec3d normal, double tolerance)
		{
			Polygon unique;
			for(Vec3d point : points) append_unique(unique, point, tolerance);
			if(unique.size() < 3) return {};
			normal = normalized(normal);
			const Vec3d seed = std::abs(normal.x) < 0.8 ? Vec3d{1,0,0} : Vec3d{0,1,0};
			const Vec3d u = normalized(cross(normal, seed)), v = cross(normal, u);
			std::sort(unique.begin(), unique.end(), [&](Vec3d a, Vec3d b)
			{
				const double au=dot(a,u),bu=dot(b,u);
				if(au!=bu)return au<bu;
				return dot(a,v)<dot(b,v);
			});
			auto turn=[&](Vec3d a,Vec3d b,Vec3d c)
			{
				return dot(cross(b-a,c-a),normal);
			};
			auto append_hull=[&](Polygon& hull,Vec3d point)
			{
				while(hull.size()>=2)
				{
					const Vec3d a=hull[hull.size()-2],b=hull.back();
					const double scale=std::sqrt(length2(b-a))+std::sqrt(length2(point-b));
					if(turn(a,b,point)>tolerance*scale)break;
					hull.pop_back();
				}
				hull.push_back(point);
			};
			Polygon lower,upper;
			for(Vec3d point:unique)append_hull(lower,point);
			for(auto iterator=unique.rbegin();iterator!=unique.rend();++iterator)
				append_hull(upper,*iterator);
			if(lower.size()<2||upper.size()<2)return {};
			lower.pop_back();upper.pop_back();
			lower.insert(lower.end(),upper.begin(),upper.end());
			if(lower.size()<3)return {};
			if(dot(cross(lower[1]-lower[0],lower[2]-lower[0]),normal)<0.0)
				std::reverse(lower.begin(),lower.end());
			return lower;
		}

		Measure polygon_measure(const Polygon& polygon, Vec3d normal)
		{
			Measure result;
			if(polygon.size() < 3) return result;
			normal = normalized(normal);
			const Vec3d reference=polygon[0];
			for(std::size_t index = 1; index+1 < polygon.size(); ++index)
			{
				const Vec3d first=polygon[index]-reference,second=polygon[index+1]-reference;
				const double area = 0.5*std::abs(dot(cross(first,second),normal));
				if(!(area > 0.0)) continue;
				result.value += area;
				result.centroid = result.centroid+(first+second)*(area/3.0);
			}
			if(result.value > 0.0) result.centroid = reference+result.centroid/result.value;
			return result;
		}

		double polygon_area_3d(const Polygon& polygon)
		{
			if(polygon.size() < 3) return 0.0;
			double area = 0.0;
			for(std::size_t index=1;index+1<polygon.size();++index)
				area+=0.5*std::sqrt(length2(cross(polygon[index]-polygon[0],polygon[index+1]-polygon[0])));
			return area;
		}

		bool point_on_convex_polygon(Vec3d point,const Polygon& polygon,Vec3d normal,double tolerance)
		{
			if(polygon.size()<3)return false;
			normal=normalized(normal);
			if(std::abs(dot(normal,point-polygon[0]))>tolerance)return false;
			bool positive=false,negative=false;
			for(std::size_t index=0;index<polygon.size();++index)
			{
				const Vec3d a=polygon[index],b=polygon[(index+1)%polygon.size()];
				const Vec3d edge=b-a;const double edge_length=std::sqrt(length2(edge));
				if(!(edge_length>0.0))continue;
				const double side=dot(cross(edge,point-a),normal);
				const double limit=tolerance*edge_length;
				positive=positive||side>limit;negative=negative||side<-limit;
				if(positive&&negative)return false;
			}
			return true;
		}

		Measure polyhedron_measure(const Polyhedron& polyhedron)
		{
			Measure result;
			std::vector<Vec3d> vertices;
			for(const auto& face : polyhedron.faces) for(Vec3d point : face.polygon) vertices.push_back(point);
			if(vertices.empty()) return result;
			const Vec3d anchor=vertices.front();Vec3d reference_offset{};
			for(Vec3d point : vertices) reference_offset = reference_offset+(point-anchor);
			const Vec3d reference=anchor+reference_offset/static_cast<double>(vertices.size());
			long double total_volume=0.0L,moment_x=0.0L,moment_y=0.0L,moment_z=0.0L;
			for(const auto& face : polyhedron.faces)
			{
				for(std::size_t index = 1; index+1 < face.polygon.size(); ++index)
				{
					const Vec3d first=face.polygon[0]-reference;
					const Vec3d second=face.polygon[index]-reference;
					const Vec3d third=face.polygon[index+1]-reference;
					const long double triple=
						static_cast<long double>(first.x)*(static_cast<long double>(second.y)*third.z-static_cast<long double>(second.z)*third.y)-
						static_cast<long double>(first.y)*(static_cast<long double>(second.x)*third.z-static_cast<long double>(second.z)*third.x)+
						static_cast<long double>(first.z)*(static_cast<long double>(second.x)*third.y-static_cast<long double>(second.y)*third.x);
					const long double volume=std::abs(triple)/6.0L;
					if(!(volume>0.0L))continue;
					total_volume+=volume;
					moment_x+=(static_cast<long double>(first.x)+second.x+third.x)*(volume/4.0L);
					moment_y+=(static_cast<long double>(first.y)+second.y+third.y)*(volume/4.0L);
					moment_z+=(static_cast<long double>(first.z)+second.z+third.z)*(volume/4.0L);
				}
			}
			result.value=static_cast<double>(total_volume);
			if(total_volume>0.0L)result.centroid=reference+Vec3d{
				static_cast<double>(moment_x/total_volume),static_cast<double>(moment_y/total_volume),
				static_cast<double>(moment_z/total_volume)};
			return result;
		}

		std::vector<Vec3d> polyhedron_vertices(const Polyhedron& polyhedron, double tolerance)
		{
			std::vector<Vec3d> vertices;
			for(const auto& face : polyhedron.faces) for(Vec3d point : face.polygon) append_unique(vertices,point,tolerance);
			return vertices;
		}

		int boundary_plane_id(int axis, bool upper) { return -(2*axis+(upper?2:1)); }

		bool decode_boundary_plane(int plane_id, int& axis, bool& upper)
		{
			if(plane_id >= 0 || plane_id < -6) return false;
			const int value = -plane_id-1;
			axis = value/2; upper = (value%2) != 0; return true;
		}

		Polyhedron cube_polyhedron(const Aabb3d& box, int plane_count)
		{
			const double x0=box.lo.x,y0=box.lo.y,z0=box.lo.z,x1=box.hi.x,y1=box.hi.y,z1=box.hi.z;
			Polyhedron result;
			result.side.assign(plane_count,0);
			result.faces = {
				{{{x0,y0,z0},{x0,y0,z1},{x0,y1,z1},{x0,y1,z0}},boundary_plane_id(0,false)},
				{{{x1,y0,z0},{x1,y1,z0},{x1,y1,z1},{x1,y0,z1}},boundary_plane_id(0,true)},
				{{{x0,y0,z0},{x1,y0,z0},{x1,y0,z1},{x0,y0,z1}},boundary_plane_id(1,false)},
				{{{x0,y1,z0},{x0,y1,z1},{x1,y1,z1},{x1,y1,z0}},boundary_plane_id(1,true)},
				{{{x0,y0,z0},{x0,y1,z0},{x1,y1,z0},{x1,y0,z0}},boundary_plane_id(2,false)},
				{{{x0,y0,z1},{x1,y0,z1},{x1,y1,z1},{x0,y1,z1}},boundary_plane_id(2,true)}};
			return result;
		}

		Polygon polyhedron_plane_cap(const Polyhedron& input,const Plane& plane,
			double tolerance,double snap_tolerance)
		{
			std::vector<Vec3d> cut_points;
			for(const auto& face:input.faces)
				for(std::size_t index=0;index<face.polygon.size();++index)
				{
					const Vec3d a=face.polygon[index],b=face.polygon[(index+1)%face.polygon.size()];
					const double raw_da=plane.distance(a),raw_db=plane.distance(b);
					const double da=snapped_distance(raw_da,snap_tolerance);
					const double db=snapped_distance(raw_db,snap_tolerance);
					if(da==0.0)append_unique(cut_points,a-plane.normal*raw_da,tolerance);
					if(da*db<0.0)
					{
						const double parameter=std::clamp(da/(da-db),0.0,1.0);
						append_unique(cut_points,a+(b-a)*parameter,tolerance);
					}
				}
			return ordered_polygon(std::move(cut_points),plane.normal,tolerance);
		}

		Polyhedron clip_polyhedron(const Polyhedron& input, const Plane& plane, int keep_side,
			int plane_id,std::uint64_t cap_lineage,const Polygon& shared_cap,double tolerance,
			double snap_tolerance,double area_tolerance)
		{
			Polyhedron output; output.side = input.side; output.side[plane_id] = static_cast<std::int8_t>(keep_side);
			for(const auto& face : input.faces)
			{
				Polygon clipped=clip_polygon(face.polygon,plane,keep_side,0.0,
					snap_tolerance,tolerance);
				if(clipped.size()>=3&&polygon_area_3d(clipped)>area_tolerance)
					output.faces.push_back({std::move(clipped),face.plane_id,face.lineage});
			}
			Polygon cap=shared_cap;
			if(keep_side>0)std::reverse(cap.begin(),cap.end());
			if(cap.size()>=3&&polygon_measure(cap,plane.normal).value>area_tolerance)
				output.faces.push_back({std::move(cap),plane_id,cap_lineage});
			return output;
		}

		Polygon clip_triangle_to_box(const LocalSurfaceTriangle& triangle,const Aabb3d& box,
			double snap_tolerance,double merge_tolerance)
		{
			Polygon polygon{triangle.a,triangle.b,triangle.c};
			const Plane planes[6]={{{-1,0,0},-box.lo.x},{{1,0,0},box.hi.x},
				{{0,-1,0},-box.lo.y},{{0,1,0},box.hi.y},{{0,0,-1},-box.lo.z},{{0,0,1},box.hi.z}};
			for(const Plane& plane:planes) polygon=clip_polygon(polygon,plane,-1,0.0,
				snap_tolerance,merge_tolerance);
			return polygon;
		}

		bool clip_segment_to_box(Vec3d a,Vec3d b,const Aabb3d& box,double tolerance,
			Vec3d& clipped_a,Vec3d& clipped_b)
		{
			const Vec3d delta=b-a;double first=0.0,last=1.0;
			for(int axis=0;axis<3;++axis)
			{
				const double lower=box.lo[axis]-tolerance,upper=box.hi[axis]+tolerance;
				if(std::abs(delta[axis])<std::numeric_limits<double>::min())
				{
					if(a[axis]<lower||a[axis]>upper)return false;
					continue;
				}
				double enter=(lower-a[axis])/delta[axis],leave=(upper-a[axis])/delta[axis];
				if(enter>leave)std::swap(enter,leave);
				first=std::max(first,enter);last=std::min(last,leave);
				if(first>last)return false;
			}
			clipped_a=a+delta*first;clipped_b=a+delta*last;
			return true;
		}

		Vec3d closest_point_on_triangle(Vec3d point,const LocalSurfaceTriangle& triangle)
		{
			const Vec3d ab=triangle.b-triangle.a,ac=triangle.c-triangle.a,ap=point-triangle.a;
			const double d1=dot(ab,ap),d2=dot(ac,ap);
			if(d1<=0.0&&d2<=0.0)return triangle.a;
			const Vec3d bp=point-triangle.b;const double d3=dot(ab,bp),d4=dot(ac,bp);
			if(d3>=0.0&&d4<=d3)return triangle.b;
			const double vc=d1*d4-d3*d2;
			if(vc<=0.0&&d1>=0.0&&d3<=0.0)return triangle.a+ab*(d1/(d1-d3));
			const Vec3d cp=point-triangle.c;const double d5=dot(ab,cp),d6=dot(ac,cp);
			if(d6>=0.0&&d5<=d6)return triangle.c;
			const double vb=d5*d2-d1*d6;
			if(vb<=0.0&&d2>=0.0&&d6<=0.0)return triangle.a+ac*(d2/(d2-d6));
			const double va=d3*d6-d5*d4;
			if(va<=0.0&&d4-d3>=0.0&&d5-d6>=0.0)
				return triangle.b+(triangle.c-triangle.b)*((d4-d3)/((d4-d3)+(d5-d6)));
			const double inverse=1.0/(va+vb+vc);
			return triangle.a+ab*(vb*inverse)+ac*(vc*inverse);
		}

		std::array<Vec3d,8> box_corners(const Aabb3d& box)
		{
			return {{{box.lo.x,box.lo.y,box.lo.z},{box.hi.x,box.lo.y,box.lo.z},
				{box.lo.x,box.hi.y,box.lo.z},{box.hi.x,box.hi.y,box.lo.z},
				{box.lo.x,box.lo.y,box.hi.z},{box.hi.x,box.lo.y,box.hi.z},
				{box.lo.x,box.hi.y,box.hi.z},{box.hi.x,box.hi.y,box.hi.z}}};
		}

		bool plane_cuts_box(const Plane& plane,const Aabb3d& box,double tolerance)
		{
			double minimum=std::numeric_limits<double>::infinity(),maximum=-minimum;
			for(Vec3d point:box_corners(box)){const double value=plane.distance(point);minimum=std::min(minimum,value);maximum=std::max(maximum,value);}
			return minimum < -tolerance && maximum > tolerance;
		}

		int coplanar_boundary_face(const LocalSurfaceTriangle& triangle,const Polygon& polygon,
			Vec3d triangle_normal,const Aabb3d& box,double distance_tolerance,
			double angular_tolerance)
		{
			for(int axis=0;axis<3;++axis) for(int side=0;side<2;++side)
			{
				Vec3d boundary_normal{};
				boundary_normal[axis]=1.0;
				const double orientation=std::abs(dot(triangle_normal,boundary_normal));
				const double angular_distance=std::sqrt(std::max(0.0,2.0*(1.0-orientation)));
				if(angular_distance>angular_tolerance)continue;
				const double coordinate=side?box.hi[axis]:box.lo[axis];
				if(std::abs(triangle.a[axis]-coordinate)>distance_tolerance)continue;
				bool all=true;
				for(Vec3d point:polygon)
					if(std::abs(point[axis]-coordinate)>distance_tolerance){all=false;break;}
				if(all)return 2*axis+side;
			}
			return -1;
		}

		// A two-component floating expansion is sufficient for the sign of the 2-D
		// orientation predicates used below.  It evaluates the determinant of the
		// input binary64 coordinates at roughly double precision without introducing a
		// third-party exact-geometry dependency.  Intersection coordinates remain
		// binary64; only the combinatorial inside/outside decision needs the wider type.
		struct DoubleDouble{double hi=0.0,lo=0.0;};

		DoubleDouble two_sum(double a,double b)
		{
			const double sum=a+b;
			const double b_virtual=sum-a;
			return {sum,(a-(sum-b_virtual))+(b-b_virtual)};
		}

		DoubleDouble two_difference(double a,double b)
		{
			const double difference=a-b;
			const double b_virtual=a-difference;
			return {difference,(a-(difference+b_virtual))+(b_virtual-b)};
		}

		DoubleDouble two_product(double a,double b)
		{
			const double product=a*b;
			return {product,std::fma(a,b,-product)};
		}

		DoubleDouble add(DoubleDouble a,DoubleDouble b)
		{
			DoubleDouble sum=two_sum(a.hi,b.hi);
			DoubleDouble correction=two_sum(sum.lo,a.lo+b.lo);
			DoubleDouble result=two_sum(sum.hi,correction.hi);
			result.lo+=correction.lo;
			return two_sum(result.hi,result.lo);
		}

		DoubleDouble subtract(DoubleDouble a,DoubleDouble b)
		{
			return add(a,{-b.hi,-b.lo});
		}

		DoubleDouble multiply(DoubleDouble a,DoubleDouble b)
		{
			DoubleDouble result=two_product(a.hi,b.hi);
			result=add(result,two_product(a.hi,b.lo));
			result=add(result,two_product(a.lo,b.hi));
			return add(result,two_product(a.lo,b.lo));
		}

		DoubleDouble orientation_2d(Vec3d a,Vec3d b,Vec3d point,int axis0,int axis1)
		{
			const DoubleDouble ab0=two_difference(b[axis0],a[axis0]);
			const DoubleDouble ab1=two_difference(b[axis1],a[axis1]);
			const DoubleDouble ap0=two_difference(point[axis0],a[axis0]);
			const DoubleDouble ap1=two_difference(point[axis1],a[axis1]);
			return subtract(multiply(ab0,ap1),multiply(ab1,ap0));
		}

		int expansion_sign(DoubleDouble value)
		{
			if(value.hi>0.0)return 1;
			if(value.hi<0.0)return -1;
			return value.lo>0.0?1:value.lo<0.0?-1:0;
		}

		double expansion_value(DoubleDouble value){return value.hi+value.lo;}

		Polygon intersect_coplanar_convex(const Polygon& subject,const Polygon& clip,Vec3d normal,double tolerance)
		{
			if(subject.size()<3||clip.size()<3||length2(normal)==0.0)return {};
			normal=normalized(normal);
			int drop=0;
			if(std::abs(normal.y)>std::abs(normal[drop]))drop=1;
			if(std::abs(normal.z)>std::abs(normal[drop]))drop=2;
			const int axis0=(drop+1)%3,axis1=(drop+2)%3;
			int clip_orientation=0;
			for(std::size_t index=0;index<clip.size()&&clip_orientation==0;++index)
				clip_orientation=expansion_sign(orientation_2d(clip[index],
					clip[(index+1)%clip.size()],clip[(index+2)%clip.size()],axis0,axis1));
			if(clip_orientation==0)return {};
			Polygon output=subject;
			for(std::size_t edge=0;edge<clip.size()&&!output.empty();++edge)
			{
				const Vec3d clip_a=clip[edge],clip_b=clip[(edge+1)%clip.size()];
				Polygon next;next.reserve(output.size()+2);
				auto append=[&](Vec3d point)
				{
					if(next.empty()||length2(next.back()-point)>tolerance*tolerance)
						next.push_back(point);
				};
				for(std::size_t index=0;index<output.size();++index)
				{
					const Vec3d a=output[index],b=output[(index+1)%output.size()];
					const DoubleDouble raw_a=orientation_2d(clip_a,clip_b,a,axis0,axis1);
					const DoubleDouble raw_b=orientation_2d(clip_a,clip_b,b,axis0,axis1);
					const int sign_a=clip_orientation*expansion_sign(raw_a);
					const int sign_b=clip_orientation*expansion_sign(raw_b);
					const bool inside_a=sign_a>=0,inside_b=sign_b>=0;
					if(inside_a)append(a);
					if(inside_a!=inside_b)
					{
						const double da=expansion_value(raw_a),db=expansion_value(raw_b);
						const double denominator=da-db;
						if(denominator!=0.0)
							append(a+(b-a)*std::clamp(da/denominator,0.0,1.0));
					}
				}
				if(next.size()>1&&length2(next.front()-next.back())<=tolerance*tolerance)
					next.pop_back();
				output=std::move(next);
			}
			return output.size()>=3?output:Polygon{};
		}

		// Split one convex polygon by the infinite line through line_a/line_b in its
		// support plane.  Both children receive the very same computed intersection
		// vertices.  That detail is important here: independently clipping the two
		// half-planes can leave roundoff-width gaps, which would turn a tessellation
		// edge into a physical vent (or close a real narrow vent after common refinement).
		std::pair<Polygon,Polygon> split_coplanar_convex_by_line(const Polygon& input,
			Vec3d line_a,Vec3d line_b,Vec3d normal,double merge_tolerance)
		{
			if(input.size()<3||length2(normal)==0.0||length2(line_b-line_a)==0.0)
				return {input,{}};
			normal=normalized(normal);
			int drop=0;
			if(std::abs(normal.y)>std::abs(normal[drop]))drop=1;
			if(std::abs(normal.z)>std::abs(normal[drop]))drop=2;
			const int axis0=(drop+1)%3,axis1=(drop+2)%3;
			Polygon negative,positive;
			negative.reserve(input.size()+2);positive.reserve(input.size()+2);
			auto append=[&](Polygon& polygon,Vec3d point)
			{
				if(polygon.empty()||length2(polygon.back()-point)>
					merge_tolerance*merge_tolerance)polygon.push_back(point);
			};
			for(std::size_t index=0;index<input.size();++index)
			{
				const Vec3d a=input[index],b=input[(index+1)%input.size()];
				const DoubleDouble raw_a=orientation_2d(line_a,line_b,a,axis0,axis1);
				const DoubleDouble raw_b=orientation_2d(line_a,line_b,b,axis0,axis1);
				const int sign_a=expansion_sign(raw_a),sign_b=expansion_sign(raw_b);
				if(sign_a<=0)append(negative,a);
				if(sign_a>=0)append(positive,a);
				if(sign_a*sign_b<0)
				{
					const double da=expansion_value(raw_a),db=expansion_value(raw_b);
					const double denominator=da-db;
					if(denominator!=0.0)
					{
						const Vec3d intersection=a+(b-a)*std::clamp(da/denominator,0.0,1.0);
						append(negative,intersection);append(positive,intersection);
					}
				}
			}
			auto finish=[&](Polygon& polygon)
			{
				if(polygon.size()>1&&length2(polygon.front()-polygon.back())<=
					merge_tolerance*merge_tolerance)polygon.pop_back();
				if(polygon.size()<3)polygon.clear();
			};
			finish(negative);finish(positive);
			return {std::move(negative),std::move(positive)};
		}

		const LocalArrangementFace* face_on_plane(const LocalArrangementAtom& atom,int plane)
		{
			for(const auto& face : atom.faces)
				if(face.plane_id == plane) return &face;
			return nullptr;
		}

		void fail(LocalSurfaceArrangement& result,std::string message)
		{
			result.valid=false;result.error=std::move(message);
		}
	}

	LocalSurfaceArrangement build_local_surface_arrangement(const Aabb3d& cell,
		const std::vector<LocalSurfaceTriangle>& triangles,const LocalSurfaceArrangementOptions& options)
	{
		LocalSurfaceArrangement result;result.cell=cell;
		if(!cell.valid()){fail(result,"local arrangement cell AABB is invalid");return result;}
		const Vec3d extent=cell.hi-cell.lo;const double h=std::max({extent.x,extent.y,extent.z});
		if(!(extent.x>0&&extent.y>0&&extent.z>0)){fail(result,"local arrangement cell has zero extent");return result;}
		const double scale=maximum_coordinate(cell),machine=std::numeric_limits<double>::epsilon();
		const double local_scale=std::max(1.0,h);
		const double clipping_roundoff=std::max(64.0*machine*local_scale,8.0*machine*scale);
		const double default_contact=std::max(128.0*machine*local_scale,16.0*machine*scale);
		const double contact=options.contact_tolerance>0?options.contact_tolerance:default_contact;
		const double ambiguity=options.ambiguity_tolerance>0?options.ambiguity_tolerance:contact;
		const double angular_contact=options.angular_contact_tolerance>0?options.angular_contact_tolerance:256.0*machine;
		const double angular_ambiguity=options.angular_ambiguity_tolerance>0?options.angular_ambiguity_tolerance:angular_contact;
		const double cartesian_face=options.cartesian_face_tolerance>0?
			options.cartesian_face_tolerance:clipping_roundoff;
		const double angular_cartesian_face=options.angular_cartesian_face_tolerance>0?
			options.angular_cartesian_face_tolerance:256.0*machine;
		const double minimum_face_area=std::min({extent.x*extent.y,extent.x*extent.z,extent.y*extent.z});
		const double minimum_extent=std::min({extent.x,extent.y,extent.z});
		const double area_tolerance=std::max(options.minimum_area_fraction*minimum_face_area,contact*contact);
		const double volume_tolerance=std::max(options.minimum_volume_fraction*extent.x*extent.y*extent.z,contact*contact*contact);
		// Plane arrangements must remain conforming on both sides of every facet.  The
		// product-level small-fragment merger handles physically tiny control volumes
		// later; dropping one half of a grazing split here makes neighbouring atoms retain
		// different facet subdivisions and creates a topological crack.  Use only an
		// arithmetic underflow guard while constructing the exact complex.
		const double arithmetic_volume_tolerance=std::max(
			std::numeric_limits<double>::denorm_min(),
			64.0*std::numeric_limits<double>::epsilon()*
			std::numeric_limits<double>::epsilon()*extent.x*extent.y*extent.z);
		if(!(contact>=0&&contact<0.25*minimum_extent&&ambiguity>=contact&&
			angular_contact>=0&&angular_ambiguity>=angular_contact&&
			cartesian_face>=0&&cartesian_face<0.25*minimum_extent&&
			angular_cartesian_face>=0)||
			options.maximum_planes<1||options.maximum_atoms<1)
		{fail(result,"local arrangement options are invalid");return result;}

		result.contact_tolerance=contact;result.ambiguity_tolerance=ambiguity;
		result.cartesian_face_tolerance=cartesian_face;
		result.angular_cartesian_face_tolerance=angular_cartesian_face;
		std::vector<LocalSurfaceTriangle> canonical_triangles=triangles;
		std::string canonicalization_error;
		if(!canonicalize_certified_attachments(canonical_triangles,contact,clipping_roundoff,
			canonicalization_error))
		{
			fail(result,std::move(canonicalization_error));return result;
		}
		std::vector<PreparedTriangle> prepared;prepared.reserve(canonical_triangles.size());
		std::vector<int> prepared_by_input(canonical_triangles.size(),-1);
		for(int input_index=0;input_index<static_cast<int>(canonical_triangles.size());++input_index)
		{
			const auto& input=canonical_triangles[input_index];const Vec3d raw_normal=cross(input.b-input.a,input.c-input.a);
			if(length2(raw_normal)<=area_tolerance*area_tolerance)continue;
			PreparedTriangle triangle;triangle.input=input;triangle.input_index=input_index;triangle.normal=normalized(raw_normal);
			triangle.polygon=clip_triangle_to_box(input,cell,clipping_roundoff,contact);const Measure measure=polygon_measure(triangle.polygon,triangle.normal);
			if(!(measure.value>area_tolerance))continue;
			triangle.area=measure.value;triangle.centroid=measure.centroid;
			triangle.boundary_face=coplanar_boundary_face(input,triangle.polygon,triangle.normal,
				cell,cartesian_face,angular_cartesian_face);
			prepared_by_input[input_index]=static_cast<int>(prepared.size());
			prepared.push_back(std::move(triangle));result.clipped_input_surface_area+=measure.value;
		}

		// Build geometric fabric components and certify every attached edge locally.
		// This prevents an FP-sized crack at a sewn seam from being interpreted as an
		// intentional opening.  Attachment provenance narrows matching when supplied,
		// but geometry must agree within the contact tolerance in every case.
		DisjointSet fabric_components(static_cast<int>(prepared.size()));
		struct TriangleEdgePolicy
		{
			std::array<bool,3> reaches_cell_interior{};
			std::array<bool,3> certified_attached{};
			std::array<int,3> provenance_peer_count{};
			std::array<double,3> worst_peer_distance{};
			bool unknown=false,uncertified_attached=false;
			int first_uncertified_edge=-1;
		};
		std::vector<TriangleEdgePolicy> triangle_edge_policy(prepared.size());
		const double contact2=contact*contact;
		for(int triangle_index=0;triangle_index<static_cast<int>(prepared.size());++triangle_index)
		{
			const Vec3d vertex[3]={prepared[triangle_index].input.a,prepared[triangle_index].input.b,
				prepared[triangle_index].input.c};
			for(int edge=0;edge<3;++edge)
			{
				Vec3d edge_a{},edge_b{};
				if(!clip_segment_to_box(vertex[edge],vertex[(edge+1)%3],cell,contact,edge_a,edge_b)||
					length2(edge_b-edge_a)<=contact2)continue;
				const Vec3d edge_midpoint=(edge_a+edge_b)*0.5;
				bool edge_reaches_cell_interior=true;
				for(int axis=0;axis<3;++axis)
					edge_reaches_cell_interior=edge_reaches_cell_interior&&
						edge_midpoint[axis]>cell.lo[axis]+contact&&edge_midpoint[axis]<cell.hi[axis]-contact;
				const FabricEdgeKind kind=prepared[triangle_index].input.edge_kind[edge];
				triangle_edge_policy[triangle_index].reaches_cell_interior[edge]=
					edge_reaches_cell_interior;
				if(kind==FabricEdgeKind::confirmed_free)
					continue;
				if(kind==FabricEdgeKind::unknown)
				{
					// An uncertain edge that enters the open cell interior can change fluid
					// connectivity and is therefore fail-closed.  An edge wholly on the
					// Cartesian boundary cannot reconnect atoms inside this cell; its neighbour
					// is checked later by exact shared-face assembly.
					triangle_edge_policy[triangle_index].unknown=
						triangle_edge_policy[triangle_index].unknown||edge_reaches_cell_interior;
					if(!edge_reaches_cell_interior)
					{
						for(int axis=0;axis<3;++axis)for(int upper_index=0;upper_index<2;++upper_index)
						{
							const bool upper=upper_index!=0;const double coordinate=upper?cell.hi[axis]:cell.lo[axis];
							if(std::abs(edge_a[axis]-coordinate)>contact||std::abs(edge_b[axis]-coordinate)>contact)continue;
							result.boundary_edge_hazards.push_back({
								prepared[triangle_index].input.source_triangle_id,
								prepared[triangle_index].input.source_face_id,
								static_cast<std::int8_t>(axis),upper,edge_a,edge_b,contact,ambiguity});
						}
					}
					continue;
				}

				const std::uint64_t attachment=prepared[triangle_index].input.edge_attachment_id[edge];
				bool has_provenance_peer=false;int provenance_peer_count=0;
				if(attachment!=no_attachment)
				{
					for(int candidate=0;candidate<static_cast<int>(canonical_triangles.size());++candidate)
						if(candidate!=prepared[triangle_index].input_index)
							for(std::uint64_t candidate_attachment:canonical_triangles[candidate].edge_attachment_id)
								if(candidate_attachment==attachment)
								{has_provenance_peer=true;++provenance_peer_count;}
				}
				bool certified=true;double worst_peer_distance=0.0;
				// Endpoint coordinates of independently tessellated incident faces need not
				// terminate at the same parameter on their common CAD edge (notably at a
				// 1:N seam or a many-face cusp).  The reciprocal attachment token already
				// identifies the candidate class; certify the positive-length interior of
				// the clipped edge and leave an ambiguity-sized endpoint neighbourhood to
				// the other incident edge sectors.  A missing or separated peer still fails
				// every interior sample.  This is dimensional and applies to any CAD seam.
				const double edge_length=std::sqrt(length2(edge_b-edge_a));
				const double endpoint_parameter=edge_length>0.0?
					std::min(0.25,16.0*ambiguity/edge_length):0.0;
				for(int sample_index=0;sample_index<3;++sample_index)
				{
					const double parameter=sample_index==0?endpoint_parameter:
						(sample_index==1?0.5:1.0-endpoint_parameter);
					const Vec3d sample=edge_a+(edge_b-edge_a)*parameter;
					bool sample_touched=false;double nearest_distance2=std::numeric_limits<double>::infinity();
					for(int candidate=0;candidate<static_cast<int>(canonical_triangles.size());++candidate)
					{
						if(candidate==prepared[triangle_index].input_index)continue;
						if(has_provenance_peer)
						{
							bool same_attachment=false;
							for(std::uint64_t candidate_attachment:canonical_triangles[candidate].edge_attachment_id)
								same_attachment=same_attachment||candidate_attachment==attachment;
							if(!same_attachment)continue;
						}
						const Vec3d candidate_normal=cross(canonical_triangles[candidate].b-canonical_triangles[candidate].a,
							canonical_triangles[candidate].c-canonical_triangles[candidate].a);
						if(length2(candidate_normal)<=area_tolerance*area_tolerance)continue;
						const Vec3d nearest=closest_point_on_triangle(sample,canonical_triangles[candidate]);
						nearest_distance2=std::min(nearest_distance2,length2(nearest-sample));
						if(length2(nearest-sample)<=contact2)
						{
							sample_touched=true;
							if(prepared_by_input[candidate]>=0)
								fabric_components.join(triangle_index,prepared_by_input[candidate]);
						}
					}
					if(std::isfinite(nearest_distance2))
						worst_peer_distance=std::max(worst_peer_distance,std::sqrt(nearest_distance2));
					else worst_peer_distance=std::numeric_limits<double>::infinity();
					certified=certified&&sample_touched;
				}
				triangle_edge_policy[triangle_index].provenance_peer_count[edge]=
					provenance_peer_count;
				triangle_edge_policy[triangle_index].worst_peer_distance[edge]=worst_peer_distance;
				triangle_edge_policy[triangle_index].uncertified_attached=
					triangle_edge_policy[triangle_index].uncertified_attached||
					!certified;
				if(!certified&&triangle_edge_policy[triangle_index].first_uncertified_edge<0)
					triangle_edge_policy[triangle_index].first_uncertified_edge=edge;
				triangle_edge_policy[triangle_index].certified_attached[edge]=certified;
			}
		}
		struct ComponentEdgePolicy
		{
			bool unknown=false,uncertified_attached=false;
			int uncertified_triangle=-1,uncertified_edge=-1;
		};
		std::map<int,ComponentEdgePolicy> component_edge_policy;
		for(int triangle_index=0;triangle_index<static_cast<int>(prepared.size());++triangle_index)
		{
			auto& policy=component_edge_policy[fabric_components.root(triangle_index)];
			policy.unknown=policy.unknown||triangle_edge_policy[triangle_index].unknown;
			if(triangle_edge_policy[triangle_index].uncertified_attached)
			{
				policy.uncertified_attached=true;
				if(policy.uncertified_triangle<0)
				{
					policy.uncertified_triangle=triangle_index;
					policy.uncertified_edge=
						triangle_edge_policy[triangle_index].first_uncertified_edge;
				}
			}
		}
		for(const auto& item:component_edge_policy)
		{
			if(item.second.unknown)
			{
				fail(result,"fabric component has an unknown edge in the cell interior; topology is unresolved");
				return result;
			}
			if(item.second.uncertified_attached)
			{
				const int triangle_index=item.second.uncertified_triangle;
				const int edge=item.second.uncertified_edge;
				const auto& triangle=prepared[triangle_index].input;
				fail(result,"attached edge in the cell interior lacks incident-geometry certification"
					" (source triangle="+std::to_string(triangle.source_triangle_id)+
					", local edge="+std::to_string(edge)+", attachment="+
					std::to_string(edge>=0?triangle.edge_attachment_id[edge]:no_attachment)+
					", local provenance peers="+std::to_string(edge>=0?
						triangle_edge_policy[triangle_index].provenance_peer_count[edge]:0)+
					", worst peer distance="+diagnostic_number(edge>=0?
						triangle_edge_policy[triangle_index].worst_peer_distance[edge]:
						std::numeric_limits<double>::infinity())+
					", cell=["+std::to_string(cell.lo.x)+","+std::to_string(cell.lo.y)+","+
					std::to_string(cell.lo.z)+"]..["+std::to_string(cell.hi.x)+","+
					std::to_string(cell.hi.y)+","+std::to_string(cell.hi.z)+"])");
				return result;
			}
		}

		std::string plane_error;
		std::vector<std::vector<int>> support_prepared_by_plane;
		const Vec3d plane_reference=(cell.lo+cell.hi)*0.5;
		auto add_plane=[&](Plane raw,int support_prepared,bool edge_limiter)->int
		{
			raw=canonical_plane(raw.normal,raw.offset);if(length2(raw.normal)==0.0)return -1;
			for(int index=0;index<static_cast<int>(result.planes.size());++index)
			{
				const double orientation=dot(result.planes[index].normal,raw.normal)>=0.0?1.0:-1.0;
				const Vec3d aligned_normal=raw.normal*orientation;
				const double aligned_offset=raw.offset*orientation;
				const double angle=std::sqrt(length2(result.planes[index].normal-aligned_normal));
				const double existing_local_offset=result.planes[index].offset-
					dot(result.planes[index].normal,plane_reference);
				const double raw_local_offset=aligned_offset-dot(aligned_normal,plane_reference);
				const double distance=std::abs(existing_local_offset-raw_local_offset);
				if(angle<=angular_contact&&distance<=contact)
				{
					if(support_prepared>=0)
					{
						result.planes[index].support_triangles.push_back(prepared[support_prepared].input_index);
						result.planes[index].support_source_triangles.push_back(
							prepared[support_prepared].input.source_triangle_id);
						support_prepared_by_plane[index].push_back(support_prepared);
					}
					result.planes[index].has_edge_limiter_role=result.planes[index].has_edge_limiter_role||edge_limiter;return index;
				}
				if(angle<=angular_ambiguity&&distance<=ambiguity)
				{
					plane_error="distinct local partition planes fall inside the geometric ambiguity band";return -1;
				}
			}
			if(static_cast<int>(result.planes.size())>=options.maximum_planes)
			{
				plane_error="local arrangement exceeds maximum_planes (limit="+
					std::to_string(options.maximum_planes)+", clipped triangles="+
					std::to_string(prepared.size())+")";
				return -1;
			}
			LocalArrangementPlane plane;plane.normal=raw.normal;plane.offset=raw.offset;plane.has_edge_limiter_role=edge_limiter;
			std::vector<int> support_prepared_indices;
			if(support_prepared>=0)
			{
				plane.support_triangles.push_back(prepared[support_prepared].input_index);
				plane.support_source_triangles.push_back(
					prepared[support_prepared].input.source_triangle_id);
				support_prepared_indices.push_back(support_prepared);
			}
			result.planes.push_back(std::move(plane));
			support_prepared_by_plane.push_back(std::move(support_prepared_indices));
			return static_cast<int>(result.planes.size())-1;
		};

		// Only physical support planes partition the 3-D cell.  A triangle edge bounds
		// fabric, not fluid volume: extending each edge to an infinite limiter plane creates
		// O(N^3) artificial atoms for an otherwise ordinary curved tessellation.  The finite
		// polygons are common-refined on each support facet below.  Neighbouring non-coplanar
		// support planes already provide the exact edge line at folds, seams, and junctions;
		// a genuinely free edge leaves positive-area open fluid on the same facet.
		for(int prepared_index=0;prepared_index<static_cast<int>(prepared.size())&&!plane_error.size();++prepared_index)
		{
			auto& triangle=prepared[prepared_index];
			if(triangle.boundary_face<0)
			{
				const Plane support{triangle.normal,dot(triangle.normal,triangle.input.a)};
				triangle.support_plane=add_plane(support,prepared_index,false);
				if(!plane_error.empty())break;
			}
		}
		if(!plane_error.empty()){fail(result,plane_error);return result;}
		for(int plane_index=0;plane_index<static_cast<int>(result.planes.size());++plane_index)
		{
			auto& public_support=result.planes[plane_index].support_triangles;
			std::sort(public_support.begin(),public_support.end());
			public_support.erase(std::unique(public_support.begin(),public_support.end()),public_support.end());
			auto& source_support=result.planes[plane_index].support_source_triangles;
			std::sort(source_support.begin(),source_support.end());
			source_support.erase(std::unique(source_support.begin(),source_support.end()),
				source_support.end());
			auto& prepared_support=support_prepared_by_plane[plane_index];
			std::sort(prepared_support.begin(),prepared_support.end());
			prepared_support.erase(std::unique(prepared_support.begin(),prepared_support.end()),prepared_support.end());
		}
		// A contact-certified set of nominally coplanar CAD triangles shares one
		// canonical support plane.  Project its finite polygons onto that same plane
		// before doing the 2-D common refinement.  Otherwise the support facet and a
		// triangle can be parallel but separated by a few ULPs: 3-D half-space clipping
		// then reports an empty overlap even when the facet centroid is inside the
		// triangle, opening a false passage at a sewn cusp.
		for(PreparedTriangle& triangle:prepared)
			if(triangle.support_plane>=0)
			{
				const LocalArrangementPlane& support=result.planes[triangle.support_plane];
				for(Vec3d& point:triangle.polygon)
					point=point-support.normal*(dot(support.normal,point)-support.offset);
				const Measure projected=polygon_measure(triangle.polygon,support.normal);
				triangle.area=projected.value;triangle.centroid=projected.centroid;
			}

		std::vector<Polyhedron> polyhedra{cube_polyhedron(cell,static_cast<int>(result.planes.size()))};
		std::uint64_t next_cap_lineage=1;
		for(int plane_index=0;plane_index<static_cast<int>(result.planes.size());++plane_index)
		{
			const Plane plane{result.planes[plane_index].normal,result.planes[plane_index].offset};std::vector<Polyhedron> next;
			for(Polyhedron& polyhedron:polyhedra)
			{
				const std::vector<Vec3d> vertices=polyhedron_vertices(polyhedron,clipping_roundoff);double minimum=std::numeric_limits<double>::infinity(),maximum=-minimum;
				for(Vec3d point:vertices){const double value=plane.distance(point);minimum=std::min(minimum,value);maximum=std::max(maximum,value);}
				if(minimum < -clipping_roundoff&&maximum > clipping_roundoff)
				{
					const Polygon shared_cap=polyhedron_plane_cap(polyhedron,plane,
						clipping_roundoff,clipping_roundoff);
					const std::uint64_t cap_lineage=next_cap_lineage++;
					Polyhedron minus=clip_polyhedron(polyhedron,plane,-1,plane_index,cap_lineage,
						shared_cap,clipping_roundoff,clipping_roundoff,area_tolerance);
					Polyhedron plus=clip_polyhedron(polyhedron,plane,1,plane_index,cap_lineage,
						shared_cap,clipping_roundoff,clipping_roundoff,area_tolerance);
					const Measure minus_measure=polyhedron_measure(minus),plus_measure=polyhedron_measure(plus);
					if(minus_measure.value>arithmetic_volume_tolerance&&
						plus_measure.value>arithmetic_volume_tolerance)
					{
						next.push_back(std::move(minus));
						next.push_back(std::move(plus));
					}
					else
					{
						// A grazing split below the representable volume is not a pressure
						// control volume. Retain the dominant atom for limiter and support
						// planes alike. This cannot silently erase fabric: a support patch needs
						// an opposite atom pair, and the strict represented-vs-input surface-area
						// invariant below rejects the arrangement if any positive-area patch was
						// actually lost.
						// The retained polyhedron still geometrically straddles this plane.  A
						// dominant-side label would be false topology: later cuts may subdivide
						// its neighbours differently, and point classification in the discarded
						// sliver must still resolve to this same control volume.  Zero is an
						// explicit wildcard for a deliberately unrepresented grazing cut.
						polyhedron.side[plane_index] = 0;
						next.push_back(std::move(polyhedron));
					}
				}
				else
				{
					const double representative=0.5*(minimum+maximum);polyhedron.side[plane_index]=representative<0?-1:1;next.push_back(std::move(polyhedron));
				}
				if(static_cast<int>(next.size())>options.maximum_atoms)
				{
					fail(result,"local arrangement exceeds maximum_atoms (limit="+
						std::to_string(options.maximum_atoms)+", plane="+
						std::to_string(plane_index+1)+"/"+
						std::to_string(result.planes.size())+")");
					return result;
				}
			}
			polyhedra=std::move(next);
		}

		result.atoms.reserve(polyhedra.size());
		for(Polyhedron& polyhedron:polyhedra)
		{
			const Measure measure=polyhedron_measure(polyhedron);
			if(!(measure.value>arithmetic_volume_tolerance))
			{
				fail(result,"local arrangement retained a numerically degenerate atom");return result;
			}
			LocalArrangementAtom atom;atom.volume=measure.value;atom.centroid=measure.centroid;atom.plane_side=std::move(polyhedron.side);atom.faces=std::move(polyhedron.faces);result.atoms.push_back(std::move(atom));
		}
		DisjointSet components(static_cast<int>(result.atoms.size()));std::vector<AtomicPatch> atomic_patches;
		struct OpenAdjacency
		{
			int a=-1,b=-1,plane=-1;double covered=0.0,area=0.0;Vec3d centroid{};
			bool bounded_by_confirmed_free_edge=false;
		};
		std::vector<OpenAdjacency> open_adjacencies;
		struct LocalEdgeEvidence{bool certified_attached=false,confirmed_free=false;};
		auto local_edge_evidence=[&](int triangle_index,const Polygon& polygon,
			int support_plane,Vec3d plane_normal)
		{
			LocalEdgeEvidence evidence;
			const PreparedTriangle& triangle=prepared[triangle_index];
			const LocalArrangementPlane& support=result.planes[support_plane];
			// A BSP intersection at angle theta amplifies coordinate roundoff by
			// O(1/sin(theta)).  The sqrt(epsilon)*h transition is the standard
			// scale at which that conditioned displacement and the retained facet
			// area have equal relative error.  This is used only to identify which
			// provenance edge bounds a sub-resolution facet; it never turns an
			// unlabelled or confirmed-free edge into fabric.
			const double line_tolerance=std::max(contact,
				64.0*std::sqrt(machine)*h);
			const Measure polygon_geometry=polygon_measure(polygon,plane_normal);
			const double provenance_area_bound=std::max({area_tolerance,64.0*contact*h,
				64.0*std::sqrt(machine)*minimum_face_area});
			for(int edge=0;edge<3;++edge)
			{
				if(!triangle_edge_policy[triangle_index].reaches_cell_interior[edge])continue;
				Vec3d edge_a=triangle_vertex(triangle.input,edge);
				Vec3d edge_b=triangle_vertex(triangle.input,(edge+1)%3);
				edge_a=edge_a-plane_normal*(dot(plane_normal,edge_a)-support.offset);
				edge_b=edge_b-plane_normal*(dot(plane_normal,edge_b)-support.offset);
				const Vec3d edge_delta=edge_b-edge_a;
				const double edge_length=std::sqrt(length2(edge_delta));
				if(!(edge_length>line_tolerance))continue;
				bool incident=false;
				for(std::size_t index=0;index<polygon.size()&&!incident;++index)
				{
					const Vec3d a=polygon[index];
					const Vec3d b=polygon[(index+1)%polygon.size()];
					const double distance_a=std::abs(dot(cross(edge_delta,a-edge_a),plane_normal))/edge_length;
					const double distance_b=std::abs(dot(cross(edge_delta,b-edge_a),plane_normal))/edge_length;
					if(distance_a>line_tolerance||distance_b>line_tolerance)continue;
					const double parameter_a=dot(a-edge_a,edge_delta)/length2(edge_delta);
					const double parameter_b=dot(b-edge_a,edge_delta)/length2(edge_delta);
					const double overlap_begin=std::max(0.0,std::min(parameter_a,parameter_b));
					const double overlap_end=std::min(1.0,std::max(parameter_a,parameter_b));
					incident=(overlap_end-overlap_begin)*edge_length>
						8.0*clipping_roundoff;
				}
				if(!incident&&polygon_geometry.value<=provenance_area_bound)
				{
					const double parameter=std::clamp(dot(polygon_geometry.centroid-edge_a,
						edge_delta)/length2(edge_delta),0.0,1.0);
					incident=std::sqrt(length2(polygon_geometry.centroid-
						(edge_a+edge_delta*parameter)))<=line_tolerance;
				}
				if(!incident)continue;
				const FabricEdgeKind kind=triangle.input.edge_kind[edge];
				evidence.confirmed_free=evidence.confirmed_free||
					kind==FabricEdgeKind::confirmed_free;
				evidence.certified_attached=evidence.certified_attached||
					(kind==FabricEdgeKind::attached&&
					triangle_edge_policy[triangle_index].certified_attached[edge]);
			}
			return evidence;
		};

		for(int plane_index=0;plane_index<static_cast<int>(result.planes.size());++plane_index)
		{
			// Intersect the actual convex facets on the two sides instead of pairing
			// atoms by their complete BSP sign vectors.  Sign-vector equality assumes
			// every other plane split both neighbours identically; that is not true when
			// a grazing, sub-volume-tolerance sliver is deliberately retained unsplit.
			// The polygon overlaps are the common refinement and therefore remain valid
			// for mismatched tessellations and any order of the partition planes.
			struct Facet
			{
				int atom=-1;
				const LocalArrangementFace* face=nullptr;
				double area=0.0;
				double matched_area=0.0;
				std::uint64_t lineage=0;
			};
			std::vector<Facet> minus_facets,plus_facets;
			const Vec3d plane_normal=result.planes[plane_index].normal;
			for(int atom=0;atom<static_cast<int>(result.atoms.size());++atom)
			{
				const LocalArrangementFace* face=face_on_plane(result.atoms[atom],plane_index);
				if(!face)continue;
				const Measure measure=polygon_measure(face->polygon,plane_normal);
				if(!(measure.value>area_tolerance))continue;
				const std::int8_t side=result.atoms[atom].plane_side[plane_index];
				if(side<0)minus_facets.push_back({atom,face,measure.value,0.0,face->lineage});
				else if(side>0)plus_facets.push_back({atom,face,measure.value,0.0,face->lineage});
				else
				{
					fail(result,"grazing-unsplit atom unexpectedly owns a facet on the skipped plane");
					return result;
				}
			}

			for(Facet& minus:minus_facets)for(Facet& plus:plus_facets)
			{
				if(minus.lineage!=plus.lineage)continue;
				// Later planes subdivide both copies of this original cap.  Only descendants
				// occupying the same side of every other represented plane can be opposite
				// sides of one facet.  A zero is the explicit wildcard for an arithmetic
				// grazing split retained unsplit; treating it as a hard mismatch would crack
				// the complex, while omitting this check connects incompatible sectors at a
				// shallow attached cusp.
				bool compatible=true;
				for(int other=0;other<static_cast<int>(result.planes.size());++other)
					if(other!=plane_index)
					{
						const std::int8_t minus_side=result.atoms[minus.atom].plane_side[other];
						const std::int8_t plus_side=result.atoms[plus.atom].plane_side[other];
						if(minus_side!=0&&plus_side!=0&&minus_side!=plus_side)
						{compatible=false;break;}
					}
				if(!compatible)continue;
				Polygon overlap=intersect_coplanar_convex(minus.face->polygon,plus.face->polygon,
					plane_normal,contact);
				const Measure overlap_measure=polygon_measure(overlap,plane_normal);
				if(!(overlap_measure.value>area_tolerance))continue;
				minus.matched_area+=overlap_measure.value;
				plus.matched_area+=overlap_measure.value;

				const double coverage_tolerance=std::max(area_tolerance,64.0*contact*h);
				std::vector<AtomicPatch> pieces;double covered=0.0;
				for(int triangle_index:support_prepared_by_plane[plane_index])
				{
					const PreparedTriangle& triangle=prepared[triangle_index];
					Polygon polygon=intersect_coplanar_convex(
						overlap,triangle.polygon,plane_normal,contact);
					const Measure measure=polygon_measure(polygon,plane_normal);
					// This is an exact-predicate convex intersection on one canonical support
					// plane.  Do not compare a physical fabric piece with sqrt(epsilon) times
					// the parent triangle/facet area: at a closed sharp cusp, later support
					// planes legitimately subdivide a fully covered facet below that scale.
					// Dropping the piece then turns a tessellation seam into an open route and
					// joins the two pressure sides.  `area_tolerance` is the declared minimum
					// representable measure and is already scale-aware; positive pieces above
					// it must participate regardless of their fraction of a larger triangle.
					if(measure.value>area_tolerance)
					{
						pieces.push_back({triangle_index,minus.atom,plus.atom,std::move(polygon),
							measure.value,measure.centroid,overlap_measure.value,0.0,plane_index});
						covered+=measure.value;
					}
				}
				// At a highly conditioned cusp, later plane cuts can leave a support
				// subfacet far below the runtime aperture scale.  Its binary64 polygon
				// intersection may collapse even though the subfacet centroid is certified
				// inside one finite fabric triangle.  Resolve only such sub-resolution,
				// unambiguously owned facets as fabric; a centroid in open fluid remains an
				// aperture, and the global surface-area invariant below catches double
				// ownership.  This is a dimensional arithmetic threshold, not a profile or
				// source-face exception.
				const double subfacet_tolerance=std::max(coverage_tolerance,
					64.0*std::sqrt(machine)*minimum_face_area);
				if(overlap_measure.value<=subfacet_tolerance)
				{
					int attached_owner=-1;bool touches_confirmed_free=false;
					for(int triangle_index:support_prepared_by_plane[plane_index])
					{
						const PreparedTriangle& triangle=prepared[triangle_index];
						const LocalEdgeEvidence evidence=local_edge_evidence(triangle_index,
							overlap,plane_index,plane_normal);
						touches_confirmed_free=touches_confirmed_free||evidence.confirmed_free;
						if(attached_owner<0&&evidence.certified_attached&&
							point_on_convex_polygon(overlap_measure.centroid,triangle.polygon,
								plane_normal,contact))attached_owner=triangle_index;
					}
					if(!touches_confirmed_free&&attached_owner>=0)
					{
						pieces.clear();
						pieces.push_back({attached_owner,minus.atom,plus.atom,overlap,
							overlap_measure.value,overlap_measure.centroid,overlap_measure.value,
							overlap_measure.value,plane_index});
						covered=overlap_measure.value;
					}
				}
				bool locally_certified_attached=false,locally_confirmed_free=false;
				for(const AtomicPatch& piece:pieces)
				{
					const LocalEdgeEvidence evidence=local_edge_evidence(piece.triangle,
						piece.polygon,piece.support_plane,plane_normal);
					locally_certified_attached=locally_certified_attached||evidence.certified_attached;
					locally_confirmed_free=locally_confirmed_free||evidence.confirmed_free;
				}
				for(int triangle_index:support_prepared_by_plane[plane_index])
				{
					const LocalEdgeEvidence evidence=local_edge_evidence(triangle_index,overlap,
						plane_index,plane_normal);
					if(evidence.certified_attached)
					{
						locally_certified_attached=true;
					}
					if(evidence.confirmed_free)
						locally_confirmed_free=true;
				}
				const double uncovered=overlap_measure.value-covered;
				if(!pieces.empty()&&!locally_confirmed_free&&locally_certified_attached&&uncovered>0.0&&
					uncovered<=subfacet_tolerance)
				{
					covered=overlap_measure.value;
				}
				for(std::size_t a=0;a<pieces.size();++a)for(std::size_t b=a+1;b<pieces.size();++b)
				{
					const Polygon fabric_overlap=intersect_coplanar_convex(pieces[a].polygon,
						pieces[b].polygon,plane_normal,contact);
					if(polygon_measure(fabric_overlap,plane_normal).value>coverage_tolerance)
					{
						fail(result,"positive-area overlapping fabric triangles share a support facet");
						return result;
					}
				}
				const double coverage_fraction_tolerance=64.0*std::sqrt(machine);
				if(pieces.empty())
				{
					open_adjacencies.push_back({minus.atom,plus.atom,plane_index,covered,
						overlap_measure.value,overlap_measure.centroid,false});
					components.join(minus.atom,plus.atom);
				}
				else
				{
					if(covered>(1.0+coverage_fraction_tolerance)*overlap_measure.value)
					{
						fail(result,"finite fabric polygons over-cover a support facet");
						return result;
					}
					// A positive open remainder connects the two convex atoms around a finite
					// edge.  Keep every covered piece as a wall patch as well, however small
					// its fraction of this support facet.  An area-fraction shortcut here used
					// to erase genuine arithmetic-scale fabric while still joining its sides.
					// The component
					// audit below accepts this only when CAD/BVH provenance confirmed a real
					// free edge; an attached seam cannot turn a numerical coverage crack into
					// fluid connectivity.
					const bool attached_seam_carrier=locally_certified_attached&&
						!locally_confirmed_free&&
						covered<=coverage_fraction_tolerance*overlap_measure.value;
					for(AtomicPatch& piece:pieces)
						piece.attached_seam_carrier=attached_seam_carrier;
					const double open_remainder_tolerance=locally_confirmed_free?
						area_tolerance:coverage_fraction_tolerance*overlap_measure.value;
					if(overlap_measure.value-covered>open_remainder_tolerance)
					{
						open_adjacencies.push_back({minus.atom,plus.atom,plane_index,covered,
							overlap_measure.value,overlap_measure.centroid,
							locally_confirmed_free});
						components.join(minus.atom,plus.atom);
					}
					for(AtomicPatch& piece:pieces)piece.parent_fabric_area=covered;
					atomic_patches.insert(atomic_patches.end(),std::make_move_iterator(pieces.begin()),
						std::make_move_iterator(pieces.end()));
				}
			}

			// Every retained split facet must participate in the common refinement.
			// This catches a true crack in the atom complex rather than silently treating
			// an unmatched numerical partition as fabric or as an opening.
			const double match_tolerance=std::max(area_tolerance,128.0*contact*h);
			for(const Facet& facet:minus_facets)
				if(std::abs(facet.matched_area-facet.area)>match_tolerance)
				{
					double plus_area=0.0,lineage_area=0.0,lineage_hausdorff=0.0;
					std::size_t lineage_count=0,lineage_vertices=0;
					for(const Facet& plus:plus_facets)
					{
						plus_area+=plus.area;
						if(plus.lineage==facet.lineage)
						{
							lineage_area+=plus.area;++lineage_count;
							lineage_vertices=plus.face->polygon.size();
							for(Vec3d a:facet.face->polygon)
							{
								double nearest=std::numeric_limits<double>::infinity();
								for(Vec3d b:plus.face->polygon)nearest=std::min(nearest,length2(a-b));
								lineage_hausdorff=std::max(lineage_hausdorff,std::sqrt(nearest));
							}
						}
					}
					fail(result,"minus-side atom facet lacks a conservative opposite-side match (plane="+
						std::to_string(plane_index)+", matched="+
						std::to_string(facet.matched_area)+", area="+
						std::to_string(facet.area)+", plus facets/area="+
						std::to_string(plus_facets.size())+"/"+std::to_string(plus_area)+
						", lineage facets/area="+std::to_string(lineage_count)+"/"+
						std::to_string(lineage_area)+", vertices="+
						std::to_string(facet.face->polygon.size())+"/"+
						std::to_string(lineage_vertices)+", hausdorff="+
						diagnostic_number(lineage_hausdorff)+
						", cell=["+std::to_string(cell.lo.x)+","+std::to_string(cell.lo.y)+","+
						std::to_string(cell.lo.z)+"]..["+std::to_string(cell.hi.x)+","+
						std::to_string(cell.hi.y)+","+std::to_string(cell.hi.z)+"])");
					return result;
				}
			for(const Facet& facet:plus_facets)
				if(std::abs(facet.matched_area-facet.area)>match_tolerance)
				{
					fail(result,"plus-side atom facet lacks a conservative opposite-side match");
					return result;
				}
		}

		std::map<int,int> fragment_for_root;
		for(int atom=0;atom<static_cast<int>(result.atoms.size());++atom)
		{
			const int root=components.root(atom);auto inserted=fragment_for_root.emplace(root,static_cast<int>(fragment_for_root.size()));const int fragment=inserted.first->second;result.atoms[atom].fragment=fragment;
			if(inserted.second) result.fragments.push_back({});
			auto& output = result.fragments[fragment];
			output.volume += result.atoms[atom].volume;
			output.centroid = output.centroid +
				(result.atoms[atom].centroid-plane_reference)*result.atoms[atom].volume;
			output.atoms.push_back(atom);
		}
		for(auto& fragment:result.fragments)if(fragment.volume>0)
			fragment.centroid=plane_reference+fragment.centroid/fragment.volume;

		for(AtomicPatch& atomic:atomic_patches)
		{
			int minus_fragment=result.atoms[atomic.canonical_minus_atom].fragment;
			int plus_fragment=result.atoms[atomic.canonical_plus_atom].fragment;
			if(minus_fragment==plus_fragment&&atomic.attached_seam_carrier)
			{
				// At a conditioned attached cusp, projecting the two incident CAD faces
				// onto their canonical support planes can leave an arithmetic-scale piece
				// of one source triangle on the adjacent sector.  Its open remainder joins
				// that local atom pair, but deleting the finite piece loses real pressure
				// area.  A resolved portion of the *same input triangle* supplies the only
				// winding-consistent pressure-side certificate: retain this polygon and map
				// it to that portion's two fluid components.  No other CAD face, profile, or
				// proximity inference may authorise the remap.
				AtomicPatch* donor=nullptr;
				for(AtomicPatch& candidate:atomic_patches)
					if(candidate.triangle==atomic.triangle&&
						result.atoms[candidate.canonical_minus_atom].fragment!=
						result.atoms[candidate.canonical_plus_atom].fragment&&
						(!donor||candidate.area>donor->area))donor=&candidate;
				if(donor)
				{
					atomic.canonical_minus_atom=donor->canonical_minus_atom;
					atomic.canonical_plus_atom=donor->canonical_plus_atom;
					minus_fragment=result.atoms[atomic.canonical_minus_atom].fragment;
					plus_fragment=result.atoms[atomic.canonical_plus_atom].fragment;
				}
			}
			if(minus_fragment!=plus_fragment)continue;
			// A same-fragment patch is valid only when the actual atom-to-atom path
			// crosses an open remainder bounded by a confirmed-free edge on that very
			// facet.  Merely having a free edge somewhere on a joined CAD component is
			// insufficient (notably at ribs and T-junctions).
			const int state_count=2*static_cast<int>(result.atoms.size());
			std::vector<int> predecessor(state_count,-1),predecessor_edge(state_count,-1);
			std::vector<int> frontier{2*atomic.canonical_minus_atom};
			predecessor[frontier.front()]=frontier.front();
			for(std::size_t cursor=0;cursor<frontier.size()&&
				predecessor[2*atomic.canonical_plus_atom+1]<0;++cursor)
			{
				const int atom=frontier[cursor]/2;
				const bool already_crossed_free=(frontier[cursor]&1)!=0;
				for(int edge=0;edge<static_cast<int>(open_adjacencies.size());++edge)
				{
					const OpenAdjacency& adjacency=open_adjacencies[edge];
					int next_atom=-1;
					if(adjacency.a==atom)next_atom=adjacency.b;
					else if(adjacency.b==atom)next_atom=adjacency.a;
					if(next_atom<0)continue;
					const bool crosses_component_free=
						adjacency.bounded_by_confirmed_free_edge;
					const int next_state=2*next_atom+
						((already_crossed_free||crosses_component_free)?1:0);
					if(predecessor[next_state]<0)
					{
						predecessor[next_state]=frontier[cursor];
						predecessor_edge[next_state]=edge;
						frontier.push_back(next_state);
					}
				}
			}
			if(predecessor[2*atomic.canonical_plus_atom+1]<0)
			{
				std::string path;
				int state=2*atomic.canonical_plus_atom;
				if(predecessor[state]<0)state=2*atomic.canonical_plus_atom+1;
				for(;state!=2*atomic.canonical_minus_atom&&state>=0;state=predecessor[state])
				{
					const int edge=predecessor_edge[state];if(edge<0)break;
					const OpenAdjacency& adjacency=open_adjacencies[edge];
					int strict_inside=0,tolerant_inside=0,strict_vertices=0,tolerant_vertices=0;
					for(int triangle:support_prepared_by_plane[adjacency.plane])
					{
						strict_inside+=point_on_convex_polygon(adjacency.centroid,
							prepared[triangle].polygon,result.planes[adjacency.plane].normal,
							clipping_roundoff)?1:0;
						tolerant_inside+=point_on_convex_polygon(adjacency.centroid,
							prepared[triangle].polygon,result.planes[adjacency.plane].normal,
							contact)?1:0;
						const LocalArrangementFace* adjacency_face=face_on_plane(
							result.atoms[adjacency.a],adjacency.plane);
						if(adjacency_face)for(Vec3d vertex:adjacency_face->polygon)
						{
							strict_vertices+=point_on_convex_polygon(vertex,
								prepared[triangle].polygon,result.planes[adjacency.plane].normal,
								clipping_roundoff)?1:0;
							tolerant_vertices+=point_on_convex_polygon(vertex,
								prepared[triangle].polygon,result.planes[adjacency.plane].normal,
								contact)?1:0;
						}
					}
					path=" p"+std::to_string(adjacency.plane)+":"+
						diagnostic_number(adjacency.covered)+"/"+
						diagnostic_number(adjacency.area)+"i"+
						std::to_string(strict_inside)+"/"+std::to_string(tolerant_inside)+"v"+
						std::to_string(strict_vertices)+"/"+std::to_string(tolerant_vertices)+"f"+
						(adjacency.bounded_by_confirmed_free_edge?"1":"0")+path;
				}
				fail(result,"fluid reconnects around fabric without a confirmed free interior edge (source triangle="+
					std::to_string(prepared[atomic.triangle].input.source_triangle_id)+", plane="+
					std::to_string(atomic.support_plane)+", fabric/facet="+
					std::to_string(atomic.parent_fabric_area)+"/"+
					std::to_string(atomic.parent_facet_area)+", cell=["+
					std::to_string(cell.lo.x)+","+std::to_string(cell.lo.y)+","+
					std::to_string(cell.lo.z)+"]..["+std::to_string(cell.hi.x)+","+
					std::to_string(cell.hi.y)+","+std::to_string(cell.hi.z)+"], supports="+
					std::to_string(result.planes.size())+", path="+path+")");
				return result;
			}
		}

		for(AtomicPatch& atomic:atomic_patches)
		{
			const PreparedTriangle& triangle=prepared[atomic.triangle];const LocalArrangementPlane& plane=result.planes[triangle.support_plane];const bool aligned=dot(triangle.normal,plane.normal)>=0;
			LocalArrangementSurfacePatch patch;patch.source_triangle_id=triangle.input.source_triangle_id;patch.source_face_id=triangle.input.source_face_id;patch.area=atomic.area;patch.centroid=atomic.centroid;patch.normal=triangle.normal;patch.polygon=std::move(atomic.polygon);
			patch.plus_fragment=result.atoms[aligned?atomic.canonical_plus_atom:atomic.canonical_minus_atom].fragment;patch.minus_fragment=result.atoms[aligned?atomic.canonical_minus_atom:atomic.canonical_plus_atom].fragment;
			result.represented_surface_area+=patch.area;result.surface_patches.push_back(std::move(patch));
		}

		for(int atom=0;atom<static_cast<int>(result.atoms.size());++atom)for(const auto& face:result.atoms[atom].faces)
		{
			int axis=0;bool upper=false;if(!decode_boundary_plane(face.plane_id,axis,upper))continue;Vec3d normal{};normal[axis]=upper?1.0:-1.0;const Measure face_measure=polygon_measure(face.polygon,normal);if(!(face_measure.value>area_tolerance))continue;
			std::vector<std::pair<int,Polygon>> pieces;double covered=0.0;
			for(int triangle_index=0;triangle_index<static_cast<int>(prepared.size());++triangle_index)
			{
				const PreparedTriangle& triangle=prepared[triangle_index];if(triangle.boundary_face!=2*axis+(upper?1:0))continue;Polygon polygon=intersect_coplanar_convex(face.polygon,triangle.polygon,normal,contact);const double area=polygon_measure(polygon,normal).value;if(area>area_tolerance){covered+=area;pieces.push_back({triangle_index,std::move(polygon)});}
			}
			const double boundary_coverage_tolerance=std::max(area_tolerance,
				64.0*clipping_roundoff*h);
			for(std::size_t a=0;a<pieces.size();++a)for(std::size_t b=a+1;b<pieces.size();++b)
			{
				const Polygon overlap=intersect_coplanar_convex(pieces[a].second,pieces[b].second,normal,contact);
				if(polygon_measure(overlap,normal).value>boundary_coverage_tolerance){fail(result,"positive-area overlapping fabric triangles lie on a Cartesian boundary face");return result;}
			}
			// A coordinate-error area bound must not become an aperture-size cutoff.
			// When a source edge is explicitly certified free and reaches the relative
			// interior of this Cartesian face, its two-dimensional subdivision is a
			// topological constraint.  Retain carriers down to the squared positional
			// uncertainty.  Closed/attached seams keep the wider arithmetic coverage
			// bound and are additionally audited by the shared-face side certificate.
			bool has_confirmed_free_interior_edge=false;
			for(const auto& piece:pieces)
			{
				const PreparedTriangle& triangle=prepared[piece.first];
				for(int edge=0;edge<3&&!has_confirmed_free_interior_edge;++edge)
				{
					if(triangle.input.edge_kind[edge]!=FabricEdgeKind::confirmed_free)continue;
					Vec3d edge_a{},edge_b{};
					if(!clip_segment_to_box(triangle_vertex(triangle.input,edge),
						triangle_vertex(triangle.input,(edge+1)%3),cell,contact,edge_a,edge_b)||
						length2(edge_b-edge_a)<=contact*contact)continue;
					const double coordinate=upper?cell.hi[axis]:cell.lo[axis];
					if(std::abs(edge_a[axis]-coordinate)>contact||
						std::abs(edge_b[axis]-coordinate)>contact)continue;
					const Vec3d midpoint=(edge_a+edge_b)*0.5;
					bool relative_interior=true;
					for(int tangent=0;tangent<3;++tangent)if(tangent!=axis)
						relative_interior=relative_interior&&
							midpoint[tangent]>cell.lo[tangent]+contact&&
							midpoint[tangent]<cell.hi[tangent]-contact;
					has_confirmed_free_interior_edge=relative_interior;
				}
			}
			// `contact` is a coordinate/provenance uncertainty, not an aperture-size
			// threshold.  In particular, using O(contact*h) here silently sealed narrow
			// physical vents.  Coverage decisions use only the arithmetic construction
			// tolerance; any positive resolved complement is represented explicitly.
			const double partition_area_tolerance=has_confirmed_free_interior_edge?
				std::max(contact*contact,64.0*clipping_roundoff*clipping_roundoff):
				boundary_coverage_tolerance;
			if(covered<=partition_area_tolerance)
			{
				result.boundary_apertures.push_back({static_cast<std::int8_t>(axis),upper,
					face_measure.value,face_measure.centroid,face.polygon,
					result.atoms[atom].fragment,contact});
			}
			else if(std::abs(covered-face_measure.value)<=partition_area_tolerance)
			{
				for(auto& piece:pieces)
				{
					const PreparedTriangle& triangle=prepared[piece.first];const Measure measure=polygon_measure(piece.second,normal);const Vec3d inward=normal*-1.0;
					LocalArrangementBoundarySurfacePatch patch;patch.source_triangle_id=triangle.input.source_triangle_id;patch.source_face_id=triangle.input.source_face_id;patch.axis=static_cast<std::int8_t>(axis);patch.upper=upper;patch.area=measure.value;patch.centroid=measure.centroid;patch.normal=triangle.normal;patch.polygon=std::move(piece.second);patch.interior_fragment=result.atoms[atom].fragment;patch.interior_is_plus=dot(triangle.normal,inward)>0;patch.contact_tolerance=contact;
					result.represented_surface_area+=patch.area;result.boundary_surface_patches.push_back(std::move(patch));
				}
			}
			else
			{
				// A finite sheet may cover only part of a Cartesian face.  Support planes
				// partition 3-D fluid atoms, but a triangle edge is not an infinite fluid
				// wall and therefore deliberately does not enter that BSP.  Common-refine
				// this boundary facet by the finite fabric-piece edges in 2-D, then label
				// each convex tile as fabric or open fluid.  This retains separate vents
				// (including sub-grid-area vents) without blocking the entire MAC face.
				std::vector<Polygon> partitions{face.polygon};
				for(const auto& piece:pieces)
					for(std::size_t edge=0;edge<piece.second.size();++edge)
					{
						std::vector<Polygon> next;next.reserve(2*partitions.size());
						const Vec3d line_a=piece.second[edge];
						const Vec3d line_b=piece.second[(edge+1)%piece.second.size()];
						for(Polygon& partition:partitions)
						{
							auto [negative,positive]=split_coplanar_convex_by_line(
								partition,line_a,line_b,normal,clipping_roundoff);
							const double negative_area=polygon_measure(negative,normal).value;
							const double positive_area=polygon_measure(positive,normal).value;
							if(negative_area>partition_area_tolerance&&
								positive_area>partition_area_tolerance)
							{
								next.push_back(std::move(negative));next.push_back(std::move(positive));
							}
							else next.push_back(std::move(partition));
						}
						partitions=std::move(next);
						if(static_cast<int>(partitions.size())>options.maximum_atoms)
						{
							fail(result,"Cartesian boundary common refinement exceeds maximum_atoms");
							return result;
						}
					}

				for(Polygon& partition:partitions)
				{
					const Measure partition_measure=polygon_measure(partition,normal);
					if(!(partition_measure.value>partition_area_tolerance))continue;
					std::vector<std::pair<int,Polygon>> fabric;double partition_covered=0.0;
					for(const auto& piece:pieces)
					{
						Polygon overlap=intersect_coplanar_convex(partition,piece.second,
							normal,clipping_roundoff);
						const double overlap_area=polygon_measure(overlap,normal).value;
						if(overlap_area>partition_area_tolerance)
						{
							partition_covered+=overlap_area;
							fabric.push_back({piece.first,std::move(overlap)});
						}
					}
					bool partition_is_open=partition_covered<=partition_area_tolerance;
					bool partition_is_fabric=
						std::abs(partition_covered-partition_measure.value)<=partition_area_tolerance;
					if(!partition_is_open&&!partition_is_fabric&&
						std::abs(partition_covered-partition_measure.value)<=boundary_coverage_tolerance)
					{
						// Intersections on a common fabric edge can leave an arithmetic-width
						// overlap on an otherwise atomic tile.  Its exact label is binary; select
						// the nearer of zero and full coverage.  Crucially, a tiny aperture bounded
						// by confirmed-free edges remains an open tile even when its absolute area
						// is below the global subtraction-error bound.
						partition_is_fabric=partition_covered>0.5*partition_measure.value;
						partition_is_open=!partition_is_fabric;
					}
					if(partition_is_open)
					{
						result.boundary_apertures.push_back({static_cast<std::int8_t>(axis),upper,
							partition_measure.value,partition_measure.centroid,std::move(partition),
							result.atoms[atom].fragment,contact});
						continue;
					}
					if(!partition_is_fabric||fabric.empty())
					{
						fail(result,"Cartesian boundary common refinement retained mixed fabric/open coverage"
							" (partition="+diagnostic_number(partition_measure.value)+
							", fabric="+diagnostic_number(partition_covered)+
							", tolerance="+diagnostic_number(partition_area_tolerance)+")");
						return result;
					}
					// The partition is one atomic fabric tile.  At an arithmetic-width
					// tessellation overlap several source triangles can claim it; retain it
					// exactly once under the largest (then lowest-ID) source carrier.
					auto owner=std::max_element(fabric.begin(),fabric.end(),[&](const auto& a,const auto& b)
					{
						const double area_a=polygon_measure(a.second,normal).value;
						const double area_b=polygon_measure(b.second,normal).value;
						if(area_a!=area_b)return area_a<area_b;
						return prepared[a.first].input.source_triangle_id>
							prepared[b.first].input.source_triangle_id;
					});
					const PreparedTriangle& triangle=prepared[owner->first];const Vec3d inward=normal*-1.0;
					LocalArrangementBoundarySurfacePatch patch;
					patch.source_triangle_id=triangle.input.source_triangle_id;
					patch.source_face_id=triangle.input.source_face_id;
					patch.axis=static_cast<std::int8_t>(axis);patch.upper=upper;
					patch.area=partition_measure.value;patch.centroid=partition_measure.centroid;
					patch.normal=triangle.normal;patch.polygon=std::move(partition);
					patch.interior_fragment=result.atoms[atom].fragment;
					patch.interior_is_plus=dot(triangle.normal,inward)>0;
					patch.contact_tolerance=contact;
					result.represented_surface_area+=patch.area;
					result.boundary_surface_patches.push_back(std::move(patch));
				}
			}
		}
		for(int axis=0;axis<3;++axis)for(int upper_index=0;upper_index<2;++upper_index)
		{
			const bool upper=upper_index!=0;const int tangent0=(axis+1)%3,tangent1=(axis+2)%3;
			const double expected_area=extent[tangent0]*extent[tangent1];Vec3d expected_centroid=(cell.lo+cell.hi)*0.5;
			expected_centroid[axis]=upper?cell.hi[axis]:cell.lo[axis];double represented_area=0.0;Vec3d centred_first_moment{};
			for(const auto& aperture:result.boundary_apertures)if(aperture.axis==axis&&aperture.upper==upper)
			{
				represented_area+=aperture.area;centred_first_moment=centred_first_moment+
					(aperture.centroid-expected_centroid)*aperture.area;
			}
			for(const auto& patch:result.boundary_surface_patches)if(patch.axis==axis&&patch.upper==upper)
			{
				represented_area+=patch.area;centred_first_moment=centred_first_moment+
					(patch.centroid-expected_centroid)*patch.area;
			}
			const double perimeter=2.0*(extent[tangent0]+extent[tangent1]);
			const double area_invariant_tolerance=std::max({
				static_cast<double>(options.maximum_atoms)*area_tolerance,
				16.0*contact*perimeter,4096.0*machine*expected_area});
			const double moment_invariant_tolerance=area_invariant_tolerance*std::max(1.0,h)+
				contact*expected_area;
			if(std::abs(represented_area-expected_area)>area_invariant_tolerance||
				std::sqrt(length2(centred_first_moment))>moment_invariant_tolerance)
			{
				fail(result,"local arrangement does not conservatively tile a Cartesian boundary face");return result;
			}
		}

		const double cell_volume=extent.x*extent.y*extent.z;double volume_sum=0.0;Vec3d first_moment{};
		const Vec3d expected_centroid=(cell.lo+cell.hi)*0.5;
		const double cell_surface_area=2.0*(extent.x*extent.y+extent.x*extent.z+extent.y*extent.z);
		// Contact is a measured coordinate-uncertainty bound (normally supplied by
		// TriangleBvh).  Moving every cell face by that amount sweeps at most
		// contact*surface_area volume to first order.  This gives a dimensional,
		// translation-independent invariant tolerance; it is not a case-specific
		// residual relaxation.
		const double invariant_tolerance=std::max({256.0*volume_tolerance,
			2048.0*machine*cell_volume,contact*cell_surface_area});
		const double first_moment_tolerance=invariant_tolerance*std::max(1.0,h)+
			contact*cell_volume;
		for(const auto& fragment:result.fragments)
		{
			volume_sum+=fragment.volume;
			first_moment=first_moment+(fragment.centroid-expected_centroid)*fragment.volume;
		}
		result.volume_conservation_error=std::abs(volume_sum-cell_volume);
		result.first_moment_conservation_error=std::sqrt(length2(first_moment));
		result.volume_conservation_tolerance=invariant_tolerance;
		result.first_moment_conservation_tolerance=first_moment_tolerance;
		if(result.volume_conservation_error>invariant_tolerance||
			result.first_moment_conservation_error>first_moment_tolerance)
		{fail(result,"local arrangement violates volume/first-moment conservation");return result;}
		const double surface_tolerance=std::max(256.0*area_tolerance,128.0*contact*h);
		if(std::abs(result.represented_surface_area-result.clipped_input_surface_area)>surface_tolerance)
		{fail(result,"local arrangement does not conserve clipped fabric area");return result;}
		result.valid=true;result.error.clear();return result;
	}

	int locate_local_arrangement_fragment(const LocalSurfaceArrangement& arrangement,Vec3d point,double tolerance)
	{
		if(!arrangement.valid)return -1;
		if(!(tolerance>0))tolerance=arrangement.contact_tolerance>0?arrangement.contact_tolerance:
			128.0*std::numeric_limits<double>::epsilon()*maximum_coordinate(arrangement.cell);
		if(point.x<arrangement.cell.lo.x-tolerance||point.x>arrangement.cell.hi.x+tolerance||point.y<arrangement.cell.lo.y-tolerance||point.y>arrangement.cell.hi.y+tolerance||point.z<arrangement.cell.lo.z-tolerance||point.z>arrangement.cell.hi.z+tolerance)return -1;
		for(const auto& patch:arrangement.surface_patches)
			if(point_on_convex_polygon(point,patch.polygon,patch.normal,tolerance))return -1;
		for(const auto& patch:arrangement.boundary_surface_patches)
			if(point_on_convex_polygon(point,patch.polygon,patch.normal,tolerance))return -1;
		int fragment=-1;
		for(const auto& atom:arrangement.atoms)
		{
			bool matches=true;
			for(int plane=0;plane<static_cast<int>(arrangement.planes.size());++plane)
			{
				const double value=dot(arrangement.planes[plane].normal,point)-arrangement.planes[plane].offset;
				if(std::abs(value) <= tolerance) continue;
				const int side = value < 0 ? -1 : 1;
				if(atom.plane_side[plane] != 0 && atom.plane_side[plane] != side)
				{
					matches = false;
					break;
				}
			}
			if(!matches) continue;
			if(fragment < 0) fragment = atom.fragment;
			else if(fragment != atom.fragment) return -1;
		}
		return fragment;
	}

	bool intersect_local_boundary_apertures(const LocalArrangementBoundaryAperture& a,
		const LocalArrangementBoundaryAperture& b,LocalArrangementBoundaryOverlap& overlap,double tolerance)
	{
		overlap={};if(a.axis!=b.axis||a.upper==b.upper||a.polygon.size()<3||b.polygon.size()<3)return false;
		if(!(tolerance>0))
		{
			const double coordinate_scale=std::max(polygon_coordinate_scale(a.polygon),
				polygon_coordinate_scale(b.polygon));
			tolerance=std::max({a.contact_tolerance,b.contact_tolerance,
				32.0*std::numeric_limits<double>::epsilon()*coordinate_scale});
		}
		const double a_coordinate=a.polygon.front()[a.axis];
		const double b_coordinate=b.polygon.front()[b.axis];
		if(std::abs(a_coordinate-b_coordinate) > tolerance) return false;
		Vec3d normal{};
		normal[a.axis] = 1.0;
		Polygon polygon=intersect_coplanar_convex(a.polygon,b.polygon,normal,tolerance);const Measure measure=polygon_measure(polygon,normal);if(!(measure.value>tolerance*tolerance))return false;
		overlap.axis=a.axis;overlap.area=measure.value;overlap.centroid=measure.centroid;overlap.polygon=std::move(polygon);overlap.fragment_a=a.fragment;overlap.fragment_b=b.fragment;return true;
	}

	LocalArrangementBoundarySurfaceAssembly assemble_local_boundary_surface_patches(
		const std::vector<LocalArrangementBoundarySurfacePatch>& a,
		const std::vector<LocalArrangementBoundarySurfacePatch>& b,double tolerance)
	{
		LocalArrangementBoundarySurfaceAssembly result;
		if(a.empty()&&b.empty())
		{
			result.valid=true;
			return result;
		}
		if(a.empty()||b.empty())
		{
			result.error="common Cartesian face has fabric on only one cell side";
			return result;
		}
		double coordinate_scale=1.0,stored_contact=0.0;
		for(const auto& patch:a)
		{
			coordinate_scale=std::max(coordinate_scale,polygon_coordinate_scale(patch.polygon));
			stored_contact=std::max(stored_contact,patch.contact_tolerance);
		}
		for(const auto& patch:b)
		{
			coordinate_scale=std::max(coordinate_scale,polygon_coordinate_scale(patch.polygon));
			stored_contact=std::max(stored_contact,patch.contact_tolerance);
		}
		if(!(tolerance>0.0))tolerance=std::max(stored_contact,
			32.0*std::numeric_limits<double>::epsilon()*coordinate_scale);
		const int axis=a.front().axis;
		const bool a_upper=a.front().upper,b_upper=b.front().upper;
		if(axis<0||axis>2)
		{
			result.error="boundary surface list has an invalid Cartesian axis";
			return result;
		}
		if(a.front().polygon.empty()||b.front().polygon.empty())
		{
			result.error="boundary surface list contains an empty polygon";
			return result;
		}
		const double a_coordinate=a.front().polygon.front()[axis];
		const double b_coordinate=b.front().polygon.front()[axis];
		if(b.front().axis!=axis||a_upper==b_upper||
			std::abs(a_coordinate-b_coordinate)>tolerance)
		{
			result.error="boundary surface lists are not opposite sides of one Cartesian face";
			return result;
		}
		auto list_is_consistent=[&](const std::vector<LocalArrangementBoundarySurfacePatch>& list,
			bool upper,double coordinate)
		{
			for(const auto& patch:list)
				if(patch.axis!=axis||patch.upper!=upper||patch.polygon.size()<3||
					std::abs(patch.polygon.front()[axis]-coordinate)>tolerance)return false;
			return true;
		};
		const double coordinate=a_coordinate+0.5*(b_coordinate-a_coordinate);
		if(!list_is_consistent(a,a_upper,coordinate)||!list_is_consistent(b,b_upper,coordinate))
		{
			result.error="boundary surface list contains a patch from another face";
			return result;
		}
		Vec3d face_normal{};face_normal[axis]=1.0;
		std::vector<double> covered_a(a.size(),0.0),covered_b(b.size(),0.0);
		const double angular_limit=2048.0*std::numeric_limits<double>::epsilon();
		for(std::size_t a_index=0;a_index<a.size();++a_index)
			for(std::size_t b_index=0;b_index<b.size();++b_index)
			{
				const auto& patch_a=a[a_index];const auto& patch_b=b[b_index];
				if(patch_a.source_triangle_id!=patch_b.source_triangle_id||
					patch_a.source_face_id!=patch_b.source_face_id)continue;
				const double normal_dot=dot(normalized(patch_a.normal),normalized(patch_b.normal));
				if(1.0-std::abs(normal_dot)>angular_limit)
				{
					result.error="provenance-matched boundary fabric has inconsistent support normals";
					return result;
				}
				Polygon polygon=intersect_coplanar_convex(patch_a.polygon,patch_b.polygon,
					face_normal,tolerance);
				const Measure measure=polygon_measure(polygon,face_normal);
				if(!(measure.value>tolerance*tolerance))continue;
				const bool b_plus_in_a_winding=normal_dot>=0.0?patch_b.interior_is_plus:
					!patch_b.interior_is_plus;
				if(patch_a.interior_is_plus==b_plus_in_a_winding)
				{
					result.error="boundary fabric does not assign opposite fluid sides across the common face";
					return result;
				}
				LocalArrangementBoundarySurfaceOverlap overlap;
				overlap.source_triangle_id=patch_a.source_triangle_id;
				overlap.source_face_id=patch_a.source_face_id;
				overlap.axis=static_cast<std::int8_t>(axis);overlap.area=measure.value;
				overlap.centroid=measure.centroid;overlap.normal=patch_a.normal;
				overlap.polygon=std::move(polygon);overlap.owner_is_a=a_upper;
				if(patch_a.interior_is_plus)
				{
					overlap.plus_fragment=patch_a.interior_fragment;
					overlap.minus_fragment=patch_b.interior_fragment;
				}
				else
				{
					overlap.minus_fragment=patch_a.interior_fragment;
					overlap.plus_fragment=patch_b.interior_fragment;
				}
				covered_a[a_index]+=measure.value;covered_b[b_index]+=measure.value;
				result.area+=measure.value;result.patches.push_back(std::move(overlap));
			}

		auto coverage_matches=[&](const std::vector<LocalArrangementBoundarySurfacePatch>& list,
			const std::vector<double>& covered)
		{
			for(std::size_t index=0;index<list.size();++index)
			{
				const double allowed=std::max({tolerance*tolerance,
					64.0*tolerance*std::sqrt(std::max(0.0,list[index].area)),
					512.0*std::numeric_limits<double>::epsilon()*list[index].area});
				if(std::abs(covered[index]-list[index].area)>allowed)return false;
			}
			return true;
		};
		if(!coverage_matches(a,covered_a)||!coverage_matches(b,covered_b))
		{
			result.patches.clear();result.area=0.0;
			result.error="boundary fabric subdivisions do not form one complete provenance-matched surface";
			return result;
		}
		result.valid=true;
		return result;
	}

	LocalArrangementSharedFaceAssembly assemble_canonical_local_shared_face(
		const std::array<Vec3d,4>& face_square,std::int8_t axis,
		const std::vector<LocalArrangementSharedFaceTrace>& traces,
		const std::vector<LocalArrangementSharedFaceBarrier>& barriers,
		const std::vector<LocalArrangementSharedFaceOwner>& owners_a,
		const std::vector<LocalArrangementSharedFaceOwner>& owners_b,
		const LocalArrangementSharedFaceOptions& options)
	{
		LocalArrangementSharedFaceAssembly result;
		auto fail_shared_face=[&](std::string message)
		{
			result.valid=false;result.error=std::move(message);
			result.apertures.clear();result.blocked_surfaces.clear();
			result.open_area=0.0;result.blocked_area=0.0;
		};
		if(axis<0||axis>2)
		{
			fail_shared_face("canonical shared face has an invalid Cartesian axis");return result;
		}
		if(options.maximum_tiles<1)
		{
			fail_shared_face("canonical shared face maximum_tiles must be positive");return result;
		}
		auto finite_point=[](Vec3d point)
		{
			return std::isfinite(point.x)&&std::isfinite(point.y)&&std::isfinite(point.z);
		};
		for(Vec3d point:face_square)if(!finite_point(point))
		{
			fail_shared_face("canonical shared face contains a non-finite corner");return result;
		}

		const int tangent0=(axis+1)%3,tangent1=(axis+2)%3;
		const Vec3d world_origin=face_square.front();
		const double face_coordinate=world_origin[axis];
		double world_scale=1.0;
		for(Vec3d point:face_square)world_scale=std::max({world_scale,std::abs(point.x),
			std::abs(point.y),std::abs(point.z)});
		const double machine=std::numeric_limits<double>::epsilon();
		const double plane_roundoff=16.0*machine*world_scale;
		for(Vec3d point:face_square)if(std::abs(point[axis]-face_coordinate)>plane_roundoff)
		{
			fail_shared_face("canonical shared face corners are not coplanar on the requested axis");
			return result;
		}
		auto to_local=[&](Vec3d point)
		{
			return Vec3d{point[tangent0]-world_origin[tangent0],
				point[tangent1]-world_origin[tangent1],0.0};
		};
		auto to_world=[&](Vec3d point)
		{
			Vec3d world=world_origin;world[axis]=face_coordinate;
			world[tangent0]=world_origin[tangent0]+point.x;
			world[tangent1]=world_origin[tangent1]+point.y;return world;
		};
		auto polygon_to_world=[&](const Polygon& polygon)
		{
			std::vector<Vec3d> world;world.reserve(polygon.size());
			for(Vec3d point:polygon)world.push_back(to_world(point));return world;
		};
		Polygon local_face;local_face.reserve(face_square.size());
		for(Vec3d point:face_square)local_face.push_back(to_local(point));
		double local_scale=1.0;
		for(Vec3d point:local_face)local_scale=std::max({local_scale,std::abs(point.x),
			std::abs(point.y)});
		const double length_roundoff=64.0*machine*local_scale;
		const long double arithmetic_area_floor=
			static_cast<long double>(length_roundoff)*length_roundoff;
		const Vec3d local_normal{0,0,1};

		struct LongMeasure
		{
			long double area=0.0L,moment_x=0.0L,moment_y=0.0L;
			Vec3d centroid{};
		};
		auto long_measure=[](const Polygon& polygon)
		{
			LongMeasure measured;if(polygon.size()<3)return measured;
			const long double reference_x=polygon.front().x;
			const long double reference_y=polygon.front().y;
			long double twice_area=0.0L,centroid_x_numerator=0.0L,
				centroid_y_numerator=0.0L;
			for(std::size_t index=0;index<polygon.size();++index)
			{
				const long double ax=static_cast<long double>(polygon[index].x)-reference_x;
				const long double ay=static_cast<long double>(polygon[index].y)-reference_y;
				const long double bx=static_cast<long double>(polygon[(index+1)%polygon.size()].x)-reference_x;
				const long double by=static_cast<long double>(polygon[(index+1)%polygon.size()].y)-reference_y;
				const long double cross_value=ax*by-ay*bx;
				twice_area+=cross_value;
				centroid_x_numerator+=(ax+bx)*cross_value;
				centroid_y_numerator+=(ay+by)*cross_value;
			}
			if(twice_area==0.0L)return measured;
			const long double signed_area=0.5L*twice_area;
			const long double centroid_x=reference_x+centroid_x_numerator/(3.0L*twice_area);
			const long double centroid_y=reference_y+centroid_y_numerator/(3.0L*twice_area);
			measured.area=std::abs(signed_area);
			measured.moment_x=measured.area*centroid_x;
			measured.moment_y=measured.area*centroid_y;
			measured.centroid={static_cast<double>(centroid_x),
				static_cast<double>(centroid_y),0.0};return measured;
		};
		const LongMeasure face_measure=long_measure(local_face);
		if(!(face_measure.area>arithmetic_area_floor))
		{
			fail_shared_face("canonical shared face square has non-positive area");return result;
		}
		result.face_area=static_cast<double>(face_measure.area);

		auto point_is_on_face=[&](Vec3d point)
		{
			return finite_point(point)&&std::abs(point[axis]-face_coordinate)<=plane_roundoff;
		};
		auto polygon_is_convex=[&](const Polygon& polygon)
		{
			if(polygon.size()<3)return false;int orientation=0;
			for(std::size_t index=0;index<polygon.size();++index)
			{
				const int sign=expansion_sign(orientation_2d(polygon[index],
					polygon[(index+1)%polygon.size()],polygon[(index+2)%polygon.size()],0,1));
				if(sign==0)continue;if(orientation==0)orientation=sign;
				else if(orientation!=sign)return false;
			}
			return orientation!=0;
		};
		struct PreparedOwner{Polygon polygon;int fragment=no_shared_face_fragment;};
		auto prepare_owners=[&](const std::vector<LocalArrangementSharedFaceOwner>& input,
			std::vector<PreparedOwner>& output,const char* side)
		{
			if(input.empty())
			{
				fail_shared_face(std::string("canonical shared face has an empty ")+side+
					" ownership cover");return false;
			}
			output.reserve(input.size());
			for(const auto& owner:input)
			{
				if(owner.fragment==no_shared_face_fragment||owner.polygon.size()<3)
				{
					fail_shared_face(std::string("canonical shared face has an invalid ")+side+
						" ownership polygon");return false;
				}
				Polygon polygon;polygon.reserve(owner.polygon.size());
				for(Vec3d point:owner.polygon)
				{
					if(!point_is_on_face(point))
					{
						fail_shared_face(std::string("canonical shared face ")+side+
							" ownership polygon is not on the face");return false;
					}
					polygon.push_back(to_local(point));
				}
				if(!polygon_is_convex(polygon)||!(long_measure(polygon).area>arithmetic_area_floor))
				{
					fail_shared_face(std::string("canonical shared face has a degenerate/non-convex ")+
						side+" ownership polygon");return false;
				}
				output.push_back({std::move(polygon),owner.fragment});
			}
			return true;
		};
		std::vector<PreparedOwner> prepared_a,prepared_b;
		if(!prepare_owners(owners_a,prepared_a,"side-A")||
			!prepare_owners(owners_b,prepared_b,"side-B"))return result;

		struct PreparedBarrier
		{
			std::uint32_t source_triangle_id=0,source_face_id=0;
			Vec3d normal{};Polygon polygon;
		};
		std::vector<PreparedBarrier> prepared_barriers;prepared_barriers.reserve(barriers.size());
		for(const auto& barrier:barriers)
		{
			if(barrier.polygon.size()<3||!finite_point(barrier.normal)||!(length2(barrier.normal)>0.0))
			{
				fail_shared_face("canonical shared face has an invalid coplanar barrier");return result;
			}
			Polygon polygon;polygon.reserve(barrier.polygon.size());
			for(Vec3d point:barrier.polygon)
			{
				if(!point_is_on_face(point))
				{
					fail_shared_face("canonical shared-face barrier is not coplanar with the face");
					return result;
				}
				polygon.push_back(to_local(point));
			}
			if(!polygon_is_convex(polygon)||!(long_measure(polygon).area>arithmetic_area_floor))
			{
				fail_shared_face("canonical shared face has a degenerate/non-convex barrier polygon");
				return result;
			}
			Vec3d axis_normal{};axis_normal[axis]=1.0;
			if(1.0-std::abs(dot(normalized(barrier.normal),axis_normal))>128.0*machine)
			{
				fail_shared_face("canonical shared-face barrier normal is not normal to the face");
				return result;
			}
			prepared_barriers.push_back({barrier.source_triangle_id,barrier.source_face_id,
				barrier.normal,std::move(polygon)});
		}

		struct CanonicalLine{double a=0.0,b=0.0,c=0.0;};
		std::vector<CanonicalLine> lines;
		auto append_line=[&](Vec3d a,Vec3d b,const char* kind)
		{
			const Vec3d local_a=to_local(a),local_b=to_local(b);
			const double dx=local_b.x-local_a.x,dy=local_b.y-local_a.y;
			const double magnitude=std::hypot(dx,dy);
			if(!(magnitude>length_roundoff))
			{
				fail_shared_face(std::string("canonical shared face has a degenerate ")+kind);
				return false;
			}
			CanonicalLine line{-dy/magnitude,dx/magnitude,0.0};
			if(line.a<0.0||(line.a==0.0&&line.b<0.0)){line.a=-line.a;line.b=-line.b;}
			line.c=-(line.a*local_a.x+line.b*local_a.y);lines.push_back(line);return true;
		};
		for(const auto& trace:traces)
		{
			if(!point_is_on_face(trace.a)||!point_is_on_face(trace.b))
			{
				fail_shared_face("canonical shared-face trace is not on the face");return result;
			}
			if(!append_line(trace.a,trace.b,"physical trace"))return result;
		}
		for(const auto& barrier:barriers)
			for(std::size_t edge=0;edge<barrier.polygon.size();++edge)
				if(!append_line(barrier.polygon[edge],
					barrier.polygon[(edge+1)%barrier.polygon.size()],"barrier edge"))return result;
		// The physical trace is the topological authority, but each local 3-D
		// arrangement may further subdivide the same open region.  Include both
		// ownership-cover boundaries in this one common refinement so every tile is
		// wholly owned by one fragment on each side.  Artificial extensions are safe:
		// they only subdivide a tile and cannot turn open area into a barrier.
		for(const auto& owner:owners_a)
			for(std::size_t edge=0;edge<owner.polygon.size();++edge)
				if(!append_line(owner.polygon[edge],owner.polygon[(edge+1)%owner.polygon.size()],
					"side-A ownership edge"))return result;
		for(const auto& owner:owners_b)
			for(std::size_t edge=0;edge<owner.polygon.size();++edge)
				if(!append_line(owner.polygon[edge],owner.polygon[(edge+1)%owner.polygon.size()],
					"side-B ownership edge"))return result;
		std::sort(lines.begin(),lines.end(),[](const CanonicalLine& a,const CanonicalLine& b)
		{
			if(a.a!=b.a)return a.a<b.a;
			if(a.b!=b.b)return a.b<b.b;
			return a.c<b.c;
		});
		const double angular_roundoff=128.0*machine;
		lines.erase(std::unique(lines.begin(),lines.end(),[&](const CanonicalLine& a,
			const CanonicalLine& b)
		{
			return std::abs(a.a-b.a)<=angular_roundoff&&std::abs(a.b-b.b)<=angular_roundoff&&
				std::abs(a.c-b.c)<=length_roundoff;
		}),lines.end());

		std::vector<Polygon> tiles{local_face};
		for(const CanonicalLine& line:lines)
		{
			const Vec3d line_a{-line.a*line.c,-line.b*line.c,0.0};
			const Vec3d line_b=line_a+Vec3d{line.b,-line.a,0.0};
			std::vector<Polygon> next;next.reserve(2*tiles.size());
			for(Polygon& tile:tiles)
			{
				auto [negative,positive]=split_coplanar_convex_by_line(tile,line_a,line_b,
					local_normal,length_roundoff);
				const LongMeasure negative_measure=long_measure(negative);
				const LongMeasure positive_measure=long_measure(positive);
				if(negative_measure.area>arithmetic_area_floor&&
					positive_measure.area>arithmetic_area_floor)
				{
					next.push_back(std::move(negative));next.push_back(std::move(positive));
				}
				else next.push_back(std::move(tile));
			}
			tiles=std::move(next);
			if(static_cast<int>(tiles.size())>options.maximum_tiles)
			{
				fail_shared_face("canonical shared-face arrangement exceeds maximum_tiles");return result;
			}
		}
		std::sort(tiles.begin(),tiles.end(),[&](const Polygon& a,const Polygon& b)
		{
			const LongMeasure ma=long_measure(a),mb=long_measure(b);
			if(ma.centroid.x!=mb.centroid.x)return ma.centroid.x<mb.centroid.x;
			if(ma.centroid.y!=mb.centroid.y)return ma.centroid.y<mb.centroid.y;
			return ma.area<mb.area;
		});

		auto tile_tolerance=[&](long double tile_area,std::size_t item_count)
		{
			return std::max(arithmetic_area_floor*(64.0L+8.0L*item_count),
				512.0L*static_cast<long double>(machine)*tile_area*(1.0L+item_count));
		};
		auto owner_contains_interior_point=[&](const PreparedOwner& owner,Vec3d point)
		{
			int orientation=0;
			for(std::size_t index=0;index<owner.polygon.size()&&orientation==0;++index)
				orientation=expansion_sign(orientation_2d(owner.polygon[index],
					owner.polygon[(index+1)%owner.polygon.size()],
					owner.polygon[(index+2)%owner.polygon.size()],0,1));
			if(orientation==0)return false;
			for(std::size_t edge=0;edge<owner.polygon.size();++edge)
				if(orientation*expansion_sign(orientation_2d(owner.polygon[edge],
					owner.polygon[(edge+1)%owner.polygon.size()],point,0,1))<0)return false;
			return true;
		};
		auto owner_for_tile=[&](const Polygon& tile,const LongMeasure& tile_measure,
			const std::vector<PreparedOwner>& owners,const char* side,int& owner_fragment)
		{
			// Every ownership boundary was already inserted into the canonical line
			// arrangement.  A positive tile therefore lies wholly inside or outside each
			// owner polygon.  Re-clipping those nearly coincident polygons here can lose an
			// arithmetic-scale corner (and silently close a real microscopic opening).
			// Classify the tile's strict interior point with expansion-sign predicates
			// instead.  Zero/multiple owners remain transactional gap/overlap failures.
			Vec3d interior{};for(Vec3d point:tile)interior=interior+point;
			interior=interior/static_cast<double>(tile.size());int matches=0;
			for(const PreparedOwner& owner:owners)
				if(owner_contains_interior_point(owner,interior))
				{
					owner_fragment=owner.fragment;++matches;
				}
			if(matches==1)return true;
			double nearest_distance=std::numeric_limits<double>::infinity();int nearest_fragment=0;
			for(const PreparedOwner& owner:owners)for(std::size_t edge=0;edge<owner.polygon.size();++edge)
			{
				const Vec3d a=owner.polygon[edge],b=owner.polygon[(edge+1)%owner.polygon.size()];
				const Vec3d delta=b-a;const double denominator=length2(delta);
				const double t=denominator>0.0?std::clamp(dot(interior-a,delta)/denominator,0.0,1.0):0.0;
				const double distance=std::sqrt(length2(interior-(a+delta*t)));
				if(distance<nearest_distance){nearest_distance=distance;nearest_fragment=owner.fragment;}
			}
			fail_shared_face(std::string("canonical shared-face ")+side+
				(matches==0?" ownership cover has a gap":" ownership cover overlaps")+
				" (tile area="+diagnostic_number(static_cast<double>(tile_measure.area))+
				", centroid="+diagnostic_number(tile_measure.centroid.x)+","+
				diagnostic_number(tile_measure.centroid.y)+", nearest="+
				diagnostic_number(nearest_distance)+", fragment="+
				std::to_string(nearest_fragment)+")");
			return false;
		};

		long double open_area=0.0L,blocked_area=0.0L,total_area=0.0L;
		long double total_moment_x=0.0L,total_moment_y=0.0L;
		std::vector<LocalArrangementSharedFaceAperture> staged_apertures;
		std::vector<LocalArrangementSharedFaceBlockedSurface> staged_surfaces;
		for(const Polygon& tile:tiles)
		{
			const LongMeasure measured=long_measure(tile);
			if(!(measured.area>arithmetic_area_floor))continue;
			int fragment_a=no_shared_face_fragment,fragment_b=no_shared_face_fragment;
			if(!owner_for_tile(tile,measured,prepared_a,"side-A",fragment_a)||
				!owner_for_tile(tile,measured,prepared_b,"side-B",fragment_b))return result;
			long double barrier_coverage=0.0L;int barrier_owner=-1;
			for(int index=0;index<static_cast<int>(prepared_barriers.size());++index)
			{
				const Polygon overlap=intersect_coplanar_convex(tile,
					prepared_barriers[index].polygon,local_normal,length_roundoff);
				const long double area=long_measure(overlap).area;
				if(!(area>arithmetic_area_floor))continue;
				barrier_coverage+=area;
				const long double allowed=tile_tolerance(measured.area,prepared_barriers.size());
				if(std::abs(area-measured.area)<=allowed)
				{
					if(barrier_owner>=0)
					{
						fail_shared_face("canonical shared-face tile is owned by overlapping barriers");
						return result;
					}
					barrier_owner=index;
				}
			}
			const long double barrier_allowed=tile_tolerance(measured.area,prepared_barriers.size());
			const bool blocked=barrier_coverage>barrier_allowed;
			if(blocked&&(barrier_owner<0||
				std::abs(barrier_coverage-measured.area)>barrier_allowed))
			{
				fail_shared_face("canonical shared-face tile has mixed/overlapping barrier coverage");
				return result;
			}
			const double area=static_cast<double>(measured.area);
			const Vec3d centroid=to_world(measured.centroid);
			std::vector<Vec3d> world_polygon=polygon_to_world(tile);
			if(!blocked)
			{
				LocalArrangementSharedFaceAperture aperture;aperture.area=area;
				aperture.centroid=centroid;aperture.polygon=std::move(world_polygon);
				aperture.fragment_a=fragment_a;aperture.fragment_b=fragment_b;
				staged_apertures.push_back(std::move(aperture));open_area+=measured.area;
			}
			else
			{
				const PreparedBarrier& barrier=prepared_barriers[barrier_owner];
				Vec3d axis_normal{};axis_normal[axis]=1.0;
				LocalArrangementSharedFaceBlockedSurface surface;
				surface.source_triangle_id=barrier.source_triangle_id;
				surface.source_face_id=barrier.source_face_id;surface.area=area;
				surface.centroid=centroid;surface.normal=barrier.normal;
				surface.polygon=std::move(world_polygon);
				if(dot(barrier.normal,axis_normal)>0.0)
				{
					surface.minus_fragment=fragment_a;surface.plus_fragment=fragment_b;
				}
				else
				{
					surface.plus_fragment=fragment_a;surface.minus_fragment=fragment_b;
				}
				staged_surfaces.push_back(std::move(surface));blocked_area+=measured.area;
			}
			total_area+=measured.area;total_moment_x+=measured.moment_x;
			total_moment_y+=measured.moment_y;
		}

		const long double operation_count=static_cast<long double>(1+lines.size()+tiles.size());
		const long double area_tolerance=std::max(arithmetic_area_floor*(64.0L+operation_count),
			64.0L*machine*operation_count*face_measure.area);
		const long double moment_scale=std::max({1.0L,
			static_cast<long double>(local_scale),
			static_cast<long double>(std::abs(face_measure.centroid.x)),
			static_cast<long double>(std::abs(face_measure.centroid.y))});
		const long double moment_tolerance=area_tolerance*moment_scale+
			arithmetic_area_floor*local_scale*operation_count;
		const long double area_error=std::abs(total_area-face_measure.area);
		const long double moment_x_error=total_moment_x-face_measure.moment_x;
		const long double moment_y_error=total_moment_y-face_measure.moment_y;
		const long double moment_error=std::sqrt(moment_x_error*moment_x_error+
			moment_y_error*moment_y_error);
		result.area_conservation_error=static_cast<double>(area_error);
		result.first_moment_conservation_error=static_cast<double>(moment_error);
		result.area_conservation_tolerance=static_cast<double>(area_tolerance);
		result.first_moment_conservation_tolerance=static_cast<double>(moment_tolerance);
		if(area_error>area_tolerance||moment_error>moment_tolerance||
			std::abs(open_area+blocked_area-face_measure.area)>area_tolerance)
		{
			fail_shared_face("canonical shared-face arrangement violates area/first-moment conservation");
			return result;
		}
		result.apertures=std::move(staged_apertures);
		result.blocked_surfaces=std::move(staged_surfaces);
		result.open_area=static_cast<double>(open_area);
		result.blocked_area=static_cast<double>(blocked_area);
		result.valid=true;result.error.clear();return result;
	}

	LocalArrangementSharedFaceAssembly assemble_canonical_local_shared_face(
		const std::array<Vec3d,4>& face_square,std::int8_t axis,
		const std::vector<LocalArrangementSharedFaceTrace>& traces,
		const std::vector<LocalArrangementSharedFaceBarrier>& barriers,
		const LocalArrangementSharedFaceSide& side_a,
		const LocalArrangementSharedFaceSide& side_b,
		const LocalArrangementSharedFaceOptions& options)
	{
		LocalArrangementSharedFaceAssembly result;
		auto fail_shared_face=[&](std::string message)
		{
			result.valid=false;result.error=std::move(message);
			result.apertures.clear();result.blocked_surfaces.clear();
			result.open_area=0.0;result.blocked_area=0.0;
		};
		if(axis<0||axis>2)
		{
			fail_shared_face("implicit shared face has an invalid Cartesian axis");return result;
		}
		if(options.maximum_tiles<1)
		{
			fail_shared_face("implicit shared face maximum_tiles must be positive");return result;
		}
		if(side_a.inward_axis_sign!=-1||side_b.inward_axis_sign!=1)
		{
			fail_shared_face("implicit shared face requires inward axis signs A=-1 and B=+1");
			return result;
		}
		auto finite_point=[](Vec3d point)
		{
			return std::isfinite(point.x)&&std::isfinite(point.y)&&std::isfinite(point.z);
		};
		for(Vec3d point:face_square)if(!finite_point(point))
		{
			fail_shared_face("implicit shared face contains a non-finite corner");return result;
		}

		const int tangent0=(axis+1)%3,tangent1=(axis+2)%3;
		const Vec3d world_origin=face_square.front();
		const double face_coordinate=world_origin[axis];
		double world_scale=1.0;
		for(Vec3d point:face_square)world_scale=std::max({world_scale,std::abs(point.x),
			std::abs(point.y),std::abs(point.z)});
		const double machine=std::numeric_limits<double>::epsilon();
		const double plane_roundoff=16.0*machine*world_scale;
		for(Vec3d point:face_square)if(std::abs(point[axis]-face_coordinate)>plane_roundoff)
		{
			fail_shared_face("implicit shared face corners are not coplanar on the requested axis");
			return result;
		}
		auto to_local=[&](Vec3d point)
		{
			return Vec3d{point[tangent0]-world_origin[tangent0],
				point[tangent1]-world_origin[tangent1],0.0};
		};
		auto to_world=[&](Vec3d point)
		{
			Vec3d world=world_origin;world[axis]=face_coordinate;
			world[tangent0]=world_origin[tangent0]+point.x;
			world[tangent1]=world_origin[tangent1]+point.y;return world;
		};
		auto polygon_to_world=[&](const Polygon& polygon)
		{
			std::vector<Vec3d> world;world.reserve(polygon.size());
			for(Vec3d point:polygon)world.push_back(to_world(point));return world;
		};
		Polygon local_face;local_face.reserve(face_square.size());
		for(Vec3d point:face_square)local_face.push_back(to_local(point));
		double local_scale=1.0;
		for(Vec3d point:local_face)local_scale=std::max({local_scale,std::abs(point.x),
			std::abs(point.y)});
		const Vec3d local_normal{0,0,1};

		struct LongMeasure
		{
			long double area=0.0L,moment_x=0.0L,moment_y=0.0L;
			Vec3d centroid{};
		};
		auto long_measure=[](const Polygon& polygon)
		{
			LongMeasure measured;if(polygon.size()<3)return measured;
			const long double reference_x=polygon.front().x;
			const long double reference_y=polygon.front().y;
			long double twice_area=0.0L,centroid_x_numerator=0.0L,
				centroid_y_numerator=0.0L;
			for(std::size_t index=0;index<polygon.size();++index)
			{
				const long double ax=static_cast<long double>(polygon[index].x)-reference_x;
				const long double ay=static_cast<long double>(polygon[index].y)-reference_y;
				const long double bx=static_cast<long double>(polygon[(index+1)%polygon.size()].x)-reference_x;
				const long double by=static_cast<long double>(polygon[(index+1)%polygon.size()].y)-reference_y;
				const long double cross_value=ax*by-ay*bx;
				twice_area+=cross_value;
				centroid_x_numerator+=(ax+bx)*cross_value;
				centroid_y_numerator+=(ay+by)*cross_value;
			}
			if(twice_area==0.0L)return measured;
			const long double signed_area=0.5L*twice_area;
			const long double centroid_x=reference_x+
				centroid_x_numerator/(3.0L*twice_area);
			const long double centroid_y=reference_y+
				centroid_y_numerator/(3.0L*twice_area);
			measured.area=std::abs(signed_area);
			measured.moment_x=measured.area*centroid_x;
			measured.moment_y=measured.area*centroid_y;
			measured.centroid={static_cast<double>(centroid_x),
				static_cast<double>(centroid_y),0.0};return measured;
		};
		const LongMeasure face_measure=long_measure(local_face);
		if(!(face_measure.area>0.0L))
		{
			fail_shared_face("implicit shared face square has non-positive area");return result;
		}
		result.face_area=static_cast<double>(face_measure.area);

		auto polygon_is_convex=[&](const Polygon& polygon)
		{
			if(polygon.size()<3)return false;int orientation=0;
			for(std::size_t index=0;index<polygon.size();++index)
			{
				const int sign=expansion_sign(orientation_2d(polygon[index],
					polygon[(index+1)%polygon.size()],polygon[(index+2)%polygon.size()],0,1));
				if(sign==0)continue;if(orientation==0)orientation=sign;
				else if(orientation!=sign)return false;
			}
			return orientation!=0;
		};
		if(!polygon_is_convex(local_face))
		{
			fail_shared_face("implicit shared face square is degenerate or non-convex");return result;
		}
		auto point_is_on_face=[&](Vec3d point)
		{
			return finite_point(point)&&std::abs(point[axis]-face_coordinate)<=plane_roundoff;
		};

		struct CanonicalLine
		{
			double a=0.0,b=0.0,c=0.0;
			Vec3d point_a{},point_b{};
			bool has_source_segment=false;
			std::uint64_t attachment_id=std::numeric_limits<std::uint64_t>::max();
		};
		auto canonical_line_from_coefficients=[](long double raw_a,long double raw_b,
			long double raw_c,CanonicalLine& line,int& orientation_to_raw)
		{
			const long double magnitude=std::hypot(raw_a,raw_b);
			if(!(magnitude>0.0L)||!std::isfinite(magnitude))return false;
			long double a=raw_a/magnitude,b=raw_b/magnitude,c=raw_c/magnitude;
			orientation_to_raw=1;
			if(a<0.0L||(a==0.0L&&b<0.0L))
			{
				a=-a;b=-b;c=-c;orientation_to_raw=-1;
			}
			line.a=static_cast<double>(a);line.b=static_cast<double>(b);
			line.c=static_cast<double>(c);
			line.point_a={-line.a*line.c,-line.b*line.c,0.0};
			line.point_b=line.point_a+Vec3d{line.b,-line.a,0.0};
			return std::isfinite(line.a)&&std::isfinite(line.b)&&std::isfinite(line.c)&&
				(line.a!=0.0||line.b!=0.0);
		};
		auto line_points=[](const CanonicalLine& line)
		{
			return std::pair<Vec3d,Vec3d>{line.point_a,line.point_b};
		};
		auto canonical_line_from_points=[&](Vec3d a,Vec3d b,CanonicalLine& line,
			int& orientation_to_input)
		{
			auto point_less=[](Vec3d p,Vec3d q)
			{
				if(p.x!=q.x)return p.x<q.x;
				if(p.y!=q.y)return p.y<q.y;
				return p.z<q.z;
			};
			const bool reversed=point_less(b,a);
			const Vec3d ordered_a=reversed?b:a,ordered_b=reversed?a:b;
			int orientation_to_ordered=0;
			if(!canonical_line_from_coefficients(
				-static_cast<long double>(ordered_b.y-ordered_a.y),
				static_cast<long double>(ordered_b.x-ordered_a.x),
				static_cast<long double>(ordered_b.y-ordered_a.y)*ordered_a.x-
					static_cast<long double>(ordered_b.x-ordered_a.x)*ordered_a.y,
				line,orientation_to_ordered))return false;
			line.point_a=orientation_to_ordered>0?ordered_a:ordered_b;
			line.point_b=orientation_to_ordered>0?ordered_b:ordered_a;
			line.has_source_segment=true;
			orientation_to_input=orientation_to_ordered*(reversed?-1:1);
			return true;
		};
		struct PreparedTrace
		{
			std::uint32_t source_triangle_id=0;
			CanonicalLine line{};
		};
		std::vector<PreparedTrace> prepared_traces;prepared_traces.reserve(traces.size());
		for(const auto& trace:traces)
		{
			if(!finite_point(trace.a)||!finite_point(trace.b)||
				!point_is_on_face(trace.a)||!point_is_on_face(trace.b))
			{
				fail_shared_face("implicit shared face has an invalid physical trace");return result;
			}
			CanonicalLine line;int orientation_to_input=0;
			if(!canonical_line_from_points(to_local(trace.a),to_local(trace.b),line,
				orientation_to_input))
			{
				fail_shared_face("implicit shared face has a degenerate physical trace");return result;
			}
			line.attachment_id=trace.attachment_id;
			prepared_traces.push_back({trace.source_triangle_id,line});
		}
		// One reciprocal mesh attachment is one exact physical segment.  Select its
		// representative independently of triangle/cell order and reuse it below for
		// every incident support plane.
		std::map<std::uint64_t,CanonicalLine> attached_trace_line;
		for(const PreparedTrace& trace:prepared_traces)
			if(trace.line.attachment_id!=std::numeric_limits<std::uint64_t>::max())
			{
				auto [found,inserted]=attached_trace_line.emplace(trace.line.attachment_id,
					trace.line);
				if(!inserted)
				{
					auto point_less=[](Vec3d p,Vec3d q)
					{
						if(p.x!=q.x)return p.x<q.x;
						if(p.y!=q.y)return p.y<q.y;
						return p.z<q.z;
					};
					if(point_less(trace.line.point_a,found->second.point_a)||
						(!point_less(found->second.point_a,trace.line.point_a)&&
						 !point_less(trace.line.point_a,found->second.point_a)&&
						 point_less(trace.line.point_b,found->second.point_b)))
						found->second=trace.line;
				}
			}
		for(PreparedTrace& trace:prepared_traces)
			if(trace.line.attachment_id!=std::numeric_limits<std::uint64_t>::max())
				trace.line=attached_trace_line.at(trace.line.attachment_id);
		std::vector<CanonicalLine> lines;

		struct PreparedPlane
		{
			CanonicalLine line{};
			bool has_line=false;
			int line_index=-1;
			std::int8_t constant_side=0;
			int orientation_to_plane=1;
		};
		struct PreparedSide
		{
			std::vector<PreparedPlane> planes;
			const std::vector<LocalArrangementSharedFaceRegion>* regions=nullptr;
			const char* name=nullptr;
		};
		auto prepare_side=[&](const LocalArrangementSharedFaceSide& input,
			PreparedSide& output,const char* name)
		{
			if(input.regions.empty())
			{
				fail_shared_face(std::string("implicit shared face has no ")+name+" regions");
				return false;
			}
			output.name=name;output.regions=&input.regions;output.planes.reserve(input.planes.size());
			for(std::size_t plane_index=0;plane_index<input.planes.size();++plane_index)
			{
				const auto& plane=input.planes[plane_index];
				if(!finite_point(plane.normal)||!std::isfinite(plane.offset)||
					!(length2(plane.normal)>0.0))
				{
					fail_shared_face(std::string("implicit shared face has an invalid ")+name+
						" partition plane");return false;
				}
				const long double local_a=plane.normal[tangent0];
				const long double local_b=plane.normal[tangent1];
				const long double local_c=
					static_cast<long double>(plane.normal.x)*world_origin.x+
					static_cast<long double>(plane.normal.y)*world_origin.y+
					static_cast<long double>(plane.normal.z)*world_origin.z-plane.offset;
				PreparedPlane prepared;
				if(local_a==0.0L&&local_b==0.0L)
				{
					if(local_c==0.0L)
					{
						const double inward_value=plane.normal[axis]*input.inward_axis_sign;
						if(inward_value==0.0)
						{
							fail_shared_face(std::string("implicit shared face cannot select the inward side of a ")+
								name+" face-coplanar plane");return false;
						}
						prepared.constant_side=inward_value<0.0?-1:1;
					}
					else prepared.constant_side=local_c<0.0L?-1:1;
				}
				else
				{
					if(!canonical_line_from_coefficients(local_a,local_b,local_c,
						prepared.line,prepared.orientation_to_plane))
					{
						fail_shared_face(std::string("implicit shared face cannot represent a ")+
							name+" partition-plane intersection");return false;
					}
					const CanonicalLine* physical_line=nullptr;bool one_physical_line=true;
					for(std::uint32_t support_triangle:plane.support_triangles)
						for(const PreparedTrace& trace:prepared_traces)
							if(support_triangle==trace.source_triangle_id)
							{
								if(!physical_line)physical_line=&trace.line;
								else one_physical_line=one_physical_line&&
									physical_line->a==trace.line.a&&
									physical_line->b==trace.line.b&&
									physical_line->c==trace.line.c;
							}
					if(physical_line&&one_physical_line)
					{
						// Both lines use the same deterministic half-plane orientation;
						// only the physical trace is allowed to replace the rounded offset.
						if(prepared.line.a*physical_line->a+
							prepared.line.b*physical_line->b<=0.0)
						{
							fail_shared_face(std::string("implicit shared face has inconsistent ")+
								name+" trace/partition orientation");return false;
						}
						prepared.line=*physical_line;
					}
					prepared.has_line=true;lines.push_back(prepared.line);
				}
				output.planes.push_back(prepared);
			}
			for(const auto& region:input.regions)
			{
				if(region.fragment==no_shared_face_fragment||
					region.plane_side.size()!=input.planes.size())
				{
					fail_shared_face(std::string("implicit shared face has an invalid ")+name+
						" region row");return false;
				}
				for(std::int8_t sign:region.plane_side)if(sign<-1||sign>1)
				{
					fail_shared_face(std::string("implicit shared face has a non-sign entry in a ")+
						name+" region row");return false;
				}
			}
			return true;
		};
		PreparedSide prepared_a,prepared_b;
		if(!prepare_side(side_a,prepared_a,"side-A")||
			!prepare_side(side_b,prepared_b,"side-B"))return result;

		struct PreparedBarrier
		{
			std::uint32_t source_triangle_id=0,source_face_id=0;
			struct Edge
			{
				CanonicalLine line{};
				int line_index=-1;
				std::int8_t interior_side=0;
			};
			Vec3d normal{};Polygon polygon;std::vector<Edge> edges;
		};
		std::vector<PreparedBarrier> prepared_barriers;prepared_barriers.reserve(barriers.size());
		for(const auto& barrier:barriers)
		{
			if(barrier.polygon.size()<3||!finite_point(barrier.normal)||
				!(length2(barrier.normal)>0.0))
			{
				fail_shared_face("implicit shared face has an invalid coplanar barrier");return result;
			}
			Polygon polygon;polygon.reserve(barrier.polygon.size());
			for(Vec3d point:barrier.polygon)
			{
				if(!point_is_on_face(point))
				{
					fail_shared_face("implicit shared-face barrier is not coplanar with the face");
					return result;
				}
				polygon.push_back(to_local(point));
			}
			if(!polygon_is_convex(polygon)||!(long_measure(polygon).area>0.0L))
			{
				fail_shared_face("implicit shared face has a degenerate/non-convex barrier polygon");
				return result;
			}
			if(barrier.normal[tangent0]!=0.0||barrier.normal[tangent1]!=0.0||
				barrier.normal[axis]==0.0)
			{
				fail_shared_face("implicit shared-face barrier normal is not normal to the face");
				return result;
			}
			int polygon_orientation=0;
			for(std::size_t edge=0;edge<polygon.size()&&polygon_orientation==0;++edge)
				polygon_orientation=expansion_sign(orientation_2d(polygon[edge],
					polygon[(edge+1)%polygon.size()],polygon[(edge+2)%polygon.size()],0,1));
			if(polygon_orientation==0)
			{
				fail_shared_face("implicit shared face has a barrier with no winding orientation");
				return result;
			}
			std::vector<PreparedBarrier::Edge> barrier_edges;barrier_edges.reserve(polygon.size());
			for(std::size_t edge=0;edge<polygon.size();++edge)
			{
				const Vec3d a=polygon[edge],b=polygon[(edge+1)%polygon.size()];
				auto point_less=[](Vec3d p,Vec3d q)
				{
					if(p.x!=q.x)return p.x<q.x;
					if(p.y!=q.y)return p.y<q.y;
					return p.z<q.z;
				};
				// Form the implicit line from a deterministic endpoint order.  The
				// algebraically equivalent expression based on the opposite endpoint
				// can round c by one ULP, which would create two partition lines along
				// one shared triangle edge and a fictitious overlap sliver.
				const bool reversed=point_less(b,a);
				const Vec3d ordered_a=reversed?b:a,ordered_b=reversed?a:b;
				CanonicalLine line;int orientation_to_ordered=0;
				if(!canonical_line_from_coefficients(
					-static_cast<long double>(ordered_b.y-ordered_a.y),
					static_cast<long double>(ordered_b.x-ordered_a.x),
					static_cast<long double>(ordered_b.y-ordered_a.y)*ordered_a.x-
						static_cast<long double>(ordered_b.x-ordered_a.x)*ordered_a.y,
					line,orientation_to_ordered))
				{
					fail_shared_face("implicit shared face has a degenerate barrier edge");return result;
				}
				// Preserve the canonical coefficient orientation.  The splitter's
				// negative/positive labels are defined by the directed representative
				// segment, so lexicographically sorting these endpoints would invert
				// some otherwise identical lines.
				line.point_a=orientation_to_ordered>0?ordered_a:ordered_b;
				line.point_b=orientation_to_ordered>0?ordered_b:ordered_a;
				line.has_source_segment=true;
				lines.push_back(line);barrier_edges.push_back({line,-1,
					static_cast<std::int8_t>(polygon_orientation*orientation_to_ordered*
						(reversed?-1:1))});
			}
			prepared_barriers.push_back({barrier.source_triangle_id,barrier.source_face_id,
				barrier.normal,std::move(polygon),std::move(barrier_edges)});
		}

		// Only exactly identical representable lines are duplicates.  In particular,
		// no contact/angle tolerance is allowed to merge near-distinct partitions: the
		// positive tile between them may be a real vent.
		std::sort(lines.begin(),lines.end(),[](const CanonicalLine& a,const CanonicalLine& b)
		{
			if(a.a!=b.a)return a.a<b.a;
			if(a.b!=b.b)return a.b<b.b;
			if(a.c!=b.c)return a.c<b.c;
			// When two inputs describe the same exact implicit line, keep an
			// original barrier edge as its geometric representative.  Rebuilding
			// that line from normalized coefficients need not pass bit-for-bit
			// through the source vertices and can manufacture vanishing overlap
			// tiles at shared endpoints.
			if(a.has_source_segment!=b.has_source_segment)
				return a.has_source_segment>b.has_source_segment;
			auto point_less=[](Vec3d p,Vec3d q)
			{
				if(p.x!=q.x)return p.x<q.x;
				if(p.y!=q.y)return p.y<q.y;
				return p.z<q.z;
			};
			if(point_less(a.point_a,b.point_a))return true;
			if(point_less(b.point_a,a.point_a))return false;
			return point_less(a.point_b,b.point_b);
		});
		lines.erase(std::unique(lines.begin(),lines.end(),[](const CanonicalLine& a,
			const CanonicalLine& b)
		{
			return a.a==b.a&&a.b==b.b&&a.c==b.c;
		}),lines.end());
		auto line_index=[&](const CanonicalLine& line)
		{
			const auto found=std::lower_bound(lines.begin(),lines.end(),line,
				[](const CanonicalLine& a,const CanonicalLine& b)
				{
					if(a.a!=b.a)return a.a<b.a;if(a.b!=b.b)return a.b<b.b;
					return a.c<b.c;
				});
			return found!=lines.end()&&found->a==line.a&&found->b==line.b&&found->c==line.c?
				static_cast<int>(found-lines.begin()):-1;
		};
		for(PreparedSide* side:{&prepared_a,&prepared_b})for(PreparedPlane& plane:side->planes)
			if(plane.has_line&&((plane.line_index=line_index(plane.line))<0))
			{
				fail_shared_face("implicit shared-face partition plane lost its canonical line identity");
				return result;
			}
		for(PreparedBarrier& barrier:prepared_barriers)for(auto& edge:barrier.edges)
			if((edge.line_index=line_index(edge.line))<0)
			{
				fail_shared_face("implicit shared-face barrier edge lost its canonical line identity");
				return result;
			}

		struct CanonicalTile
		{
			Polygon polygon;
			std::vector<std::int8_t> line_side;
		};
		std::vector<CanonicalTile> tiles{{local_face,
			std::vector<std::int8_t>(lines.size(),0)}};
		for(int current_line=0;current_line<static_cast<int>(lines.size());++current_line)
		{
			const CanonicalLine& line=lines[current_line];
			const auto [line_a,line_b]=line_points(line);
			std::vector<CanonicalTile> next;next.reserve(2*tiles.size());
			for(CanonicalTile& tile:tiles)
			{
				auto [negative,positive]=split_coplanar_convex_by_line(tile.polygon,line_a,line_b,
					local_normal,0.0);
				const LongMeasure negative_measure=long_measure(negative);
				const LongMeasure positive_measure=long_measure(positive);
				if(negative_measure.area>0.0L&&positive_measure.area>0.0L)
				{
					CanonicalTile negative_tile{std::move(negative),tile.line_side};
					negative_tile.line_side[current_line]=-1;
					CanonicalTile positive_tile{std::move(positive),std::move(tile.line_side)};
					positive_tile.line_side[current_line]=1;
					next.push_back(std::move(negative_tile));next.push_back(std::move(positive_tile));
				}
				else if(negative_measure.area>0.0L)
				{
					tile.polygon=std::move(negative);tile.line_side[current_line]=-1;
					next.push_back(std::move(tile));
				}
				else if(positive_measure.area>0.0L)
				{
					tile.polygon=std::move(positive);tile.line_side[current_line]=1;
					next.push_back(std::move(tile));
				}
				else
				{
					fail_shared_face("implicit shared-face partition erased a positive parent tile");
					return result;
				}
			}
			tiles=std::move(next);
			if(static_cast<int>(tiles.size())>options.maximum_tiles)
			{
				fail_shared_face("implicit shared-face arrangement exceeds maximum_tiles");return result;
			}
		}
		std::sort(tiles.begin(),tiles.end(),[&](const CanonicalTile& a,const CanonicalTile& b)
		{
			const LongMeasure ma=long_measure(a.polygon),mb=long_measure(b.polygon);
			if(ma.centroid.x!=mb.centroid.x)return ma.centroid.x<mb.centroid.x;
			if(ma.centroid.y!=mb.centroid.y)return ma.centroid.y<mb.centroid.y;
			return ma.area<mb.area;
		});

		auto fragment_for_tile=[&](const PreparedSide& side,const CanonicalTile& tile,
			const LongMeasure& tile_measure,int& fragment)
		{
			std::vector<std::int8_t> signs;signs.reserve(side.planes.size());
			for(const PreparedPlane& plane:side.planes)
			{
				if(!plane.has_line)
				{
					signs.push_back(plane.constant_side);continue;
				}
				const int canonical_sign=plane.line_index>=0&&
					plane.line_index<static_cast<int>(tile.line_side.size())?
					tile.line_side[plane.line_index]:0;
				if(canonical_sign==0)
				{
					fail_shared_face(std::string("implicit shared-face positive tile lies on a ")+
						side.name+" partition plane");return false;
				}
				signs.push_back(static_cast<std::int8_t>(canonical_sign*
					plane.orientation_to_plane));
			}
			int matches=0;
			for(const auto& region:*side.regions)
			{
				bool matches_region=true;
				for(std::size_t plane=0;plane<signs.size();++plane)
					if(region.plane_side[plane]!=0&&region.plane_side[plane]!=signs[plane])
					{matches_region=false;break;}
				if(matches_region){fragment=region.fragment;++matches;}
			}
			if(matches==1)return true;
			std::ostringstream sign_text;
			for(std::int8_t sign:signs)sign_text<<(sign<0?'-':'+');
			fail_shared_face(std::string("implicit shared-face tile matches ")+
				(matches==0?"no ":"multiple ")+side.name+" regions (signs="+
				sign_text.str()+", centroid="+diagnostic_number(tile_measure.centroid.x)+","+
				diagnostic_number(tile_measure.centroid.y)+")");return false;
		};
		auto tile_in_barrier=[](const PreparedBarrier& barrier,const CanonicalTile& tile)
		{
			for(const auto& edge:barrier.edges)
			{
				if(edge.line_index<0||edge.line_index>=static_cast<int>(tile.line_side.size())||
					tile.line_side[edge.line_index]==0)return -1;
				if(tile.line_side[edge.line_index]!=edge.interior_side)return 0;
			}
			return 1;
		};

		long double open_area=0.0L,blocked_area=0.0L,total_area=0.0L;
		long double total_moment_x=0.0L,total_moment_y=0.0L;
		std::vector<LocalArrangementSharedFaceAperture> staged_apertures;
		std::vector<LocalArrangementSharedFaceBlockedSurface> staged_surfaces;
		for(const CanonicalTile& canonical_tile:tiles)
		{
			const Polygon& tile=canonical_tile.polygon;
			const LongMeasure measured=long_measure(tile);
			if(!(measured.area>0.0L))
			{
				fail_shared_face("implicit shared-face arrangement retained a zero-area tile");
				return result;
			}
			int fragment_a=no_shared_face_fragment,fragment_b=no_shared_face_fragment;
			if(!fragment_for_tile(prepared_a,canonical_tile,measured,fragment_a)||
				!fragment_for_tile(prepared_b,canonical_tile,measured,fragment_b))return result;
			int barrier_owner=-1;
			for(int barrier=0;barrier<static_cast<int>(prepared_barriers.size());++barrier)
			{
				const int membership=tile_in_barrier(prepared_barriers[barrier],canonical_tile);
				if(membership<0)
				{
					fail_shared_face("implicit shared-face tile straddles a barrier edge despite canonical subdivision");
					return result;
				}
				if(membership==0)continue;
				if(barrier_owner>=0)
				{
					std::ostringstream message;message<<
						"implicit shared-face tile is owned by multiple barriers (triangles "<<
						prepared_barriers[barrier_owner].source_triangle_id<<" and "<<
						prepared_barriers[barrier].source_triangle_id<<", area "<<
						std::setprecision(17)<<static_cast<double>(measured.area)<<", centroid "<<
						static_cast<double>(measured.centroid.x)<<","<<
						static_cast<double>(measured.centroid.y);
					auto append_edges=[&](int owner)
					{
						message<<"; tri "<<prepared_barriers[owner].source_triangle_id<<" edges";
						for(const auto& edge:prepared_barriers[owner].edges)
							message<<" ["<<edge.line_index<<":"<<edge.line.a<<","<<
								edge.line.b<<","<<edge.line.c<<"/"<<
								static_cast<int>(edge.interior_side)<<"]";
					};
					append_edges(barrier_owner);append_edges(barrier);message<<")";
					fail_shared_face(message.str());
					return result;
				}
				barrier_owner=barrier;
			}
			const double area=static_cast<double>(measured.area);
			const Vec3d centroid=to_world(measured.centroid);
			std::vector<Vec3d> world_polygon=polygon_to_world(tile);
			if(barrier_owner<0)
			{
				LocalArrangementSharedFaceAperture aperture;aperture.area=area;
				aperture.centroid=centroid;aperture.polygon=std::move(world_polygon);
				aperture.fragment_a=fragment_a;aperture.fragment_b=fragment_b;
				staged_apertures.push_back(std::move(aperture));open_area+=measured.area;
			}
			else
			{
				const PreparedBarrier& barrier=prepared_barriers[barrier_owner];
				LocalArrangementSharedFaceBlockedSurface surface;
				surface.source_triangle_id=barrier.source_triangle_id;
				surface.source_face_id=barrier.source_face_id;surface.area=area;
				surface.centroid=centroid;surface.normal=barrier.normal;
				surface.polygon=std::move(world_polygon);
				if(barrier.normal[axis]>0.0)
				{
					surface.minus_fragment=fragment_a;surface.plus_fragment=fragment_b;
				}
				else
				{
					surface.plus_fragment=fragment_a;surface.minus_fragment=fragment_b;
				}
				staged_surfaces.push_back(std::move(surface));blocked_area+=measured.area;
			}
			total_area+=measured.area;total_moment_x+=measured.moment_x;
			total_moment_y+=measured.moment_y;
		}

		const long double operation_count=static_cast<long double>(1+lines.size()+tiles.size());
		const long double area_tolerance=64.0L*machine*operation_count*face_measure.area;
		const long double moment_scale=std::max({1.0L,
			static_cast<long double>(local_scale),
			static_cast<long double>(std::abs(face_measure.centroid.x)),
			static_cast<long double>(std::abs(face_measure.centroid.y))});
		const long double moment_tolerance=area_tolerance*moment_scale;
		const long double area_error=std::abs(total_area-face_measure.area);
		const long double moment_x_error=total_moment_x-face_measure.moment_x;
		const long double moment_y_error=total_moment_y-face_measure.moment_y;
		const long double moment_error=std::sqrt(moment_x_error*moment_x_error+
			moment_y_error*moment_y_error);
		result.area_conservation_error=static_cast<double>(area_error);
		result.first_moment_conservation_error=static_cast<double>(moment_error);
		result.area_conservation_tolerance=static_cast<double>(area_tolerance);
		result.first_moment_conservation_tolerance=static_cast<double>(moment_tolerance);
		if(area_error>area_tolerance||moment_error>moment_tolerance||
			std::abs(open_area+blocked_area-face_measure.area)>area_tolerance)
		{
			fail_shared_face("implicit shared-face arrangement violates area/first-moment conservation");
			return result;
		}
		result.apertures=std::move(staged_apertures);
		result.blocked_surfaces=std::move(staged_surfaces);
		result.open_area=static_cast<double>(open_area);
		result.blocked_area=static_cast<double>(blocked_area);
		result.valid=true;result.error.clear();return result;
	}

	LocalArrangementSharedFaceAssembly assemble_canonical_local_shared_face(
		const std::array<Vec3d,4>& face_square,std::int8_t axis,
		const std::vector<LocalArrangementSharedFaceBarrier>& barriers,
		const LocalArrangementSharedFaceSide& side_a,
		const LocalArrangementSharedFaceSide& side_b,
		const LocalArrangementSharedFaceOptions& options)
	{
		static const std::vector<LocalArrangementSharedFaceTrace> no_traces;
		return assemble_canonical_local_shared_face(face_square,axis,no_traces,barriers,
			side_a,side_b,options);
	}
}
