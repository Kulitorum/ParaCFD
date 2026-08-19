#include "core/geometry/embedded_boundary.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace paracfd::core
{
	Aabb3d UniformEbGrid::cell_box(int i,int j,int k) const { Vec3d lo=origin+Vec3d{i*h,j*h,k*h}; return {lo,lo+Vec3d{h,h,h}}; }
	Vec3d UniformEbGrid::cell_centroid(int c) const { auto q=cell_coord(c); return origin+Vec3d{(q[0]+0.5)*h,(q[1]+0.5)*h,(q[2]+0.5)*h}; }

	namespace
	{
		using Polygon=std::vector<Vec3d>;
		struct Plane { Vec3d n; double d=0; double signed_distance(Vec3d p)const{return dot(n,p)-d;} };
		struct PolyMeasure { double volume=0; Vec3d centroid{}; Polygon cap; };
		struct ClippedTriangle { std::uint32_t id=0; Polygon polygon; Vec3d normal{}; Vec3d centroid{}; double area=0; };

		Vec3d exact_point(const ExactCellPoint& point)
		{
			return {point[0],point[1],point[2]};
		}

		ExactCellPoint exact_point(Vec3d point)
		{
			return {point.x,point.y,point.z};
		}

		std::pair<double,Vec3d> triangle_measure(Vec3d a,Vec3d b,Vec3d c)
		{
			const double area=0.5*std::sqrt(length2(cross(b-a,c-a)));
			return {area,(a+b+c)/3.0};
		}

		Vec3d closest_point_on_triangle(Vec3d p,const EbOrientedBoundaryTriangle& triangle)
		{
			const Vec3d a=triangle.a,b=triangle.b,c=triangle.c;
			const Vec3d ab=b-a,ac=c-a,ap=p-a;const double d1=dot(ab,ap),d2=dot(ac,ap);
			if(d1<=0&&d2<=0)return a;const Vec3d bp=p-b;const double d3=dot(ab,bp),d4=dot(ac,bp);
			if(d3>=0&&d4<=d3)return b;const double vc=d1*d4-d3*d2;if(vc<=0&&d1>=0&&d3<=0){const double v=d1/(d1-d3);return a+ab*v;}
			const Vec3d cp=p-c;const double d5=dot(ab,cp),d6=dot(ac,cp);if(d6>=0&&d5<=d6)return c;
			const double vb=d5*d2-d1*d6;if(vb<=0&&d2>=0&&d6<=0){const double w=d2/(d2-d6);return a+ac*w;}
			const double va=d3*d6-d5*d4;if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0){const double w=(d4-d3)/((d4-d3)+(d5-d6));return b+(c-b)*w;}
			const double denominator=1.0/(va+vb+vc),v=vb*denominator,w=vc*denominator;return a+ab*v+ac*w;
		}

		double oriented_solid_angle(const EbOrientedBoundaryTriangle& triangle,Vec3d query)
		{
			const Vec3d a=triangle.a-query,b=triangle.b-query,c=triangle.c-query;
			const double la=std::sqrt(length2(a)),lb=std::sqrt(length2(b)),lc=std::sqrt(length2(c));
			const double numerator=dot(a,cross(b,c));
			const double denominator=la*lb*lc+dot(a,b)*lc+dot(b,c)*la+dot(c,a)*lb;
			return 2.0*std::atan2(numerator,denominator);
		}

		int locate_exact_fragment(const ExactEbCellLocator& locator,Vec3d point,double requested_tolerance)
		{
			const double tolerance=std::max(requested_tolerance,locator.ambiguity_tolerance),tolerance2=tolerance*tolerance;
			int containing=-1;constexpr double winding_tolerance=1.0e-6;
			for(int fragment=0;fragment<static_cast<int>(locator.fragments.size());++fragment)
			{
				double sum=0.0,compensation=0.0;
				for(const EbOrientedBoundaryTriangle& triangle:locator.fragments[fragment].triangles)
				{
					if(length2(point-closest_point_on_triangle(point,triangle))<=tolerance2)return -1;
					const double value=oriented_solid_angle(triangle,point)-compensation;
					const double next=sum+value;compensation=(next-sum)-value;sum=next;
				}
				const double winding=sum/(4.0*std::acos(-1.0)),magnitude=std::abs(winding);
				if(magnitude<=winding_tolerance)continue;
				if(std::abs(magnitude-1.0)>winding_tolerance||containing>=0)return -1;
				containing=fragment;
			}
			return containing;
		}

		Polygon clip_polygon_plane(const Polygon& in,const Plane& p,int keep_side,double eps=1e-12)
		{
			Polygon out;if(in.empty())return out;
			// side -1 is the geometric minus half-space (signed distance <= 0), side +1
			// is the plus half-space (signed distance >= 0).
			auto inside=[&](Vec3d v){const double d=p.signed_distance(v);return keep_side<0?d<=eps:d>=-eps;};
			for(std::size_t i=0;i<in.size();++i)
			{
				Vec3d a=in[i],b=in[(i+1)%in.size()];double da=p.signed_distance(a),db=p.signed_distance(b);bool ia=inside(a),ib=inside(b);
				if(ia)out.push_back(a);if(ia!=ib){double t=da/(da-db);out.push_back(a+(b-a)*t);}
			}
			return out;
		}

		Polygon clip_triangle_box(const BvhTriangle& t,const Aabb3d& b)
		{
			Polygon p{t.a,t.b,t.c};
			const Plane planes[6]={{{-1,0,0},-b.lo.x},{{1,0,0},b.hi.x},{{0,-1,0},-b.lo.y},{{0,1,0},b.hi.y},{{0,0,-1},-b.lo.z},{{0,0,1},b.hi.z}};
			for(const Plane& q:planes)p=clip_polygon_plane(p,q,-1,1e-12);return p;
		}

		std::pair<double,Vec3d> polygon_measure(const Polygon& p,Vec3d normal)
		{
			if(p.size()<3)return {};
			double area=0;Vec3d c{};for(std::size_t i=1;i+1<p.size();++i){Vec3d cr=cross(p[i]-p[0],p[i+1]-p[0]);double a=0.5*std::abs(dot(cr,normal));area+=a;c=c+(p[0]+p[i]+p[i+1])*(a/3.0);}if(area>0)c=c/area;return {area,c};
		}

		std::vector<Polygon> cube_faces(const Aabb3d& b)
		{
			const double x0=b.lo.x,y0=b.lo.y,z0=b.lo.z,x1=b.hi.x,y1=b.hi.y,z1=b.hi.z;
			return {
				{{x0,y0,z0},{x0,y0,z1},{x0,y1,z1},{x0,y1,z0}},
				{{x1,y0,z0},{x1,y1,z0},{x1,y1,z1},{x1,y0,z1}},
				{{x0,y0,z0},{x1,y0,z0},{x1,y0,z1},{x0,y0,z1}},
				{{x0,y1,z0},{x0,y1,z1},{x1,y1,z1},{x1,y1,z0}},
				{{x0,y0,z0},{x0,y1,z0},{x1,y1,z0},{x1,y0,z0}},
				{{x0,y0,z1},{x1,y0,z1},{x1,y1,z1},{x0,y1,z1}}};
		}

		Polygon ordered_cap(std::vector<Vec3d> points,Vec3d desired_normal,double tol)
		{
			Polygon unique;for(Vec3d p:points){bool seen=false;for(Vec3d q:unique)if(length2(p-q)<=tol*tol){seen=true;break;}if(!seen)unique.push_back(p);}if(unique.size()<3)return {};
			Vec3d c{};for(Vec3d p:unique)c=c+p;c=c/static_cast<double>(unique.size());Vec3d n=normalized(desired_normal);Vec3d seed=std::abs(n.x)<0.8?Vec3d{1,0,0}:Vec3d{0,1,0};Vec3d u=normalized(cross(n,seed)),v=cross(n,u);
			std::sort(unique.begin(),unique.end(),[&](Vec3d a,Vec3d b){return std::atan2(dot(a-c,v),dot(a-c,u))<std::atan2(dot(b-c,v),dot(b-c,u));});
			if(dot(cross(unique[1]-unique[0],unique[2]-unique[0]),n)<0)std::reverse(unique.begin(),unique.end());return unique;
		}

		PolyMeasure clip_box_halfspace(const Aabb3d& b,const Plane& plane,int keep_side)
		{
			std::vector<Polygon> faces;std::vector<Vec3d> cuts;const double tol=1e-10*std::max(1.0,b.hi.x-b.lo.x);
			for(const Polygon& f:cube_faces(b))
			{
				for(std::size_t i=0;i<f.size();++i){double da=plane.signed_distance(f[i]),db=plane.signed_distance(f[(i+1)%f.size()]);if((da<0&&db>0)||(da>0&&db<0))cuts.push_back(f[i]+(f[(i+1)%f.size()]-f[i])*(da/(da-db)));else if(std::abs(da)<=tol)cuts.push_back(f[i]);}
				Polygon q=clip_polygon_plane(f,plane,keep_side,tol);if(q.size()>=3)faces.push_back(std::move(q));
			}
			Polygon cap=ordered_cap(cuts,plane.n*(keep_side<0?1.0:-1.0),tol);if(cap.size()>=3)faces.push_back(cap);
			std::vector<Vec3d> verts;for(const auto& f:faces)for(Vec3d p:f)verts.push_back(p);PolyMeasure out;out.cap=cap;if(verts.empty())return out;Vec3d ref{};for(Vec3d p:verts)ref=ref+p;ref=ref/static_cast<double>(verts.size());
			for(const Polygon& f:faces)for(std::size_t i=1;i+1<f.size();++i){double vol=std::abs(dot(f[0]-ref,cross(f[i]-ref,f[i+1]-ref)))/6.0;if(vol>0){out.volume+=vol;out.centroid=out.centroid+(ref+f[0]+f[i]+f[i+1])*(vol/4.0);}}
			if(out.volume>0)out.centroid=out.centroid/out.volume;return out;
		}

		Polygon face_square(const Aabb3d& a,int axis,bool upper)
		{
			const double x=upper?a.hi.x:a.lo.x,y=upper?a.hi.y:a.lo.y,z=upper?a.hi.z:a.lo.z;
			if(axis==0)return {{x,a.lo.y,a.lo.z},{x,a.hi.y,a.lo.z},{x,a.hi.y,a.hi.z},{x,a.lo.y,a.hi.z}};
			if(axis==1)return {{a.lo.x,y,a.lo.z},{a.lo.x,y,a.hi.z},{a.hi.x,y,a.hi.z},{a.hi.x,y,a.lo.z}};
			return {{a.lo.x,a.lo.y,z},{a.hi.x,a.lo.y,z},{a.hi.x,a.hi.y,z},{a.lo.x,a.hi.y,z}};
		}

		std::pair<double,Vec3d> planar_polygon_measure(const Polygon& p,int axis)
		{
			Vec3d n{};n[axis]=1;return polygon_measure(p,n);
		}

		Polygon intersect_coplanar_convex(const Polygon& a,const Polygon& b,Vec3d normal,double tolerance)
		{
			if(a.size()<3||b.size()<3||length2(normal)==0)return {};normal=normalized(normal);double orientation=0.0;
			for(std::size_t i=1;i+1<b.size();++i)orientation+=dot(cross(b[i]-b[0],b[i+1]-b[0]),normal);
			if(std::abs(orientation)<=tolerance*tolerance)return {};Polygon result=a;
			for(std::size_t edge=0;edge<b.size()&&!result.empty();++edge)
			{
				const Vec3d edge_vector=b[(edge+1)%b.size()]-b[edge];if(length2(edge_vector)<=tolerance*tolerance)continue;
				Vec3d inward=normalized(cross(normal,edge_vector));if(orientation<0)inward=inward*-1.0;
				result=clip_polygon_plane(result,{inward,dot(inward,b[edge])},1,tolerance);
			}
			return result;
		}

		double point_segment_distance2(Vec3d point,Vec3d a,Vec3d b)
		{
			const Vec3d delta=b-a;const double denominator=length2(delta);if(!(denominator>0))return length2(point-a);const double t=std::clamp(dot(point-a,delta)/denominator,0.0,1.0);return length2(point-(a+delta*t));
		}

		bool clip_segment_box(Vec3d a,Vec3d b,const Aabb3d& box,Vec3d& clipped_a,Vec3d& clipped_b)
		{
			if(!box.valid())return false;const Vec3d delta=b-a;double t0=0,t1=1;
			for(int axis=0;axis<3;++axis)
			{
				if(std::abs(delta[axis])<1e-300){if(a[axis]<box.lo[axis]||a[axis]>box.hi[axis])return false;continue;}
				double q0=(box.lo[axis]-a[axis])/delta[axis],q1=(box.hi[axis]-a[axis])/delta[axis];if(q0>q1)std::swap(q0,q1);t0=std::max(t0,q0);t1=std::min(t1,q1);if(t0>t1)return false;
			}
			clipped_a=a+delta*t0;clipped_b=a+delta*t1;return length2(clipped_b-clipped_a)>0;
		}

		std::vector<FragmentRef> cell_fragments(const EmbeddedBoundary& eb,int cell)
		{
			const auto& c=eb.cells[cell];if(c.state==EbCellState::regular)return {regular_fragment(cell)};std::vector<FragmentRef> r;for(int i=0;i<c.fragment_count;++i)r.push_back(irregular_fragment(c.first_fragment+i));return r;
		}

		std::uint8_t fragment_surface_side_mask(const EmbeddedBoundary& eb,FragmentRef ref,std::uint32_t source_face_id)
		{
			if(fragment_is_regular(ref)||ref==invalid_fragment)return 0;const FluidFragment& fragment=eb.fragments[irregular_fragment_index(ref)];
			for(int q=0;q<fragment.surface_side_count;++q){const FragmentSurfaceSide& side=eb.fragment_surface_sides[fragment.surface_side_offset+q];if(side.source_face_id==source_face_id)return side.side_mask;}return 0;
		}

		bool fragments_are_side_compatible(const EmbeddedBoundary& eb,const TriangleBvh& bvh,Vec3d face_point,int cell_a,FragmentRef a,int cell_b,FragmentRef b,std::uint32_t* conflicting_source_face=nullptr)
		{
			if(conflicting_source_face)*conflicting_source_face=~std::uint32_t{0};
			if(a==invalid_fragment||b==invalid_fragment)return true;const EbCellTopology& ta=eb.cells[cell_a];const EbCellTopology& tb=eb.cells[cell_b];
			// A sampled fragment at a junction can touch two or more CAD faces.  If it is paired
			// with an ordinary neighbour, verify that the neighbour lies on a side represented by
			// that fragment for every nearby face.  This closes a sharp, zero-area trailing-edge
			// cusp that subcell-centre rasterisation otherwise lets percolate into a regular cell.
			// A fragment touching only one face is deliberately not constrained: fluid must remain
			// free to pass around a real fabric edge/opening.
			auto junction_matches_regular=[&](FragmentRef sampled,int sampled_cell,int regular_cell)
			{
				if(fragment_is_regular(sampled)||!eb.cells[sampled_cell].sampled_resolution)return true;
				const FluidFragment& fragment=eb.fragments[irregular_fragment_index(sampled)];if(fragment.surface_side_count<2)return true;
				const Vec3d regular_probe=face_point+(eb.grid.cell_centroid(regular_cell)-face_point)*0.25;
				int nearby_faces=0;bool mismatch=false;std::uint32_t mismatch_face=~std::uint32_t{0};
				for(int q=0;q<fragment.surface_side_count;++q)
				{
					const FragmentSurfaceSide& side=eb.fragment_surface_sides[fragment.surface_side_offset+q];
					const NearestSurfacePoint nearest=bvh.nearest_on_face(regular_probe,side.source_face_id,2.0*eb.grid.h);if(!nearest.found)continue;
					++nearby_faces;const std::uint8_t regular_side=dot(regular_probe-nearest.point,nearest.geometric_normal)>=0?2:1;
					if(!(side.side_mask&regular_side)){mismatch=true;mismatch_face=side.source_face_id;}
				}
				const bool compatible=nearby_faces<2||!mismatch;if(!compatible&&conflicting_source_face)*conflicting_source_face=mismatch_face;return compatible;
			};
			if(fragment_is_regular(a)&&!fragment_is_regular(b))return junction_matches_regular(b,cell_b,cell_a);
			if(!fragment_is_regular(a)&&fragment_is_regular(b))return junction_matches_regular(a,cell_a,cell_b);
			if(fragment_is_regular(a)||fragment_is_regular(b))return true;
			if(!ta.sampled_resolution&&!tb.sampled_resolution&&ta.source_face_id!=~std::uint32_t{0}&&ta.source_face_id==tb.source_face_id){const int side_a=eb.fragments[irregular_fragment_index(a)].side,side_b=eb.fragments[irregular_fragment_index(b)].side;if(side_a*side_b*dot(ta.plane_normal,tb.plane_normal)<=0){if(conflicting_source_face)*conflicting_source_face=ta.source_face_id;return false;}}
			auto opposite_special_side=[&](const EbCellTopology& analytic,FragmentRef analytic_ref,FragmentRef special_ref){if(analytic.source_face_id==~std::uint32_t{0}||analytic.sampled_resolution||analytic.arrangement_index>=0)return false;const int side=eb.fragments[irregular_fragment_index(analytic_ref)].side;const std::uint8_t mask=fragment_surface_side_mask(eb,special_ref,analytic.source_face_id);if(!mask)return false;const NearestSurfacePoint nearest=bvh.nearest_on_face(face_point,analytic.source_face_id,eb.grid.h);if(!nearest.found)return false;std::uint8_t expected=side>0?2:1;if(dot(analytic.plane_normal,nearest.geometric_normal)<0)expected=expected==2?1:2;if(mask!=3)return !(mask&expected);const std::uint8_t local_side=dot(face_point-nearest.point,nearest.geometric_normal)>=0?2:1;return local_side!=expected;};
			const bool special_a=ta.sampled_resolution||ta.arrangement_index>=0,special_b=tb.sampled_resolution||tb.arrangement_index>=0;
			if(!special_a&&special_b&&opposite_special_side(ta,a,b)){if(conflicting_source_face)*conflicting_source_face=ta.source_face_id;return false;}if(special_a&&!special_b&&opposite_special_side(tb,b,a)){if(conflicting_source_face)*conflicting_source_face=tb.source_face_id;return false;}
			if(special_a&&special_b){const FluidFragment& fa=eb.fragments[irregular_fragment_index(a)];for(int sa=0;sa<fa.surface_side_count;++sa){const FragmentSurfaceSide& side_a=eb.fragment_surface_sides[fa.surface_side_offset+sa];const std::uint8_t side_b=fragment_surface_side_mask(eb,b,side_a.source_face_id);if(side_b&&!(side_a.side_mask&side_b)){if(conflicting_source_face)*conflicting_source_face=side_a.source_face_id;return false;}}}
			return true;
		}
	}

	FragmentRef EmbeddedBoundary::fragment_for_side(int cell,int side)const
	{
		const auto& c=cells[cell];if(c.state==EbCellState::regular)return regular_fragment(cell);if(c.state!=EbCellState::split)return invalid_fragment;for(int i=0;i<c.fragment_count;++i)if(fragments[c.first_fragment+i].side==side)return irregular_fragment(c.first_fragment+i);return invalid_fragment;
	}
	FragmentRef EmbeddedBoundary::fragment_containing_point(int cell,Vec3d point,double tolerance)const
	{
		if(cell<0||cell>=static_cast<int>(cells.size()))return invalid_fragment;
		const EbCellTopology& topology=cells[cell];
		if(topology.state==EbCellState::regular)return regular_fragment(cell);
		if(topology.state!=EbCellState::split)return invalid_fragment;
		if(topology.exact_locator_index>=0)
		{
			if(topology.exact_locator_index>=static_cast<int>(exact_locators.size()))return invalid_fragment;
			const int local=locate_exact_fragment(exact_locators[topology.exact_locator_index],point,tolerance);
			return local>=0&&local<topology.fragment_count?
				irregular_fragment(topology.first_fragment+local):invalid_fragment;
		}
		if(topology.arrangement_index>=0)
		{
			if(topology.arrangement_index>=static_cast<int>(arrangements.size()))return invalid_fragment;
			const int local=locate_local_arrangement_fragment(arrangements[topology.arrangement_index],
				point,tolerance);
			return local>=0&&local<topology.fragment_count?
				irregular_fragment(topology.first_fragment+local):invalid_fragment;
		}
		if(topology.sampled_resolution&&topology.sampled_voxel_offset>=0)
		{
			const int r=topology.sampled_resolution;const auto coord=grid.cell_coord(cell);
			const Aabb3d box=grid.cell_box(coord[0],coord[1],coord[2]);int q[3];
			for(int axis=0;axis<3;++axis)q[axis]=std::clamp(static_cast<int>(
				std::floor((point[axis]-box.lo[axis])/grid.h*r)),0,r-1);
			const int sample=topology.sampled_voxel_offset+(q[2]*r+q[1])*r+q[0];
			return sample>=0&&sample<static_cast<int>(sampled_voxel_fragments.size())?
				sampled_voxel_fragments[sample]:invalid_fragment;
		}
		return fragment_for_side(cell,dot(topology.plane_normal,point)-topology.plane_offset>=0?1:-1);
	}
	Vec3d EmbeddedBoundary::fragment_centroid(FragmentRef r)const{return fragment_is_regular(r)?grid.cell_centroid(regular_fragment_cell(r)):fragments[irregular_fragment_index(r)].centroid;}
	double EmbeddedBoundary::fragment_volume(FragmentRef r)const{return fragment_is_regular(r)?grid.h*grid.h*grid.h:fragments[irregular_fragment_index(r)].volume;}

	EmbeddedBoundary build_embedded_boundary(const TriMesh& mesh,const TriangleBvh& bvh,const UniformEbGrid& grid,const EmbeddedBoundaryBuildOptions& opt)
	{
		if(!(grid.h>0)||grid.nx<=0||grid.ny<=0||grid.nz<=0||!(opt.min_volume_fraction>0&&opt.min_volume_fraction<0.5)||!(opt.min_aperture_area_fraction>=0&&opt.min_aperture_area_fraction<0.5))
			throw std::invalid_argument("invalid embedded-boundary grid or stabilization fraction");
		EmbeddedBoundary eb;eb.grid=grid;eb.min_volume_fraction=opt.min_volume_fraction;eb.min_aperture_area_fraction=opt.min_aperture_area_fraction;eb.cells.resize(grid.cell_count());eb.cut_face_mask.assign(grid.cell_count(),0);const double cell_volume=grid.h*grid.h*grid.h;
		// A surface whose support is inside the explicit coplanar tolerance of a
		// Cartesian face has one canonical topology: it belongs to that face.  Query
		// the tolerance shell on both cells and snap only a genuinely axis-aligned
		// triangle.  This makes ownership independent of which closed BVH/AABB test
		// first discovers coordinates rounded a few ULPs to either side of the face.
		const double aligned_snap_tolerance=std::max(opt.coplanar_distance_tolerance,1e-9*grid.h);
		const double shared_face_length_tolerance=std::max(bvh.edge_contact_tolerance(),64.0*std::numeric_limits<double>::epsilon()*grid.h);
		const double shared_face_area_tolerance=std::max(128.0*shared_face_length_tolerance*grid.h,4096.0*std::numeric_limits<double>::epsilon()*grid.h*grid.h);
		const double topology_roundoff_area_bound=std::max(shared_face_area_tolerance,64.0*std::sqrt(std::numeric_limits<double>::epsilon())*grid.h*grid.h);
		auto expanded_query=[&](const Aabb3d& box)
		{
			Aabb3d query=box;const Vec3d margin{aligned_snap_tolerance,aligned_snap_tolerance,aligned_snap_tolerance};
			query.lo=query.lo-margin;query.hi=query.hi+margin;return bvh.query_aabb(query);
		};
		auto canonical_triangle=[&](std::uint32_t id,const Aabb3d& box)
		{
			BvhTriangle triangle=bvh.triangle(id);const Vec3d raw_normal=normalized(cross(triangle.b-triangle.a,triangle.c-triangle.a));
			Vec3d* vertex[3]={&triangle.a,&triangle.b,&triangle.c};
			for(int axis=0;axis<3;++axis)
			{
				if(std::abs(std::abs(raw_normal[axis])-1.0)>opt.coplanar_angle_tolerance)continue;
				for(int upper_index=0;upper_index<2;++upper_index)
				{
					const double coordinate=upper_index?box.hi[axis]:box.lo[axis];bool all=true;
					for(Vec3d* point:vertex)all=all&&std::abs((*point)[axis]-coordinate)<=aligned_snap_tolerance;
					if(!all)continue;for(Vec3d* point:vertex)(*point)[axis]=coordinate;return triangle;
				}
			}
			return triangle;
		};
		for(int k=0;k<grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<grid.nx;++i)
		{
			const int cell=grid.cell_index(i,j,k);const Aabb3d box=grid.cell_box(i,j,k);std::vector<std::uint32_t> ids=expanded_query(box);if(ids.empty())continue;
			std::vector<ClippedTriangle> clipped;Vec3d reference_normal{};
			for(std::uint32_t id:ids)
			{
				const BvhTriangle t=canonical_triangle(id,box);Polygon p=clip_triangle_box(t,box);Vec3d n=normalized(cross(t.b-t.a,t.c-t.a));auto [area,cent]=polygon_measure(p,n);if(area<=1e-16*grid.h*grid.h)continue;
				if(clipped.empty())reference_normal=n;
				clipped.push_back({id,std::move(p),n,cent,area});
			}
			if(clipped.empty())continue;
			// Fabric coincident with a Cartesian face does not divide the interior of
			// either adjacent cell. Peel off only a complete, non-overlapping face here,
			// stage it once on the canonical lower owner, and let any remaining interior
			// sheet use its normal analytic/exact classifier. Incomplete faces stay in the
			// arrangement path so vents and separate apertures are represented exactly.
			std::array<std::vector<int>,6> boundary_groups;
			for(int q=0;q<static_cast<int>(clipped.size());++q)
			{
				for(int axis=0;axis<3;++axis)
				{
					if(std::abs(std::abs(clipped[q].normal[axis])-1.0)>opt.coplanar_angle_tolerance)continue;
					for(int upper_index=0;upper_index<2;++upper_index)
					{
						const double coordinate=upper_index?box.hi[axis]:box.lo[axis];bool all=true;
						for(Vec3d point:clipped[q].polygon)all=all&&std::abs(point[axis]-coordinate)<=aligned_snap_tolerance;
						if(all){boundary_groups[2*axis+upper_index].push_back(q);break;}
					}
				}
			}
			std::vector<std::uint8_t> peeled(clipped.size(),0);bool boundary_failure=false;
			const double boundary_length_tolerance=std::max(bvh.edge_contact_tolerance(),64.0*std::numeric_limits<double>::epsilon()*grid.h);
			const double boundary_area_tolerance=std::max(128.0*boundary_length_tolerance*grid.h,4096.0*std::numeric_limits<double>::epsilon()*grid.h*grid.h);
			for(int face=0;face<6&&!boundary_failure;++face)
			{
				const std::vector<int>& group=boundary_groups[face];if(group.empty())continue;
				const int axis=face/2;const bool upper=(face%2)!=0;double covered=0.0;bool overlap=false,finite_edge=false;
				const int coordinate_index=axis==0?i:(axis==1?j:k),coordinate_count=axis==0?grid.nx:(axis==1?grid.ny:grid.nz);
				if((!upper&&coordinate_index==0)||(upper&&coordinate_index+1==coordinate_count))
				{
					eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,
						"fabric snaps to the physical domain boundary; increase the aerodynamic domain margin"});boundary_failure=true;break;
				}
				for(int q:group)covered+=clipped[q].area;
				Vec3d face_normal{};face_normal[axis]=1.0;
				for(std::size_t a=0;a<group.size()&&!overlap;++a)for(std::size_t b=a+1;b<group.size();++b)
				{
					const Polygon common=intersect_coplanar_convex(clipped[group[a]].polygon,clipped[group[b]].polygon,face_normal,boundary_length_tolerance);
					overlap=polygon_measure(common,face_normal).first>boundary_area_tolerance;
				}
				for(int q:group)
				{
					const BvhTriangle triangle=canonical_triangle(clipped[q].id,box);const Vec3d vertex[3]={triangle.a,triangle.b,triangle.c};
					for(int edge=0;edge<3&&!finite_edge;++edge)
					{
						const FabricEdgeCertificate certificate=bvh.edge_certificate(clipped[q].id,edge);
						if(certificate.kind==FabricEdgeKind::attached&&
							certificate.cad_contact_id!=TriMesh::kNoCadContactId&&
							certificate.incident_fan_degree>=2)
						{
							Vec3d contact_a{},contact_b{};
							if(clip_segment_box(vertex[edge],vertex[(edge+1)%3],box,contact_a,contact_b)&&
								length2(contact_b-contact_a)>boundary_length_tolerance*boundary_length_tolerance)
							{
								// Keep certified contacts, including contacts on a Cartesian
								// edge, in the exact decomposition.  Peeling this face would
								// discard the fan-incidence audit constraint.
								finite_edge=true;continue;
							}
						}
						if(certificate.kind==FabricEdgeKind::attached&&certificate.role!=FabricEdgeRole::junction&&
							certificate.role!=FabricEdgeRole::unknown)continue;
						Vec3d a{},b{};if(!clip_segment_box(vertex[edge],vertex[(edge+1)%3],box,a,b)||length2(b-a)<=boundary_length_tolerance*boundary_length_tolerance)continue;
						const Vec3d midpoint=(a+b)*0.5;bool relative_interior=true;
						for(int tangent=0;tangent<3;++tangent)if(tangent!=axis)relative_interior=relative_interior&&
							midpoint[tangent]>box.lo[tangent]+boundary_length_tolerance&&midpoint[tangent]<box.hi[tangent]-boundary_length_tolerance;
						finite_edge=relative_interior;
					}
				}
				if(overlap||finite_edge||std::abs(covered-grid.h*grid.h)>boundary_area_tolerance)continue;

				int owner_i=i,owner_j=j,owner_k=k,neighbour_i=i,neighbour_j=j,neighbour_k=k;
				if(!upper){owner_i-=(axis==0);owner_j-=(axis==1);owner_k-=(axis==2);}
				else{neighbour_i+=(axis==0);neighbour_j+=(axis==1);neighbour_k+=(axis==2);}
				if(owner_i<0||owner_j<0||owner_k<0||neighbour_i>=grid.nx||neighbour_j>=grid.ny||neighbour_k>=grid.nz)
				{
					eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,
						"fabric snaps to the physical domain boundary; increase the aerodynamic domain margin"});boundary_failure=true;break;
				}
				const int owner=grid.cell_index(owner_i,owner_j,owner_k),neighbour=grid.cell_index(neighbour_i,neighbour_j,neighbour_k);
				if(!(eb.cut_face_mask[owner]&static_cast<std::uint8_t>(1u<<axis)))
				{
					eb.cut_face_mask[owner]|=static_cast<std::uint8_t>(1u<<axis);const double face_coordinate=grid.cell_box(owner_i,owner_j,owner_k).hi[axis];
					for(int q:group)
					{
						const BvhTriangle& triangle=bvh.triangle(clipped[q].id);Vec3d centroid=clipped[q].centroid;centroid[axis]=face_coordinate;
						const bool neighbour_is_plus=dot(clipped[q].normal,grid.cell_centroid(neighbour)-centroid)>0;
						SurfacePatch patch;patch.source_triangle_id=clipped[q].id;patch.source_face_id=triangle.source_face_id;patch.area=clipped[q].area;patch.centroid=centroid;patch.normal=clipped[q].normal;patch.plus_fragment=regular_fragment(neighbour_is_plus?neighbour:owner);patch.minus_fragment=regular_fragment(neighbour_is_plus?owner:neighbour);patch.provisional_face_cell=owner;patch.provisional_face_axis=static_cast<std::int8_t>(axis);eb.patches.push_back(patch);
					}
				}
				for(int q:group)peeled[q]=1;
			}
			if(boundary_failure)continue;
			if(std::any_of(peeled.begin(),peeled.end(),[](std::uint8_t value){return value!=0;}))
			{
				std::vector<ClippedTriangle> interior;interior.reserve(clipped.size());
				for(int q=0;q<static_cast<int>(clipped.size());++q)if(!peeled[q])interior.push_back(std::move(clipped[q]));
				clipped=std::move(interior);if(clipped.empty())continue;reference_normal=clipped.front().normal;
			}
			auto exact_occt_fallback=[&]()->bool
			{
				if(!opt.exact_cell_decomposer)return false;
				auto reject=[&](std::string reason)
				{
					EbCellTopology& topology=eb.cells[cell];topology=EbCellTopology{};
					topology.state=EbCellState::unresolved;
					std::vector<std::uint32_t> source;source.reserve(clipped.size());
					for(const ClippedTriangle& triangle:clipped)source.push_back(triangle.id);
					std::sort(source.begin(),source.end());source.erase(std::unique(source.begin(),source.end()),source.end());
					// The exact decomposer rejected the operation as a whole.  Every local
					// triangle was an input, but OCCT did not identify which one caused the
					// warning.  Preserve these as search context, not causal lineage.
					eb.unresolved.push_back({cell,std::move(source),std::move(reason),true});
				};
				ExactCellInput input;input.cell_min=exact_point(box.lo);input.cell_max=exact_point(box.hi);
				std::unordered_map<std::uint32_t,std::uint32_t> vertex_map;
				const double coordinate_tolerance=std::max(256.0*std::numeric_limits<double>::epsilon()*grid.h,
					opt.coplanar_distance_tolerance);
				for(const ClippedTriangle& clipped_triangle:clipped)
				{
					const std::uint32_t id=clipped_triangle.id;
					if(3*static_cast<std::size_t>(id)+2>=mesh.indices.size())
					{
						reject("exact cell decomposition input has an invalid source triangle index");return true;
					}
					const BvhTriangle triangle=canonical_triangle(id,box);const Vec3d point[3]={triangle.a,triangle.b,triangle.c};
					ExactCellTriangle exact_triangle;exact_triangle.triangle_id=id;exact_triangle.source_face_id=triangle.source_face_id;
					for(int edge=0;edge<3;++edge)
					{
						const FabricEdgeCertificate certificate=bvh.edge_certificate(id,edge);
						if(certificate.kind==FabricEdgeKind::attached&&
							certificate.cad_contact_id!=TriMesh::kNoCadContactId&&
							certificate.incident_fan_degree>=2)
						{
							exact_triangle.cad_contact_ids[edge]=certificate.cad_contact_id;
							exact_triangle.certified_fan_degrees[edge]=certificate.incident_fan_degree;
						}
					}
					for(int corner=0;corner<3;++corner)
					{
						const std::uint32_t global=mesh.indices[3*static_cast<std::size_t>(id)+corner];
						auto found=vertex_map.find(global);
						if(found==vertex_map.end())
						{
							const std::uint32_t local=static_cast<std::uint32_t>(input.vertices.size());
							input.vertices.push_back(exact_point(point[corner]));vertex_map.emplace(global,local);
							exact_triangle.vertices[corner]=local;
						}
						else
						{
							const Vec3d retained=exact_point(input.vertices[found->second]);
							if(length2(retained-point[corner])>coordinate_tolerance*coordinate_tolerance)
							{
								reject("exact cell decomposition found inconsistent coordinates for one indexed mesh vertex");return true;
							}
							exact_triangle.vertices[corner]=found->second;
						}
					}
					input.triangles.push_back(exact_triangle);
				}
				ExactCellDecomposition exact;
				try{exact=opt.exact_cell_decomposer(input);}
				catch(const std::exception& failure){reject(std::string("exact OCCT cell decomposition threw: ")+failure.what());return true;}
				if(!exact.valid()||!exact.warnings.empty())
				{
					std::ostringstream reason;reason<<"exact OCCT cell decomposition rejected the cell";
					int shown=0;for(const std::string& error:exact.errors){reason<<(shown++?"; ":": ")<<error;if(shown==3)break;}
					for(const std::string& warning:exact.warnings){reason<<(shown++?"; ":": ")<<warning;if(shown==3)break;}
					reject(reason.str());return true;
				}
				if(exact.fragments.empty()||exact.fragments.size()>static_cast<std::size_t>(std::numeric_limits<std::uint16_t>::max()))
				{
					reject("exact OCCT cell decomposition produced an invalid fragment count");return true;
				}
				auto valid_local_fragment=[&](int local){return local>=0&&local<static_cast<int>(exact.fragments.size());};
				auto polygon_from_piece=[](const std::array<ExactCellPoint,3>& piece)
				{
					return Polygon{exact_point(piece[0]),exact_point(piece[1]),exact_point(piece[2])};
				};
				LocalSurfaceArrangement arrangement;arrangement.cell=box;arrangement.valid=true;
				arrangement.contact_tolerance=std::max(bvh.edge_contact_tolerance(),exact.effective_fuzzy_tolerance);
				arrangement.ambiguity_tolerance=arrangement.contact_tolerance;
				arrangement.represented_surface_area=exact.history_surface_area_sum;
				arrangement.clipped_input_surface_area=exact.input_clipped_surface_area_sum;
				arrangement.volume_conservation_error=std::abs(exact.volume_conservation_residual);
				arrangement.first_moment_conservation_error=exact.first_moment_conservation_residual;
				arrangement.fragments.reserve(exact.fragments.size());
				ExactEbCellLocator locator;locator.ambiguity_tolerance=std::max(
					8.0*exact.effective_fuzzy_tolerance,256.0*std::numeric_limits<double>::epsilon()*grid.h);
				locator.fragments.reserve(exact.fragments.size());
				for(int local=0;local<static_cast<int>(exact.fragments.size());++local)
				{
					const ExactCellFragment& source=exact.fragments[local];
					if(source.id!=local||!(source.volume>0.0)||source.boundary_triangles.empty())
					{
						reject("exact OCCT cell decomposition returned an invalid fragment record");return true;
					}
					arrangement.fragments.push_back({source.volume,exact_point(source.centroid),{}});
					ExactEbFragmentBoundary boundary;boundary.interior_witness=exact_point(source.interior_witness);
					boundary.triangles.reserve(source.boundary_triangles.size());
					for(const ExactCellOrientedTriangle& triangle:source.boundary_triangles)
						boundary.triangles.push_back({exact_point(triangle.vertices[0]),exact_point(triangle.vertices[1]),exact_point(triangle.vertices[2])});
					locator.fragments.push_back(std::move(boundary));
				}
				for(int local=0;local<static_cast<int>(locator.fragments.size());++local)
					if(locate_exact_fragment(locator,locator.fragments[local].interior_witness,0.0)!=local)
					{
						reject("retained generalized-winding locator does not uniquely reach every exact fragment");return true;
					}
				for(const ExactCellSurfacePatch& source:exact.surface_patches)
				{
					if(source.source_triangle_id>std::numeric_limits<std::uint32_t>::max()||source.source_face_id>std::numeric_limits<std::uint32_t>::max()||
						source.region.convex_pieces.empty())
					{
						reject("exact OCCT cell decomposition returned invalid surface provenance or an empty planar region");return true;
					}
					const Vec3d normal=exact_point(source.normal);
					if(source.boundary_axis>=0)
					{
						const int axis=source.boundary_axis,coordinate=axis==0?i:(axis==1?j:k),count_axis=axis==0?grid.nx:(axis==1?grid.ny:grid.nz);
						if((!source.boundary_upper&&coordinate==0)||(source.boundary_upper&&coordinate+1==count_axis))
						{
							reject("exact fabric lies on the physical domain boundary; increase the aerodynamic domain margin");return true;
						}
						const bool plus_valid=valid_local_fragment(source.plus_fragment),minus_valid=valid_local_fragment(source.minus_fragment);
						if(plus_valid==minus_valid)
						{
							reject("exact Cartesian-face fabric does not have exactly one local interior fragment");return true;
						}
						const int interior=plus_valid?source.plus_fragment:source.minus_fragment;
						for(const auto& piece:source.region.convex_pieces)
						{
							Polygon polygon=polygon_from_piece(piece);const auto [area,centroid]=triangle_measure(polygon[0],polygon[1],polygon[2]);
							LocalArrangementBoundarySurfacePatch patch;patch.source_triangle_id=static_cast<std::uint32_t>(source.source_triangle_id);patch.source_face_id=static_cast<std::uint32_t>(source.source_face_id);patch.axis=source.boundary_axis;patch.upper=source.boundary_upper;patch.area=area;patch.centroid=centroid;patch.normal=normal;patch.polygon=std::move(polygon);patch.interior_fragment=interior;patch.interior_is_plus=plus_valid;patch.contact_tolerance=arrangement.contact_tolerance;arrangement.boundary_surface_patches.push_back(std::move(patch));
						}
					}
					else
					{
						if(!valid_local_fragment(source.plus_fragment)||!valid_local_fragment(source.minus_fragment))
						{
							reject("exact interior fabric patch has an invalid pressure-side fragment");return true;
						}
						for(const auto& piece:source.region.convex_pieces)
						{
							Polygon polygon=polygon_from_piece(piece);const auto [area,centroid]=triangle_measure(polygon[0],polygon[1],polygon[2]);
							LocalArrangementSurfacePatch patch;patch.source_triangle_id=static_cast<std::uint32_t>(source.source_triangle_id);patch.source_face_id=static_cast<std::uint32_t>(source.source_face_id);patch.area=area;patch.centroid=centroid;patch.normal=normal;patch.polygon=std::move(polygon);patch.plus_fragment=source.plus_fragment;patch.minus_fragment=source.minus_fragment;arrangement.surface_patches.push_back(std::move(patch));
						}
					}
				}
				for(const ExactCellBoxAperture& source:exact.box_apertures)
				{
					if(source.axis<0||source.axis>2||!valid_local_fragment(source.fragment)||source.region.convex_pieces.empty())
					{
						reject("exact OCCT cell decomposition returned an invalid Cartesian-face aperture");return true;
					}
					for(const auto& piece:source.region.convex_pieces)
					{
						Polygon polygon=polygon_from_piece(piece);const auto [area,centroid]=triangle_measure(polygon[0],polygon[1],polygon[2]);
						LocalArrangementBoundaryAperture aperture;aperture.axis=source.axis;aperture.upper=source.upper;aperture.area=area;aperture.centroid=centroid;aperture.polygon=std::move(polygon);aperture.fragment=source.fragment;aperture.contact_tolerance=arrangement.contact_tolerance;arrangement.boundary_apertures.push_back(std::move(aperture));
					}
				}

				const int fragment_begin=static_cast<int>(eb.fragments.size()),fragment_count=static_cast<int>(arrangement.fragments.size());
				std::vector<std::map<std::uint32_t,std::uint8_t>> side_masks(fragment_count);
				std::vector<SurfacePatch> staged_patches;staged_patches.reserve(arrangement.surface_patches.size());
				std::size_t accepted_same_fragment_count=0;double accepted_same_fragment_area=0.0;
				for(const LocalArrangementSurfacePatch& local_patch:arrangement.surface_patches)
				{
					SurfacePatch patch;patch.source_triangle_id=local_patch.source_triangle_id;patch.source_face_id=local_patch.source_face_id;patch.area=local_patch.area;patch.centroid=local_patch.centroid;patch.normal=local_patch.normal;patch.plus_fragment=irregular_fragment(fragment_begin+local_patch.plus_fragment);patch.minus_fragment=irregular_fragment(fragment_begin+local_patch.minus_fragment);side_masks[local_patch.plus_fragment][local_patch.source_face_id]|=2;side_masks[local_patch.minus_fragment][local_patch.source_face_id]|=1;if(local_patch.plus_fragment==local_patch.minus_fragment){++accepted_same_fragment_count;accepted_same_fragment_area+=local_patch.area;}staged_patches.push_back(patch);
				}
				for(const LocalArrangementBoundarySurfacePatch& patch:arrangement.boundary_surface_patches)
					side_masks[patch.interior_fragment][patch.source_face_id]|=patch.interior_is_plus?2:1;
				std::vector<FragmentSurfaceSide> staged_sides;std::vector<FluidFragment> staged_fragments;staged_fragments.reserve(fragment_count);
				for(int local=0;local<fragment_count;++local)
				{
					const LocalArrangementFragment& source=arrangement.fragments[local];FluidFragment fragment;fragment.parent_cell=cell;fragment.volume=source.volume;fragment.centroid=source.centroid;fragment.side=0;fragment.merge_target=irregular_fragment(fragment_begin+local);fragment.pressure_dof=grid.cell_count()+fragment_begin+local;fragment.surface_side_offset=static_cast<int>(eb.fragment_surface_sides.size()+staged_sides.size());for(const auto& item:side_masks[local])staged_sides.push_back({item.first,item.second});fragment.surface_side_count=static_cast<int>(eb.fragment_surface_sides.size()+staged_sides.size())-fragment.surface_side_offset;staged_fragments.push_back(fragment);
				}
				EbCellTopology topology;topology.state=EbCellState::split;topology.first_fragment=fragment_begin;topology.fragment_count=static_cast<std::uint16_t>(fragment_count);topology.arrangement_index=static_cast<int>(eb.arrangements.size());topology.exact_locator_index=static_cast<int>(eb.exact_locators.size());
				eb.fragments.insert(eb.fragments.end(),staged_fragments.begin(),staged_fragments.end());eb.fragment_surface_sides.insert(eb.fragment_surface_sides.end(),staged_sides.begin(),staged_sides.end());eb.patches.insert(eb.patches.end(),staged_patches.begin(),staged_patches.end());eb.arrangements.push_back(std::move(arrangement));eb.exact_locators.push_back(std::move(locator));eb.cells[cell]=topology;eb.irregular_cells.push_back(cell);eb.accepted_free_edge_same_fragment_patches+=accepted_same_fragment_count;eb.accepted_free_edge_same_fragment_area+=accepted_same_fragment_area;++eb.exact_decomposition_cells;eb.exact_decomposition_milliseconds+=exact.general_fuse_milliseconds;
				return true;
			};
			if(exact_occt_fallback())continue;
			auto exact_arrangement_fallback=[&]()->bool
			{
				// The subdivision option remains the product-level switch for complex-cell
				// reconstruction, but the enabled path is geometric: it has no raster and no
				// assumptions about a particular CAD model, seam layout, or fragment count.
				if(opt.complex_subdivisions<2||
					opt.allow_unverified_same_fragment_patches_for_diagnostics)return false;
				std::vector<LocalSurfaceTriangle> input;input.reserve(ids.size());
				double arrangement_contact=bvh.edge_contact_tolerance(),arrangement_clearance=bvh.edge_clearance_tolerance();
				for(std::uint32_t id:ids)
				{
					const BvhTriangle triangle=canonical_triangle(id,box);LocalSurfaceTriangle local;
					local.a=triangle.a;local.b=triangle.b;local.c=triangle.c;
					local.source_triangle_id=id;local.source_face_id=triangle.source_face_id;
					for(int edge=0;edge<3;++edge)
					{
						const FabricEdgeCertificate certificate=bvh.edge_certificate(id,edge);local.edge_kind[edge]=certificate.kind;
						if(certificate.attachment_id!=FabricEdgeCertificate::no_attachment)
							local.edge_attachment_id[edge]=certificate.attachment_id;
						arrangement_contact=std::max(arrangement_contact,certificate.contact_tolerance);
						arrangement_clearance=std::max(arrangement_clearance,certificate.clearance_tolerance);
					}
					input.push_back(local);
				}
				LocalSurfaceArrangementOptions arrangement_options;
				arrangement_options.contact_tolerance=arrangement_contact;
				arrangement_options.ambiguity_tolerance=std::max(
					arrangement_options.contact_tolerance,arrangement_clearance);
				arrangement_options.cartesian_face_tolerance=aligned_snap_tolerance;
				arrangement_options.angular_cartesian_face_tolerance=
					opt.coplanar_angle_tolerance;
				arrangement_options.angular_contact_tolerance=std::max(
					1024.0*std::numeric_limits<double>::epsilon(),
					arrangement_options.contact_tolerance/grid.h);
				arrangement_options.angular_ambiguity_tolerance=std::max(
					arrangement_options.angular_contact_tolerance,
					arrangement_options.ambiguity_tolerance/grid.h);
				LocalSurfaceArrangement arrangement=build_local_surface_arrangement(
					box,input,arrangement_options);
				auto reject=[&](std::string reason)
				{
					EbCellTopology& topology=eb.cells[cell];topology=EbCellTopology{};
					topology.state=EbCellState::unresolved;
					eb.unresolved.push_back({cell,ids,std::move(reason)});
				};
				if(!arrangement.valid)
				{
					reject("exact finite-triangle local arrangement rejected the complex cell: "+
						arrangement.error+"; refine or repair the source geometry");return true;
				}
				if(arrangement.fragments.empty()||arrangement.fragments.size()>
					static_cast<std::size_t>(std::numeric_limits<std::uint16_t>::max()))
				{
					reject("exact finite-triangle local arrangement produced an invalid fragment count; refine the cell");return true;
				}

				// Stage every globally referenced record before touching EmbeddedBoundary.  A
				// validation failure therefore cannot leave orphan fragments or patches behind.
				const int fragment_begin=static_cast<int>(eb.fragments.size());
				const int fragment_count=static_cast<int>(arrangement.fragments.size());
				auto valid_local_fragment=[&](int local){return local>=0&&local<fragment_count;};
				std::vector<std::map<std::uint32_t,std::uint8_t>> side_masks(fragment_count);
				std::vector<SurfacePatch> staged_patches;staged_patches.reserve(arrangement.surface_patches.size());
				std::size_t accepted_same_fragment_count=0;double accepted_same_fragment_area=0.0;
				for(const LocalArrangementSurfacePatch& local_patch:arrangement.surface_patches)
				{
					if(!valid_local_fragment(local_patch.plus_fragment)||
						!valid_local_fragment(local_patch.minus_fragment))
					{
						reject("exact finite-triangle local arrangement returned a surface patch with an invalid fragment reference");return true;
					}
					SurfacePatch patch;patch.source_triangle_id=local_patch.source_triangle_id;
					patch.source_face_id=local_patch.source_face_id;patch.area=local_patch.area;
					patch.centroid=local_patch.centroid;patch.normal=local_patch.normal;
					patch.plus_fragment=irregular_fragment(fragment_begin+local_patch.plus_fragment);
					patch.minus_fragment=irregular_fragment(fragment_begin+local_patch.minus_fragment);
					side_masks[local_patch.plus_fragment][local_patch.source_face_id]|=2;
					side_masks[local_patch.minus_fragment][local_patch.source_face_id]|=1;
					if(local_patch.plus_fragment==local_patch.minus_fragment)
					{
						++accepted_same_fragment_count;accepted_same_fragment_area+=local_patch.area;
					}
					staged_patches.push_back(patch);
				}
				for(const LocalArrangementBoundaryAperture& aperture:arrangement.boundary_apertures)
					if(!valid_local_fragment(aperture.fragment))
					{
						reject("exact finite-triangle local arrangement returned a boundary aperture with an invalid fragment reference");return true;
					}
				for(const LocalArrangementBoundarySurfacePatch& patch:arrangement.boundary_surface_patches)
				{
					if(!valid_local_fragment(patch.interior_fragment))
					{
						reject("exact finite-triangle local arrangement returned a boundary surface with an invalid fragment reference");return true;
					}
					side_masks[patch.interior_fragment][patch.source_face_id]|=
						patch.interior_is_plus?2:1;
				}
				std::vector<FragmentSurfaceSide> staged_sides;
				std::vector<FluidFragment> staged_fragments;staged_fragments.reserve(fragment_count);
				for(int local=0;local<fragment_count;++local)
				{
					const LocalArrangementFragment& source=arrangement.fragments[local];
					if(!(source.volume>0.0))
					{
						reject("exact finite-triangle local arrangement returned a non-positive fluid volume");return true;
					}
					FluidFragment fragment;fragment.parent_cell=cell;fragment.volume=source.volume;
					fragment.centroid=source.centroid;fragment.side=0;
					fragment.merge_target=irregular_fragment(fragment_begin+local);
					fragment.pressure_dof=grid.cell_count()+fragment_begin+local;
					fragment.surface_side_offset=static_cast<int>(eb.fragment_surface_sides.size()+staged_sides.size());
					for(const auto& item:side_masks[local])staged_sides.push_back({item.first,item.second});
					fragment.surface_side_count=static_cast<int>(eb.fragment_surface_sides.size()+staged_sides.size())-
						fragment.surface_side_offset;staged_fragments.push_back(fragment);
				}

				EbCellTopology topology;topology.state=EbCellState::split;
				topology.first_fragment=fragment_begin;
				topology.fragment_count=static_cast<std::uint16_t>(fragment_count);
				topology.arrangement_index=static_cast<int>(eb.arrangements.size());
				eb.fragments.insert(eb.fragments.end(),staged_fragments.begin(),staged_fragments.end());
				eb.fragment_surface_sides.insert(eb.fragment_surface_sides.end(),staged_sides.begin(),staged_sides.end());
				eb.patches.insert(eb.patches.end(),staged_patches.begin(),staged_patches.end());
				eb.arrangements.push_back(std::move(arrangement));eb.cells[cell]=topology;
				eb.irregular_cells.push_back(cell);
				eb.accepted_free_edge_same_fragment_patches+=accepted_same_fragment_count;
				eb.accepted_free_edge_same_fragment_area+=accepted_same_fragment_area;
				return true;
			};
			auto sampled_fallback=[&]()->bool
			{
				const int n=opt.complex_subdivisions;if(n<2)return false;const int count=n*n*n;const double micro_h=grid.h/n,micro_volume=micro_h*micro_h*micro_h;
				const std::size_t fragment_begin=eb.fragments.size(),sample_begin=eb.sampled_voxel_fragments.size();
				auto micro_index=[&](int x,int y,int z){return(z*n+y)*n+x;};auto micro_center=[&](int x,int y,int z){return box.lo+Vec3d{(x+0.5)*micro_h,(y+0.5)*micro_h,(z+0.5)*micro_h};};
				const double topology_tolerance=1e-9*grid.h;
				std::vector<int> parent(count);std::iota(parent.begin(),parent.end(),0);auto root=[&](int q){while(parent[q]!=q){parent[q]=parent[parent[q]];q=parent[q];}return q;};auto join=[&](int a,int b){a=root(a);b=root(b);if(a!=b)parent[std::max(a,b)]=std::min(a,b);};
				for(int z=0;z<n;++z)for(int y=0;y<n;++y)for(int x=0;x<n;++x)for(int axis=0;axis<3;++axis)
				{
					int q[3]={x,y,z};if(++q[axis]>=n)continue;const Vec3d a=micro_center(x,y,z),b=micro_center(q[0],q[1],q[2]);
					if(!bvh.segment_touches_surface(a,b,topology_tolerance))join(micro_index(x,y,z),micro_index(q[0],q[1],q[2]));
				}
				std::map<int,int> component;std::vector<double> component_volume;std::vector<Vec3d> component_centroid;
				for(int z=0;z<n;++z)for(int y=0;y<n;++y)for(int x=0;x<n;++x){const int r=root(micro_index(x,y,z));auto inserted=component.emplace(r,static_cast<int>(component.size()));if(inserted.second){component_volume.push_back(0);component_centroid.push_back({});}const int c=inserted.first->second;component_volume[c]+=micro_volume;component_centroid[c]=component_centroid[c]+micro_center(x,y,z)*micro_volume;}
				EbCellTopology& topology=eb.cells[cell];topology.state=EbCellState::split;topology.first_fragment=static_cast<int>(eb.fragments.size());topology.fragment_count=static_cast<std::uint16_t>(component.size());topology.sampled_voxel_offset=static_cast<int>(eb.sampled_voxel_fragments.size());topology.sampled_resolution=static_cast<std::uint8_t>(n);eb.irregular_cells.push_back(cell);eb.sampled_voxel_fragments.resize(eb.sampled_voxel_fragments.size()+count);
				for(int c=0;c<static_cast<int>(component.size());++c){component_centroid[c]=component_centroid[c]/component_volume[c];const int fragment=static_cast<int>(eb.fragments.size());eb.fragments.push_back({cell,component_volume[c],component_centroid[c],0,irregular_fragment(fragment),grid.cell_count()+fragment,0,0});}
				for(int z=0;z<n;++z)for(int y=0;y<n;++y)for(int x=0;x<n;++x){const int local=micro_index(x,y,z),c=component[root(local)];eb.sampled_voxel_fragments[topology.sampled_voxel_offset+local]=irregular_fragment(topology.first_fragment+c);}
				// Edge authorization is scoped by CAD-face provenance. A real vent on a
				// neighbouring rib or skin must not legitimize a raster leak through this
				// patch's closed seam merely because both happen to share one sampled fluid
				// component.
				std::vector<std::map<std::uint32_t,std::uint8_t>> component_edge_provenance(component.size());
				const double interior_inset=0.5*micro_h+bvh.edge_clearance_tolerance(),incident_radius=0.5*std::sqrt(3.0)*micro_h+bvh.edge_clearance_tolerance();const Vec3d inset{interior_inset,interior_inset,interior_inset};const Aabb3d interior_box{box.lo+inset,box.hi-inset};
				auto mark_edge_components=[&](Vec3d a,Vec3d b,FabricEdgeKind kind,std::uint32_t source_face_id)
				{
					if(kind==FabricEdgeKind::attached)return;Vec3d clipped_a{},clipped_b{};if(!clip_segment_box(a,b,interior_box,clipped_a,clipped_b))return;const double radius2=incident_radius*incident_radius;
					for(int z=0;z<n;++z)for(int y=0;y<n;++y)for(int x=0;x<n;++x)if(point_segment_distance2(micro_center(x,y,z),clipped_a,clipped_b)<=radius2){const int c=component[root(micro_index(x,y,z))];component_edge_provenance[c][source_face_id]|=kind==FabricEdgeKind::confirmed_free?1:2;}
				};
				for(const ClippedTriangle& cp:clipped)
				{
					const BvhTriangle& triangle=bvh.triangle(cp.id);const Vec3d point[3]={triangle.a,triangle.b,triangle.c};for(int edge=0;edge<3;++edge)mark_edge_components(point[edge],point[(edge+1)%3],bvh.edge_kind(cp.id,edge),triangle.source_face_id);
				}
				auto nearest_side_fragment=[&](Vec3d point,Vec3d normal,int side)
				{
					double best=std::numeric_limits<double>::infinity();FragmentRef result=invalid_fragment;for(int z=0;z<n;++z)for(int y=0;y<n;++y)for(int x=0;x<n;++x){Vec3d centre=micro_center(x,y,z);if(side*dot(centre-point,normal)<=1e-12)continue;double distance=length2(centre-point);if(distance<best){best=distance;result=eb.sampled_voxel_fragments[topology.sampled_voxel_offset+micro_index(x,y,z)];}}return result;
				};
				std::vector<std::map<std::uint32_t,std::uint8_t>> surface_side_masks(component.size());std::vector<SurfacePatch> sampled_patches;sampled_patches.reserve(clipped.size());std::vector<std::uint32_t> collapsed_triangles;bool ambiguous_collapse=false;std::size_t accepted_count=0,rejected_count=0,diagnostic_unverified_count=0;double accepted_area=0,rejected_area=0,diagnostic_unverified_area=0;
				for(const ClippedTriangle& cp:clipped)
				{
					const BvhTriangle& t=bvh.triangle(cp.id);SurfacePatch patch;patch.source_triangle_id=cp.id;patch.source_face_id=t.source_face_id;patch.area=cp.area;patch.centroid=cp.centroid;patch.normal=cp.normal;patch.plus_fragment=nearest_side_fragment(cp.centroid,cp.normal,1);patch.minus_fragment=nearest_side_fragment(cp.centroid,cp.normal,-1);
					if(patch.plus_fragment!=invalid_fragment&&patch.plus_fragment==patch.minus_fragment)
					{
						const int c=irregular_fragment_index(patch.plus_fragment)-topology.first_fragment;std::uint8_t edge_mask=0;if(c>=0&&c<static_cast<int>(component.size())){const auto found=component_edge_provenance[c].find(patch.source_face_id);if(found!=component_edge_provenance[c].end())edge_mask=found->second;}if(edge_mask&1){++accepted_count;accepted_area+=patch.area;}
						else if(opt.allow_unverified_same_fragment_patches_for_diagnostics){++diagnostic_unverified_count;diagnostic_unverified_area+=patch.area;}
						else{++rejected_count;rejected_area+=patch.area;collapsed_triangles.push_back(cp.id);if(edge_mask&2)ambiguous_collapse=true;}
					}
					if(patch.plus_fragment!=invalid_fragment)surface_side_masks[irregular_fragment_index(patch.plus_fragment)-topology.first_fragment][t.source_face_id]|=2;if(patch.minus_fragment!=invalid_fragment)surface_side_masks[irregular_fragment_index(patch.minus_fragment)-topology.first_fragment][t.source_face_id]|=1;sampled_patches.push_back(patch);
				}
				if(rejected_count)
				{
					eb.fragments.resize(fragment_begin);eb.sampled_voxel_fragments.resize(sample_begin);if(!eb.irregular_cells.empty()&&eb.irregular_cells.back()==cell)eb.irregular_cells.pop_back();topology=EbCellTopology{};topology.state=EbCellState::unresolved;std::sort(collapsed_triangles.begin(),collapsed_triangles.end());collapsed_triangles.erase(std::unique(collapsed_triangles.begin(),collapsed_triangles.end()),collapsed_triangles.end());eb.rejected_same_fragment_patches+=rejected_count;eb.rejected_same_fragment_area+=rejected_area;if(ambiguous_collapse)++eb.ambiguous_edge_collapse_cells;eb.unresolved.push_back({cell,std::move(collapsed_triangles),ambiguous_collapse?"sampled topology joined both fabric sides near a geometrically ambiguous edge; repair the STEP edge or increase geometry precision":"sampled topology joined both fabric sides without a confirmed free interior edge (closed seam, junction, or under-resolved topology); refine and rebuild"});return true;
				}
				for(int c=0;c<static_cast<int>(component.size());++c){FluidFragment& fragment=eb.fragments[topology.first_fragment+c];fragment.surface_side_offset=static_cast<int>(eb.fragment_surface_sides.size());for(const auto& item:surface_side_masks[c])eb.fragment_surface_sides.push_back({item.first,item.second});fragment.surface_side_count=static_cast<int>(eb.fragment_surface_sides.size())-fragment.surface_side_offset;}
				eb.patches.insert(eb.patches.end(),sampled_patches.begin(),sampled_patches.end());eb.accepted_free_edge_same_fragment_patches+=accepted_count;eb.accepted_free_edge_same_fragment_area+=accepted_area;eb.diagnostic_unverified_same_fragment_patches+=diagnostic_unverified_count;eb.diagnostic_unverified_same_fragment_area+=diagnostic_unverified_area;
				return true;
			};
			// A local plane is valid only for a complete smooth sheet. Surface-area coverage
			// alone can miss a small but topologically real opening, so a finite free or
			// ambiguous edge entering the cell (or lying in the relative interior of one
			// Cartesian face) selects the exact arrangement. Attached tessellation edges do
			// not: CAD face provenance is not topology, and the angle/distance/overlap checks
			// below already distinguish a smooth continuation from a crease, junction, or
			// second sheet. This keeps ordinary multi-face tessellation on the structured path.
			bool finite_edge_preflight=false;
			for(const ClippedTriangle& cp:clipped)
			{
				const std::uint32_t id=cp.id;
				const BvhTriangle& triangle=bvh.triangle(id);const Vec3d vertex[3]={triangle.a,triangle.b,triangle.c};
				for(int edge=0;edge<3;++edge)
				{
					Vec3d edge_a{},edge_b{};if(!clip_segment_box(vertex[edge],vertex[(edge+1)%3],box,edge_a,edge_b))continue;
					const Vec3d midpoint=(edge_a+edge_b)*0.5;bool enters_interior=true;
					const double inset=std::max(1e-12*grid.h,bvh.edge_contact_tolerance());
					for(int axis=0;axis<3;++axis)enters_interior=enters_interior&&
						midpoint[axis]>box.lo[axis]+inset&&midpoint[axis]<box.hi[axis]-inset;
					bool enters_boundary_face=false;
					if(length2(edge_b-edge_a)>inset*inset)
					{
						for(int normal_axis=0;normal_axis<3&&!enters_boundary_face;++normal_axis)
							for(int upper_index=0;upper_index<2&&!enters_boundary_face;++upper_index)
							{
								const int coordinate_index=normal_axis==0?i:(normal_axis==1?j:k);
								const int coordinate_count=normal_axis==0?grid.nx:(normal_axis==1?grid.ny:grid.nz);
								const bool has_neighbour=upper_index?coordinate_index+1<coordinate_count:coordinate_index>0;
								// A finite edge on a physical-domain face cannot reconnect the two
								// sides inside this cell. Only an internal Cartesian face needs exact
								// shared-face topology with its neighbouring cell.
								if(!has_neighbour)continue;
								const double coordinate=upper_index?box.hi[normal_axis]:box.lo[normal_axis];
								if(std::abs(edge_a[normal_axis]-coordinate)>aligned_snap_tolerance||
									std::abs(edge_b[normal_axis]-coordinate)>aligned_snap_tolerance)continue;
								bool tangent_interior=true;
								for(int tangent=0;tangent<3;++tangent)if(tangent!=normal_axis)
									tangent_interior=tangent_interior&&midpoint[tangent]>box.lo[tangent]+inset&&
										midpoint[tangent]<box.hi[tangent]-inset;
								enters_boundary_face=tangent_interior;
							}
					}
					const FabricEdgeCertificate certificate=bvh.edge_certificate(id,edge);const FabricEdgeKind kind=certificate.kind;
					finite_edge_preflight=finite_edge_preflight||
						(enters_interior||enters_boundary_face)&&
						(kind==FabricEdgeKind::confirmed_free||kind==FabricEdgeKind::unknown||
							certificate.role==FabricEdgeRole::junction||certificate.role==FabricEdgeRole::unknown||
							certificate.incident_fan_degree>2);
				}
			}
			bool positive_area_overlap=false;const double overlap_length_tolerance=std::max(
				bvh.edge_contact_tolerance(),64.0*std::numeric_limits<double>::epsilon()*grid.h);
			const double overlap_area_tolerance=std::max(64.0*overlap_length_tolerance*grid.h,
				4096.0*std::numeric_limits<double>::epsilon()*grid.h*grid.h);
			for(std::size_t a=0;a<clipped.size()&&!positive_area_overlap;++a)for(std::size_t b=a+1;b<clipped.size();++b)
			{
				const double alignment=dot(clipped[a].normal,clipped[b].normal);if(1.0-std::abs(alignment)>opt.coplanar_angle_tolerance)continue;
				const Vec3d b_normal=alignment>=0?clipped[b].normal:clipped[b].normal*-1.0;
				if(std::abs(dot(clipped[a].normal,clipped[b].centroid-clipped[a].centroid))>overlap_length_tolerance||
					std::abs(dot(b_normal,clipped[a].centroid-clipped[b].centroid))>overlap_length_tolerance)continue;
				const Polygon overlap=intersect_coplanar_convex(clipped[a].polygon,clipped[b].polygon,
					clipped[a].normal,overlap_length_tolerance);
				if(polygon_measure(overlap,clipped[a].normal).first>overlap_area_tolerance){positive_area_overlap=true;break;}
			}
			if(finite_edge_preflight||positive_area_overlap)
			{
				if(exact_arrangement_fallback())continue;
				if(sampled_fallback())continue;
				eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,
					positive_area_overlap?"positive-area coplanar fabric overlap requires exact complex-cell reconstruction":"finite or ambiguous fabric edge enters a cell without exact complex-cell reconstruction; enable complex_subdivisions or refine"});continue;
			}
			// Fit an unoriented local plane. Align each triangle normal to the first one
			// before averaging so triangle winding swaps plus/minus but cannot cancel the
			// fit. Projected-area weights make the result insensitive to tessellation size.
			Vec3d normal_sum{},centroid_sum{};double weight_sum=0;
			for(const ClippedTriangle& cp:clipped){Vec3d n=dot(cp.normal,reference_normal)<0?cp.normal*-1.0:cp.normal;normal_sum=normal_sum+n*cp.area;centroid_sum=centroid_sum+cp.centroid*cp.area;weight_sum+=cp.area;}
			Vec3d repn=normalized(normal_sum);const Vec3d fitted_centroid=weight_sum>0?centroid_sum/weight_sum:grid.cell_centroid(cell);double repd=dot(repn,fitted_centroid);bool angle_violation=length2(repn)==0,distance_violation=false;
			const double cos_limit=std::cos(std::max(0.0,opt.smooth_sheet_angle_tolerance));const double distance_limit=std::max(opt.coplanar_distance_tolerance,opt.smooth_sheet_distance_fraction*grid.h);
			for(const ClippedTriangle& cp:clipped)
			{
				if(std::abs(dot(cp.normal,repn))<cos_limit)angle_violation=true;
				for(Vec3d vertex:cp.polygon)if(std::abs(dot(repn,vertex)-repd)>distance_limit){distance_violation=true;break;}
			}
			if(angle_violation||distance_violation)
			{
				if(exact_arrangement_fallback())continue;
				if(sampled_fallback())continue;
				const char* reason=angle_violation&&distance_violation?"local fabric exceeds both smooth-sheet angle and distance bounds; refine and rebuild":(angle_violation?"local fabric exceeds smooth-sheet angle bound; refine and rebuild":"local fabric exceeds smooth-sheet distance bound; refine and rebuild");
				eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,reason});continue;
			}
			Plane plane{repn,repd};PolyMeasure minus=clip_box_halfspace(box,plane,-1),plus=clip_box_halfspace(box,plane,1);auto [section_area,section_c]=polygon_measure(minus.cap,repn);double covered=0;for(const auto& cp:clipped)covered+=polygon_measure(cp.polygon,repn).first;

			// A membrane exactly on a Cartesian face blocks that face but does not split either
			// adjacent cell. Assign ownership to the lower-index cell's positive face so the
			// tessellation is emitted once even though closed AABB queries see it from both cells.
			int aligned_axis=-1;double aligned_coordinate=0;const double aligned_tol=aligned_snap_tolerance;
			for(int axis=0;axis<3;++axis)if(std::abs(std::abs(repn[axis])-1.0)<=opt.coplanar_angle_tolerance){aligned_axis=axis;aligned_coordinate=repd/repn[axis];break;}
			if(aligned_axis>=0)
			{
				const bool on_lo=std::abs(aligned_coordinate-box.lo[aligned_axis])<=aligned_tol,on_hi=std::abs(aligned_coordinate-box.hi[aligned_axis])<=aligned_tol;
				if(on_lo||on_hi)
				{
					// Canonical owner is always the lower-coordinate cell's positive face,
					// including when only the upper cell first discovered a near-face triangle.
					int owner_i=i,owner_j=j,owner_k=k,neighbour_i=i,neighbour_j=j,neighbour_k=k;
					if(on_lo)
					{
						owner_i-=(aligned_axis==0);owner_j-=(aligned_axis==1);owner_k-=(aligned_axis==2);
					}
					else
					{
						neighbour_i+=(aligned_axis==0);neighbour_j+=(aligned_axis==1);neighbour_k+=(aligned_axis==2);
					}
					if(owner_i<0||owner_j<0||owner_k<0||neighbour_i>=grid.nx||neighbour_j>=grid.ny||neighbour_k>=grid.nz)
					{
						eb.cells[cell].state=EbCellState::unresolved;
						eb.unresolved.push_back({cell,ids,"fabric snaps to the physical domain boundary; increase the aerodynamic domain margin"});continue;
					}
					const int owner=grid.cell_index(owner_i,owner_j,owner_k);
					const int neighbour=grid.cell_index(neighbour_i,neighbour_j,neighbour_k);
					bool aligned_edge_preflight=false;
					for(const ClippedTriangle& cp:clipped)
					{
						const BvhTriangle& triangle=bvh.triangle(cp.id);const Vec3d vertex[3]={triangle.a,triangle.b,triangle.c};
						for(int edge=0;edge<3;++edge)
						{
							Vec3d edge_a{},edge_b{};if(!clip_segment_box(vertex[edge],vertex[(edge+1)%3],box,edge_a,edge_b))continue;
							if(std::abs(edge_a[aligned_axis]-aligned_coordinate)>aligned_tol||std::abs(edge_b[aligned_axis]-aligned_coordinate)>aligned_tol)continue;
							const Vec3d midpoint=(edge_a+edge_b)*0.5;bool relative_interior=true;
							for(int tangent=0;tangent<3;++tangent)if(tangent!=aligned_axis)relative_interior=relative_interior&&midpoint[tangent]>box.lo[tangent]+overlap_length_tolerance&&midpoint[tangent]<box.hi[tangent]-overlap_length_tolerance;
							const FabricEdgeKind kind=bvh.edge_kind(cp.id,edge);aligned_edge_preflight=aligned_edge_preflight||relative_interior&&(kind==FabricEdgeKind::confirmed_free||kind==FabricEdgeKind::unknown);
						}
					}
					const double exact_area_tolerance=std::max(128.0*bvh.edge_contact_tolerance()*grid.h,4096.0*std::numeric_limits<double>::epsilon()*grid.h*grid.h);
					const bool exact_area_mismatch=!(section_area>0)||std::abs(covered-section_area)>exact_area_tolerance;
					if(aligned_edge_preflight||exact_area_mismatch)
					{
						if(exact_arrangement_fallback())continue;
						eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,
							"Cartesian-aligned fabric face is not exactly complete; enable complex_subdivisions for exact opening/overlap topology"});continue;
					}
					if(eb.cut_face_mask[owner]&static_cast<std::uint8_t>(1u<<aligned_axis))continue;
					eb.cut_face_mask[owner]|=static_cast<std::uint8_t>(1u<<aligned_axis);
					for(const auto& cp:clipped)
					{
						const auto& t=bvh.triangle(cp.id);Vec3d n=cp.normal;auto [patch_area,patch_centroid]=polygon_measure(cp.polygon,n);
						patch_centroid[aligned_axis]=grid.cell_box(owner_i,owner_j,owner_k).hi[aligned_axis];
						const bool neighbour_is_plus=dot(n,grid.cell_centroid(neighbour)-patch_centroid)>0;
						SurfacePatch patch;patch.source_triangle_id=cp.id;patch.source_face_id=t.source_face_id;patch.area=patch_area;patch.centroid=patch_centroid;patch.normal=n;patch.plus_fragment=regular_fragment(neighbour_is_plus?neighbour:owner);patch.minus_fragment=regular_fragment(neighbour_is_plus?owner:neighbour);patch.provisional_face_cell=owner;patch.provisional_face_axis=static_cast<std::int8_t>(aligned_axis);eb.patches.push_back(patch);
					}
					continue;
				}
			}
			if(!(section_area>0)||covered<section_area*(1.0-opt.surface_coverage_tolerance)||covered>section_area*(1.0+opt.surface_coverage_tolerance))
			{
				if(exact_arrangement_fallback())continue;
				if(sampled_fallback())continue;
				eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,"fabric coverage is not one complete local sheet (termination, overlap, or sub-cell opening); refine rather than merge opposite sides"});continue;
			}
			EbCellTopology& ct=eb.cells[cell];ct.state=EbCellState::split;ct.first_fragment=static_cast<int>(eb.fragments.size());ct.fragment_count=2;ct.plane_normal=repn;ct.plane_offset=repd;ct.source_face_id=bvh.triangle(clipped.front().id).source_face_id;for(const ClippedTriangle& cp:clipped)if(bvh.triangle(cp.id).source_face_id!=ct.source_face_id){ct.source_face_id=~std::uint32_t{0};break;}
			ct.plane_support_offset=static_cast<int>(eb.analytic_plane_support_triangles.size());
			for(const ClippedTriangle& cp:clipped)
				eb.analytic_plane_support_triangles.push_back(cp.id);
			std::sort(eb.analytic_plane_support_triangles.begin()+ct.plane_support_offset,
				eb.analytic_plane_support_triangles.end());
			eb.analytic_plane_support_triangles.erase(std::unique(
				eb.analytic_plane_support_triangles.begin()+ct.plane_support_offset,
				eb.analytic_plane_support_triangles.end()),
				eb.analytic_plane_support_triangles.end());
			ct.plane_support_count=static_cast<std::uint32_t>(
				eb.analytic_plane_support_triangles.size()-ct.plane_support_offset);
			eb.irregular_cells.push_back(cell);
			eb.fragments.push_back({cell,minus.volume,minus.centroid,-1,irregular_fragment(ct.first_fragment),grid.cell_count()+ct.first_fragment,0,0});
			eb.fragments.push_back({cell,plus.volume,plus.centroid,1,irregular_fragment(ct.first_fragment+1),grid.cell_count()+ct.first_fragment+1,0,0});
			for(const auto& cp:clipped)
			{
				const auto& t=bvh.triangle(cp.id);Vec3d n=cp.normal;auto [area,cent]=polygon_measure(cp.polygon,n);const bool same=dot(n,repn)>=0;
				SurfacePatch patch;patch.source_triangle_id=cp.id;patch.source_face_id=t.source_face_id;patch.area=area;patch.centroid=cent;patch.normal=n;patch.plus_fragment=eb.fragment_for_side(cell,same?1:-1);patch.minus_fragment=eb.fragment_for_side(cell,same?-1:1);eb.patches.push_back(patch);
			}
		}

		// A sheet side can become thinner than one fallback microcell near a sharp
		// trailing edge. The surface still exists and must not disappear from pressure
		// traction or wall momentum. Map only that unresolved pressure sample to the
		// nearest represented side of the same CAD face. Normal reversal swaps the
		// candidate's plus/minus label, preserving winding invariance, and the bounded
		// distance prevents a folded/distant part of one CAD face from becoming a donor.
		for(std::size_t patch_index=0;patch_index<eb.patches.size();++patch_index)
		{
			SurfacePatch& patch=eb.patches[patch_index];
			// This recovery maps only a pressure sampling side; it creates no volume or flux
			// connection.  Use CAD-face provenance plus locality, rather than the much stricter
			// planar-cell angle bound: a valid smooth airfoil face turns rapidly at its nose.
			constexpr double minimum_recovery_alignment=1e-3;
			auto recover=[&](FragmentRef& target,bool plus)
			{
				if(target!=invalid_fragment)return 0.0;double best=std::numeric_limits<double>::infinity();FragmentRef donor=invalid_fragment;
				for(std::size_t candidate_index=0;candidate_index<eb.patches.size();++candidate_index)
				{
					if(candidate_index==patch_index)continue;const SurfacePatch& candidate=eb.patches[candidate_index];if(candidate.source_face_id!=patch.source_face_id)continue;const double alignment=dot(candidate.normal,patch.normal);if(std::abs(alignment)<minimum_recovery_alignment)continue;const FragmentRef side=alignment>=0?(plus?candidate.plus_fragment:candidate.minus_fragment):(plus?candidate.minus_fragment:candidate.plus_fragment);if(side==invalid_fragment)continue;const double distance2=length2(candidate.centroid-patch.centroid);if(distance2<best){best=distance2;donor=side;}
				}
				const double distance=std::sqrt(best);if(donor==invalid_fragment||distance>2.5*grid.h)return distance;target=donor;++eb.recovered_subgrid_surface_sides;eb.maximum_surface_side_recovery_distance=std::max(eb.maximum_surface_side_recovery_distance,distance);return distance;
			};
			const double plus_recovery=recover(patch.plus_fragment,true),minus_recovery=recover(patch.minus_fragment,false);
			if(patch.plus_fragment==invalid_fragment||patch.minus_fragment==invalid_fragment)
			{
				if(opt.allow_unverified_same_fragment_patches_for_diagnostics)
				{
					++eb.diagnostic_unmapped_surface_patches;
					eb.diagnostic_unmapped_surface_area+=patch.area;
				}
				else
				{
					const int i=std::clamp(static_cast<int>(std::floor((patch.centroid.x-grid.origin.x)/grid.h)),0,grid.nx-1),j=std::clamp(static_cast<int>(std::floor((patch.centroid.y-grid.origin.y)/grid.h)),0,grid.ny-1),k=std::clamp(static_cast<int>(std::floor((patch.centroid.z-grid.origin.z)/grid.h)),0,grid.nz-1),cell=grid.cell_index(i,j,k);const double missing_distance=patch.plus_fragment==invalid_fragment?plus_recovery:minus_recovery;const std::string detail=std::isfinite(missing_distance)?("nearest same-face donor is "+std::to_string(missing_distance/grid.h)+" h away") : "no same-face donor exists";eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,{patch.source_triangle_id},"sampled sharp-edge surface side has no local same-face donor ("+detail+"); refine and rebuild"});
				}
			}
		}
		if(opt.allow_unverified_same_fragment_patches_for_diagnostics)
			std::erase_if(eb.patches,[](const SurfacePatch& patch)
			{
				return patch.plus_fragment==invalid_fragment||patch.minus_fragment==invalid_fragment;
			});

		// Cartesian-aligned patches are discovered while cells are classified in index
		// order.  If a neighbour later becomes split, the canonical shared-face
		// transaction below replaces the provisional patch and assigns both pressure
		// sides from the completed implicit region partitions.  Do not epsilon-probe an
		// exact locator here: the requested point can lie inside its ambiguity band.
		auto connection_leg_crosses_fabric=[&](FragmentRef a,FragmentRef b,Vec3d face_centroid)
		{
			// The finite-volume connection is a pair of straight centroid-to-aperture
			// legs. If either leg reaches fabric (including at the nominal aperture
			// centroid), the clipped polygon is a false overlap at a sheet seam/cusp,
			// not an open fluid face. This catches a sharp closed trailing edge without
			// closing a real off-centre opening in a partially covered Cartesian face.
			return bvh.intersect_segment(eb.fragment_centroid(a),face_centroid,1e-9,1.0).hit||
				bvh.intersect_segment(eb.fragment_centroid(b),face_centroid,1e-9,1.0).hit;
		};
		auto reject_cross_fabric_aperture=[&](double area){++eb.rejected_cross_fabric_apertures;eb.rejected_cross_fabric_aperture_area+=area;};
		auto sampled_junction_fragment=[&](FragmentRef ref,int cell)
		{
			if(fragment_is_regular(ref))return false;const EbCellTopology& topology=eb.cells[cell];
			// An analytic fit spanning more than one CAD face is itself a junction cell.
			// `source_face_id` is deliberately invalidated during classification when the
			// contributing clipped triangles do not share one provenance ID.
			if(!topology.sampled_resolution)return topology.state==EbCellState::split&&topology.source_face_id==~std::uint32_t{0};
			// A component on one side of a closed CAD seam may itself touch only one of
			// the incident faces. Classify the sampled CELL as a junction from the union
			// of face provenance on all of its fragments, otherwise two neighbouring
			// sampled cells can each expose one half of a closed leading-edge seam.
			std::uint32_t first_face=~std::uint32_t{0};
			for(int local=0;local<topology.fragment_count;++local)
			{
				const FluidFragment& fragment=eb.fragments[topology.first_fragment+local];
				for(int q=0;q<fragment.surface_side_count;++q)
				{
					const std::uint32_t face=eb.fragment_surface_sides[fragment.surface_side_offset+q].source_face_id;
					if(first_face==~std::uint32_t{0})first_face=face;else if(face!=first_face)return true;
				}
			}
			return false;
		};
		// A side conflict may close an apparent overlap only when the EXACT source CAD
		// face is a closed two-sector sheet continuation.  Certify that globally from
		// its complete boundary topology: every non-tessellation half-edge must be an
		// attached manifold seam with fan degree two.  A real vent/free edge, unknown
		// provenance, or 3+ sector junction therefore defeats closure automatically.
		struct SourceFaceClosure{bool closed=true;std::size_t boundary_half_edges=0;};
		std::map<std::uint32_t,SourceFaceClosure> source_face_closure;
		for(std::uint32_t triangle_id=0;triangle_id<bvh.triangle_count();++triangle_id)
		{
			SourceFaceClosure& closure=source_face_closure[bvh.triangle(triangle_id).source_face_id];
			for(int edge=0;edge<3;++edge)
			{
				const FabricEdgeCertificate certificate=bvh.edge_certificate(triangle_id,edge);if(certificate.role==FabricEdgeRole::tessellation_interior)continue;++closure.boundary_half_edges;
				closure.closed=closure.closed&&certificate.kind==FabricEdgeKind::attached&&certificate.role==FabricEdgeRole::manifold_seam&&certificate.incident_fan_degree==2;
			}
		}
		auto source_face_is_closed=[&](std::uint32_t source_face_id){const auto found=source_face_closure.find(source_face_id);return source_face_id!=~std::uint32_t{0}&&found!=source_face_closure.end()&&found->second.closed;};
		auto subresolution_closed_corner_conflict=[&](FragmentRef corner,FragmentRef other,Vec3d face_point)
		{
			// A roundoff-scale carrier at the intersection of two or more closed CAD
			// faces is not evidence of an opening.  Prove the local regions differ by
			// evaluating the other carrier just across the Cartesian face against every
			// winding-aware side certificate on the corner fragment.  Any free/unknown/
			// junction boundary makes its source face non-closed above and defeats this
			// certificate, so a real microscopic vent remains connected.
			if(corner==invalid_fragment||other==invalid_fragment||fragment_is_regular(corner))return false;
			const FluidFragment& fragment=eb.fragments[irregular_fragment_index(corner)];if(fragment.surface_side_count<2)return false;
			const Vec3d probe=face_point+(eb.fragment_centroid(other)-face_point)*0.25;std::set<std::uint32_t> distinct_faces;int nearby_faces=0,mismatches=0;
			for(int q=0;q<fragment.surface_side_count;++q)
			{
				const FragmentSurfaceSide& side=eb.fragment_surface_sides[fragment.surface_side_offset+q];if(side.side_mask==0||side.side_mask==3||!source_face_is_closed(side.source_face_id))return false;distinct_faces.insert(side.source_face_id);
				const NearestSurfacePoint nearest=bvh.nearest_on_face(probe,side.source_face_id,2.0*grid.h);if(!nearest.found)continue;++nearby_faces;const std::uint8_t probe_side=dot(probe-nearest.point,nearest.geometric_normal)>=0?2:1;if(!(side.side_mask&probe_side))++mismatches;
			}
			return distinct_faces.size()>=2&&nearby_faces>=2&&mismatches>0;
		};
		struct ExactBoundaryPiece{LocalArrangementBoundaryAperture aperture;};
		auto arrangement_and_map=[&](int cell,const LocalSurfaceArrangement*& arrangement,const std::vector<FragmentRef>*& map)->bool
		{
			const EbCellTopology& topology=eb.cells[cell];if(topology.arrangement_index>=0)
			{
				if(topology.arrangement_index>=static_cast<int>(eb.arrangements.size()))return false;arrangement=&eb.arrangements[topology.arrangement_index];map=nullptr;return true;
			}
			// A regular or already-certified analytic cell has no hidden finite topology.
			// Its exact face partition is synthesized directly from the stored fitted plane
			// below; rebuilding all source triangles would reintroduce tessellation-dependent
			// support/limiter planes at an otherwise ordinary shared face.
			if(topology.state==EbCellState::regular||
				(topology.state==EbCellState::split&&!topology.sampled_resolution))
			{
				arrangement=nullptr;map=nullptr;return true;
			}
			return false;
		};
		auto map_local_fragment=[&](int cell,int local,const std::vector<FragmentRef>* map)
		{
			if(map)return local>=0&&local<static_cast<int>(map->size())?(*map)[local]:invalid_fragment;const EbCellTopology& topology=eb.cells[cell];return local>=0&&local<topology.fragment_count?irregular_fragment(topology.first_fragment+local):invalid_fragment;
		};
		auto exact_boundary_pieces=[&](int cell,int axis,bool upper)
		{
			std::vector<ExactBoundaryPiece> result;const LocalSurfaceArrangement* arrangement=nullptr;const std::vector<FragmentRef>* map=nullptr;if(!arrangement_and_map(cell,arrangement,map))return result;
			if(!arrangement)
			{
				const auto coordinate=grid.cell_coord(cell);const Aabb3d box=grid.cell_box(coordinate[0],coordinate[1],coordinate[2]);
				int owner=cell;if(!upper){int oi=coordinate[0]-(axis==0),oj=coordinate[1]-(axis==1),ok=coordinate[2]-(axis==2);if(oi>=0&&oj>=0&&ok>=0)owner=grid.cell_index(oi,oj,ok);}
				if((eb.cut_face_mask[owner]&(1u<<axis))!=0)return result;
				const Polygon face=face_square(box,axis,upper);const EbCellTopology& topology=eb.cells[cell];
				auto append_piece=[&](Polygon polygon,FragmentRef fragment)
				{
					const auto [area,centroid]=planar_polygon_measure(polygon,axis);if(!(area>0.0)||fragment==invalid_fragment)return;
					LocalArrangementBoundaryAperture aperture;aperture.axis=static_cast<std::int8_t>(axis);aperture.upper=upper;aperture.area=area;aperture.centroid=centroid;aperture.polygon=std::move(polygon);aperture.fragment=fragment;aperture.contact_tolerance=bvh.edge_contact_tolerance();result.push_back({std::move(aperture)});
				};
				if(topology.state==EbCellState::regular)append_piece(face,regular_fragment(cell));
				else
				{
					append_piece(clip_polygon_plane(face,{topology.plane_normal,topology.plane_offset},-1,1e-12*grid.h),eb.fragment_for_side(cell,-1));
					append_piece(clip_polygon_plane(face,{topology.plane_normal,topology.plane_offset},1,1e-12*grid.h),eb.fragment_for_side(cell,1));
				}
				return result;
			}
			for(const LocalArrangementBoundaryAperture& source:arrangement->boundary_apertures)if(source.axis==axis&&source.upper==upper){ExactBoundaryPiece piece;piece.aperture=source;piece.aperture.fragment=map_local_fragment(cell,source.fragment,map);result.push_back(std::move(piece));}return result;
		};
		auto exact_boundary_surfaces=[&](int cell,int axis,bool upper)
		{
			std::vector<LocalArrangementBoundarySurfacePatch> result;const LocalSurfaceArrangement* arrangement=nullptr;const std::vector<FragmentRef>* map=nullptr;if(!arrangement_and_map(cell,arrangement,map))return result;
			if(arrangement)
			{
				for(const LocalArrangementBoundarySurfacePatch& source:arrangement->boundary_surface_patches)if(source.axis==axis&&source.upper==upper){LocalArrangementBoundarySurfacePatch patch=source;patch.interior_fragment=map_local_fragment(cell,source.interior_fragment,map);result.push_back(std::move(patch));}return result;
			}
			const auto coordinate=grid.cell_coord(cell);const Aabb3d box=grid.cell_box(coordinate[0],coordinate[1],coordinate[2]);
			int owner=cell;if(!upper){int oi=coordinate[0]-(axis==0),oj=coordinate[1]-(axis==1),ok=coordinate[2]-(axis==2);if(oi<0||oj<0||ok<0)return result;owner=grid.cell_index(oi,oj,ok);}
			for(const SurfacePatch& staged:eb.patches)
			{
				if(staged.provisional_face_cell!=owner||staged.provisional_face_axis!=axis)continue;
				const BvhTriangle triangle=canonical_triangle(staged.source_triangle_id,box);Polygon polygon=clip_triangle_box(triangle,box);if(polygon.size()<3)continue;
				const double coordinate_value=upper?box.hi[axis]:box.lo[axis];bool on_face=true;for(Vec3d point:polygon)on_face=on_face&&std::abs(point[axis]-coordinate_value)<=aligned_snap_tolerance;if(!on_face)continue;
				const auto [area,centroid]=planar_polygon_measure(polygon,axis);if(!(area>0.0))continue;
				Vec3d probe=centroid;probe[axis]+=(upper?-1.0:1.0)*std::max(1e-7*grid.h,2.0*aligned_snap_tolerance);
				const FragmentRef interior=eb.fragment_containing_point(cell,probe);if(interior==invalid_fragment)continue;
				double patch_contact=bvh.edge_contact_tolerance();for(int edge=0;edge<3;++edge)patch_contact=std::max(patch_contact,bvh.edge_contact_tolerance(staged.source_triangle_id,edge));
				LocalArrangementBoundarySurfacePatch patch;patch.source_triangle_id=staged.source_triangle_id;patch.source_face_id=staged.source_face_id;patch.axis=static_cast<std::int8_t>(axis);patch.upper=upper;patch.area=area;patch.centroid=centroid;patch.normal=staged.normal;patch.polygon=std::move(polygon);patch.interior_fragment=interior;patch.interior_is_plus=dot(staged.normal,probe-staged.centroid)>0;patch.contact_tolerance=patch_contact;result.push_back(std::move(patch));
			}
			return result;
		};
		auto exact_boundary_hazards=[&](int cell,int axis,bool upper)
		{
			std::vector<LocalArrangementBoundaryEdgeHazard> result;const LocalSurfaceArrangement* arrangement=nullptr;const std::vector<FragmentRef>* map=nullptr;if(!arrangement_and_map(cell,arrangement,map)||!arrangement)return result;
			for(const LocalArrangementBoundaryEdgeHazard& hazard:arrangement->boundary_edge_hazards)if(hazard.axis==axis&&hazard.upper==upper)result.push_back(hazard);return result;
		};
		auto exact_partition_side=[&](int cell,std::int8_t inward_axis_sign,
			LocalArrangementSharedFaceSide& side,std::string& error)
		{
			side={};side.inward_axis_sign=inward_axis_sign;const LocalSurfaceArrangement* arrangement=nullptr;
			const std::vector<FragmentRef>* map=nullptr;if(!arrangement_and_map(cell,arrangement,map))
			{error="shared-face structural adapter cannot access the cell topology";return false;}
			const EbCellTopology& topology=eb.cells[cell];
			if(!arrangement)
			{
				if(topology.state==EbCellState::regular)
				{
					side.regions.push_back({{},regular_fragment(cell)});return true;
				}
				if(topology.state!=EbCellState::split||topology.sampled_resolution)
				{error="shared-face structural adapter received a non-analytic cell";return false;}
				LocalArrangementSharedFacePartitionPlane plane;
				plane.normal=topology.plane_normal;plane.offset=topology.plane_offset;
				if(topology.plane_support_offset>=0&&topology.plane_support_count>0&&
					static_cast<std::size_t>(topology.plane_support_offset)+
						topology.plane_support_count<=eb.analytic_plane_support_triangles.size())
				{
					const auto begin=eb.analytic_plane_support_triangles.begin()+
						topology.plane_support_offset;
					plane.support_triangles.assign(begin,begin+topology.plane_support_count);
				}
				side.planes.push_back(std::move(plane));
				side.regions.push_back({{-1},eb.fragment_for_side(cell,-1)});
				side.regions.push_back({{1},eb.fragment_for_side(cell,1)});return true;
			}
			// The OCCT per-cell diagnostic represents regions by generalized winding
			// boundaries rather than one shared plane/sign arrangement.  Keep it on the
			// polygon-cover diagnostic overload; production finite-triangle arrangements
			// always retain their structural atoms and support planes here.
			if(arrangement->atoms.empty())return false;
			for(const LocalArrangementPlane& plane:arrangement->planes)
				side.planes.push_back({plane.normal,plane.offset,
					plane.support_source_triangles});
			std::map<std::vector<std::int8_t>,FragmentRef> region_by_sign;
			for(const LocalArrangementAtom& atom:arrangement->atoms)
			{
				if(atom.plane_side.size()!=side.planes.size())
				{error="local arrangement atom has an incomplete structural sign vector";return false;}
				const FragmentRef fragment=map_local_fragment(cell,atom.fragment,map);
				if(fragment==invalid_fragment)
				{error="local arrangement atom has no global pressure fragment";return false;}
				auto [found,inserted]=region_by_sign.emplace(atom.plane_side,fragment);
				if(!inserted&&found->second!=fragment)
				{error="one local structural sign vector maps to multiple pressure fragments";return false;}
			}
			for(const auto& item:region_by_sign)side.regions.push_back({item.first,item.second});
			return !side.regions.empty();
		};
		struct CanonicalSharedFaceGeometry
		{
			bool valid=false;std::string error;
			std::vector<LocalArrangementSharedFaceTrace> traces;
			std::vector<LocalArrangementSharedFaceBarrier> barriers;
			std::vector<CoincidentSurfaceProvenance> coincident_provenance;
		};
		auto canonical_shared_face_geometry=[&](const Aabb3d& lower_cell,int axis,
			int parent_face_cell)
		{
			CanonicalSharedFaceGeometry result;const double coordinate=lower_cell.hi[axis];
			Aabb3d face_box=lower_cell;face_box.lo[axis]=coordinate;face_box.hi[axis]=coordinate;
			Aabb3d query=face_box;const Vec3d query_margin{aligned_snap_tolerance,
				aligned_snap_tolerance,aligned_snap_tolerance};query.lo=query.lo-query_margin;
			query.hi=query.hi+query_margin;std::vector<std::uint32_t> ids=bvh.query_aabb(query);
			std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
			// Plane membership must acknowledge the representable spacing of the world
			// coordinate itself.  Tangential trace merging must not: scaling that test by
			// |world origin| turns a translated face into a micron-scale topology weld.
			// Its bound is therefore expressed solely in the local face length scale.
			const double machine=std::numeric_limits<double>::epsilon();
			const double face_length_roundoff=std::max(
				128.0*machine*grid.h,std::numeric_limits<double>::denorm_min());
			const double coordinate_ulp=std::max(
				std::abs(std::nextafter(coordinate,std::numeric_limits<double>::infinity())-
					coordinate),
				std::abs(coordinate-std::nextafter(coordinate,
					-std::numeric_limits<double>::infinity())));
			const double plane_roundoff=std::max(face_length_roundoff,4.0*coordinate_ulp);
			auto polygons_are_exact_duplicates=[](const Polygon& a,const Polygon& b)
			{
				if(a.size()!=b.size()||a.empty())return false;
				auto same=[](Vec3d p,Vec3d q){return p.x==q.x&&p.y==q.y&&p.z==q.z;};
				for(std::size_t start=0;start<b.size();++start)if(same(a.front(),b[start]))
					for(int direction:{-1,1})
					{
						bool equal=true;for(std::size_t q=0;q<a.size();++q)
						{
							const std::size_t index=static_cast<std::size_t>(
								(static_cast<long long>(start)+direction*static_cast<long long>(q)+
								static_cast<long long>(b.size()))%static_cast<long long>(b.size()));
							equal=equal&&same(a[q],b[index]);
						}
						if(equal)return true;
					}
				return false;
			};
			auto append_unique_point=[&](std::vector<Vec3d>& points,Vec3d point)
			{
				point[axis]=coordinate;for(Vec3d retained:points)
					if(length2(point-retained)<=face_length_roundoff*face_length_roundoff)return;
				points.push_back(point);
			};
			for(std::uint32_t id:ids)
			{
				BvhTriangle triangle=canonical_triangle(id,lower_cell);
				const Vec3d vertex[3]={triangle.a,triangle.b,triangle.c};
				for(int edge=0;edge<3;++edge)
				{
					const FabricEdgeCertificate certificate=bvh.edge_certificate(id,edge);
					if(certificate.kind!=FabricEdgeKind::unknown)continue;
					const double edge_plane_bound=std::max(plane_roundoff,
						certificate.contact_tolerance);Vec3d edge_a=vertex[edge],
						edge_b=vertex[(edge+1)%3];
					if(std::abs(edge_a[axis]-coordinate)>edge_plane_bound||
						std::abs(edge_b[axis]-coordinate)>edge_plane_bound)continue;
					edge_a[axis]=edge_b[axis]=coordinate;Vec3d clipped_a{},clipped_b{};
					if(clip_segment_box(edge_a,edge_b,face_box,clipped_a,clipped_b)&&
						length2(clipped_b-clipped_a)>face_length_roundoff*face_length_roundoff)
					{
						result.error="unknown finite fabric edge lies on the shared Cartesian face; refine or repair CAD topology";
						return result;
					}
				}
				const Vec3d raw_normal=cross(triangle.b-triangle.a,triangle.c-triangle.a);
				if(!(length2(raw_normal)>0.0))continue;const Vec3d normal=normalized(raw_normal);
				bool coplanar=std::abs(std::abs(normal[axis])-1.0)<=opt.coplanar_angle_tolerance;
				for(Vec3d point:vertex)coplanar=coplanar&&
					std::abs(point[axis]-coordinate)<=aligned_snap_tolerance;
				if(coplanar)
				{
					triangle.a[axis]=triangle.b[axis]=triangle.c[axis]=coordinate;
					Polygon polygon=clip_triangle_box(triangle,face_box);
					for(Vec3d& point:polygon)point[axis]=coordinate;
					const auto [area,centroid]=planar_polygon_measure(polygon,axis);
					if(!(area>0.0))continue;
					LocalArrangementSharedFaceBarrier candidate{id,triangle.source_face_id,
						normalized(normal),std::move(polygon)};int duplicate_index=-1;
					for(const auto& retained:result.barriers)
						if(polygons_are_exact_duplicates(candidate.polygon,retained.polygon))
						{duplicate_index=static_cast<int>(&retained-result.barriers.data());break;}
					if(duplicate_index<0)result.barriers.push_back(std::move(candidate));
					else
					{
						const auto& canonical=result.barriers[duplicate_index];
						auto provenance=std::find_if(result.coincident_provenance.begin(),
							result.coincident_provenance.end(),[&](const auto& item)
							{return item.canonical_source_triangle_id==
								canonical.source_triangle_id;});
						if(provenance==result.coincident_provenance.end())
						{
							CoincidentSurfaceProvenance item;item.parent_face_cell=parent_face_cell;
							item.axis=static_cast<std::int8_t>(axis);
							item.canonical_source_triangle_id=canonical.source_triangle_id;
							item.coincident_source_triangle_ids.push_back(
								canonical.source_triangle_id);
							item.coincident_source_face_ids.push_back(canonical.source_face_id);
							result.coincident_provenance.push_back(std::move(item));
							provenance=std::prev(result.coincident_provenance.end());
						}
						provenance->coincident_source_triangle_ids.push_back(id);
						provenance->coincident_source_face_ids.push_back(triangle.source_face_id);
					}
					continue;
				}

				std::vector<Vec3d> intersections;
				for(int edge=0;edge<3;++edge)
				{
					const Vec3d a=vertex[edge],b=vertex[(edge+1)%3];
					const double da=a[axis]-coordinate,db=b[axis]-coordinate;
					if(std::abs(da)<=plane_roundoff)append_unique_point(intersections,a);
					if((da< -plane_roundoff&&db>plane_roundoff)||
						(da>plane_roundoff&&db< -plane_roundoff))
						append_unique_point(intersections,a+(b-a)*(da/(da-db)));
				}
				if(intersections.size()<2)continue;std::size_t far_a=0,far_b=1;
				double far_distance=length2(intersections[1]-intersections[0]);
				for(std::size_t a=0;a<intersections.size();++a)
					for(std::size_t b=a+1;b<intersections.size();++b)
						if(length2(intersections[b]-intersections[a])>far_distance)
						{far_a=a;far_b=b;far_distance=length2(intersections[b]-intersections[a]);}
				Vec3d clipped_a{},clipped_b{};
				if(!clip_segment_box(intersections[far_a],intersections[far_b],face_box,
					clipped_a,clipped_b)||length2(clipped_b-clipped_a)<=
					face_length_roundoff*face_length_roundoff)continue;
				clipped_a[axis]=clipped_b[axis]=coordinate;
				// A certified trace can meet a Cartesian face edge a few representable
				// roundoff units away when its CAD vertex came from trigonometric/NURBS
				// evaluation.  Canonicalize that endpoint coordinate to the face edge so
				// neighbouring face tiles do not acquire a positive arithmetic carrier.
				// This is a coordinate-ULP bound, not an aperture-area cutoff: resolved
				// openings remain untouched.
				for(int tangent_axis=0;tangent_axis<3;++tangent_axis)if(tangent_axis!=axis)
					for(double boundary:{face_box.lo[tangent_axis],face_box.hi[tangent_axis]})
					{
						if(std::abs(clipped_a[tangent_axis]-boundary)<=face_length_roundoff)
							clipped_a[tangent_axis]=boundary;
						if(std::abs(clipped_b[tangent_axis]-boundary)<=face_length_roundoff)
							clipped_b[tangent_axis]=boundary;
					}
				std::uint64_t attachment_id=FabricEdgeCertificate::no_attachment;
				for(int edge=0;edge<3;++edge)
				{
					const FabricEdgeCertificate certificate=bvh.edge_certificate(id,edge);
					if(certificate.kind!=FabricEdgeKind::attached||
						certificate.attachment_id==FabricEdgeCertificate::no_attachment)continue;
					Vec3d edge_a=vertex[edge],edge_b=vertex[(edge+1)%3];
					if(std::abs(edge_a[axis]-coordinate)>plane_roundoff||
						std::abs(edge_b[axis]-coordinate)>plane_roundoff)continue;
					edge_a[axis]=edge_b[axis]=coordinate;Vec3d attached_a{},attached_b{};
					if(!clip_segment_box(edge_a,edge_b,face_box,attached_a,attached_b)||
						length2(attached_b-attached_a)<=
						face_length_roundoff*face_length_roundoff)continue;
					attachment_id=std::min(attachment_id,certificate.attachment_id);
				}
				result.traces.push_back({id,triangle.source_face_id,clipped_a,clipped_b,normal,
					attachment_id});
			}
			result.valid=true;return result;
		};
		auto boundary_surface_union_matches=[&](const Aabb3d& lower_box,int axis,
			const CanonicalSharedFaceGeometry& physical,
			const std::vector<LocalArrangementBoundarySurfacePatch>& local,
			const char* side,std::string& error)
		{
			// This is an audit only.  It neither grows nor erodes polygons and therefore
			// cannot turn a microscopic opening into fabric.  The independently built
			// cell-boundary representation must describe the same positive-area union as
			// the one physical BVH barrier set used by the canonical transaction.
			auto exact_duplicate=[](const Polygon& a,const Polygon& b)
			{
				if(a.size()!=b.size()||a.empty())return false;
				auto same=[](Vec3d p,Vec3d q)
				{return p.x==q.x&&p.y==q.y&&p.z==q.z;};
				for(std::size_t start=0;start<b.size();++start)if(same(a.front(),b[start]))
					for(int direction:{-1,1})
					{
						bool equal=true;for(std::size_t q=0;q<a.size();++q)
						{
							const auto index=static_cast<std::size_t>((static_cast<long long>(start)+
								direction*static_cast<long long>(q)+static_cast<long long>(b.size()))%
								static_cast<long long>(b.size()));
							equal=equal&&same(a[q],b[index]);
						}
						if(equal)return true;
					}
				return false;
			};
			std::map<std::uint32_t,std::uint32_t> canonical_source;
			for(const auto& barrier:physical.barriers)
				canonical_source[barrier.source_triangle_id]=barrier.source_triangle_id;
			for(const auto& group:physical.coincident_provenance)
				for(std::uint32_t source:group.coincident_source_triangle_ids)
					canonical_source[source]=group.canonical_source_triangle_id;
			struct LocalPolygon{std::uint32_t canonical=0;const Polygon* polygon=nullptr;};
			std::vector<LocalPolygon> unique_local;unique_local.reserve(local.size());
			for(const auto& patch:local)
			{
				if(patch.polygon.size()<3)continue;
				const auto found=canonical_source.find(patch.source_triangle_id);
				if(found==canonical_source.end())
				{
					error=std::string(side)+
						" boundary surface has no once-derived physical barrier provenance";
					return false;
				}
				bool duplicate=false;for(const auto& retained:unique_local)
					if(retained.canonical==found->second&&
						exact_duplicate(*retained.polygon,patch.polygon))
					{duplicate=true;break;}
				if(!duplicate)unique_local.push_back({found->second,&patch.polygon});
			}
			struct Accumulation
			{
				long double area=0,moment_x=0,moment_y=0,moment_z=0;
			};
			const Vec3d face_center=(lower_box.lo+lower_box.hi)*0.5+
				Vec3d{axis==0?0.5*grid.h:0.0,axis==1?0.5*grid.h:0.0,
					axis==2?0.5*grid.h:0.0};
			auto append=[&](Accumulation& sum,const Polygon& polygon)
			{
				const auto [area,centroid]=planar_polygon_measure(polygon,axis);
				if(!(area>0.0))return;const long double a=area;sum.area+=a;
				sum.moment_x+=a*static_cast<long double>(centroid.x-face_center.x);
				sum.moment_y+=a*static_cast<long double>(centroid.y-face_center.y);
				sum.moment_z+=a*static_cast<long double>(centroid.z-face_center.z);
			};
			Accumulation physical_sum,local_sum,overlap_sum;
			Vec3d face_normal{};face_normal[axis]=1.0;
			for(const auto& barrier:physical.barriers)append(physical_sum,barrier.polygon);
			for(const auto& patch:unique_local)append(local_sum,*patch.polygon);
			for(const auto& barrier:physical.barriers)for(const auto& patch:unique_local)
				append(overlap_sum,intersect_coplanar_convex(barrier.polygon,
					*patch.polygon,face_normal,0.0));
			const long double operation_count=static_cast<long double>(1+
				physical.barriers.size()+unique_local.size()+
				physical.barriers.size()*unique_local.size());
			const long double face_area=static_cast<long double>(grid.h)*grid.h;
			const long double area_scale=std::max({face_area,physical_sum.area,
				local_sum.area,overlap_sum.area,std::numeric_limits<long double>::min()});
			const long double area_tolerance=256.0L*
				std::numeric_limits<double>::epsilon()*operation_count*area_scale;
			const long double moment_tolerance=area_tolerance*grid.h;
			auto moment_difference=[](const Accumulation& a,const Accumulation& b)
			{
				const long double x=a.moment_x-b.moment_x,y=a.moment_y-b.moment_y,
					z=a.moment_z-b.moment_z;return std::sqrt(x*x+y*y+z*z);
			};
			const bool conserved=
				std::abs(local_sum.area-physical_sum.area)<=area_tolerance&&
				std::abs(overlap_sum.area-physical_sum.area)<=area_tolerance&&
				std::abs(overlap_sum.area-local_sum.area)<=area_tolerance&&
				moment_difference(local_sum,physical_sum)<=moment_tolerance&&
				moment_difference(overlap_sum,physical_sum)<=moment_tolerance&&
				moment_difference(overlap_sum,local_sum)<=moment_tolerance;
			if(!conserved)
			{
				std::ostringstream message;message<<side<<
					" boundary-surface union differs from the once-derived physical barrier union"
					<<" (physical/local/intersection area "<<static_cast<double>(physical_sum.area)
					<<" / "<<static_cast<double>(local_sum.area)<<" / "
					<<static_cast<double>(overlap_sum.area)<<")";error=message.str();return false;
			}
			return true;
		};
		auto mark_exact_face_unresolved=[&](int cell,const std::vector<std::uint32_t>& triangles,const std::string& reason)
		{
			if(cell<0||cell>=static_cast<int>(eb.cells.size())||eb.cells[cell].state==EbCellState::unresolved)return;
			// This list is a BVH query spanning both adjacent cells. It locates the
			// failed transaction but does not identify a causal CAD triangle.
			eb.cells[cell].state=EbCellState::unresolved;
			eb.unresolved.push_back({cell,triangles,reason,true});
		};

		// Split only those Cartesian faces that touch an irregular cell. Pairwise half-space clipping
		// produces separate apertures when a membrane divides a MAC face.
		for(int k=0;k<grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<grid.nx;++i)
		{
			const int ca=grid.cell_index(i,j,k);if(eb.cells[ca].state==EbCellState::unresolved)continue;
			for(int axis=0;axis<3;++axis)
			{
				int ni=i+(axis==0),nj=j+(axis==1),nk=k+(axis==2);if(ni>=grid.nx||nj>=grid.ny||nk>=grid.nz)continue;int cb=grid.cell_index(ni,nj,nk);if(eb.cells[cb].state==EbCellState::unresolved)continue;const bool explicitly_cut=(eb.cut_face_mask[ca]&(1u<<axis))!=0;const bool has_exact_arrangement=eb.cells[ca].arrangement_index>=0||eb.cells[cb].arrangement_index>=0;const bool has_non_sampled_split=(eb.cells[ca].state==EbCellState::split&&!eb.cells[ca].sampled_resolution)||(eb.cells[cb].state==EbCellState::split&&!eb.cells[cb].sampled_resolution);const bool exact_face_required=has_exact_arrangement||has_non_sampled_split||(explicitly_cut&&(eb.cells[ca].state!=EbCellState::regular||eb.cells[cb].state!=EbCellState::regular));if(explicitly_cut&&!exact_face_required)continue;if(eb.cells[ca].state==EbCellState::regular&&eb.cells[cb].state==EbCellState::regular&&!explicitly_cut)continue;
				if(exact_face_required&&!eb.cells[ca].sampled_resolution&&!eb.cells[cb].sampled_resolution)
				{
					auto adapter_valid=[&](int owner)
					{
						const EbCellTopology& topology=eb.cells[owner];
						if(topology.arrangement_index>=0)
							return topology.arrangement_index<static_cast<int>(eb.arrangements.size());
						return topology.state==EbCellState::regular||
							(topology.state==EbCellState::split&&!topology.sampled_resolution);
					};
					bool face_valid=adapter_valid(ca)&&adapter_valid(cb);std::string face_error;
					if(!face_valid)face_error="shared-face exact adapter is invalid";
					std::vector<ExactBoundaryPiece> pieces_a,pieces_b;
					std::vector<LocalArrangementBoundarySurfacePatch> surfaces_a,surfaces_b;
					std::vector<LocalArrangementBoundaryEdgeHazard> hazards_a,hazards_b;
					if(face_valid)
					{
						pieces_a=exact_boundary_pieces(ca,axis,true);
						pieces_b=exact_boundary_pieces(cb,axis,false);
						surfaces_a=exact_boundary_surfaces(ca,axis,true);
						surfaces_b=exact_boundary_surfaces(cb,axis,false);
						hazards_a=exact_boundary_hazards(ca,axis,true);
						hazards_b=exact_boundary_hazards(cb,axis,false);
						if(!hazards_a.empty()||!hazards_b.empty())
						{
							face_valid=false;
							face_error="ambiguous fabric edge lies on a shared Cartesian face; refine or repair geometry rather than infer an opening/attachment";
						}
					}
					LocalArrangementSharedFaceSide structural_a,structural_b;bool use_structural=false;
					if(face_valid)
					{
						std::string structural_error_a,structural_error_b;
						const bool have_a=exact_partition_side(ca,-1,structural_a,structural_error_a);
						const bool have_b=exact_partition_side(cb,1,structural_b,structural_error_b);
						if(!structural_error_a.empty()||!structural_error_b.empty())
						{
							face_valid=false;face_error=!structural_error_a.empty()?
								structural_error_a:structural_error_b;
						}
						else use_structural=have_a&&have_b;
					}
					std::vector<LocalArrangementSharedFaceOwner> owners_a,owners_b;
					auto append_owners=[&](const std::vector<ExactBoundaryPiece>& pieces,
						const std::vector<LocalArrangementBoundarySurfacePatch>& surfaces,
						std::vector<LocalArrangementSharedFaceOwner>& owners,const char* side)
					{
						for(const ExactBoundaryPiece& piece:pieces)
						{
							if(piece.aperture.fragment==invalid_fragment||piece.aperture.polygon.size()<3)
							{
								face_valid=false;face_error=std::string(side)+
									" exact face piece has no pressure-fragment owner";return;
							}
							owners.push_back({piece.aperture.polygon,piece.aperture.fragment});
						}
						for(const auto& surface:surfaces)
						{
							if(surface.interior_fragment==invalid_fragment||surface.polygon.size()<3)
							{
								face_valid=false;face_error=std::string(side)+
									" boundary fabric has no pressure-fragment owner";return;
							}
							owners.push_back({surface.polygon,surface.interior_fragment});
						}
					};
					if(face_valid&&!use_structural)
					{
						append_owners(pieces_a,surfaces_a,owners_a,"lower");
						if(face_valid)append_owners(pieces_b,surfaces_b,owners_b,"upper");
					}
					CanonicalSharedFaceGeometry physical_geometry;
					LocalArrangementSharedFaceAssembly assembly;
					if(face_valid)
					{
						const Aabb3d lower_box=grid.cell_box(i,j,k);
						physical_geometry=canonical_shared_face_geometry(lower_box,axis,ca);
						if(!physical_geometry.valid)
						{
							face_valid=false;face_error=physical_geometry.error;
						}
						else
						{
							if(!boundary_surface_union_matches(lower_box,axis,
								physical_geometry,surfaces_a,"lower",face_error)||
								!boundary_surface_union_matches(lower_box,axis,
								physical_geometry,surfaces_b,"upper",face_error))
							{
								face_valid=false;
							}
							const Polygon square=face_square(lower_box,axis,true);
							if(face_valid&&square.size()!=4)
							{
								face_valid=false;face_error="canonical shared face is not a quadrilateral";
							}
							else if(face_valid)
							{
								const std::array<Vec3d,4> face{square[0],square[1],square[2],square[3]};
								if(use_structural)
									assembly=assemble_canonical_local_shared_face(face,
										static_cast<std::int8_t>(axis),physical_geometry.traces,
										physical_geometry.barriers,
										structural_a,structural_b);
								else assembly=assemble_canonical_local_shared_face(face,
									static_cast<std::int8_t>(axis),physical_geometry.traces,
									physical_geometry.barriers,owners_a,owners_b);
								if(!assembly.valid){face_valid=false;face_error=assembly.error;}
							}
						}
					}
					std::vector<FaceAperture> staged_apertures;
					std::vector<FragmentConnection> staged_connections;
					std::vector<SurfacePatch> staged_boundary_patches;
					if(face_valid)
					{
						staged_apertures.reserve(assembly.apertures.size());
						staged_connections.reserve(assembly.apertures.size());
						staged_boundary_patches.reserve(assembly.blocked_surfaces.size());
						auto arrangement_mask_conflict=[&](FragmentRef a,FragmentRef b,
							std::uint32_t& source_face)
						{
							if(fragment_is_regular(a)||fragment_is_regular(b))return false;
							const EbCellTopology& topology_a=eb.cells[ca];
							const EbCellTopology& topology_b=eb.cells[cb];
							if(topology_a.arrangement_index<0||topology_b.arrangement_index<0)return false;
							const FluidFragment& fragment_a=eb.fragments[irregular_fragment_index(a)];
							for(int q=0;q<fragment_a.surface_side_count;++q)
							{
								const FragmentSurfaceSide& side_a=
									eb.fragment_surface_sides[fragment_a.surface_side_offset+q];
								const std::uint8_t side_b=fragment_surface_side_mask(eb,b,
									side_a.source_face_id);
								if(side_b&&!(side_a.side_mask&side_b))
								{source_face=side_a.source_face_id;return true;}
							}
							return false;
						};
						auto tile_overlaps_physical_trace=[&](const std::vector<Vec3d>& polygon,
							std::uint32_t source_face)
						{
							const int tangent0=(axis+1)%3,tangent1=(axis+2)%3;
							const double roundoff=128.0*std::numeric_limits<double>::epsilon()*grid.h;
							for(const auto& trace:physical_geometry.traces)
							{
								if(trace.source_face_id!=source_face)continue;
								const double dx=trace.b[tangent0]-trace.a[tangent0];
								const double dy=trace.b[tangent1]-trace.a[tangent1];
								const double trace_length=std::hypot(dx,dy);
								if(!(trace_length>0.0))continue;
								for(std::size_t edge=0;edge<polygon.size();++edge)
								{
									const Vec3d p=polygon[edge],q=polygon[(edge+1)%polygon.size()];
									auto cross_to_trace=[&](Vec3d point)
									{
										return dx*(point[tangent1]-trace.a[tangent1])-
											dy*(point[tangent0]-trace.a[tangent0]);
									};
									if(std::abs(cross_to_trace(p))>roundoff*trace_length||
										std::abs(cross_to_trace(q))>roundoff*trace_length)continue;
									const int projection_axis=std::abs(dx)>=std::abs(dy)?tangent0:tangent1;
									const double trace_lo=std::min(trace.a[projection_axis],trace.b[projection_axis]);
									const double trace_hi=std::max(trace.a[projection_axis],trace.b[projection_axis]);
									const double edge_lo=std::min(p[projection_axis],q[projection_axis]);
									const double edge_hi=std::max(p[projection_axis],q[projection_axis]);
									if(std::min(trace_hi,edge_hi)-std::max(trace_lo,edge_lo)>roundoff)
										return true;
								}
							}
							return false;
						};
						for(const auto& source:assembly.apertures)
						{
							const FragmentRef a=static_cast<FragmentRef>(source.fragment_a);
							const FragmentRef b=static_cast<FragmentRef>(source.fragment_b);
							if(a==invalid_fragment||b==invalid_fragment||!(source.area>0.0))
							{
								face_valid=false;face_error=
									"canonical shared-face aperture has an invalid pressure owner or area";break;
							}
							std::uint32_t conflicting_source_face=~std::uint32_t{0};
							if(arrangement_mask_conflict(a,b,conflicting_source_face)&&
								tile_overlaps_physical_trace(source.polygon,conflicting_source_face))
							{
								face_valid=false;face_error=
									"canonical shared-face tile crosses a continuous transverse fabric trace";
								break;
							}
							staged_apertures.push_back({ca,static_cast<std::int8_t>(axis),
								source.area,source.centroid,a,b});
							staged_connections.push_back({a,b,source.area,source.centroid,
								std::max(1e-12,std::sqrt(length2(eb.fragment_centroid(a)-
									eb.fragment_centroid(b)))),static_cast<std::int8_t>(axis)});
						}
						for(const auto& source:assembly.blocked_surfaces)if(face_valid)
						{
							const FragmentRef plus=static_cast<FragmentRef>(source.plus_fragment);
							const FragmentRef minus=static_cast<FragmentRef>(source.minus_fragment);
							if(plus==invalid_fragment||minus==invalid_fragment||plus==minus||
								!(source.area>0.0))
							{
								face_valid=false;face_error=
									"canonical shared-face barrier has invalid independent pressure sides";break;
							}
							SurfacePatch patch;patch.source_triangle_id=source.source_triangle_id;
							patch.source_face_id=source.source_face_id;patch.area=source.area;
							patch.centroid=source.centroid;patch.normal=normalized(source.normal);
							patch.plus_fragment=plus;patch.minus_fragment=minus;
							patch.provisional_face_cell=-1;patch.provisional_face_axis=-1;
							staged_boundary_patches.push_back(patch);
						}
					}
					if(!face_valid)
					{
						if(explicitly_cut)std::erase_if(eb.patches,[&](const SurfacePatch& patch)
						{
							return patch.provisional_face_cell==ca&&patch.provisional_face_axis==axis;
						});
						std::vector<std::uint32_t> face_triangles=bvh.query_aabb(grid.cell_box(i,j,k));
						const auto neighbour_coordinate=grid.cell_coord(cb);
						const std::vector<std::uint32_t> neighbour_triangles=bvh.query_aabb(
							grid.cell_box(neighbour_coordinate[0],neighbour_coordinate[1],
								neighbour_coordinate[2]));
						face_triangles.insert(face_triangles.end(),neighbour_triangles.begin(),
							neighbour_triangles.end());std::sort(face_triangles.begin(),face_triangles.end());
						face_triangles.erase(std::unique(face_triangles.begin(),face_triangles.end()),
							face_triangles.end());
						const std::string reason="exact shared-face transaction rejected: "+face_error;
						mark_exact_face_unresolved(ca,face_triangles,reason);
						mark_exact_face_unresolved(cb,face_triangles,reason);break;
					}
					if(explicitly_cut)std::erase_if(eb.patches,[&](const SurfacePatch& patch)
					{
						return patch.provisional_face_cell==ca&&patch.provisional_face_axis==axis;
					});
					eb.apertures.insert(eb.apertures.end(),staged_apertures.begin(),
						staged_apertures.end());
					eb.connections.insert(eb.connections.end(),staged_connections.begin(),
						staged_connections.end());
					eb.patches.insert(eb.patches.end(),staged_boundary_patches.begin(),
						staged_boundary_patches.end());
					eb.coincident_surface_provenance.insert(
						eb.coincident_surface_provenance.end(),
						physical_geometry.coincident_provenance.begin(),
						physical_geometry.coincident_provenance.end());continue;
				}
				const int sampled_resolution=std::max(static_cast<int>(eb.cells[ca].sampled_resolution),static_cast<int>(eb.cells[cb].sampled_resolution));
				if(sampled_resolution>0)
				{
					if(explicitly_cut)
					{
						std::erase_if(eb.patches,[&](const SurfacePatch& patch){return patch.provisional_face_cell==ca&&patch.provisional_face_axis==axis;});
						const std::vector<std::uint32_t> face_triangles=expanded_query(grid.cell_box(i,j,k));
						const std::string reason="sampled diagnostic topology cannot map an aligned fabric face transactionally; use the exact arrangement path";
						mark_exact_face_unresolved(ca,face_triangles,reason);mark_exact_face_unresolved(cb,face_triangles,reason);break;
					}
					struct ApertureAccumulation{FragmentRef a=invalid_fragment,b=invalid_fragment;double area=0;Vec3d weighted_centroid{};};const int tile_count=sampled_resolution*sampled_resolution;std::vector<int> tile_parent(tile_count,-1);std::vector<FragmentRef> tile_a(tile_count,invalid_fragment),tile_b(tile_count,invalid_fragment);std::vector<Vec3d> tile_centres(tile_count);const double tile_h=grid.h/sampled_resolution,tile_area=tile_h*tile_h;const Aabb3d cell_box=grid.cell_box(i,j,k);
					auto tile_root=[&](int q){while(tile_parent[q]!=q){tile_parent[q]=tile_parent[tile_parent[q]];q=tile_parent[q];}return q;};
					auto join_tiles=[&](int a,int b){a=tile_root(a);b=tile_root(b);if(a!=b)tile_parent[std::max(a,b)]=std::min(a,b);};
					auto face_fragment=[&](int owner,bool upper,Vec3d point,int u,int v)->FragmentRef
					{
						const EbCellTopology& topology=eb.cells[owner];if(topology.state==EbCellState::regular)return regular_fragment(owner);if(topology.arrangement_index>=0){Vec3d probe=point;probe[axis]+=(upper?-1.0:1.0)*(1e-7*grid.h);return eb.fragment_containing_point(owner,probe);}if(topology.sampled_resolution)
						{
							const int r=topology.sampled_resolution;const Aabb3d owner_box=grid.cell_box(grid.cell_coord(owner)[0],grid.cell_coord(owner)[1],grid.cell_coord(owner)[2]);int q[3];for(int d=0;d<3;++d)q[d]=std::clamp(static_cast<int>((point[d]-owner_box.lo[d])/grid.h*r),0,r-1);q[axis]=upper?r-1:0;const int local=(q[2]*r+q[1])*r+q[0];return eb.sampled_voxel_fragments[topology.sampled_voxel_offset+local];
						}
						const int side=dot(topology.plane_normal,point)-topology.plane_offset>=0?1:-1;return eb.fragment_for_side(owner,side);
					};
					for(int v=0;v<sampled_resolution;++v)for(int u=0;u<sampled_resolution;++u)
					{
						const int tile=v*sampled_resolution+u;Vec3d centre{};centre[axis]=cell_box.hi[axis];const int tangent0=(axis+1)%3,tangent1=(axis+2)%3;centre[tangent0]=cell_box.lo[tangent0]+(u+0.5)*tile_h;centre[tangent1]=cell_box.lo[tangent1]+(v+0.5)*tile_h;Vec3d left=centre,right=centre;left[axis]-=0.5*tile_h;right[axis]+=0.5*tile_h;
						const double topology_tolerance=1e-9*grid.h;
						if(bvh.segment_touches_surface(left,right,topology_tolerance))continue;const FragmentRef a=face_fragment(ca,true,centre,u,v),b=face_fragment(cb,false,centre,u,v);if(a==invalid_fragment||b==invalid_fragment||!fragments_are_side_compatible(eb,bvh,centre,ca,a,cb,b))continue;tile_parent[tile]=tile;tile_a[tile]=a;tile_b[tile]=b;tile_centres[tile]=centre;
					}
					for(int v=0;v<sampled_resolution;++v)for(int u=0;u<sampled_resolution;++u){const int tile=v*sampled_resolution+u;if(tile_parent[tile]<0)continue;if(u>0){const int other=tile-1;if(tile_parent[other]>=0&&tile_a[other]==tile_a[tile]&&tile_b[other]==tile_b[tile])join_tiles(tile,other);}if(v>0){const int other=tile-sampled_resolution;if(tile_parent[other]>=0&&tile_a[other]==tile_a[tile]&&tile_b[other]==tile_b[tile])join_tiles(tile,other);}}
					std::map<int,ApertureAccumulation> groups;for(int tile=0;tile<tile_count;++tile)if(tile_parent[tile]>=0){const int root=tile_root(tile);auto& accumulation=groups[root];accumulation.a=tile_a[tile];accumulation.b=tile_b[tile];accumulation.area+=tile_area;accumulation.weighted_centroid=accumulation.weighted_centroid+tile_centres[tile]*tile_area;}
					// A sampled fragment normally has a valid route around one free sheet edge, so
					// a blanket centroid-visibility test would close real vents. At a junction
					// touching two or more CAD faces, however, a centroid-to-aperture leg that
					// intersects fabric is a subcell percolation error. Apply this to sampled-to-
					// sampled junctions as well as sampled-to-regular transitions: a rotated sharp
					// leading/trailing edge can put both adjacent cells into fallback topology.
					for(const auto& item:groups){const FragmentRef a=item.second.a,b=item.second.b;const double area=item.second.area;const Vec3d centroid=item.second.weighted_centroid/area;const bool touches_sampled_junction=sampled_junction_fragment(a,ca)||sampled_junction_fragment(b,cb);const bool crosses_junction=touches_sampled_junction&&(connection_leg_crosses_fabric(a,b,centroid)||bvh.segment_touches_surface(eb.fragment_centroid(a),eb.fragment_centroid(b),1e-9*grid.h));if(crosses_junction){reject_cross_fabric_aperture(area);continue;}eb.apertures.push_back({ca,static_cast<std::int8_t>(axis),area,centroid,a,b});eb.connections.push_back({a,b,area,centroid,std::max(1e-12,std::sqrt(length2(eb.fragment_centroid(a)-eb.fragment_centroid(b)))),static_cast<std::int8_t>(axis)});}
					continue;
				}
				for(FragmentRef a:cell_fragments(eb,ca))for(FragmentRef b:cell_fragments(eb,cb))
				{
					// Two tangent planes fitted independently to the same curved CAD face can
					// overlap in a narrow strip. Pairing every half-space combination would then
					// create a spurious aperture from the minus side of the membrane to its plus
					// side. CAD provenance gives an unambiguous local correspondence without any
					// global inside/outside assumption: reversing the face winding reverses both
					// its normal and side labels, leaving this test invariant.
					Polygon poly=face_square(grid.cell_box(i,j,k),axis,true);
					if(!fragment_is_regular(a)){const FluidFragment& f=eb.fragments[irregular_fragment_index(a)];const auto& c=eb.cells[ca];poly=clip_polygon_plane(poly,{c.plane_normal,c.plane_offset},f.side,1e-11*grid.h);}
					if(!fragment_is_regular(b)){const FluidFragment& f=eb.fragments[irregular_fragment_index(b)];const auto& c=eb.cells[cb];poly=clip_polygon_plane(poly,{c.plane_normal,c.plane_offset},f.side,1e-11*grid.h);}
					auto [area,cent]=planar_polygon_measure(poly,axis);if(area<=1e-14*grid.h*grid.h||!fragments_are_side_compatible(eb,bvh,cent,ca,a,cb,b))continue;const bool closed_corner_roundoff=area<=topology_roundoff_area_bound&&(subresolution_closed_corner_conflict(a,b,cent)||subresolution_closed_corner_conflict(b,a,cent));if(closed_corner_roundoff||connection_leg_crosses_fabric(a,b,cent)){reject_cross_fabric_aperture(area);continue;}FaceAperture ap{ca,static_cast<std::int8_t>(axis),area,cent,a,b};eb.apertures.push_back(ap);eb.connections.push_back({a,b,area,cent,std::max(1e-12,length2(eb.fragment_centroid(a)-eb.fragment_centroid(b))>0?std::sqrt(length2(eb.fragment_centroid(a)-eb.fragment_centroid(b))):grid.h),static_cast<std::int8_t>(axis)});
				}
			}
		}

		// Cell construction and shared Cartesian-face assembly form one topology
		// transaction.  A face can reject after an earlier face has already committed
		// apertures/connections, and a locally rejected cell can still have provisional
		// patches from a peeled aligned face.  Selectively compacting those records here
		// would require remapping every FragmentRef, arrangement/locator index and sampled
		// voxel reference.  There is no useful partial topology once even one cell is
		// unresolved, so roll the whole attempted topology back to diagnostics-only state
		// before stabilization can consume an orphan.  The unresolved records and their
		// source-triangle evidence remain intact for refinement/error reporting.
		if(!eb.unresolved.empty())
		{
			for(EbCellTopology& topology:eb.cells)
			{
				const bool unresolved=topology.state==EbCellState::unresolved;
				topology=EbCellTopology{};
				if(unresolved)topology.state=EbCellState::unresolved;
			}
			std::fill(eb.cut_face_mask.begin(),eb.cut_face_mask.end(),std::uint8_t{0});
			eb.irregular_cells.clear();
			eb.fragments.clear();
			eb.fragment_surface_sides.clear();
			eb.connections.clear();
			eb.apertures.clear();
			eb.patches.clear();
			eb.sampled_voxel_fragments.clear();
			eb.arrangements.clear();
			eb.exact_locators.clear();
			eb.coincident_surface_provenance.clear();
			return eb;
		}

		// Audit, but never delete, small positive-area apertures.  Aperture measure is
		// not a topology certificate: one physical vent may be partitioned into many
		// clipping/common-refinement pieces, each below an arbitrary per-piece cutoff.
		// Removing those pieces silently seals the vent and changes the pressure graph.
		//
		// A future optimization may coalesce pieces only after proving that they are
		// subdivision carriers for the same physical connection.  Such a replacement
		// must preserve both sum(A) and sum(A*x_face).  Until that proof is represented
		// explicitly, retaining every aperture is the conservative finite-volume policy.
		const double minimum_aperture_area=opt.min_aperture_area_fraction*grid.h*grid.h;
		if(minimum_aperture_area>0)
		{
			for(const FaceAperture& aperture:eb.apertures)if(aperture.area<minimum_aperture_area){++eb.retained_subgrid_apertures;eb.retained_subgrid_aperture_area+=aperture.area;}
		}

		// Conservative small-fragment stabilization. Aggregate only local neighbours and
		// stop as soon as the aggregate reaches the requested minimum volume. Uniting an
		// entire connected component of tiny fragments lets a sliver chain percolate along
		// a fabric sheet, creating one non-local control volume spanning many cells.
		const double minimum_volume=opt.min_volume_fraction*cell_volume,maximum_merge_span=2*grid.h;const int fragment_count=static_cast<int>(eb.fragments.size());std::vector<int> small_parent(fragment_count);std::vector<double> aggregate_volume(fragment_count);std::vector<Vec3d> aggregate_lo(fragment_count),aggregate_hi(fragment_count);std::iota(small_parent.begin(),small_parent.end(),0);for(int fragment=0;fragment<fragment_count;++fragment){aggregate_volume[fragment]=eb.fragments[fragment].volume;aggregate_lo[fragment]=aggregate_hi[fragment]=eb.fragments[fragment].centroid;}auto small_root=[&](int q){while(small_parent[q]!=q){small_parent[q]=small_parent[small_parent[q]];q=small_parent[q];}return q;};auto is_small=[&](int fragment){return eb.fragments[fragment].volume<minimum_volume;};auto join_small=[&](int a,int b)
		{
			a=small_root(a);b=small_root(b);if(a==b||aggregate_volume[a]>=minimum_volume||aggregate_volume[b]>=minimum_volume)return;const Vec3d lo{std::min(aggregate_lo[a].x,aggregate_lo[b].x),std::min(aggregate_lo[a].y,aggregate_lo[b].y),std::min(aggregate_lo[a].z,aggregate_lo[b].z)},hi{std::max(aggregate_hi[a].x,aggregate_hi[b].x),std::max(aggregate_hi[a].y,aggregate_hi[b].y),std::max(aggregate_hi[a].z,aggregate_hi[b].z)};if(std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z})>maximum_merge_span+1e-12*grid.h)return;const int keep=std::min(a,b),drop=std::max(a,b);small_parent[drop]=keep;aggregate_volume[keep]+=aggregate_volume[drop];aggregate_lo[keep]=lo;aggregate_hi[keep]=hi;
		};
		for(const FaceAperture& aperture:eb.apertures)if(!fragment_is_regular(aperture.fragment_a)&&!fragment_is_regular(aperture.fragment_b)){const int a=irregular_fragment_index(aperture.fragment_a),b=irregular_fragment_index(aperture.fragment_b);if(is_small(a)&&is_small(b))join_small(a,b);}
		std::map<int,std::vector<int>> small_groups;for(int fragment=0;fragment<fragment_count;++fragment)if(is_small(fragment))small_groups[small_root(fragment)].push_back(fragment);
		for(const auto& item:small_groups)
		{
			const std::vector<int>& members=item.second;double total_volume=0;int aggregate=members.front();for(int fragment:members){total_volume+=eb.fragments[fragment].volume;if(eb.fragments[fragment].volume>eb.fragments[aggregate].volume)aggregate=fragment;}
			if(total_volume>=minimum_volume){for(int fragment:members)if(fragment!=aggregate){eb.fragments[fragment].merge_target=irregular_fragment(aggregate);eb.fragments[fragment].pressure_dof=eb.fragments[aggregate].pressure_dof;}continue;}
			FragmentRef external=invalid_fragment;double external_volume=-1,external_area=-1;for(const FaceAperture& aperture:eb.apertures)
			{
				auto consider=[&](FragmentRef member_ref,FragmentRef other){if(fragment_is_regular(member_ref)||small_root(irregular_fragment_index(member_ref))!=item.first)return;double volume=eb.fragment_volume(other);if(!fragment_is_regular(other)){const int other_fragment=irregular_fragment_index(other),other_root=small_root(other_fragment);if(is_small(other_fragment)&&(other_root==item.first||aggregate_volume[other_root]<minimum_volume))return;if(is_small(other_fragment))volume=aggregate_volume[other_root];}if(volume>external_volume||((std::abs(volume-external_volume)<=1e-14)&&aperture.area>external_area)){external=other;external_volume=volume;external_area=aperture.area;}};consider(aperture.fragment_a,aperture.fragment_b);consider(aperture.fragment_b,aperture.fragment_a);
			}
			const FragmentRef target=external!=invalid_fragment?external:irregular_fragment(aggregate);for(int fragment:members)if(fragment!=aggregate||external!=invalid_fragment){eb.fragments[fragment].merge_target=target;eb.fragments[fragment].pressure_dof=fragment_is_regular(target)?regular_fragment_cell(target):eb.fragments[irregular_fragment_index(target)].pressure_dof;}
			if(external==invalid_fragment)eb.fragments[aggregate].pressure_static=true;
		}

		// Recompute pressure-static status from the FINAL merged aperture graph. The
		// group-local decision above cannot see a later group choosing it as a merge
		// target; retaining that provisional flag can leave an inactive root connected
		// by a real aperture to active fluid. Preserve a static designation only when the
		// final component contains static seeds exclusively and reaches no ordinary cell.
		// A resolved closed cavity contains ordinary-sized roots, so it remains active and
		// receives the pressure gauge added by the composite solver.
		std::vector<unsigned char> provisional_static(fragment_count,0);for(int fragment=0;fragment<fragment_count;++fragment){provisional_static[fragment]=eb.fragments[fragment].pressure_static?1:0;eb.fragments[fragment].pressure_static=false;}
		std::vector<int> component_parent(fragment_count),component_rank(fragment_count,0);
		std::iota(component_parent.begin(),component_parent.end(),0);
		auto component_root=[&](int q){while(component_parent[q]!=q){component_parent[q]=component_parent[component_parent[q]];q=component_parent[q];}return q;};
		auto component_join=[&](int a,int b){a=component_root(a);b=component_root(b);if(a==b)return;if(component_rank[a]<component_rank[b])std::swap(a,b);component_parent[b]=a;if(component_rank[a]==component_rank[b])++component_rank[a];};
		std::vector<unsigned char> resolving_final(fragment_count,0);
		std::function<FragmentRef(FragmentRef)> final_ref=[&](FragmentRef ref)->FragmentRef
		{
			if(fragment_is_regular(ref))return ref;const int fragment=irregular_fragment_index(ref);
			if(fragment<0||fragment>=fragment_count)throw std::runtime_error("invalid final EB merge reference");
			const FragmentRef target=eb.fragments[fragment].merge_target;if(target==ref)return ref;
			if(resolving_final[fragment])throw std::runtime_error("cyclic final EB merge reference");
			resolving_final[fragment]=1;const FragmentRef result=final_ref(target);resolving_final[fragment]=0;return result;
		};
		for(const FaceAperture& aperture:eb.apertures)
		{
			const FragmentRef a=final_ref(aperture.fragment_a),b=final_ref(aperture.fragment_b);
			if(!fragment_is_regular(a)&&!fragment_is_regular(b))component_join(irregular_fragment_index(a),irregular_fragment_index(b));
		}
		std::vector<unsigned char> touches_regular(fragment_count,0),has_static_seed(fragment_count,0),has_active_root(fragment_count,0);
		for(const FaceAperture& aperture:eb.apertures)
		{
			const FragmentRef a=final_ref(aperture.fragment_a),b=final_ref(aperture.fragment_b);
			if(fragment_is_regular(a)&&!fragment_is_regular(b))touches_regular[component_root(irregular_fragment_index(b))]=1;
			if(fragment_is_regular(b)&&!fragment_is_regular(a))touches_regular[component_root(irregular_fragment_index(a))]=1;
		}
		for(int fragment=0;fragment<fragment_count;++fragment)
		{
			const FragmentRef root=final_ref(irregular_fragment(fragment));
			if(!fragment_is_regular(root)&&root==irregular_fragment(fragment)){const int component=component_root(fragment);if(provisional_static[fragment])has_static_seed[component]=1;else has_active_root[component]=1;}
		}
		for(int fragment=0;fragment<fragment_count;++fragment){const FragmentRef root=final_ref(irregular_fragment(fragment));if(!fragment_is_regular(root)&&root==irregular_fragment(fragment)){const int component=component_root(fragment);if(!touches_regular[component]&&has_static_seed[component]&&!has_active_root[component])eb.fragments[fragment].pressure_static=true;}}

		// Compact per-fragment connection ranges, ordered by owning irregular endpoint. Connections
		// remain unique face fluxes; the offsets are a work-list convenience, not duplicate physics.
		std::vector<int> counts(eb.fragments.size(),0);for(const auto& c:eb.connections){if(!fragment_is_regular(c.fragment_a))++counts[irregular_fragment_index(c.fragment_a)];if(!fragment_is_regular(c.fragment_b))++counts[irregular_fragment_index(c.fragment_b)];}
		int off=0;for(std::size_t i=0;i<eb.fragments.size();++i){eb.fragments[i].connection_offset=off;eb.fragments[i].connection_count=counts[i];off+=counts[i];}
		return eb;
	}
}
