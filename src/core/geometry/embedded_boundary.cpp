#include "core/geometry/embedded_boundary.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

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

		std::vector<FragmentRef> cell_fragments(const EmbeddedBoundary& eb,int cell)
		{
			const auto& c=eb.cells[cell];if(c.state==EbCellState::regular)return {regular_fragment(cell)};std::vector<FragmentRef> r;for(int i=0;i<c.fragment_count;++i)r.push_back(irregular_fragment(c.first_fragment+i));return r;
		}
	}

	FragmentRef EmbeddedBoundary::fragment_for_side(int cell,int side)const
	{
		const auto& c=cells[cell];if(c.state==EbCellState::regular)return regular_fragment(cell);if(c.state!=EbCellState::split)return invalid_fragment;for(int i=0;i<c.fragment_count;++i)if(fragments[c.first_fragment+i].side==side)return irregular_fragment(c.first_fragment+i);return invalid_fragment;
	}
	Vec3d EmbeddedBoundary::fragment_centroid(FragmentRef r)const{return fragment_is_regular(r)?grid.cell_centroid(regular_fragment_cell(r)):fragments[irregular_fragment_index(r)].centroid;}
	double EmbeddedBoundary::fragment_volume(FragmentRef r)const{return fragment_is_regular(r)?grid.h*grid.h*grid.h:fragments[irregular_fragment_index(r)].volume;}

	EmbeddedBoundary build_embedded_boundary(const TriMesh& mesh,const TriangleBvh& bvh,const UniformEbGrid& grid,const EmbeddedBoundaryBuildOptions& opt)
	{
		if(!(grid.h>0)||grid.nx<=0||grid.ny<=0||grid.nz<=0||!(opt.min_volume_fraction>0&&opt.min_volume_fraction<0.5)||!(opt.min_aperture_area_fraction>=0&&opt.min_aperture_area_fraction<0.5))
			throw std::invalid_argument("invalid embedded-boundary grid or stabilization fraction");
		EmbeddedBoundary eb;eb.grid=grid;eb.min_volume_fraction=opt.min_volume_fraction;eb.min_aperture_area_fraction=opt.min_aperture_area_fraction;eb.cells.resize(grid.cell_count());eb.cut_face_mask.assign(grid.cell_count(),0);const double cell_volume=grid.h*grid.h*grid.h;
		for(int k=0;k<grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<grid.nx;++i)
		{
			const int cell=grid.cell_index(i,j,k);const Aabb3d box=grid.cell_box(i,j,k);std::vector<std::uint32_t> ids=bvh.query_aabb(box);if(ids.empty())continue;
			std::vector<ClippedTriangle> clipped;Vec3d reference_normal{};
			for(std::uint32_t id:ids)
			{
				const auto& t=bvh.triangle(id);Polygon p=clip_triangle_box(t,box);Vec3d n=normalized(cross(t.b-t.a,t.c-t.a));auto [area,cent]=polygon_measure(p,n);if(area<=1e-16*grid.h*grid.h)continue;
				if(clipped.empty())reference_normal=n;
				clipped.push_back({id,std::move(p),n,cent,area});
			}
			if(clipped.empty())continue;
			auto sampled_fallback=[&]()->bool
			{
				const int n=opt.complex_subdivisions;if(n<2)return false;const int count=n*n*n;const double micro_h=grid.h/n,micro_volume=micro_h*micro_h*micro_h;
				auto micro_index=[&](int x,int y,int z){return(z*n+y)*n+x;};auto micro_center=[&](int x,int y,int z){return box.lo+Vec3d{(x+0.5)*micro_h,(y+0.5)*micro_h,(z+0.5)*micro_h};};
				std::vector<int> parent(count);std::iota(parent.begin(),parent.end(),0);auto root=[&](int q){while(parent[q]!=q){parent[q]=parent[parent[q]];q=parent[q];}return q;};auto join=[&](int a,int b){a=root(a);b=root(b);if(a!=b)parent[std::max(a,b)]=std::min(a,b);};
				for(int z=0;z<n;++z)for(int y=0;y<n;++y)for(int x=0;x<n;++x)for(int axis=0;axis<3;++axis)
				{
					int q[3]={x,y,z};if(++q[axis]>=n)continue;const Vec3d a=micro_center(x,y,z),b=micro_center(q[0],q[1],q[2]);if(!bvh.intersect_segment(a,b,1e-9,1.0-1e-9).hit)join(micro_index(x,y,z),micro_index(q[0],q[1],q[2]));
				}
				std::map<int,int> component;std::vector<double> component_volume;std::vector<Vec3d> component_centroid;
				for(int z=0;z<n;++z)for(int y=0;y<n;++y)for(int x=0;x<n;++x){const int r=root(micro_index(x,y,z));auto inserted=component.emplace(r,static_cast<int>(component.size()));if(inserted.second){component_volume.push_back(0);component_centroid.push_back({});}const int c=inserted.first->second;component_volume[c]+=micro_volume;component_centroid[c]=component_centroid[c]+micro_center(x,y,z)*micro_volume;}
				EbCellTopology& topology=eb.cells[cell];topology.state=EbCellState::split;topology.first_fragment=static_cast<int>(eb.fragments.size());topology.fragment_count=static_cast<std::uint16_t>(component.size());topology.sampled_voxel_offset=static_cast<int>(eb.sampled_voxel_fragments.size());topology.sampled_resolution=static_cast<std::uint8_t>(n);eb.irregular_cells.push_back(cell);eb.sampled_voxel_fragments.resize(eb.sampled_voxel_fragments.size()+count);
				for(int c=0;c<static_cast<int>(component.size());++c){component_centroid[c]=component_centroid[c]/component_volume[c];const int fragment=static_cast<int>(eb.fragments.size());eb.fragments.push_back({cell,component_volume[c],component_centroid[c],0,irregular_fragment(fragment),grid.cell_count()+fragment,0,0});}
				for(int z=0;z<n;++z)for(int y=0;y<n;++y)for(int x=0;x<n;++x){const int local=micro_index(x,y,z),c=component[root(local)];eb.sampled_voxel_fragments[topology.sampled_voxel_offset+local]=irregular_fragment(topology.first_fragment+c);}
				auto nearest_side_fragment=[&](Vec3d point,Vec3d normal,int side)
				{
					double best=std::numeric_limits<double>::infinity();FragmentRef result=invalid_fragment;for(int z=0;z<n;++z)for(int y=0;y<n;++y)for(int x=0;x<n;++x){Vec3d centre=micro_center(x,y,z);if(side*dot(centre-point,normal)<=1e-12)continue;double distance=length2(centre-point);if(distance<best){best=distance;result=eb.sampled_voxel_fragments[topology.sampled_voxel_offset+micro_index(x,y,z)];}}return result;
				};
				for(const ClippedTriangle& cp:clipped){const BvhTriangle& t=bvh.triangle(cp.id);SurfacePatch patch;patch.source_triangle_id=cp.id;patch.source_face_id=t.source_face_id;patch.area=cp.area;patch.centroid=cp.centroid;patch.normal=cp.normal;patch.plus_fragment=nearest_side_fragment(cp.centroid,cp.normal,1);patch.minus_fragment=nearest_side_fragment(cp.centroid,cp.normal,-1);if(patch.plus_fragment!=invalid_fragment&&patch.minus_fragment!=invalid_fragment)eb.patches.push_back(patch);}
				return true;
			};
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
				if(sampled_fallback())continue;
				const char* reason=angle_violation&&distance_violation?"local fabric exceeds both smooth-sheet angle and distance bounds; refine and rebuild":(angle_violation?"local fabric exceeds smooth-sheet angle bound; refine and rebuild":"local fabric exceeds smooth-sheet distance bound; refine and rebuild");
				eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,reason});continue;
			}
			Plane plane{repn,repd};PolyMeasure minus=clip_box_halfspace(box,plane,-1),plus=clip_box_halfspace(box,plane,1);auto [section_area,section_c]=polygon_measure(minus.cap,repn);double covered=0;for(const auto& cp:clipped)covered+=polygon_measure(cp.polygon,repn).first;

			// A membrane exactly on a Cartesian face blocks that face but does not split either
			// adjacent cell. Assign ownership to the lower-index cell's positive face so the
			// tessellation is emitted once even though closed AABB queries see it from both cells.
			int aligned_axis=-1;double aligned_coordinate=0;const double aligned_tol=std::max(opt.coplanar_distance_tolerance,1e-9*grid.h);
			for(int axis=0;axis<3;++axis)if(std::abs(std::abs(repn[axis])-1.0)<=opt.coplanar_angle_tolerance){aligned_axis=axis;aligned_coordinate=repd/repn[axis];break;}
			if(aligned_axis>=0)
			{
				const bool on_lo=std::abs(aligned_coordinate-box.lo[aligned_axis])<=aligned_tol,on_hi=std::abs(aligned_coordinate-box.hi[aligned_axis])<=aligned_tol;
				if(on_lo||on_hi)
				{
					if(on_lo)continue; // positive-face owner is the adjacent lower cell
					int ni=i+(aligned_axis==0),nj=j+(aligned_axis==1),nk=k+(aligned_axis==2);
					if(ni>=grid.nx||nj>=grid.ny||nk>=grid.nz){eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,"fabric touches the physical domain boundary; increase the aerodynamic domain margin"});continue;}
					const int neighbour=grid.cell_index(ni,nj,nk);eb.cut_face_mask[cell]|=static_cast<std::uint8_t>(1u<<aligned_axis);
					if(!(section_area>0)||covered>section_area*(1.0+opt.surface_coverage_tolerance))
					{
						eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,"overlapping Cartesian-aligned fabric exceeds one face; refine before aperture construction"});continue;
					}
					if(covered<section_area*(1.0-opt.surface_coverage_tolerance))
					{
						// A terminating sheet can cover only part of an otherwise regular Cartesian
						// face. Suppress its implicit full-face flux, then reconstruct each connected
						// opening as a separate sampled aperture. This remains fluid on both sides;
						// only exact BVH segment hits block a face tile.
						const int n=opt.complex_subdivisions;
						if(n<2){eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,"partially covered Cartesian face needs complex_subdivisions >= 2"});continue;}
						const int tile_count=n*n;std::vector<int> parent(tile_count,-1);std::vector<Vec3d> centre(tile_count);const double tile_h=grid.h/n,tile_area=tile_h*tile_h;const int tangent0=(aligned_axis+1)%3,tangent1=(aligned_axis+2)%3;
						auto root=[&](int q){while(parent[q]!=q){parent[q]=parent[parent[q]];q=parent[q];}return q;};auto join=[&](int a,int b){a=root(a);b=root(b);if(a!=b)parent[std::max(a,b)]=std::min(a,b);};
						for(int v=0;v<n;++v)for(int u=0;u<n;++u){const int q=v*n+u;Vec3d c{};c[aligned_axis]=box.hi[aligned_axis];c[tangent0]=box.lo[tangent0]+(u+0.5)*tile_h;c[tangent1]=box.lo[tangent1]+(v+0.5)*tile_h;Vec3d a=c,b=c;a[aligned_axis]-=0.5*tile_h;b[aligned_axis]+=0.5*tile_h;if(!bvh.intersect_segment(a,b,1e-9,1.0-1e-9).hit){parent[q]=q;centre[q]=c;}}
						for(int v=0;v<n;++v)for(int u=0;u<n;++u){const int q=v*n+u;if(parent[q]<0)continue;if(u+1<n&&parent[q+1]>=0)join(q,q+1);if(v+1<n&&parent[q+n]>=0)join(q,q+n);}
						struct Opening{double area=0;Vec3d weighted_centroid{};};std::map<int,Opening> openings;for(int q=0;q<tile_count;++q)if(parent[q]>=0){Opening& opening=openings[root(q)];opening.area+=tile_area;opening.weighted_centroid=opening.weighted_centroid+centre[q]*tile_area;}
						for(const auto& item:openings){const double area=item.second.area;const Vec3d centroid=item.second.weighted_centroid/area;const FragmentRef a=regular_fragment(cell),b=regular_fragment(neighbour);eb.apertures.push_back({cell,static_cast<std::int8_t>(aligned_axis),area,centroid,a,b});eb.connections.push_back({a,b,area,centroid,grid.h,static_cast<std::int8_t>(aligned_axis)});}
					}
					for(const auto& cp:clipped)
					{
						const auto& t=bvh.triangle(cp.id);Vec3d n=cp.normal;auto [patch_area,patch_centroid]=polygon_measure(cp.polygon,n);
						const bool neighbour_is_plus=dot(n,grid.cell_centroid(neighbour)-t.a)>0;
						SurfacePatch patch;patch.source_triangle_id=cp.id;patch.source_face_id=t.source_face_id;patch.area=patch_area;patch.centroid=patch_centroid;patch.normal=n;patch.plus_fragment=regular_fragment(neighbour_is_plus?neighbour:cell);patch.minus_fragment=regular_fragment(neighbour_is_plus?cell:neighbour);eb.patches.push_back(patch);
					}
					continue;
				}
			}
			if(!(section_area>0)||covered<section_area*(1.0-opt.surface_coverage_tolerance)||covered>section_area*(1.0+opt.surface_coverage_tolerance))
			{
				if(sampled_fallback())continue;
				eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,"fabric coverage is not one complete local sheet (termination, overlap, or sub-cell opening); refine rather than merge opposite sides"});continue;
			}
			EbCellTopology& ct=eb.cells[cell];ct.state=EbCellState::split;ct.first_fragment=static_cast<int>(eb.fragments.size());ct.fragment_count=2;ct.plane_normal=repn;ct.plane_offset=repd;eb.irregular_cells.push_back(cell);
			eb.fragments.push_back({cell,minus.volume,minus.centroid,-1,irregular_fragment(ct.first_fragment),grid.cell_count()+ct.first_fragment,0,0});
			eb.fragments.push_back({cell,plus.volume,plus.centroid,1,irregular_fragment(ct.first_fragment+1),grid.cell_count()+ct.first_fragment+1,0,0});
			for(const auto& cp:clipped)
			{
				const auto& t=bvh.triangle(cp.id);Vec3d n=cp.normal;auto [area,cent]=polygon_measure(cp.polygon,n);const bool same=dot(n,repn)>=0;
				SurfacePatch patch;patch.source_triangle_id=cp.id;patch.source_face_id=t.source_face_id;patch.area=area;patch.centroid=cent;patch.normal=n;patch.plus_fragment=eb.fragment_for_side(cell,same?1:-1);patch.minus_fragment=eb.fragment_for_side(cell,same?-1:1);eb.patches.push_back(patch);
			}
		}

		// Split only those Cartesian faces that touch an irregular cell. Pairwise half-space clipping
		// produces separate apertures when a membrane divides a MAC face.
		for(int k=0;k<grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<grid.nx;++i)
		{
			const int ca=grid.cell_index(i,j,k);if(eb.cells[ca].state==EbCellState::unresolved)continue;
			for(int axis=0;axis<3;++axis)
			{
				int ni=i+(axis==0),nj=j+(axis==1),nk=k+(axis==2);if(ni>=grid.nx||nj>=grid.ny||nk>=grid.nz)continue;int cb=grid.cell_index(ni,nj,nk);if(eb.cells[cb].state==EbCellState::unresolved)continue;if(eb.cells[ca].state==EbCellState::regular&&eb.cells[cb].state==EbCellState::regular)continue;
				const int sampled_resolution=std::max<int>(eb.cells[ca].sampled_resolution,eb.cells[cb].sampled_resolution);
				if(sampled_resolution>0)
				{
					struct ApertureAccumulation{double area=0;Vec3d weighted_centroid{};};std::map<std::pair<FragmentRef,FragmentRef>,ApertureAccumulation> groups;const double tile_h=grid.h/sampled_resolution,tile_area=tile_h*tile_h;const Aabb3d cell_box=grid.cell_box(i,j,k);
					auto face_fragment=[&](int owner,bool upper,Vec3d point,int u,int v)->FragmentRef
					{
						const EbCellTopology& topology=eb.cells[owner];if(topology.state==EbCellState::regular)return regular_fragment(owner);if(topology.sampled_resolution)
						{
							const int r=topology.sampled_resolution;const Aabb3d owner_box=grid.cell_box(grid.cell_coord(owner)[0],grid.cell_coord(owner)[1],grid.cell_coord(owner)[2]);int q[3];for(int d=0;d<3;++d)q[d]=std::clamp(static_cast<int>((point[d]-owner_box.lo[d])/grid.h*r),0,r-1);q[axis]=upper?r-1:0;const int local=(q[2]*r+q[1])*r+q[0];return eb.sampled_voxel_fragments[topology.sampled_voxel_offset+local];
						}
						const int side=dot(topology.plane_normal,point)-topology.plane_offset>=0?1:-1;return eb.fragment_for_side(owner,side);
					};
					for(int v=0;v<sampled_resolution;++v)for(int u=0;u<sampled_resolution;++u)
					{
						Vec3d centre{};centre[axis]=cell_box.hi[axis];const int tangent0=(axis+1)%3,tangent1=(axis+2)%3;centre[tangent0]=cell_box.lo[tangent0]+(u+0.5)*tile_h;centre[tangent1]=cell_box.lo[tangent1]+(v+0.5)*tile_h;Vec3d left=centre,right=centre;left[axis]-=0.5*tile_h;right[axis]+=0.5*tile_h;if(bvh.intersect_segment(left,right,1e-9,1.0-1e-9).hit)continue;const FragmentRef a=face_fragment(ca,true,centre,u,v),b=face_fragment(cb,false,centre,u,v);if(a==invalid_fragment||b==invalid_fragment)continue;auto& accumulation=groups[{a,b}];accumulation.area+=tile_area;accumulation.weighted_centroid=accumulation.weighted_centroid+centre*tile_area;
					}
					for(const auto& item:groups){const FragmentRef a=item.first.first,b=item.first.second;const double area=item.second.area;const Vec3d centroid=item.second.weighted_centroid/area;eb.apertures.push_back({ca,static_cast<std::int8_t>(axis),area,centroid,a,b});eb.connections.push_back({a,b,area,centroid,std::max(1e-12,std::sqrt(length2(eb.fragment_centroid(a)-eb.fragment_centroid(b)))),static_cast<std::int8_t>(axis)});}
					continue;
				}
				for(FragmentRef a:cell_fragments(eb,ca))for(FragmentRef b:cell_fragments(eb,cb))
				{
					Polygon poly=face_square(grid.cell_box(i,j,k),axis,true);
					if(!fragment_is_regular(a)){const FluidFragment& f=eb.fragments[irregular_fragment_index(a)];const auto& c=eb.cells[ca];poly=clip_polygon_plane(poly,{c.plane_normal,c.plane_offset},f.side,1e-11*grid.h);}
					if(!fragment_is_regular(b)){const FluidFragment& f=eb.fragments[irregular_fragment_index(b)];const auto& c=eb.cells[cb];poly=clip_polygon_plane(poly,{c.plane_normal,c.plane_offset},f.side,1e-11*grid.h);}
					auto [area,cent]=planar_polygon_measure(poly,axis);if(area<=1e-14*grid.h*grid.h)continue;FaceAperture ap{ca,static_cast<std::int8_t>(axis),area,cent,a,b};eb.apertures.push_back(ap);eb.connections.push_back({a,b,area,cent,std::max(1e-12,length2(eb.fragment_centroid(a)-eb.fragment_centroid(b))>0?std::sqrt(length2(eb.fragment_centroid(a)-eb.fragment_centroid(b))):grid.h),static_cast<std::int8_t>(axis)});
				}
			}
		}

		// Finite-resolution aperture policy. Analytic clipping occasionally produces
		// positive-area slivers many orders below a face that have negligible mass flux,
		// but their independent velocity DOF can dominate CFL. The aperture and its
		// one-to-one pressure connection are removed together. This does not close a
		// resolved STEP opening: the threshold is a configurable fraction of h^2 and the
		// discarded count/area remain visible in preprocessing telemetry.
		const double minimum_aperture_area=opt.min_aperture_area_fraction*grid.h*grid.h;
		if(minimum_aperture_area>0)
		{
			for(const FaceAperture& aperture:eb.apertures)if(aperture.area<minimum_aperture_area){++eb.discarded_subgrid_apertures;eb.discarded_subgrid_aperture_area+=aperture.area;}
			std::erase_if(eb.apertures,[&](const FaceAperture& aperture){return aperture.area<minimum_aperture_area;});
			std::erase_if(eb.connections,[&](const FragmentConnection& connection){return connection.open_area<minimum_aperture_area;});
		}

		// Conservative small-fragment stabilization. First union adjacent tiny fragments
		// through actual open apertures, then merge the group into its best adjacent large
		// control volume. If none exists, aggregate into the group's largest member. Thus
		// equal-sized fragments across cells do not dead-end, and no union can cross fabric.
		const double minimum_volume=opt.min_volume_fraction*cell_volume;const int fragment_count=static_cast<int>(eb.fragments.size());std::vector<int> small_parent(fragment_count);std::iota(small_parent.begin(),small_parent.end(),0);auto small_root=[&](int q){while(small_parent[q]!=q){small_parent[q]=small_parent[small_parent[q]];q=small_parent[q];}return q;};auto join_small=[&](int a,int b){a=small_root(a);b=small_root(b);if(a!=b)small_parent[std::max(a,b)]=std::min(a,b);};auto is_small=[&](int fragment){return eb.fragments[fragment].volume<minimum_volume;};
		for(const FaceAperture& aperture:eb.apertures)if(!fragment_is_regular(aperture.fragment_a)&&!fragment_is_regular(aperture.fragment_b)){const int a=irregular_fragment_index(aperture.fragment_a),b=irregular_fragment_index(aperture.fragment_b);if(is_small(a)&&is_small(b))join_small(a,b);}
		std::map<int,std::vector<int>> small_groups;for(int fragment=0;fragment<fragment_count;++fragment)if(is_small(fragment))small_groups[small_root(fragment)].push_back(fragment);
		for(const auto& item:small_groups)
		{
			const std::vector<int>& members=item.second;double total_volume=0;int aggregate=members.front();for(int fragment:members){total_volume+=eb.fragments[fragment].volume;if(eb.fragments[fragment].volume>eb.fragments[aggregate].volume)aggregate=fragment;}
			FragmentRef external=invalid_fragment;double external_volume=-1,external_area=-1;for(const FaceAperture& aperture:eb.apertures)
			{
				auto consider=[&](FragmentRef member_ref,FragmentRef other){if(fragment_is_regular(member_ref)||small_root(irregular_fragment_index(member_ref))!=item.first)return;if(!fragment_is_regular(other)&&is_small(irregular_fragment_index(other))&&small_root(irregular_fragment_index(other))==item.first)return;const double volume=eb.fragment_volume(other);if(volume>external_volume||((std::abs(volume-external_volume)<=1e-14)&&aperture.area>external_area)){external=other;external_volume=volume;external_area=aperture.area;}};consider(aperture.fragment_a,aperture.fragment_b);consider(aperture.fragment_b,aperture.fragment_a);
			}
			const FragmentRef target=external!=invalid_fragment?external:irregular_fragment(aggregate);for(int fragment:members)if(fragment!=aggregate||external!=invalid_fragment){eb.fragments[fragment].merge_target=target;eb.fragments[fragment].pressure_dof=fragment_is_regular(target)?regular_fragment_cell(target):eb.fragments[irregular_fragment_index(target)].pressure_dof;}
			if(external==invalid_fragment)eb.fragments[aggregate].pressure_static=true;
		}

		// Compact per-fragment connection ranges, ordered by owning irregular endpoint. Connections
		// remain unique face fluxes; the offsets are a work-list convenience, not duplicate physics.
		std::vector<int> counts(eb.fragments.size(),0);for(const auto& c:eb.connections){if(!fragment_is_regular(c.fragment_a))++counts[irregular_fragment_index(c.fragment_a)];if(!fragment_is_regular(c.fragment_b))++counts[irregular_fragment_index(c.fragment_b)];}
		int off=0;for(std::size_t i=0;i<eb.fragments.size();++i){eb.fragments[i].connection_offset=off;eb.fragments[i].connection_count=counts[i];off+=counts[i];}
		return eb;
	}
}
