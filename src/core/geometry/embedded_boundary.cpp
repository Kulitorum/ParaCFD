#include "core/geometry/embedded_boundary.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>

namespace paracfd::core
{
	Aabb3d UniformEbGrid::cell_box(int i,int j,int k) const { Vec3d lo=origin+Vec3d{i*h,j*h,k*h}; return {lo,lo+Vec3d{h,h,h}}; }
	Vec3d UniformEbGrid::cell_centroid(int c) const { auto q=cell_coord(c); return origin+Vec3d{(q[0]+0.5)*h,(q[1]+0.5)*h,(q[2]+0.5)*h}; }

	namespace
	{
		using Polygon=std::vector<Vec3d>;
		struct Plane { Vec3d n; double d=0; double signed_distance(Vec3d p)const{return dot(n,p)-d;} };
		struct PolyMeasure { double volume=0; Vec3d centroid{}; Polygon cap; };

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
		EmbeddedBoundary eb;eb.grid=grid;eb.min_volume_fraction=opt.min_volume_fraction;eb.cells.resize(grid.cell_count());eb.cut_face_mask.assign(grid.cell_count(),0);const double cell_volume=grid.h*grid.h*grid.h;
		for(int k=0;k<grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<grid.nx;++i)
		{
			const int cell=grid.cell_index(i,j,k);const Aabb3d box=grid.cell_box(i,j,k);std::vector<std::uint32_t> ids=bvh.query_aabb(box);if(ids.empty())continue;
			std::vector<std::pair<std::uint32_t,Polygon>> clipped;Vec3d repn{};double repd=0;bool complex=false;
			for(std::uint32_t id:ids)
			{
				const auto& t=bvh.triangle(id);Polygon p=clip_triangle_box(t,box);Vec3d n=normalized(cross(t.b-t.a,t.c-t.a));auto [area,cent]=polygon_measure(p,n);if(area<=1e-16*grid.h*grid.h)continue;
				if(clipped.empty()){repn=n;repd=dot(n,t.a);}else{double align=dot(repn,n);if(std::abs(std::abs(align)-1.0)>opt.coplanar_angle_tolerance||std::abs(dot(repn,t.a)-repd)>(opt.coplanar_distance_tolerance+1e-8*grid.h))complex=true;}
				clipped.push_back({id,std::move(p)});
			}
			if(clipped.empty())continue;
			if(complex)
			{
				eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,"multiple non-coplanar fabric sheets intersect one cell; refine and rebuild"});continue;
			}
			Plane plane{repn,repd};PolyMeasure minus=clip_box_halfspace(box,plane,-1),plus=clip_box_halfspace(box,plane,1);auto [section_area,section_c]=polygon_measure(minus.cap,repn);double covered=0;for(auto& cp:clipped){Vec3d n=normalized(cross(bvh.triangle(cp.first).b-bvh.triangle(cp.first).a,bvh.triangle(cp.first).c-bvh.triangle(cp.first).a));covered+=polygon_measure(cp.second,n).first;}

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
					if(!(section_area>0)||covered<section_area*(1.0-opt.surface_coverage_tolerance)){eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,"partially covered Cartesian-aligned face requires refinement before aperture construction"});continue;}
					const int neighbour=grid.cell_index(ni,nj,nk);eb.cut_face_mask[cell]|=static_cast<std::uint8_t>(1u<<aligned_axis);
					for(auto& cp:clipped)
					{
						const auto& t=bvh.triangle(cp.first);Vec3d n=normalized(cross(t.b-t.a,t.c-t.a));auto [patch_area,patch_centroid]=polygon_measure(cp.second,n);
						const bool neighbour_is_plus=dot(n,grid.cell_centroid(neighbour)-t.a)>0;
						SurfacePatch patch;patch.source_triangle_id=cp.first;patch.source_face_id=t.source_face_id;patch.area=patch_area;patch.centroid=patch_centroid;patch.normal=n;patch.plus_fragment=regular_fragment(neighbour_is_plus?neighbour:cell);patch.minus_fragment=regular_fragment(neighbour_is_plus?cell:neighbour);eb.patches.push_back(patch);
					}
					continue;
				}
			}
			if(!(section_area>0)||covered<section_area*(1.0-opt.surface_coverage_tolerance))
			{
				eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,ids,"fabric terminates or has a sub-cell opening inside the cell; refine rather than merge opposite sides"});continue;
			}
			EbCellTopology& ct=eb.cells[cell];ct.state=EbCellState::split;ct.first_fragment=static_cast<int>(eb.fragments.size());ct.fragment_count=2;ct.plane_normal=repn;ct.plane_offset=repd;eb.irregular_cells.push_back(cell);
			eb.fragments.push_back({cell,minus.volume,minus.centroid,-1,irregular_fragment(ct.first_fragment),grid.cell_count()+ct.first_fragment,0,0});
			eb.fragments.push_back({cell,plus.volume,plus.centroid,1,irregular_fragment(ct.first_fragment+1),grid.cell_count()+ct.first_fragment+1,0,0});
			for(auto& cp:clipped)
			{
				const auto& t=bvh.triangle(cp.first);Vec3d n=normalized(cross(t.b-t.a,t.c-t.a));auto [area,cent]=polygon_measure(cp.second,n);const bool same=dot(n,repn)>=0;
				SurfacePatch patch;patch.source_triangle_id=cp.first;patch.source_face_id=t.source_face_id;patch.area=area;patch.centroid=cent;patch.normal=n;patch.plus_fragment=eb.fragment_for_side(cell,same?1:-1);patch.minus_fragment=eb.fragment_for_side(cell,same?-1:1);eb.patches.push_back(patch);
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
				for(FragmentRef a:cell_fragments(eb,ca))for(FragmentRef b:cell_fragments(eb,cb))
				{
					Polygon poly=face_square(grid.cell_box(i,j,k),axis,true);
					if(!fragment_is_regular(a)){const FluidFragment& f=eb.fragments[irregular_fragment_index(a)];const auto& c=eb.cells[ca];poly=clip_polygon_plane(poly,{c.plane_normal,c.plane_offset},f.side,1e-11*grid.h);}
					if(!fragment_is_regular(b)){const FluidFragment& f=eb.fragments[irregular_fragment_index(b)];const auto& c=eb.cells[cb];poly=clip_polygon_plane(poly,{c.plane_normal,c.plane_offset},f.side,1e-11*grid.h);}
					auto [area,cent]=planar_polygon_measure(poly,axis);if(area<=1e-14*grid.h*grid.h)continue;FaceAperture ap{ca,static_cast<std::int8_t>(axis),area,cent,a,b};eb.apertures.push_back(ap);eb.connections.push_back({a,b,area,cent,std::max(1e-12,length2(eb.fragment_centroid(a)-eb.fragment_centroid(b))>0?std::sqrt(length2(eb.fragment_centroid(a)-eb.fragment_centroid(b))):grid.h),static_cast<std::int8_t>(axis)});
				}
			}
		}

		// Conservative small-fragment merge target: choose the largest open aperture to a strictly
		// larger control volume. Every candidate is an aperture connection, so crossing fabric is
		// impossible by construction. If no safe target exists, report the cell as unresolved.
		for(int fi=0;fi<static_cast<int>(eb.fragments.size());++fi)
		{
			FluidFragment& f=eb.fragments[fi];if(f.volume>=opt.min_volume_fraction*cell_volume)continue;FragmentRef self=irregular_fragment(fi),best=invalid_fragment;double best_area=-1;
			for(const FaceAperture& a:eb.apertures)
			{
				FragmentRef other=a.fragment_a==self?a.fragment_b:(a.fragment_b==self?a.fragment_a:invalid_fragment);if(other==invalid_fragment||eb.fragment_volume(other)<=f.volume*(1.0+1e-10))continue;if(a.area>best_area){best_area=a.area;best=other;}
			}
			if(best!=invalid_fragment){f.merge_target=best;f.pressure_dof=fragment_is_regular(best)?regular_fragment_cell(best):eb.fragments[irregular_fragment_index(best)].pressure_dof;}
			else eb.unresolved.push_back({f.parent_cell,{},"small fluid fragment has no larger same-side aperture neighbour at the finest permitted level"});
		}

		// Compact per-fragment connection ranges, ordered by owning irregular endpoint. Connections
		// remain unique face fluxes; the offsets are a work-list convenience, not duplicate physics.
		std::vector<int> counts(eb.fragments.size(),0);for(const auto& c:eb.connections){if(!fragment_is_regular(c.fragment_a))++counts[irregular_fragment_index(c.fragment_a)];if(!fragment_is_regular(c.fragment_b))++counts[irregular_fragment_index(c.fragment_b)];}
		int off=0;for(std::size_t i=0;i<eb.fragments.size();++i){eb.fragments[i].connection_offset=off;eb.fragments[i].connection_count=counts[i];off+=counts[i];}
		return eb;
	}
}
