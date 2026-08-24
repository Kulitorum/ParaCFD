#include "core/geometry/embedded_boundary.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <exception>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>

namespace paracfd::core
{
	namespace
	{
		using Polygon=std::vector<Vec3d>;

		struct Measure
		{
			double value=0.0;
			Vec3d centroid{};
		};

		struct Plane
		{
			Vec3d normal{};
			double offset=0.0;
			double distance(Vec3d point)const{return dot(normal,point)-offset;}
		};

		struct RawPatch
		{
			std::uint32_t triangle=0;
			std::uint32_t face=0;
			Polygon polygon;
			Vec3d normal{};
			double area=0.0;
			Vec3d centroid{};
			bool output_owner=true;
			FragmentRef plus_fragment=invalid_fragment;
		};

		struct PlaneGroup
		{
			Plane plane;
			Vec3d outward_normal{};
			double outward_offset=0.0;
			std::uint32_t source_face=0;
			std::vector<int> patches;
		};

		struct AtomFace
		{
			Polygon polygon;
			Vec3d normal{}; // outward from this convex atom
			int plane=-1;
			int box_face=-1; // 2*axis+upper
		};

		struct Atom
		{
			std::vector<AtomFace> faces;
			Measure measure;
		};

		Measure polyhedron_measure(const Atom& atom);

		enum class AtomSplitResult
		{
			unchanged,
			split,
			suppressed_subresolution,
			invalid
		};

		struct AtomAdjacency
		{
			int a=-1,b=-1,plane=-1;
			Polygon polygon; // follows atom a's outward winding
			bool source_surface=false;
			bool blocked=false;
		};

		struct StagedAperture
		{
			int axis=0;
			bool upper=false;
			double area=0.0;
			Vec3d centroid{};
			FragmentRef fragment=invalid_fragment;
			std::vector<SolidFaceTriangle> triangles;
		};

		struct StagedSplitCell
		{
			std::vector<FluidFragment> fragments;
			ExactEbCellLocator locator;
		};

		struct DisjointSet
		{
			explicit DisjointSet(int count):parent(count),rank(count,0)
			{
				std::iota(parent.begin(),parent.end(),0);
			}
			int root(int value)
			{
				while(parent[value]!=value)
				{
					parent[value]=parent[parent[value]];
					value=parent[value];
				}
				return value;
			}
			void join(int a,int b)
			{
				a=root(a);b=root(b);if(a==b)return;
				if(rank[a]<rank[b])std::swap(a,b);
				parent[b]=a;if(rank[a]==rank[b])++rank[a];
			}
			std::vector<int> parent,rank;
		};

		bool finite(Vec3d value)
		{
			return std::isfinite(value.x)&&std::isfinite(value.y)&&std::isfinite(value.z);
		}

		Vec3d polygon_vector_area(const Polygon& polygon)
		{
			Vec3d twice{};
			for(std::size_t q=0;q<polygon.size();++q)
				twice=twice+cross(polygon[q],polygon[(q+1)%polygon.size()]);
			return twice*0.5;
		}

		Measure polygon_measure(const Polygon& polygon,Vec3d normal)
		{
			Measure result;if(polygon.size()<3)return result;
			normal=normalized(normal);const Vec3d origin=polygon.front();
			for(std::size_t q=1;q+1<polygon.size();++q)
			{
				const double signed_area=0.5*dot(cross(polygon[q]-origin,
					polygon[q+1]-origin),normal);
				const double area=std::abs(signed_area);
				if(!(area>0.0))continue;
				result.value+=area;
				result.centroid=result.centroid+(origin+polygon[q]+polygon[q+1])*(area/3.0);
			}
			if(result.value>0.0)result.centroid=result.centroid/result.value;
			return result;
		}

		Polygon clip_polygon_axis(const Polygon& input,int axis,double coordinate,
			bool keep_greater,double tolerance)
		{
			Polygon output;if(input.empty())return output;
			auto distance=[&](Vec3d point){return keep_greater?
				point[axis]-coordinate:coordinate-point[axis];};
			Vec3d previous=input.back();double previous_distance=distance(previous);
			bool previous_inside=previous_distance>=-tolerance;
			for(Vec3d current:input)
			{
				const double current_distance=distance(current);
				const bool current_inside=current_distance>=-tolerance;
				if(current_inside!=previous_inside)
				{
					const double denominator=previous_distance-current_distance;
					if(std::abs(denominator)>std::numeric_limits<double>::min())
					{
						Vec3d crossing=previous+(current-previous)*
							(previous_distance/denominator);
						crossing[axis]=coordinate;output.push_back(crossing);
					}
				}
				if(current_inside)
				{
					if(std::abs(current[axis]-coordinate)<=tolerance)current[axis]=coordinate;
					output.push_back(current);
				}
				previous=current;previous_distance=current_distance;
				previous_inside=current_inside;
			}
			return output;
		}

		Polygon clip_triangle_box(const BvhTriangle& triangle,const Aabb3d& box,
			double tolerance)
		{
			Polygon polygon{triangle.a,triangle.b,triangle.c};
			for(int axis=0;axis<3&&!polygon.empty();++axis)
			{
				polygon=clip_polygon_axis(polygon,axis,box.lo[axis],true,tolerance);
				polygon=clip_polygon_axis(polygon,axis,box.hi[axis],false,tolerance);
			}
			return polygon;
		}

		Polygon clip_polygon_halfspace(const Polygon& input,const Plane& plane,
			bool keep_minus,double tolerance)
		{
			Polygon output;if(input.empty())return output;
			auto signed_distance=[&](Vec3d point)
			{
				const double value=plane.distance(point);return keep_minus?-value:value;
			};
			Vec3d previous=input.back();double previous_distance=signed_distance(previous);
			bool previous_inside=previous_distance>=-tolerance;
			for(Vec3d current:input)
			{
				const double current_distance=signed_distance(current);
				const bool current_inside=current_distance>=-tolerance;
				if(current_inside!=previous_inside)
				{
					const double denominator=previous_distance-current_distance;
					if(std::abs(denominator)>std::numeric_limits<double>::min())
					{
						Vec3d crossing=previous+(current-previous)*
							(previous_distance/denominator);
						crossing=crossing-plane.normal*plane.distance(crossing);
						output.push_back(crossing);
					}
				}
				if(current_inside)
				{
					if(std::abs(plane.distance(current))<=tolerance)
						current=current-plane.normal*plane.distance(current);
					output.push_back(current);
				}
				previous=current;previous_distance=current_distance;
				previous_inside=current_inside;
			}
			return output;
		}

		void append_unique(std::vector<Vec3d>& points,Vec3d point,double tolerance)
		{
			const double tolerance2=tolerance*tolerance;
			for(Vec3d existing:points)if(length2(existing-point)<=tolerance2)return;
			points.push_back(point);
		}

		Polygon sorted_cap(std::vector<Vec3d> points,const Plane& plane,Vec3d desired_normal,
			double tolerance)
		{
			std::vector<Vec3d> unique;for(Vec3d point:points)append_unique(unique,
				point-plane.normal*plane.distance(point),tolerance);
			if(unique.size()<3)return {};
			Vec3d centre{};for(Vec3d point:unique)centre=centre+point;
			centre=centre/static_cast<double>(unique.size());
			Vec3d reference=std::abs(desired_normal.x)<0.8?Vec3d{1,0,0}:Vec3d{0,1,0};
			const Vec3d u=normalized(cross(reference,desired_normal));
			const Vec3d v=cross(desired_normal,u);
			std::sort(unique.begin(),unique.end(),[&](Vec3d a,Vec3d b)
			{
				const Vec3d da=a-centre,db=b-centre;
				return std::atan2(dot(da,v),dot(da,u))<std::atan2(dot(db,v),dot(db,u));
			});
			if(dot(polygon_vector_area(unique),desired_normal)<0.0)
				std::reverse(unique.begin(),unique.end());
			return unique;
		}

		AtomSplitResult split_atom(const Atom& source,const Plane& plane,int plane_index,
			double length_tolerance,double area_tolerance,double volume_tolerance,
			Atom& minus,Atom& plus,Polygon& common_cap,double& suppressed_volume)
		{
			suppressed_volume=0.0;
			double minimum=std::numeric_limits<double>::infinity();
			double maximum=-std::numeric_limits<double>::infinity();
			for(const AtomFace& face:source.faces)for(Vec3d point:face.polygon)
			{
				const double value=plane.distance(point);
				minimum=std::min(minimum,value);maximum=std::max(maximum,value);
			}
			if(!(minimum<-length_tolerance&&maximum>length_tolerance))
				return AtomSplitResult::unchanged;

			std::vector<Vec3d> crossings;
			auto build_side=[&](bool keep_minus,Atom& output)
			{
				for(const AtomFace& face:source.faces)
				{
					for(std::size_t q=0;q<face.polygon.size();++q)
					{
						const Vec3d a=face.polygon[q],b=face.polygon[(q+1)%face.polygon.size()];
						const double da=plane.distance(a),db=plane.distance(b);
						if(std::abs(da)<=length_tolerance)append_unique(crossings,a,length_tolerance);
						if(da*db<0.0)
						{
							const Vec3d crossing=a+(b-a)*(da/(da-db));
							append_unique(crossings,crossing,length_tolerance);
						}
					}
					Polygon clipped=clip_polygon_halfspace(face.polygon,plane,keep_minus,
						length_tolerance);
					const Measure measure=polygon_measure(clipped,face.normal);
					if(measure.value>area_tolerance)
						output.faces.push_back({std::move(clipped),face.normal,face.plane,
							face.box_face});
				}
			};
			build_side(true,minus);build_side(false,plus);
			Polygon minus_cap=sorted_cap(crossings,plane,plane.normal,length_tolerance);
			if(polygon_measure(minus_cap,plane.normal).value<=area_tolerance)
				return AtomSplitResult::unchanged;
			common_cap=minus_cap;
			Polygon plus_cap=minus_cap;std::reverse(plus_cap.begin(),plus_cap.end());
			minus.faces.push_back({std::move(minus_cap),plane.normal,plane_index,-1});
			plus.faces.push_back({std::move(plus_cap),plane.normal*-1.0,plane_index,-1});
			minus.measure=polyhedron_measure(minus);plus.measure=polyhedron_measure(plus);
			const bool minus_valid=minus.measure.value>volume_tolerance&&
				finite(minus.measure.centroid);
			const bool plus_valid=plus.measure.value>volume_tolerance&&
				finite(plus.measure.centroid);
			if(minus_valid&&plus_valid)return AtomSplitResult::split;
			if(minus_valid==plus_valid)return AtomSplitResult::invalid;
			const Measure parent=polyhedron_measure(source);
			const double retained=minus_valid?minus.measure.value:plus.measure.value;
			const double discarded=minus_valid?plus.measure.value:minus.measure.value;
			const double closure=std::abs(parent.value-retained-std::max(0.0,discarded));
			if(parent.value>volume_tolerance&&finite(parent.centroid)&&
				std::isfinite(discarded)&&discarded>=0.0&&discarded<=volume_tolerance&&
				closure<=std::max(64.0*volume_tolerance,1.0e-10*parent.value))
			{
				suppressed_volume=std::max(0.0,discarded);
				return AtomSplitResult::suppressed_subresolution;
			}
			return AtomSplitResult::invalid;
		}

		Measure polyhedron_measure(const Atom& atom)
		{
			Measure result;std::vector<Vec3d> vertices;
			for(const AtomFace& face:atom.faces)for(Vec3d point:face.polygon)
				append_unique(vertices,point,1.0e-14);
			if(vertices.size()<4)return result;
			Vec3d reference{};for(Vec3d point:vertices)reference=reference+point;
			reference=reference/static_cast<double>(vertices.size());
			for(const AtomFace& face:atom.faces)
				for(std::size_t q=1;q+1<face.polygon.size();++q)
				{
					const Vec3d a=face.polygon[0],b=face.polygon[q],c=face.polygon[q+1];
					const double volume=std::abs(dot(a-reference,cross(b-reference,c-reference)))/6.0;
					if(!(volume>0.0))continue;
					result.value+=volume;
					result.centroid=result.centroid+(reference+a+b+c)*(volume/4.0);
				}
			if(result.value>0.0)result.centroid=result.centroid/result.value;
			return result;
		}

		Atom cube_atom(const Aabb3d& box)
		{
			const double x0=box.lo.x,x1=box.hi.x,y0=box.lo.y,y1=box.hi.y,
				z0=box.lo.z,z1=box.hi.z;
			Atom atom;
			atom.faces={
				{{{x0,y0,z0},{x0,y0,z1},{x0,y1,z1},{x0,y1,z0}},{-1,0,0},-1,0},
				{{{x1,y0,z0},{x1,y1,z0},{x1,y1,z1},{x1,y0,z1}},{1,0,0},-1,1},
				{{{x0,y0,z0},{x1,y0,z0},{x1,y0,z1},{x0,y0,z1}},{0,-1,0},-1,2},
				{{{x0,y1,z0},{x0,y1,z1},{x1,y1,z1},{x1,y1,z0}},{0,1,0},-1,3},
				{{{x0,y0,z0},{x0,y1,z0},{x1,y1,z0},{x1,y0,z0}},{0,0,-1},-1,4},
				{{{x0,y0,z1},{x1,y0,z1},{x1,y1,z1},{x0,y1,z1}},{0,0,1},-1,5}};
			return atom;
		}

		Polygon intersect_coplanar_convex(const Polygon& subject,const Polygon& clip,
			Vec3d normal,double tolerance)
		{
			if(subject.size()<3||clip.size()<3)return {};
			normal=normalized(normal);Polygon result=subject;
			const double orientation=dot(polygon_vector_area(clip),normal)>=0.0?1.0:-1.0;
			for(std::size_t edge=0;edge<clip.size()&&!result.empty();++edge)
			{
				const Vec3d a=clip[edge],b=clip[(edge+1)%clip.size()];
				Vec3d inward=cross(normal,b-a)*orientation;
				const double magnitude=std::sqrt(length2(inward));if(!(magnitude>0.0))continue;
				inward=inward/magnitude;
				auto side=[&](Vec3d point){return dot(inward,point-a);};
				Polygon output;Vec3d previous=result.back();double previous_side=side(previous);
				bool previous_inside=previous_side>=-tolerance;
				for(Vec3d current:result)
				{
					const double current_side=side(current);
					const bool current_inside=current_side>=-tolerance;
					if(current_inside!=previous_inside)
					{
						const double denominator=previous_side-current_side;
						if(std::abs(denominator)>std::numeric_limits<double>::min())
							output.push_back(previous+(current-previous)*
								(previous_side/denominator));
					}
					if(current_inside)output.push_back(current);
					previous=current;previous_side=current_side;
					previous_inside=current_inside;
				}
				result=std::move(output);
			}
			return result;
		}

		Plane canonical_plane(Vec3d normal,double offset)
		{
			int dominant=0;if(std::abs(normal.y)>std::abs(normal.x))dominant=1;
			if(std::abs(normal.z)>std::abs(normal[dominant]))dominant=2;
			if(normal[dominant]<0.0){normal=normal*-1.0;offset=-offset;}
			return {normal,offset};
		}

		std::vector<SolidFaceTriangle> triangulate(const Polygon& polygon,Vec3d normal)
		{
			std::vector<SolidFaceTriangle> result;if(polygon.size()<3)return result;
			Polygon oriented=polygon;
			if(dot(polygon_vector_area(oriented),normal)<0.0)
				std::reverse(oriented.begin(),oriented.end());
			for(std::size_t q=1;q+1<oriented.size();++q)
				result.push_back({oriented[0],oriented[q],oriented[q+1]});
			return result;
		}

		int neighbour_cell(const UniformEbGrid& grid,int cell,int axis,bool upper)
		{
			auto coordinate=grid.cell_coord(cell);coordinate[axis]+=upper?1:-1;
			if(coordinate[0]<0||coordinate[0]>=grid.nx||coordinate[1]<0||
				coordinate[1]>=grid.ny||coordinate[2]<0||coordinate[2]>=grid.nz)return -1;
			return grid.cell_index(coordinate[0],coordinate[1],coordinate[2]);
		}

		int containing_cell(const UniformEbGrid& grid,Vec3d point)
		{
			int coordinate[3];const int count[3]={grid.nx,grid.ny,grid.nz};
			for(int axis=0;axis<3;++axis)
			{
				coordinate[axis]=static_cast<int>(std::floor((point[axis]-grid.origin[axis])/grid.h));
				if(coordinate[axis]<0||coordinate[axis]>=count[axis])return -1;
			}
			return grid.cell_index(coordinate[0],coordinate[1],coordinate[2]);
		}

		template<class Work>
		unsigned parallel_for_indices(int count,Work&& work)
		{
			if(count<=0)return 0;
			const unsigned hardware=std::max(1u,std::thread::hardware_concurrency());
			unsigned requested=hardware;
			if(const char* configured=std::getenv("PARACFD_EB_THREADS"))
			{
				char* end=nullptr;const unsigned long parsed=std::strtoul(configured,&end,10);
				if(end!=configured&&*end=='\0'&&parsed>0)
					requested=static_cast<unsigned>(std::min<unsigned long>(parsed,hardware));
			}
			const unsigned workers=std::min<unsigned>(requested,
				std::max(1,(count+255)/256));
			if(workers==1)
			{
				for(int index=0;index<count;++index)work(index);
				return 1;
			}
			constexpr int grain=8;
			std::atomic<int> next{0};
			std::atomic<bool> stopped{false};
			std::exception_ptr failure;
			std::mutex failure_mutex;
			std::vector<std::thread> threads;threads.reserve(workers);
			for(unsigned worker=0;worker<workers;++worker)threads.emplace_back([&]
			{
				try
				{
					while(!stopped.load(std::memory_order_relaxed))
					{
						const int begin=next.fetch_add(grain,std::memory_order_relaxed);
						if(begin>=count)break;
						for(int index=begin;index<std::min(count,begin+grain);++index)
							work(index);
					}
				}
				catch(...)
				{
					stopped.store(true,std::memory_order_relaxed);
					const std::lock_guard<std::mutex> lock(failure_mutex);
					if(!failure)failure=std::current_exception();
				}
			});
			for(std::thread& thread:threads)thread.join();
			if(failure)std::rethrow_exception(failure);
			return workers;
		}

		template<class Work>
		unsigned parallel_for_ranges(int count,Work&& work)
		{
			if(count<=0)return 0;
			const unsigned hardware=std::max(1u,std::thread::hardware_concurrency());
			unsigned workers=hardware;
			if(const char* configured=std::getenv("PARACFD_EB_THREADS"))
			{
				char* end=nullptr;const unsigned long parsed=std::strtoul(configured,&end,10);
				if(end!=configured&&*end=='\0'&&parsed>0)
					workers=static_cast<unsigned>(std::min<unsigned long>(parsed,hardware));
			}
			workers=std::min(workers,static_cast<unsigned>(count));
			std::atomic<bool> stopped{false};std::exception_ptr failure;std::mutex failure_mutex;
			auto run=[&](unsigned worker)
			{
				try
				{
					const int begin=static_cast<int>((static_cast<long long>(count)*worker)/workers);
					const int end=static_cast<int>((static_cast<long long>(count)*(worker+1))/workers);
					if(!stopped.load(std::memory_order_relaxed))work(begin,end,worker);
				}
				catch(...)
				{
					stopped.store(true,std::memory_order_relaxed);
					const std::lock_guard<std::mutex> lock(failure_mutex);if(!failure)failure=std::current_exception();
				}
			};
			std::vector<std::thread> threads;threads.reserve(workers-1);
			for(unsigned worker=1;worker<workers;++worker)threads.emplace_back(run,worker);
			run(0);for(std::thread& thread:threads)thread.join();if(failure)std::rethrow_exception(failure);
			return workers;
		}
	}

	EmbeddedBoundary build_closed_solid_embedded_boundary(const TriMesh& display_mesh,
		const TriangleBvh& broad_phase,const ClosedSolidGeometry& solid,
		const UniformEbGrid& grid,const EmbeddedBoundaryBuildOptions& options)
	{
		if(display_mesh.empty()||broad_phase.empty())
			throw std::invalid_argument("closed-solid EB requires a closed triangle mesh");
		if(grid.nx<=0||grid.ny<=0||grid.nz<=0||!(grid.h>0.0))
			throw std::invalid_argument("closed-solid EB requires a valid Cartesian grid");
		if(!solid.bounds().valid())
			throw std::invalid_argument("closed-solid EB requires validated solid bounds");

		EmbeddedBoundary eb;eb.grid=grid;
		eb.cells.resize(static_cast<std::size_t>(grid.cell_count()));
		eb.cut_face_mask.assign(static_cast<std::size_t>(grid.cell_count()),0);
		eb.min_volume_fraction=options.min_volume_fraction;
		eb.min_aperture_area_fraction=options.min_aperture_area_fraction;
		std::vector<std::vector<StagedAperture>> cell_apertures(
			static_cast<std::size_t>(grid.cell_count()));
		std::vector<std::vector<RawPatch>> cell_patches(
			static_cast<std::size_t>(grid.cell_count()));
		std::vector<std::unique_ptr<StagedSplitCell>> split_cells(
			static_cast<std::size_t>(grid.cell_count()));
		std::vector<unsigned char> surface_cell(static_cast<std::size_t>(grid.cell_count()),0);
		std::vector<unsigned char> known_exterior(static_cast<std::size_t>(grid.cell_count()),0);
		std::vector<unsigned char> reported(static_cast<std::size_t>(grid.cell_count()),0);
		std::mutex report_mutex,reconciliation_mutex;

		auto report=[&](int cell,std::string reason,std::vector<std::uint32_t> candidates={})
		{
			if(cell<0||cell>=grid.cell_count()||reported[cell])return;
			reported[cell]=1;eb.cells[cell].state=EbCellState::unresolved;
			const std::lock_guard<std::mutex> lock(report_mutex);
			const auto coordinate=grid.cell_coord(cell);
			std::fprintf(stderr,"[closed-solid-eb] unresolved cell [%d,%d,%d]: %s\n",
				coordinate[0],coordinate[1],coordinate[2],reason.c_str());
			eb.unresolved.push_back({cell,std::move(candidates),std::move(reason),true});
		};

		const double length_tolerance=std::max(1.0e-11*grid.h,
			128.0*std::numeric_limits<double>::epsilon()*std::max(1.0,grid.h));
		const double area_tolerance=std::max(1.0e-12*grid.h*grid.h,
			length_tolerance*length_tolerance);
		const double volume_tolerance=std::max(1.0e-14*grid.h*grid.h*grid.h,
			length_tolerance*length_tolerance*length_tolerance);

		const auto cell_build_begin=std::chrono::steady_clock::now();
		const unsigned cell_workers=parallel_for_indices(grid.cell_count(),[&](int cell)
		{
			const auto coordinate=grid.cell_coord(cell);
			const Aabb3d box=grid.cell_box(coordinate[0],coordinate[1],coordinate[2]);
			Aabb3d query=box;const Vec3d padding{length_tolerance,length_tolerance,length_tolerance};
			query.lo=query.lo-padding;query.hi=query.hi+padding;
			std::vector<std::uint32_t> candidates=broad_phase.query_aabb(query);
			if(candidates.empty())return;

			std::vector<RawPatch>& patches=cell_patches[cell];
			for(std::uint32_t triangle_id:candidates)
			{
				const BvhTriangle& triangle=broad_phase.triangle(triangle_id);
				const Vec3d vector_normal=cross(triangle.b-triangle.a,triangle.c-triangle.a);
				if(!(length2(vector_normal)>0.0))continue;
				const Vec3d normal=normalized(vector_normal);
				Polygon polygon=clip_triangle_box(triangle,box,length_tolerance);
				const Measure measure=polygon_measure(polygon,normal);
				if(!(measure.value>area_tolerance))continue;
				bool output_owner=true;
				for(int axis=0;axis<3;++axis)
				{
					bool on_lower=true;for(Vec3d point:polygon)
						on_lower=on_lower&&std::abs(point[axis]-box.lo[axis])<=length_tolerance;
					if(on_lower&&grid.cell_coord(cell)[axis]>0)output_owner=false;
				}
				patches.push_back({triangle_id,triangle.source_face_id,std::move(polygon),normal,
					measure.value,measure.centroid,output_owner});
			}
			if(patches.empty())return;
			surface_cell[cell]=1;

			// OCCT's curvature-driven tessellation can put hundreds of slightly
			// different triangle planes from one smooth NURBS face in a single CFD
			// cell.  Their infinite-plane arrangement is both unnecessary and
			// combinatorial. Represent each locally intersecting CAD face by its
			// area-weighted tangent plane; the clipped source triangles remain the
			// authoritative load surface.
			std::map<std::uint32_t,std::vector<int>> patches_by_face;
			for(int patch=0;patch<static_cast<int>(patches.size());++patch)
				patches_by_face[patches[patch].face].push_back(patch);
			std::vector<PlaneGroup> planes;planes.reserve(patches_by_face.size());
			for(const auto& [source_face,face_patches]:patches_by_face)
			{
				double area=0.0;Vec3d normal_sum{},centroid_sum{};
				for(int patch:face_patches)
				{
					area+=patches[patch].area;
					normal_sum=normal_sum+patches[patch].normal*patches[patch].area;
					centroid_sum=centroid_sum+patches[patch].centroid*patches[patch].area;
				}
				if(!(area>area_tolerance)||!(length2(normal_sum)>0.0))continue;
				const Vec3d outward=normalized(normal_sum);
				const Vec3d centroid=centroid_sum/area;
				const double outward_offset=dot(outward,centroid);
				planes.push_back({canonical_plane(outward,outward_offset),outward,
					outward_offset,source_face,face_patches});
			}
			if(planes.size()>16)
			{
				report(cell,"closed triangle cell intersects more than 16 distinct CAD faces",
					std::move(candidates));return;
			}

			std::vector<Atom> atoms{cube_atom(box)};
			std::vector<AtomAdjacency> adjacencies;
			bool atom_overflow=false;
			for(int plane=0;plane<static_cast<int>(planes.size())&&!atom_overflow;++plane)
			{
				std::vector<Atom> next;next.reserve(atoms.size()*2);
				std::vector<std::array<int,2>> children(atoms.size(),{-1,-1});
				std::vector<AtomAdjacency> new_caps;
				for(int atom=0;atom<static_cast<int>(atoms.size());++atom)
				{
					Atom minus,plus;Polygon common_cap;double suppressed_volume=0.0;
					const AtomSplitResult split=split_atom(atoms[atom],planes[plane].plane,
						plane,length_tolerance,area_tolerance,volume_tolerance,minus,plus,
						common_cap,suppressed_volume);
					if(split==AtomSplitResult::split)
					{
						const int minus_index=static_cast<int>(next.size());
						next.push_back(std::move(minus));
						const int plus_index=static_cast<int>(next.size());
						next.push_back(std::move(plus));
						children[atom]={minus_index,plus_index};
						new_caps.push_back({minus_index,plus_index,plane,
							std::move(common_cap),false,false});
					}
					else
					{
						if(split==AtomSplitResult::invalid)
						{
							report(cell,"convex atom plane split failed its transactional volume audit",
								candidates);atom_overflow=true;break;
						}
						if(split==AtomSplitResult::suppressed_subresolution)
						{
							const std::lock_guard<std::mutex> lock(reconciliation_mutex);
							++eb.suppressed_subresolution_atom_splits;
							eb.suppressed_subresolution_atom_volume+=suppressed_volume;
							eb.maximum_suppressed_subresolution_atom_volume=std::max(
								eb.maximum_suppressed_subresolution_atom_volume,suppressed_volume);
						}
						const int unchanged=static_cast<int>(next.size());
						next.push_back(atoms[atom]);children[atom]={unchanged,unchanged};
					}
					if(next.size()>512){atom_overflow=true;break;}
				}
				if(atom_overflow)break;

				// Adjacencies are born as one polygon shared by the two children of a
				// split. Later planes partition that one polygon once and give the exact
				// same pieces to both sides. This makes internal interfaces conservative
				// by construction; no tolerance-based re-matching of independently
				// clipped atom faces is needed.
				std::vector<AtomAdjacency> next_adjacencies;next_adjacencies.reserve(
					adjacencies.size()*2+new_caps.size());
				for(const AtomAdjacency& adjacency:adjacencies)
				{
					const bool a_split=children[adjacency.a][0]!=children[adjacency.a][1];
					const bool b_split=children[adjacency.b][0]!=children[adjacency.b][1];
					if(!a_split&&!b_split)
					{
						AtomAdjacency mapped=adjacency;
						mapped.a=children[adjacency.a][0];mapped.b=children[adjacency.b][0];
						next_adjacencies.push_back(std::move(mapped));continue;
					}
					double partitioned_area=0.0;
					for(int side=0;side<2;++side)
					{
						Polygon piece=clip_polygon_halfspace(adjacency.polygon,
							planes[plane].plane,side==0,length_tolerance);
						const Measure measure=polygon_measure(piece,
							planes[adjacency.plane].plane.normal);
						if(!(measure.value>area_tolerance))continue;
						partitioned_area+=measure.value;
						AtomAdjacency mapped=adjacency;
						mapped.a=children[adjacency.a][side];
						mapped.b=children[adjacency.b][side];
						mapped.polygon=std::move(piece);
						next_adjacencies.push_back(std::move(mapped));
					}
					const double original_area=polygon_measure(adjacency.polygon,
						planes[adjacency.plane].plane.normal).value;
					const double residual=std::abs(partitioned_area-original_area);
					if(residual>8.0*area_tolerance)
					{
						const std::lock_guard<std::mutex> lock(reconciliation_mutex);
						++eb.reconciled_facet_pairs;
						eb.reconciled_facet_area_residual+=residual;
						eb.maximum_reconciled_facet_area_residual=std::max(
							eb.maximum_reconciled_facet_area_residual,residual);
					}
				}
				for(AtomAdjacency& cap:new_caps)next_adjacencies.push_back(std::move(cap));
				atoms=std::move(next);
				adjacencies=std::move(next_adjacencies);
			}
			if(reported[cell])return;
			if(atom_overflow)
			{
				report(cell,"closed triangle cell partition exceeds 512 convex atoms",
					std::move(candidates));return;
			}
			for(Atom& atom:atoms)atom.measure=polyhedron_measure(atom);
			if(std::any_of(atoms.begin(),atoms.end(),[&](const Atom& atom)
				{return !(atom.measure.value>volume_tolerance)||!finite(atom.measure.centroid);}))
			{
				report(cell,"closed triangle cell partition produced a degenerate convex atom",
					std::move(candidates));return;
			}

			// Mark the two sides of the fitted face from the source mesh itself.  A
			// centroid-only parity test can miss a very thin exterior wedge when the
			// curved tessellation and its local tangent plane differ slightly.  The
			// solid's oriented triangles give an unambiguous local constraint: the
			// atom on the normal side is exterior and the opposite atom is interior.
			std::vector<double> oriented_evidence(atoms.size(),0.0);
			for(AtomAdjacency& adjacency:adjacencies)
			{
				const PlaneGroup& group=planes[adjacency.plane];
				double covered_area=0.0;
				for(int patch_index:group.patches)
				{
					Polygon projected=patches[patch_index].polygon;
					for(Vec3d& point:projected)
						point=point-group.plane.normal*group.plane.distance(point);
					const Polygon overlap=intersect_coplanar_convex(adjacency.polygon,
						projected,group.plane.normal,length_tolerance);
					covered_area+=polygon_measure(overlap,group.plane.normal).value;
				}
				if(!(covered_area>area_tolerance))continue;
				adjacency.source_surface=true;
				const double side_a=dot(group.outward_normal,atoms[adjacency.a].measure.centroid)
					-group.outward_offset;
				const double side_b=dot(group.outward_normal,atoms[adjacency.b].measure.centroid)
					-group.outward_offset;
				const int outside_atom=side_a>side_b?adjacency.a:adjacency.b;
				const int inside_atom=outside_atom==adjacency.a?adjacency.b:adjacency.a;
				oriented_evidence[outside_atom]+=covered_area;
				oriented_evidence[inside_atom]-=covered_area;
			}

			std::vector<Vec3d> atom_witnesses;atom_witnesses.reserve(atoms.size());
			std::vector<int> classified_atoms;classified_atoms.reserve(atoms.size());
			std::vector<ClosedSolidPointLocation> solid_locations(atoms.size(),
				ClosedSolidPointLocation::unknown);
			for(int atom=0;atom<static_cast<int>(atoms.size());++atom)
				if(std::abs(oriented_evidence[atom])<=area_tolerance)
				{
					classified_atoms.push_back(atom);
					atom_witnesses.push_back(atoms[atom].measure.centroid);
				}
			std::vector<ClosedSolidPointLocation> classified_locations(classified_atoms.size(),
				ClosedSolidPointLocation::unknown);
			if(!classified_atoms.empty())solid.classify_points(atom_witnesses.data(),
				atom_witnesses.size(),length_tolerance,classified_locations.data());
			for(std::size_t index=0;index<classified_atoms.size();++index)
				solid_locations[classified_atoms[index]]=classified_locations[index];
			std::vector<unsigned char> atom_outside(atoms.size(),0);
			for(int atom=0;atom<static_cast<int>(atoms.size());++atom)
			{
				// A fitted interface represents an exact oriented CAD face locally. Its
				// signed sides therefore define the material of adjacent approximation
				// atoms. OCCT classifies every atom that is not incident to such an
				// interface; it is never replaced by a triangle-ray parity guess.
				if(oriented_evidence[atom]>area_tolerance)
				{
					atom_outside[atom]=1;continue;
				}
				if(oriented_evidence[atom]<-area_tolerance)continue;
				const ClosedSolidPointLocation location=solid_locations[atom];
				if(location==ClosedSolidPointLocation::unknown)
				{
					report(cell,"OCCT could not classify a positive-volume convex atom",candidates);
					break;
				}
				if(location==ClosedSolidPointLocation::boundary)
				{
					// A very thin atom can fall within OCCT's boundary tolerance. Its source
					// surface orientation is then the only signed local evidence available.
					if(std::abs(oriented_evidence[atom])<=area_tolerance)
					{
						report(cell,"positive-volume convex atom lies on the OCCT boundary without signed surface evidence",candidates);
						break;
					}
					atom_outside[atom]=oriented_evidence[atom]>0.0;continue;
				}
				atom_outside[atom]=location==ClosedSolidPointLocation::outside;
			}
			if(reported[cell])return;
			DisjointSet components(static_cast<int>(atoms.size()));
			for(AtomAdjacency& adjacency:adjacencies)
			{
				adjacency.blocked=adjacency.source_surface||
					atom_outside[adjacency.a]!=atom_outside[adjacency.b];
				if(!adjacency.blocked)components.join(adjacency.a,adjacency.b);
			}
			std::map<int,bool> outside_root;
			for(int atom=0;atom<static_cast<int>(atoms.size());++atom)
				outside_root[components.root(atom)]=atom_outside[atom]!=0;

			std::map<int,int> local_fragment;
			for(int atom=0;atom<static_cast<int>(atoms.size());++atom)
			{
				const int root=components.root(atom);
				if(outside_root[root]&&!local_fragment.contains(root))
					local_fragment[root]=static_cast<int>(local_fragment.size());
			}
			if(local_fragment.empty())
			{
				eb.cells[cell].state=EbCellState::solid;return;
			}
			if(local_fragment.size()==outside_root.size())
			{
				known_exterior[cell]=1;
				for(RawPatch& patch:patches)patch.plus_fragment=regular_fragment(cell);
				return;
			}

			EbCellTopology& topology=eb.cells[cell];topology.state=EbCellState::split;
			topology.fragment_count=static_cast<std::uint16_t>(local_fragment.size());
			auto staged=std::make_unique<StagedSplitCell>();
			staged->fragments.reserve(local_fragment.size());
			std::vector<double> volume(local_fragment.size(),0.0);
			std::vector<Vec3d> moment(local_fragment.size());
			std::vector<int> witness_atom(local_fragment.size(),-1);
			for(int atom=0;atom<static_cast<int>(atoms.size());++atom)
			{
				const int root=components.root(atom);auto found=local_fragment.find(root);
				if(found==local_fragment.end())continue;const int fragment=found->second;
				volume[fragment]+=atoms[atom].measure.value;
				moment[fragment]=moment[fragment]+atoms[atom].measure.centroid*atoms[atom].measure.value;
				if(witness_atom[fragment]<0)witness_atom[fragment]=atom;
			}
			for(int fragment=0;fragment<static_cast<int>(local_fragment.size());++fragment)
			{
				FluidFragment output;output.parent_cell=cell;output.volume=volume[fragment];
				output.centroid=moment[fragment]/volume[fragment];
				output.merge_target=irregular_fragment(fragment);
				output.pressure_dof=cell;staged->fragments.push_back(output);
			}
			// Attach the exact clipped load triangles through the fitted CAD-face
			// partition that represents them.  A tiny normal offset from the exact
			// triangle can lie on the opposite side of the fitted plane on a curved
			// face, so plane-facet ownership is the authoritative local association.
			for(int plane=0;plane<static_cast<int>(planes.size());++plane)
				for(int patch_index:planes[plane].patches)
				{
					double best_distance=std::numeric_limits<double>::infinity();
					FragmentRef best=invalid_fragment;
					const Vec3d projected=patches[patch_index].centroid-
						planes[plane].plane.normal*planes[plane].plane.distance(
							patches[patch_index].centroid);
					for(const AtomAdjacency& adjacency:adjacencies)
					{
						if(!adjacency.blocked||adjacency.plane!=plane)continue;
						const int outside_atom=atom_outside[adjacency.a]?adjacency.a:adjacency.b;
						const int root=components.root(outside_atom);
						auto found=local_fragment.find(root);if(found==local_fragment.end())continue;
						const Measure measure=polygon_measure(adjacency.polygon,
							planes[plane].plane.normal);
						const double distance=length2(projected-measure.centroid);
						if(distance<best_distance)
						{
							best_distance=distance;
							best=irregular_fragment(found->second);
						}
					}
					patches[patch_index].plus_fragment=best;
				}
			for(RawPatch& patch:patches)if(patch.plus_fragment==invalid_fragment)
			{
				const Vec3d probe=patch.centroid+patch.normal*(1.0e-4*grid.h);
				double best_score=std::numeric_limits<double>::infinity();
				for(int atom=0;atom<static_cast<int>(atoms.size());++atom)
				{
					if(!atom_outside[atom])continue;
					double violation=0.0;
					for(const AtomFace& face:atoms[atom].faces)
						violation=std::max(violation,dot(face.normal,probe)-
							dot(face.normal,face.polygon.front()));
					const double score=violation+1.0e-9*std::sqrt(length2(
						probe-atoms[atom].measure.centroid));
					if(score>=best_score)continue;best_score=score;
					const int root=components.root(atom);auto found=local_fragment.find(root);
					if(found!=local_fragment.end())patch.plus_fragment=
						irregular_fragment(found->second);
				}
			}

			ExactEbCellLocator locator;locator.ambiguity_tolerance=4.0*length_tolerance;
			locator.winding_tolerance=1.0e-4;locator.fragments.resize(local_fragment.size());
			for(int fragment=0;fragment<static_cast<int>(local_fragment.size());++fragment)
				locator.fragments[fragment].interior_witness=atoms[witness_atom[fragment]].measure.centroid;
			for(int atom=0;atom<static_cast<int>(atoms.size());++atom)
			{
				const int root=components.root(atom);auto found=local_fragment.find(root);
				if(found==local_fragment.end())continue;
				EbConvexRegion region;region.halfspaces.reserve(atoms[atom].faces.size());
				for(const AtomFace& face:atoms[atom].faces)
					region.halfspaces.push_back({face.normal,dot(face.normal,face.polygon.front())});
				locator.fragments[found->second].convex_regions.push_back(std::move(region));
			}
			auto append_boundary=[&](int atom,const Polygon& polygon,Vec3d normal)
			{
				const int root=components.root(atom);auto found=local_fragment.find(root);
				if(found==local_fragment.end())return;
				for(const SolidFaceTriangle& triangle:triangulate(polygon,normal))
					locator.fragments[found->second].triangles.push_back(
						{triangle.a,triangle.b,triangle.c});
			};
			for(int atom=0;atom<static_cast<int>(atoms.size());++atom)
				if(local_fragment.contains(components.root(atom)))
					for(const AtomFace& face:atoms[atom].faces)if(face.box_face>=0)
						append_boundary(atom,face.polygon,face.normal);
			for(const AtomAdjacency& adjacency:adjacencies)if(adjacency.blocked)
			{
				const int outside_atom=outside_root[components.root(adjacency.a)]?
					adjacency.a:adjacency.b;
				const Vec3d normal=outside_atom==adjacency.a?
					planes[adjacency.plane].plane.normal:planes[adjacency.plane].plane.normal*-1.0;
				Polygon polygon=adjacency.polygon;
				if(outside_atom==adjacency.b)std::reverse(polygon.begin(),polygon.end());
				append_boundary(outside_atom,polygon,normal);
			}
			for(int atom=0;atom<static_cast<int>(atoms.size());++atom)
			{
				const int root=components.root(atom);auto found=local_fragment.find(root);
				if(found==local_fragment.end())continue;
				const FragmentRef fragment=irregular_fragment(found->second);
				for(const AtomFace& face:atoms[atom].faces)if(face.box_face>=0)
				{
					const Measure measure=polygon_measure(face.polygon,face.normal);
					if(!(measure.value>area_tolerance))continue;
					const int axis=face.box_face/2;const bool upper=(face.box_face&1)!=0;
					cell_apertures[cell].push_back({axis,upper,measure.value,
						measure.centroid,fragment,triangulate(face.polygon,face.normal)});
				}
			}
			staged->locator=std::move(locator);split_cells[cell]=std::move(staged);
		});
		for(int cell=0;cell<grid.cell_count();++cell)if(split_cells[cell])
		{
			EbCellTopology& topology=eb.cells[cell];
			const int fragment_offset=static_cast<int>(eb.fragments.size());
			topology.first_fragment=fragment_offset;
			topology.exact_locator_index=static_cast<int>(eb.exact_locators.size());
			auto remap=[&](FragmentRef ref)
			{
				return ref==invalid_fragment||fragment_is_regular(ref)?ref:
					irregular_fragment(fragment_offset+irregular_fragment_index(ref));
			};
			for(RawPatch& patch:cell_patches[cell])patch.plus_fragment=remap(patch.plus_fragment);
			for(StagedAperture& aperture:cell_apertures[cell])aperture.fragment=remap(aperture.fragment);
			for(FluidFragment& fragment:split_cells[cell]->fragments)
				fragment.merge_target=remap(fragment.merge_target);
			eb.fragments.insert(eb.fragments.end(),
				std::make_move_iterator(split_cells[cell]->fragments.begin()),
				std::make_move_iterator(split_cells[cell]->fragments.end()));
			eb.exact_locators.push_back(std::move(split_cells[cell]->locator));
			eb.irregular_cells.push_back(cell);
		}
		const double cell_build_ms=std::chrono::duration<double,std::milli>(
			std::chrono::steady_clock::now()-cell_build_begin).count();
		std::fprintf(stderr,"[closed-solid-eb] %d cells, %u CPU workers, %.1f ms cell BSP\n",
			grid.cell_count(),cell_workers,cell_build_ms);
		const auto flood_begin=std::chrono::steady_clock::now();

		// The cut shell is already resolved. Label the remaining regular Cartesian
		// cells by a cheap graph flood, seeding the atlas boundary and every open face
		// of a certified exterior cut fragment. No CAD or ray query occurs here.
		std::deque<int> queue;std::vector<unsigned char> flooded(
			static_cast<std::size_t>(grid.cell_count()),0);
		auto seed=[&](int cell)
		{
			if(cell<0||cell>=grid.cell_count()||eb.cells[cell].state!=EbCellState::regular||
				flooded[cell])return;
			flooded[cell]=1;queue.push_back(cell);
		};
		for(int k=0;k<grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<grid.nx;++i)
			if(i==0||j==0||k==0||i+1==grid.nx||j+1==grid.ny||k+1==grid.nz)
				seed(grid.cell_index(i,j,k));
		for(int cell=0;cell<grid.cell_count();++cell)
		{
			if(surface_cell[cell]&&known_exterior[cell])seed(cell);
			for(const StagedAperture& aperture:cell_apertures[cell])
				seed(neighbour_cell(grid,cell,aperture.axis,aperture.upper));
		}
		while(!queue.empty())
		{
			const int cell=queue.front();queue.pop_front();
			const auto coordinate=grid.cell_coord(cell);
			if(coordinate[0]>0)seed(cell-1);if(coordinate[0]+1<grid.nx)seed(cell+1);
			if(coordinate[1]>0)seed(cell-grid.nx);if(coordinate[1]+1<grid.ny)seed(cell+grid.nx);
			const int slab=grid.nx*grid.ny;
			if(coordinate[2]>0)seed(cell-slab);if(coordinate[2]+1<grid.nz)seed(cell+slab);
		}
		for(int cell=0;cell<grid.cell_count();++cell)
			if(eb.cells[cell].state==EbCellState::regular&&!flooded[cell])
				eb.cells[cell].state=EbCellState::solid;
		const auto flood_end=std::chrono::steady_clock::now();

		// The tessellated solid is the load surface. Each triangle is clipped to the
		// Cartesian cells once; its outward side is attached to the exterior fragment.
		for(int cell=0;cell<grid.cell_count();++cell)for(const RawPatch& source:cell_patches[cell])
		{
			if(!source.output_owner)continue;
			FragmentRef plus=source.plus_fragment;
			for(double scale:{1.0e-6,1.0e-5,1.0e-4,1.0e-3})if(plus==invalid_fragment)
			{
				const Vec3d probe=source.centroid+source.normal*(scale*grid.h);
				const int probe_cell=containing_cell(grid,probe);if(probe_cell<0)continue;
				plus=eb.fragment_containing_point(probe_cell,probe,4.0*length_tolerance);
				if(plus!=invalid_fragment)break;
			}
			if(plus==invalid_fragment)
			{
				const Vec3d probe=source.centroid+source.normal*(1.0e-4*grid.h);
				const int probe_cell=containing_cell(grid,probe);
				std::ostringstream message;message.precision(17);
				message<<"could not attach an oriented triangle patch to exterior fluid; normal=["
					<<source.normal.x<<","<<source.normal.y<<","<<source.normal.z
					<<"] probe=["<<probe.x<<","<<probe.y<<","<<probe.z<<"] probe-cell="
					<<probe_cell;
				if(probe_cell>=0)message<<" state="<<static_cast<int>(eb.cells[probe_cell].state)
					<<" mesh-location="<<static_cast<int>(broad_phase.classify_closed_point(
						probe,length_tolerance));
				report(cell,message.str(),{source.triangle});continue;
			}
			eb.patches.push_back({source.triangle,source.face,source.area,source.centroid,
				source.normal,plus,invalid_fragment});
		}
		const auto patches_end=std::chrono::steady_clock::now();

		struct StagedFaceConnection
		{
			int lower_cell=0,axis=0,ordinal=0;double area=0;Vec3d centroid{};
			FragmentRef lower=invalid_fragment,upper=invalid_fragment;
		};
		std::vector<std::vector<StagedFaceConnection>> worker_face_connections;
		const unsigned face_workers=std::min(std::max(1u,std::thread::hardware_concurrency()),
			static_cast<unsigned>(std::max(1,grid.cell_count())));
		worker_face_connections.resize(face_workers);
		auto append_connection=[&](int lower_cell,int axis,double area,Vec3d centroid,
			FragmentRef lower,FragmentRef upper)
		{
			if(!(area>0.0)||lower==invalid_fragment||upper==invalid_fragment||lower==upper)return;
			eb.apertures.push_back({lower_cell,static_cast<std::int8_t>(axis),area,
				centroid,lower,upper});
			const double distance=std::max(1.0e-12,
				std::sqrt(length2(eb.fragment_centroid(upper)-eb.fragment_centroid(lower))));
			eb.connections.push_back({lower,upper,area,centroid,distance,
				static_cast<std::int8_t>(axis)});
			const double reporting_area=options.min_aperture_area_fraction*grid.h*grid.h;
			if(area<reporting_area)
			{
				++eb.retained_subgrid_apertures;eb.retained_subgrid_aperture_area+=area;
			}
		};
		auto face_apertures=[&](int cell,int axis,bool upper)
		{
			std::vector<const StagedAperture*> result;
			for(const StagedAperture& aperture:cell_apertures[cell])
				if(aperture.axis==axis&&aperture.upper==upper)result.push_back(&aperture);
			return result;
		};
		const double area_roundoff=std::max(area_tolerance,
			1024.0*std::numeric_limits<double>::epsilon()*grid.h*grid.h);
		const unsigned active_face_workers=parallel_for_ranges(grid.cell_count(),[&](int begin,int end,unsigned worker)
		{
			auto& staged=worker_face_connections[worker];
			for(int lower=begin;lower<end;++lower)for(int axis=0;axis<3;++axis)
			{
				const int upper=neighbour_cell(grid,lower,axis,true);if(upper<0)continue;
				const EbCellState lower_state=eb.cells[lower].state,
					upper_state=eb.cells[upper].state;
				if(lower_state==EbCellState::regular&&upper_state==EbCellState::regular)continue;
				eb.cut_face_mask[lower]|=static_cast<std::uint8_t>(1u<<axis);
				if(lower_state==EbCellState::unresolved||upper_state==EbCellState::unresolved||
					lower_state==EbCellState::solid||upper_state==EbCellState::solid)continue;
				const auto lower_open=face_apertures(lower,axis,true);
				const auto upper_open=face_apertures(upper,axis,false);
				int ordinal=0;
				auto stage=[&](double area,Vec3d centroid,FragmentRef a,FragmentRef b)
				{
					if(area>0.0&&a!=invalid_fragment&&b!=invalid_fragment&&a!=b)
						staged.push_back({lower,axis,ordinal++,area,centroid,a,b});
				};
				if(lower_state==EbCellState::regular)
				{
					for(const StagedAperture* aperture:upper_open)
						stage(aperture->area,aperture->centroid,
							regular_fragment(lower),aperture->fragment);
					continue;
				}
				if(upper_state==EbCellState::regular)
				{
					for(const StagedAperture* aperture:lower_open)
						stage(aperture->area,aperture->centroid,
							aperture->fragment,regular_fragment(upper));
					continue;
				}
				struct CommonRegion
				{
					double area=0.0;
					Vec3d moment{};
				};
				std::map<std::pair<FragmentRef,FragmentRef>,CommonRegion> common_regions;
				double common_area=0.0;
				for(const StagedAperture* a:lower_open)for(const StagedAperture* b:upper_open)
					for(const SolidFaceTriangle& ta:a->triangles)
						for(const SolidFaceTriangle& tb:b->triangles)
						{
							const Polygon pa{ta.a,ta.b,ta.c},pb{tb.a,tb.b,tb.c};
							Vec3d face_normal{};face_normal[axis]=1.0;
							const Polygon common=intersect_coplanar_convex(pa,pb,face_normal,
								length_tolerance);
							const Measure measure=polygon_measure(common,face_normal);
							if(!(measure.value>area_roundoff))continue;
							CommonRegion& region=common_regions[{a->fragment,b->fragment}];
							region.area+=measure.value;
							region.moment=region.moment+measure.centroid*measure.value;
							common_area+=measure.value;
						}
				// Build this Cartesian face once. The common refinement determines the
				// conservative open area and fragment pairing, but is coalesced so its
				// tessellation slivers never become independent pressure edges.
				for(const auto& [pair,region]:common_regions)
					stage(region.area,region.moment/region.area,pair.first,pair.second);
			}
		});
		std::vector<StagedFaceConnection> staged_connections;
		std::size_t staged_count=0;for(const auto& worker:worker_face_connections)staged_count+=worker.size();
		staged_connections.reserve(staged_count);
		for(auto& worker:worker_face_connections)staged_connections.insert(staged_connections.end(),
			std::make_move_iterator(worker.begin()),std::make_move_iterator(worker.end()));
		std::sort(staged_connections.begin(),staged_connections.end(),[](const auto& a,const auto& b)
		{
			if(a.lower_cell!=b.lower_cell)return a.lower_cell<b.lower_cell;
			if(a.axis!=b.axis)return a.axis<b.axis;return a.ordinal<b.ordinal;
		});
		for(const auto& connection:staged_connections)
			append_connection(connection.lower_cell,connection.axis,connection.area,
				connection.centroid,connection.lower,connection.upper);
		const auto connections_end=std::chrono::steady_clock::now();

		std::vector<int> counts(eb.fragments.size(),0);
		for(const FragmentConnection& connection:eb.connections)
		{
			if(!fragment_is_regular(connection.fragment_a))
				++counts[irregular_fragment_index(connection.fragment_a)];
			if(!fragment_is_regular(connection.fragment_b))
				++counts[irregular_fragment_index(connection.fragment_b)];
		}
		int offset=0;for(std::size_t fragment=0;fragment<eb.fragments.size();++fragment)
		{
			eb.fragments[fragment].connection_offset=offset;
			eb.fragments[fragment].connection_count=counts[fragment];offset+=counts[fragment];
		}
		std::sort(eb.unresolved.begin(),eb.unresolved.end(),
			[](const UnresolvedEbCell& a,const UnresolvedEbCell& b)
			{return a.parent_cell<b.parent_cell;});
		std::fprintf(stderr,"[closed-solid-eb] flood %.1f ms, patches %.1f ms, faces %.1f ms (%u CPU workers)\n",
			std::chrono::duration<double,std::milli>(flood_end-flood_begin).count(),
			std::chrono::duration<double,std::milli>(patches_end-flood_end).count(),
			std::chrono::duration<double,std::milli>(connections_end-patches_end).count(),active_face_workers);
		std::fflush(stderr);
		return eb;
	}
}
