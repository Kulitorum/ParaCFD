#include "core/fluid/amr_fields.h"
#include "core/fluid/amr_grid.h"
#include "core/geometry/embedded_boundary.h"
#include "core/geometry/step_import.h"
#include "core/geometry/triangle_bvh.h"
#include "core/paraglider_config.h"

#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using namespace paracfd::core;

int main(int argc,char**argv)
{
	std::string config_path="configs/paraglider.json",step_override;double tessellation_override=0;
	for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--config"&&i+1<argc)config_path=argv[++i];else if(a=="--step"&&i+1<argc)step_override=argv[++i];else if(a=="--tessellation-mm"&&i+1<argc)tessellation_override=std::atof(argv[++i]);else{std::fprintf(stderr,"usage: paraglider_probe [--config file] [--step wing.step] [--tessellation-mm value]\n");return 2;}}
	ParagliderConfig cfg;std::string error;if(!load_paraglider_config(config_path,cfg,&error)){std::fprintf(stderr,"[paraglider] %s\n",error.c_str());return 2;}if(!step_override.empty())cfg.step_path=step_override;if(tessellation_override>0)cfg.tessellation_deflection_mm=tessellation_override;if(cfg.step_path.empty()){std::fprintf(stderr,"[paraglider] no STEP path; set step.path or pass --step\n");return 2;}
	auto t0=std::chrono::steady_clock::now();TriMesh source=load_step_mesh(cfg.step_path,cfg.tessellation_deflection_mm,&error);if(source.empty()){std::fprintf(stderr,"[paraglider] STEP load failed: %s\n",error.c_str());return 2;}TriMesh wing=placed_mesh(source,cfg.placement);auto t1=std::chrono::steady_clock::now();TriangleBvh bvh(wing);auto t2=std::chrono::steady_clock::now();Aabb3d domain=automatic_flow_domain(wing,cfg.domain);AmrHierarchy amr=AmrHierarchy::build_static(domain,wing,bvh,cfg.amr);auto t3=std::chrono::steady_clock::now();
	struct FaceStats{Aabb3d box;double area=0;std::size_t triangles=0;};std::uint32_t max_face=0;for(std::uint32_t id:wing.source_face_ids)max_face=std::max(max_face,id);std::vector<FaceStats> face_stats(wing.source_face_ids.empty()?0:static_cast<std::size_t>(max_face)+1);
	for(std::size_t triangle=0;triangle<wing.triangle_count()&&triangle<wing.source_face_ids.size();++triangle){FaceStats& s=face_stats[wing.source_face_ids[triangle]];const std::uint32_t ia=wing.indices[3*triangle],ib=wing.indices[3*triangle+1],ic=wing.indices[3*triangle+2];Vec3d a{wing.positions[3*ia],wing.positions[3*ia+1],wing.positions[3*ia+2]},b{wing.positions[3*ib],wing.positions[3*ib+1],wing.positions[3*ib+2]},c{wing.positions[3*ic],wing.positions[3*ic+1],wing.positions[3*ic+2]};s.box.expand(a);s.box.expand(b);s.box.expand(c);s.area+=0.5*std::sqrt(length2(cross(b-a,c-a)));++s.triangles;}
	std::size_t fragments=0,patches=0,apertures=0,unresolved=0,eb_bricks=0;
	EmbeddedBoundaryBuildOptions eo;eo.min_volume_fraction=cfg.amr.min_volume_fraction;
	for(const AmrLevel& level:amr.levels())for(const BrickMetadata& brick:level.bricks)if(brick.active()&&brick.embedded_boundary())
	{
		++eb_bricks;UniformEbGrid grid{brick.origin,amr.brick_size(),amr.brick_size(),amr.brick_size(),brick.h};EmbeddedBoundary eb=build_embedded_boundary(wing,bvh,grid,eo);fragments+=eb.fragments.size();patches+=eb.patches.size();apertures+=eb.apertures.size();unresolved+=eb.unresolved.size();
	}
	auto t4=std::chrono::steady_clock::now();AmrHostFields fields(amr);
	auto ms=[](auto a,auto b){return std::chrono::duration<double,std::milli>(b-a).count();};
	std::printf("ParaCFD paraglider preprocessing\n");
	std::printf("  STEP: %zu triangles, provenance=%s, UV=%s, tessellation=%.4g mm (%.2f ms)\n",wing.triangle_count(),wing.has_face_provenance()?"yes":"no",wing.has_uv()?"yes":"no",cfg.tessellation_deflection_mm,ms(t0,t1));
	if(!face_stats.empty()){std::printf("  tessellated CAD faces: %zu (standalone STEP wires/curve sets are excluded)\n",face_stats.size());for(int axis=0;axis<3;++axis){int lo_face=-1,hi_face=-1;double lo=std::numeric_limits<double>::infinity(),hi=-lo;for(int f=0;f<(int)face_stats.size();++f)if(face_stats[f].triangles){if(face_stats[f].box.lo[axis]<lo){lo=face_stats[f].box.lo[axis];lo_face=f;}if(face_stats[f].box.hi[axis]>hi){hi=face_stats[f].box.hi[axis];hi_face=f;}}std::printf("    bbox axis %c: min %.4f by face %d (A=%.4g m^2), max %.4f by face %d (A=%.4g m^2)\n",'X'+axis,lo,lo_face,face_stats[lo_face].area,hi,hi_face,face_stats[hi_face].area);}}
	std::printf("  BVH: %zu triangles (%.2f ms)\n",bvh.triangle_count(),ms(t1,t2));
	std::printf("  domain: [%.3f %.3f %.3f] .. [%.3f %.3f %.3f] m\n",amr.domain().lo.x,amr.domain().lo.y,amr.domain().lo.z,amr.domain().hi.x,amr.domain().hi.y,amr.domain().hi.z);
	std::printf("  AMR: %zu active bricks, %zu active cells, finest h=%.6g m, balanced=%s (%.2f ms)\n",amr.active_brick_count(),amr.active_cell_count(),amr.finest_cell_size(),amr.is_two_to_one_balanced()?"yes":"NO",ms(t2,t3));
	for(const AmrLevel& l:amr.levels())std::printf("    level %d: h=%.6g m, %zu active / %zu allocated bricks\n",l.level,l.h,l.active_bricks,l.bricks.size());
	std::printf("  local EB: %zu bricks, %zu fragments, %zu apertures, %zu patches, %zu unresolved cells (%.2f ms)\n",eb_bricks,fragments,apertures,patches,unresolved,ms(t3,t4));
	std::printf("  pooled FP32 field estimate: %.2f MiB (u/v/w/p/nut/temp including halos)\n",fields.bytes()/(1024.0*1024.0));
	if(unresolved)std::fprintf(stderr,"[paraglider] NOT FLOW-READY: unresolved finest-level EB topology remains; increase max_levels/refine geometry. No sides were merged.\n");
	std::printf("  limitation: composite cross-brick EB connectivity/pressure solve is not assembled by this preprocessing probe yet.\n");
	return unresolved?3:0;
}
