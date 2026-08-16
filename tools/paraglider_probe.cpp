#include "core/aero_loads.h"
#include "core/fluid/amr_advection.h"
#include "core/fluid/amr_fields.h"
#include "core/fluid/amr_eb.h"
#include "core/fluid/amr_grid.h"
#include "core/fluid/amr_pressure.h"
#include "core/geometry/embedded_boundary.h"
#include "core/geometry/step_import.h"
#include "core/geometry/triangle_bvh.h"
#include "core/paraglider_config.h"

#include <chrono>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include <cuda_runtime.h>

using namespace paracfd::core;

namespace
{
	struct DivergenceMetrics
	{
		double maximum = 0.0;
		double maximum_volume = 0.0;
		double maximum_integrated = 0.0;
		double absolute_integrated = 0.0;
		double net_integrated = 0.0;
		double volume_weighted_rms = 0.0;
	};

	DivergenceMetrics divergence_metrics(const CompositeAmrPressureSystem& system,
		const std::vector<Real>& divergence)
	{
		DivergenceMetrics out;double weighted_square=0.0,total_volume=0.0;
		for(std::size_t q=0;q<divergence.size()&&q<system.volume.size();++q)if(system.active[q]&&system.volume[q]>0)
		{
			const double value=static_cast<double>(divergence[q]),volume=system.volume[q],absolute=std::abs(value),integrated=value*volume;
			if(absolute>out.maximum){out.maximum=absolute;out.maximum_volume=volume;}
			out.maximum_integrated=std::max(out.maximum_integrated,std::abs(integrated));out.absolute_integrated+=std::abs(integrated);out.net_integrated+=integrated;weighted_square+=value*value*volume;total_volume+=volume;
		}
		out.volume_weighted_rms=total_volume>0?std::sqrt(weighted_square/total_volume):0.0;return out;
	}
}

int main(int argc,char**argv)
{
	std::string config_path="configs/paraglider.json",step_override;double tessellation_override=0,sheet_angle_degrees=0,sheet_distance_fraction=0,projection_tolerance_override=0;int max_levels_override=0,complex_subdivisions=0,projection_max_override=0;bool project_initial=false,run_step=false;
	for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--config"&&i+1<argc)config_path=argv[++i];else if(a=="--step"&&i+1<argc)step_override=argv[++i];else if(a=="--tessellation-mm"&&i+1<argc)tessellation_override=std::atof(argv[++i]);else if(a=="--max-levels"&&i+1<argc)max_levels_override=std::atoi(argv[++i]);else if(a=="--sheet-angle-deg"&&i+1<argc)sheet_angle_degrees=std::atof(argv[++i]);else if(a=="--sheet-distance-frac"&&i+1<argc)sheet_distance_fraction=std::atof(argv[++i]);else if(a=="--complex-subdivisions"&&i+1<argc)complex_subdivisions=std::atoi(argv[++i]);else if(a=="--project-initial")project_initial=true;else if(a=="--run-step"){project_initial=true;run_step=true;}else if(a=="--projection-max"&&i+1<argc)projection_max_override=std::atoi(argv[++i]);else if(a=="--projection-tolerance"&&i+1<argc)projection_tolerance_override=std::atof(argv[++i]);else{std::fprintf(stderr,"usage: paraglider_probe [--config file] [--step wing.step] [--tessellation-mm value] [--max-levels count] [--sheet-angle-deg value] [--sheet-distance-frac value] [--complex-subdivisions N] [--project-initial|--run-step] [--projection-max N] [--projection-tolerance value]\n");return 2;}}
	ParagliderConfig cfg;std::string error;if(!load_paraglider_config(config_path,cfg,&error)){std::fprintf(stderr,"[paraglider] %s\n",error.c_str());return 2;}if(!step_override.empty())cfg.step_path=step_override;if(tessellation_override>0)cfg.tessellation_deflection_mm=tessellation_override;if(max_levels_override>0)cfg.amr.max_levels=max_levels_override;if(projection_max_override>0)cfg.solver.projection_max_iterations=projection_max_override;if(projection_tolerance_override>0)cfg.solver.projection_tolerance=projection_tolerance_override;if(cfg.step_path.empty()){std::fprintf(stderr,"[paraglider] no STEP path; set step.path or pass --step\n");return 2;}
	auto t0=std::chrono::steady_clock::now();TriMesh source=load_step_mesh(cfg.step_path,cfg.tessellation_deflection_mm,&error);if(source.empty()){std::fprintf(stderr,"[paraglider] STEP load failed: %s\n",error.c_str());return 2;}cfg.placement=frame_wing_for_external_domain(source,cfg.placement,cfg.domain.upstream_margin,cfg.domain.lateral_margin,cfg.domain.vertical_margin,cfg.amr.base_cell_size*cfg.amr.brick_size);TriMesh wing=placed_mesh(source,cfg.placement);auto t1=std::chrono::steady_clock::now();TriangleBvh bvh(wing);auto t2=std::chrono::steady_clock::now();Aabb3d domain=automatic_flow_domain(wing,cfg.domain);AmrHierarchy amr=AmrHierarchy::build_static(domain,wing,bvh,cfg.amr);auto t3=std::chrono::steady_clock::now();
	struct FaceStats{Aabb3d box;double area=0;std::size_t triangles=0;};std::uint32_t max_face=0;for(std::uint32_t id:wing.source_face_ids)max_face=std::max(max_face,id);std::vector<FaceStats> face_stats(wing.source_face_ids.empty()?0:static_cast<std::size_t>(max_face)+1);
	for(std::size_t triangle=0;triangle<wing.triangle_count()&&triangle<wing.source_face_ids.size();++triangle){FaceStats& s=face_stats[wing.source_face_ids[triangle]];const std::uint32_t ia=wing.indices[3*triangle],ib=wing.indices[3*triangle+1],ic=wing.indices[3*triangle+2];Vec3d a{wing.positions[3*ia],wing.positions[3*ia+1],wing.positions[3*ia+2]},b{wing.positions[3*ib],wing.positions[3*ib+1],wing.positions[3*ib+2]},c{wing.positions[3*ic],wing.positions[3*ic+1],wing.positions[3*ic+2]};s.box.expand(a);s.box.expand(b);s.box.expand(c);s.area+=0.5*std::sqrt(length2(cross(b-a,c-a)));++s.triangles;}
	std::size_t fragments=0,patches=0,apertures=0,unresolved=0,eb_bricks=0,static_pockets=0,discarded_apertures=0;double discarded_aperture_area=0;std::map<std::string,std::size_t> unresolved_reasons;std::map<std::string,std::array<std::size_t,2>> unresolved_face_multiplicity;
	EmbeddedBoundaryBuildOptions eo;eo.min_volume_fraction=cfg.amr.min_volume_fraction;eo.min_aperture_area_fraction=cfg.amr.min_aperture_area_fraction;eo.complex_subdivisions=cfg.amr.complex_subdivisions;if(sheet_angle_degrees>0)eo.smooth_sheet_angle_tolerance=sheet_angle_degrees*3.14159265358979323846/180.0;if(sheet_distance_fraction>0)eo.smooth_sheet_distance_fraction=sheet_distance_fraction;if(complex_subdivisions>0)eo.complex_subdivisions=complex_subdivisions;
	for(const AmrLevel& level:amr.levels())for(const BrickMetadata& brick:level.bricks)if(brick.active()&&brick.embedded_boundary())++eb_bricks;AmrEmbeddedBoundaryAtlas eb_atlas=build_amr_embedded_boundary_atlas(amr,wing,bvh,eo);
	for(const AmrEbLevelAtlas& atlas:eb_atlas.levels){const EmbeddedBoundary& eb=atlas.topology;discarded_apertures+=eb.discarded_subgrid_apertures;discarded_aperture_area+=eb.discarded_subgrid_aperture_area;for(const FluidFragment& fragment:eb.fragments)if(atlas.owned_cell[fragment.parent_cell]){++fragments;if(fragment.pressure_static)++static_pockets;}for(const FaceAperture& aperture:eb.apertures)if(atlas.owned_cell[aperture.parent_face_cell])++apertures;for(const SurfacePatch& patch:eb.patches){const BrickLocation owner=amr.locate_finest(patch.centroid);if(owner.found()&&owner.level==atlas.level)++patches;}for(const UnresolvedEbCell& problem:eb.unresolved)if(atlas.owned_cell[problem.parent_cell]){++unresolved;++unresolved_reasons[problem.reason];std::vector<std::uint32_t> faces;for(std::uint32_t triangle:problem.source_triangles)if(triangle<wing.source_face_ids.size())faces.push_back(wing.source_face_ids[triangle]);std::sort(faces.begin(),faces.end());faces.erase(std::unique(faces.begin(),faces.end()),faces.end());++unresolved_face_multiplicity[problem.reason][faces.size()>1];}}
	auto t4=std::chrono::steady_clock::now();std::size_t composite_dofs=0,coarse_fine_connections=0,embedded_connections=0,pressure_gauges=0;bool composite_ready=false;std::optional<CompositeAmrPressureSystem> pressure_system;std::string composite_error;if(!unresolved)try{pressure_system=build_composite_amr_pressure_system(amr,eb_atlas,true);composite_dofs=pressure_system->storage_size;coarse_fine_connections=pressure_system->coarse_fine.size();embedded_connections=pressure_system->embedded.size();pressure_gauges=pressure_system->gauges.size();composite_ready=true;}catch(const std::exception& exception){composite_error=exception.what();}AmrHostFields fields(amr);
	auto ms=[](auto a,auto b){return std::chrono::duration<double,std::milli>(b-a).count();};
	std::printf("ParaCFD paraglider preprocessing\n");
	std::printf("  STEP: %zu triangles, provenance=%s, UV=%s, tessellation=%.4g mm (%.2f ms)\n",wing.triangle_count(),wing.has_face_provenance()?"yes":"no",wing.has_uv()?"yes":"no",cfg.tessellation_deflection_mm,ms(t0,t1));
	if(!face_stats.empty()){std::printf("  tessellated CAD faces: %zu (standalone STEP wires/curve sets are excluded)\n",face_stats.size());for(int axis=0;axis<3;++axis){int lo_face=-1,hi_face=-1;double lo=std::numeric_limits<double>::infinity(),hi=-lo;for(int f=0;f<(int)face_stats.size();++f)if(face_stats[f].triangles){if(face_stats[f].box.lo[axis]<lo){lo=face_stats[f].box.lo[axis];lo_face=f;}if(face_stats[f].box.hi[axis]>hi){hi=face_stats[f].box.hi[axis];hi_face=f;}}std::printf("    bbox axis %c: min %.4f by face %d (A=%.4g m^2), max %.4f by face %d (A=%.4g m^2)\n",'X'+axis,lo,lo_face,face_stats[lo_face].area,hi,hi_face,face_stats[hi_face].area);}}
	std::printf("  BVH: %zu triangles (%.2f ms)\n",bvh.triangle_count(),ms(t1,t2));
	std::printf("  domain: [%.3f %.3f %.3f] .. [%.3f %.3f %.3f] m\n",amr.domain().lo.x,amr.domain().lo.y,amr.domain().lo.z,amr.domain().hi.x,amr.domain().hi.y,amr.domain().hi.z);
	std::printf("  AMR: %zu active bricks, %zu active cells, finest h=%.6g m, balanced=%s (%.2f ms)\n",amr.active_brick_count(),amr.active_cell_count(),amr.finest_cell_size(),amr.is_two_to_one_balanced()?"yes":"NO",ms(t2,t3));
	for(const AmrLevel& l:amr.levels())std::printf("    level %d: h=%.6g m, %zu active / %zu allocated bricks\n",l.level,l.h,l.active_bricks,l.bricks.size());
	std::printf("  sparse cross-brick EB atlas: %zu source bricks, %zu fragments, %zu apertures, %zu patches, %zu owned unresolved cells (%.2f ms)\n",eb_bricks,fragments,apertures,patches,unresolved,ms(t3,t4));
	std::printf("    sub-grid aperture filter: %zu discarded, %.6g m^2 total (threshold %.6g h^2)\n",discarded_apertures,discarded_aperture_area,cfg.amr.min_aperture_area_fraction);
	std::printf("    pressure-static isolated pockets: %zu (retained state/volume, excluded from projection)\n",static_pockets);
	for(const auto& reason:unresolved_reasons){const auto multiplicity=unresolved_face_multiplicity[reason.first];std::printf("    unresolved: %zu × %s [single CAD face %zu, multiple faces %zu]\n",reason.second,reason.first.c_str(),multiplicity[0],multiplicity[1]);}
	std::printf("  pooled FP32 field estimate: %.2f MiB (u/v/w/p/nut/temp including halos)\n",fields.bytes()/(1024.0*1024.0));
	if(unresolved)std::fprintf(stderr,"[paraglider] NOT FLOW-READY: unresolved finest-level EB topology remains; increase max_levels/refine geometry. No sides were merged.\n");
	if(composite_ready)std::printf("  composite pressure topology: %zu DOFs, %zu coarse/fine + %zu EB connections, %zu disconnected-component gauges (GPU projection ready)\n",composite_dofs,coarse_fine_connections,embedded_connections,pressure_gauges);else if(!unresolved)std::fprintf(stderr,"[paraglider] composite pressure topology rejected: %s\n",composite_error.c_str());
	bool projection_ok=true;if(project_initial&&pressure_system)try{DeviceAmrFields device_fields(amr);device_fields.initialize_freestream(static_cast<Real>(cfg.freestream.speed));device_fields.apply_external_aero_boundaries(static_cast<Real>(cfg.freestream.speed));DeviceCompositeAmrProjection projection(*pressure_system,device_fields);projection.initialize_special_freestream(static_cast<Real>(cfg.freestream.speed));projection.sync_coarse_fine_from_fields();projection.compute_divergence();std::vector<Real> before,after;projection.download_divergence(before);const DivergenceMetrics d0=divergence_metrics(*pressure_system,before);const Real dt=static_cast<Real>(cfg.solver.cfl*amr.finest_cell_size()/std::max(1e-9,cfg.freestream.speed));auto projection_start=std::chrono::steady_clock::now();AmrGpuSolveResult result=projection.project(static_cast<Real>(cfg.freestream.rho),dt,cfg.solver.projection_tolerance,cfg.solver.projection_max_iterations,false);auto projection_end=std::chrono::steady_clock::now();projection.download_divergence(after);const DivergenceMetrics d1=divergence_metrics(*pressure_system,after);std::size_t free_bytes=0,total_bytes=0;cudaMemGetInfo(&free_bytes,&total_bytes);std::printf("  initial freestream projection: %s, %d iterations, residual %.3e, max div %.3e -> %.3e, %.2f ms\n",result.converged?"converged":"NOT CONVERGED",result.iterations,result.relative_residual,d0.maximum,d1.maximum,ms(projection_start,projection_end));std::printf("    conservation: volume-weighted RMS %.3e, worst fragment V=%.3e m^3, max |div*V|=%.3e m^3/s, sum |div*V|=%.3e m^3/s, net=%.3e m^3/s\n",d1.volume_weighted_rms,d1.maximum_volume,d1.maximum_integrated,d1.absolute_integrated,d1.net_integrated);std::printf("    persistent GPU estimate: fields %.2f + projection %.2f = %.2f MiB; device used %.1f / %.1f MiB\n",device_fields.bytes()/(1024.0*1024.0),projection.bytes()/(1024.0*1024.0),(device_fields.bytes()+projection.bytes())/(1024.0*1024.0),(total_bytes-free_bytes)/(1024.0*1024.0),total_bytes/(1024.0*1024.0));projection_ok=result.converged;if(run_step&&projection_ok){const auto advection_build_begin=std::chrono::steady_clock::now();DeviceAmrAdvection advection(amr,bvh);const auto advection_build_end=std::chrono::steady_clock::now(),step_begin=advection_build_end;advection.advect(device_fields,dt);const auto advection_end=std::chrono::steady_clock::now();advection.diffuse_smagorinsky(device_fields,static_cast<Real>(cfg.freestream.nu),static_cast<Real>(cfg.solver.smagorinsky_cs),dt);const auto turbulence_end=std::chrono::steady_clock::now();device_fields.apply_external_aero_boundaries(static_cast<Real>(cfg.freestream.speed));projection.transport_embedded_apertures(dt,static_cast<Real>(cfg.freestream.nu));const auto embedded_transport_end=std::chrono::steady_clock::now();projection.sync_coarse_fine_from_fields();AmrGpuSolveResult step_result=projection.project(static_cast<Real>(cfg.freestream.rho),dt,cfg.solver.projection_tolerance,cfg.solver.projection_max_iterations,true);const auto step_end=std::chrono::steady_clock::now();std::vector<Real> pressure_real;projection.download_pressure(pressure_real);std::vector<double> pressure(pressure_real.begin(),pressure_real.end());AerodynamicLoads loads=compute_pressure_loads(*pressure_system,pressure,wing.triangle_count(),cfg.freestream,cfg.reference);std::printf("  side-safe advection preprocessing: %.2f ms, %zu / %zu protected faces, %.2f MiB\n",ms(advection_build_begin,advection_build_end),advection.protected_face_count(),advection.active_face_count(),advection.bytes()/(1024.0*1024.0));std::printf("  first external-aero step: %s, advection %.2f ms, LES/diffusion %.2f ms, EB transport %.2f ms, projection %.2f ms (%d iterations), total %.2f ms\n",step_result.converged?"converged":"NOT CONVERGED",ms(step_begin,advection_end),ms(advection_end,turbulence_end),ms(turbulence_end,embedded_transport_end),ms(embedded_transport_end,step_end),step_result.iterations,ms(step_begin,step_end));std::printf("    pressure-only force [N]: [%.6g %.6g %.6g], represented patches=%zu; total persistent estimate %.2f MiB\n",loads.pressure_force.x,loads.pressure_force.y,loads.pressure_force.z,pressure_system->surface_patches.size(),(device_fields.bytes()+projection.bytes()+advection.bytes())/(1024.0*1024.0));projection_ok=step_result.converged;}}catch(const std::exception& exception){std::fprintf(stderr,"[paraglider] initial GPU projection failed: %s\n",exception.what());projection_ok=false;}
	std::printf("  limitation: fabric-band fallback and compact EB aperture transport are first-order; EB graph transport uses molecular viscosity but not the Smagorinsky term.\n");
	return unresolved?3:(!composite_ready?4:(projection_ok?0:5));
}
