#include "core/fluid/external_aero_core.h"
#include "core/geometry/step_import.h"
#include "core/geometry/triangle_bvh.h"
#include "core/paraglider_config.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>

using namespace paracfd::core;

namespace
{
	struct FieldMaximum{double value=0;int level=-1,brick=-1,component=-1,i=0,j=0,k=0;Vec3d point{};};
	FieldMaximum locate_regular_maximum(const ExternalAeroCore& core,AmrHostFields& fields)
	{
		core.download_fields(fields);FieldMaximum out;const AmrHierarchy& hierarchy=core.hierarchy();const int bs=hierarchy.brick_size();
		for(int level=0;level<(int)hierarchy.levels().size();++level){const AmrLevel& metadata=hierarchy.levels()[level];const AmrHostLevelFields& values=fields.levels()[level];for(int brick=0;brick<(int)metadata.bricks.size();++brick){const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;for(int component=0;component<3;++component){const int ni=component==0?bs+1:bs,nj=component==1?bs+1:bs,nk=component==2?bs+1:bs;for(int k=0;k<nk;++k)for(int j=0;j<nj;++j)for(int i=0;i<ni;++i){const std::size_t index=component==0?values.layout.u_index(brick,i,j,k):(component==1?values.layout.v_index(brick,i,j,k):values.layout.w_index(brick,i,j,k));const double value=std::abs(static_cast<double>(component==0?values.u[index]:(component==1?values.v[index]:values.w[index])));if(value<=out.value)continue;out={value,level,brick,component,i,j,k,{record.origin.x+(i+(component==0?0.0:0.5))*record.h,record.origin.y+(j+(component==1?0.0:0.5))*record.h,record.origin.z+(k+(component==2?0.0:0.5))*record.h}};}}}}
		return out;
	}

	double elapsed_ms(std::chrono::steady_clock::time_point begin,std::chrono::steady_clock::time_point end)
	{
		return std::chrono::duration<double,std::milli>(end-begin).count();
	}

	void report(int step,const ExternalAeroCore& core,const ExternalAeroStepStats& stats,double rho,AmrHostFields* diagnostic_fields,const TriangleBvh* bvh)
	{
		const AerodynamicLoads loads=core.pressure_loads();const ExternalAeroConservationStats conservation=core.conservation_stats();
		std::printf("[paraglider-case] step=%d t=%.6f s dt=%.3e max|u_i|=%.5g [regular=%.5g special=%.5g CFL-source=%.5g] CFL=%.3f Fp=[%.8g %.8g %.8g] N div[max=%.3e rmsV=%.3e] flux-error[max=%.3e abs=%.3e net=%.3e] PCG=%d/%.3e step=%.2f ms projection=%.2f ms EB=%.2f ms\n",step,core.physical_time(),stats.dt,stats.max_abs_velocity,stats.max_abs_regular_velocity,stats.max_abs_special_velocity,stats.cfl_velocity,stats.effective_cfl,loads.pressure_force.x,loads.pressure_force.y,loads.pressure_force.z,conservation.max_abs_divergence,conservation.volume_weighted_rms_divergence,conservation.max_integrated_flux_error,conservation.absolute_integrated_flux_error,conservation.net_integrated_flux_error,stats.pressure.iterations,stats.pressure.relative_residual,stats.gpu_step_ms,stats.projection_ms,stats.embedded_transport_ms);
		CompositeAmrFluxes fluxes;core.download_special_fluxes(fluxes);double coarse_fine_max=0,embedded_max=0;std::size_t embedded_edge=0;for(double value:fluxes.coarse_fine_velocity)coarse_fine_max=std::max(coarse_fine_max,std::abs(value));for(std::size_t edge=0;edge<fluxes.embedded_velocity.size();++edge)if(std::abs(fluxes.embedded_velocity[edge])>embedded_max){embedded_max=std::abs(fluxes.embedded_velocity[edge]);embedded_edge=edge;}const CompositeAmrPressureSystem& system=core.pressure_system();if(embedded_edge<system.embedded.size()){std::vector<double> pressure;core.download_pressure(pressure);const CoarseFinePressureConnection& connection=system.embedded[embedded_edge];const int lower=connection.direction>0?connection.coarse_dof:connection.fine_dof,upper=connection.direction>0?connection.fine_dof:connection.coarse_dof;const double delta_p=pressure[upper]-pressure[lower],correction=-(stats.dt/rho)*delta_p/connection.centre_distance;std::printf("[paraglider-case]   special maxima: coarse/fine=%.5g EB=%.5g; EB edge=%zu axis=%d A=%.3e d=%.3e Va=%.3e Vb=%.3e dp=%.5g Pa du_p=%.5g m/s centroid=[%.5g %.5g %.5g]\n",coarse_fine_max,embedded_max,embedded_edge,static_cast<int>(connection.axis),connection.open_area,connection.centre_distance,system.volume[connection.coarse_dof],system.volume[connection.fine_dof],delta_p,correction,connection.face_centroid.x,connection.face_centroid.y,connection.face_centroid.z);
			const int bs=system.brick_size,last_level=static_cast<int>(system.hierarchy->levels().size())-1;const int structured_end=system.level_offset[last_level]+static_cast<int>(system.hierarchy->levels()[last_level].bricks.size())*bs*bs*bs;
			auto describe_node=[&](int node,const char* label)
			{
				int degree=0;double area_sum=0,absolute_flux=0,net_flux=0,wall_area=0;for(std::size_t edge=0;edge<system.embedded.size();++edge){const auto& incident=system.embedded[edge];const int lo=incident.direction>0?incident.coarse_dof:incident.fine_dof,hi=incident.direction>0?incident.fine_dof:incident.coarse_dof;if(node!=lo&&node!=hi)continue;const double flux=incident.open_area*fluxes.embedded_velocity[edge];++degree;area_sum+=incident.open_area;absolute_flux+=std::abs(flux);net_flux+=node==lo?flux:-flux;}for(const auto& patch:system.surface_patches)if(patch.plus_dof==node||patch.minus_dof==node)wall_area+=patch.area;const NearestSurfacePoint nearest=bvh?bvh->nearest(system.centroid[node]):NearestSurfacePoint{};std::printf("[paraglider-case]     %s node=%d %s degree=%d A-sum=%.3e |Q|=%.3e Qnet=%.3e wall-A=%.3e centre=[%.5g %.5g %.5g] fabric-distance=%.3e triangle=%u face=%u\n",label,node,node>=structured_end?"fragment":"regular",degree,area_sum,absolute_flux,net_flux,wall_area,system.centroid[node].x,system.centroid[node].y,system.centroid[node].z,nearest.found?nearest.distance:-1.0,nearest.triangle_id,nearest.source_face_id);
			};describe_node(lower,"lower");describe_node(upper,"upper");}
		if(diagnostic_fields&&bvh){const FieldMaximum maximum=locate_regular_maximum(core,*diagnostic_fields);const NearestSurfacePoint nearest=bvh->nearest(maximum.point);std::printf("[paraglider-case]   active regular maximum: %.6g m/s component=%c level=%d brick=%d ijk=[%d %d %d] point=[%.6g %.6g %.6g] fabric-distance=%.6g m triangle=%u\n",maximum.value,maximum.component>=0?"uvw"[maximum.component]:'?',maximum.level,maximum.brick,maximum.i,maximum.j,maximum.k,maximum.point.x,maximum.point.y,maximum.point.z,nearest.found?nearest.distance:-1.0,nearest.triangle_id);}
	}
}

int main(int argc,char** argv)
{
	// Long imported-wing studies are commonly redirected to a background log. Keep
	// each diagnostic sample visible while the case is running instead of buffering
	// the entire multi-minute history until process exit.
	std::setvbuf(stdout,nullptr,_IONBF,0);
	std::string config_path="configs/planb_parakite.json";int steps=50,sample_every=5,max_levels=0;double projection_tolerance=0,min_volume_fraction=0;bool locate_regular_max=false;
	for(int i=1;i<argc;++i)
	{
		const std::string argument=argv[i];
		if(argument=="--config"&&i+1<argc)config_path=argv[++i];
		else if(argument=="--steps"&&i+1<argc)steps=std::atoi(argv[++i]);
		else if(argument=="--sample-every"&&i+1<argc)sample_every=std::atoi(argv[++i]);
		else if(argument=="--projection-tolerance"&&i+1<argc)projection_tolerance=std::atof(argv[++i]);
		else if(argument=="--min-volume-fraction"&&i+1<argc)min_volume_fraction=std::atof(argv[++i]);
		else if(argument=="--max-levels"&&i+1<argc)max_levels=std::atoi(argv[++i]);
		else if(argument=="--locate-regular-max")locate_regular_max=true;
		else{std::fprintf(stderr,"usage: paraglider_case_probe [--config file] [--steps N] [--sample-every N] [--projection-tolerance value] [--min-volume-fraction value] [--max-levels N] [--locate-regular-max]\n");return 2;}
	}
	if(steps<0||sample_every<1){std::fprintf(stderr,"steps must be nonnegative and sample-every must be positive\n");return 2;}
	ParagliderConfig config;std::string error;if(!load_paraglider_config(config_path,config,&error)){std::fprintf(stderr,"[paraglider-case] %s\n",error.c_str());return 2;}if(projection_tolerance>0)config.solver.projection_tolerance=projection_tolerance;if(min_volume_fraction>0)config.amr.min_volume_fraction=min_volume_fraction;if(max_levels>0)config.amr.max_levels=max_levels;
	const auto load_begin=std::chrono::steady_clock::now();TriMesh source=load_step_mesh(config.step_path,config.tessellation_deflection_mm,&error);if(source.empty()){std::fprintf(stderr,"[paraglider-case] STEP load failed: %s\n",error.c_str());return 2;}config.placement=frame_wing_for_external_domain(source,config.placement,config.domain.upstream_margin,config.domain.lateral_margin,config.domain.vertical_margin,config.amr.base_cell_size*config.amr.brick_size);std::fprintf(stderr,"[paraglider-placement] t=[%.17g %.17g %.17g] M=[%.17g %.17g %.17g; %.17g %.17g %.17g; %.17g %.17g %.17g]\n",config.placement.tx,config.placement.ty,config.placement.tz,config.placement.m[0],config.placement.m[1],config.placement.m[2],config.placement.m[3],config.placement.m[4],config.placement.m[5],config.placement.m[6],config.placement.m[7],config.placement.m[8]);TriMesh wing=placed_mesh(source,config.placement);TriangleBvh bvh(wing);const auto load_end=std::chrono::steady_clock::now();
	try
	{
		const auto build_begin=std::chrono::steady_clock::now();ExternalAeroCore core(wing,bvh,config);const auto build_end=std::chrono::steady_clock::now();ExternalAeroStepStats stats=core.initialize();
		std::size_t discarded_apertures=0;double discarded_area=0;for(const AmrEbLevelAtlas& level:core.embedded_boundary().levels){discarded_apertures+=level.topology.discarded_subgrid_apertures;discarded_area+=level.topology.discarded_subgrid_aperture_area;}
		std::printf("[paraglider-case] triangles=%zu bricks=%zu DOFs=%d EB-edges=%zu patches=%zu discarded-apertures=%zu/%.6g-m2 load=%.2f ms preprocess=%.2f ms GPU=%.2f MiB tolerance=%.3e\n",wing.triangle_count(),core.hierarchy().active_brick_count(),core.pressure_system().storage_size,core.pressure_system().embedded.size(),core.pressure_system().surface_patches.size(),discarded_apertures,discarded_area,elapsed_ms(load_begin,load_end),elapsed_ms(build_begin,build_end),core.gpu_bytes()/(1024.0*1024.0),config.solver.projection_tolerance);
		std::unique_ptr<AmrHostFields> diagnostic_fields=locate_regular_max?std::make_unique<AmrHostFields>(core.hierarchy()):nullptr;
		if(!stats.pressure.converged){std::fprintf(stderr,"[paraglider-case] initialization did not converge: iterations=%d residual=%.3e\n",stats.pressure.iterations,stats.pressure.relative_residual);return 3;}report(0,core,stats,config.freestream.rho,diagnostic_fields.get(),&bvh);
		for(int step=1;step<=steps;++step)
		{
			stats=core.step();if(!stats.pressure.converged){std::fprintf(stderr,"[paraglider-case] step %d did not converge: iterations=%d residual=%.3e\n",step,stats.pressure.iterations,stats.pressure.relative_residual);return 4;}
			if(step==1||step==steps||step%sample_every==0)report(step,core,stats,config.freestream.rho,diagnostic_fields.get(),&bvh);
		}
	}
	catch(const std::exception& exception){std::fprintf(stderr,"[paraglider-case] failed: %s\n",exception.what());return 5;}
	return 0;
}
