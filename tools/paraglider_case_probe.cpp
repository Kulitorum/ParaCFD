#include "core/fluid/external_aero_core.h"
#include "core/geometry/mesh_clip.h"
#include "core/geometry/step_import.h"
#include "core/geometry/triangle_bvh.h"
#include "core/paraglider_config.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <deque>
#include <string>
#include <vector>

using namespace paracfd::core;

namespace
{
	ModelPlacement pitch_placement(ModelPlacement source,double degrees)
	{
		const double radians=degrees*3.14159265358979323846/180.0;
		const double c=std::cos(radians),s=std::sin(radians);
		const double rotation[9]={c,0,s,0,1,0,-s,0,c};
		ModelPlacement out=source;
		for(int row=0;row<3;++row)for(int column=0;column<3;++column)
			out.m[3*row+column]=rotation[3*row]*source.m[column]+rotation[3*row+1]*source.m[3+column]+rotation[3*row+2]*source.m[6+column];
		out.tx=out.ty=out.tz=0;
		return out;
	}

	struct FieldMaximum{double value=0;int level=-1,brick=-1,component=-1,i=0,j=0,k=0;Vec3d point{};};
	FieldMaximum locate_regular_maximum(const ExternalAeroCore& core,AmrHostFields& fields)
	{
		core.download_fields(fields);FieldMaximum out;const AmrHierarchy& hierarchy=core.hierarchy();const int bs=hierarchy.brick_size();
		for(int level=0;level<(int)hierarchy.levels().size();++level){const AmrLevel& metadata=hierarchy.levels()[level];const AmrHostLevelFields& values=fields.levels()[level];for(int brick=0;brick<(int)metadata.bricks.size();++brick){const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;for(int component=0;component<3;++component){const int ni=component==0?bs+1:bs,nj=component==1?bs+1:bs,nk=component==2?bs+1:bs;for(int k=0;k<nk;++k)for(int j=0;j<nj;++j)for(int i=0;i<ni;++i){const std::size_t index=component==0?values.layout.u_index(brick,i,j,k):(component==1?values.layout.v_index(brick,i,j,k):values.layout.w_index(brick,i,j,k));const double value=std::abs(static_cast<double>(component==0?values.u[index]:(component==1?values.v[index]:values.w[index])));if(value<=out.value)continue;out={value,level,brick,component,i,j,k,{record.origin.x+(i+(component==0?0.0:0.5))*record.h,record.origin.y+(j+(component==1?0.0:0.5))*record.h,record.origin.z+(k+(component==2?0.0:0.5))*record.h}};}}}}
		return out;
	}
	struct DirectionAccumulator
	{
		double volume=0,sum_u=0,sum_v=0,sum_w=0,sum_speed2=0,sum_transverse2=0,reverse_volume=0,cross_dominant_volume=0;
		double min_u=std::numeric_limits<double>::infinity(),max_u=-std::numeric_limits<double>::infinity();
		double max_abs_v=0,max_abs_w=0;
		void add(double u,double v,double w,double cell_volume)
		{
			if(!std::isfinite(u)||!std::isfinite(v)||!std::isfinite(w))return;const double transverse2=v*v+w*w;volume+=cell_volume;sum_u+=cell_volume*u;sum_v+=cell_volume*v;sum_w+=cell_volume*w;sum_speed2+=cell_volume*(u*u+transverse2);sum_transverse2+=cell_volume*transverse2;if(u<0)reverse_volume+=cell_volume;if(transverse2>u*u)cross_dominant_volume+=cell_volume;min_u=std::min(min_u,u);max_u=std::max(max_u,u);max_abs_v=std::max(max_abs_v,std::abs(v));max_abs_w=std::max(max_abs_w,std::abs(w));
		}
	};
	struct DirectionalStats{DirectionAccumulator all,background;};
	DirectionalStats directional_stats(const ExternalAeroCore& core,const AmrHostFields& fields,const Aabb3d& wing_box)
	{
		DirectionalStats out;const AmrHierarchy& hierarchy=core.hierarchy();const int bs=hierarchy.brick_size();const double pad=2.0;
		for(int level=0;level<(int)hierarchy.levels().size();++level){const AmrLevel& metadata=hierarchy.levels()[level];const AmrHostLevelFields& values=fields.levels()[level];const double cell_volume=metadata.h*metadata.h*metadata.h;for(int brick=0;brick<(int)metadata.bricks.size();++brick){const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const double u=0.5*(static_cast<double>(values.u[values.layout.u_index(brick,i,j,k)])+static_cast<double>(values.u[values.layout.u_index(brick,i+1,j,k)])),v=0.5*(static_cast<double>(values.v[values.layout.v_index(brick,i,j,k)])+static_cast<double>(values.v[values.layout.v_index(brick,i,j+1,k)])),w=0.5*(static_cast<double>(values.w[values.layout.w_index(brick,i,j,k)])+static_cast<double>(values.w[values.layout.w_index(brick,i,j,k+1)]));const Vec3d point=record.origin+Vec3d{(i+0.5)*metadata.h,(j+0.5)*metadata.h,(k+0.5)*metadata.h};out.all.add(u,v,w,cell_volume);const bool background=point.x<wing_box.lo.x-pad||point.y<wing_box.lo.y-pad||point.y>wing_box.hi.y+pad||point.z<wing_box.lo.z-pad||point.z>wing_box.hi.z+pad;if(background)out.background.add(u,v,w,cell_volume);}}}
		return out;
	}
	void print_direction_region(const char* label,const DirectionAccumulator& value)
	{
		if(!(value.volume>0)){std::printf("[paraglider-case]   %s direction: no represented volume\n",label);return;}const double inv=1.0/value.volume;std::printf("[paraglider-case]   %s direction: mean=[%.5g %.5g %.5g] m/s rms|u|=%.5g rms-trans=%.5g u-range=[%.5g %.5g] max|v,w|=[%.5g %.5g] reverse=%.4g%% transverse-dominant=%.4g%%\n",label,value.sum_u*inv,value.sum_v*inv,value.sum_w*inv,std::sqrt(value.sum_speed2*inv),std::sqrt(value.sum_transverse2*inv),value.min_u,value.max_u,value.max_abs_v,value.max_abs_w,100.0*value.reverse_volume*inv,100.0*value.cross_dominant_volume*inv);
	}

	struct SurfaceTemporalDiagnostic
	{
		struct Frame{std::vector<double> cp,area;double range=1;};
		std::deque<Frame> history;

		void sample(int step,double time,const TriMesh& wing,const AerodynamicLoads& loads)
		{
			Frame current;current.cp.resize(loads.triangles.size(),std::numeric_limits<double>::quiet_NaN());current.area.resize(loads.triangles.size());current.range=1e-12;
			for(std::size_t triangle=0;triangle<loads.triangles.size();++triangle)
			{
				const auto& value=loads.triangles[triangle];current.cp[triangle]=value.delta_cp;current.area[triangle]=value.represented_area;if(std::isfinite(value.delta_cp))current.range=std::max(current.range,std::abs(value.delta_cp));
			}
			auto compare=[&](const Frame& previous,std::size_t* worst_triangle=nullptr,double* worst_previous=nullptr,double* worst_current=nullptr)
			{
				double weighted_square=0,weight_sum=0,worst=0;
				for(std::size_t triangle=0;triangle<current.cp.size()&&triangle<previous.cp.size();++triangle)
				{
					if(!std::isfinite(current.cp[triangle])||!std::isfinite(previous.cp[triangle]))continue;const double weight=std::max(current.area[triangle],previous.area[triangle]);if(!(weight>0))continue;const double difference=current.cp[triangle]/current.range-previous.cp[triangle]/previous.range;weighted_square+=weight*difference*difference;weight_sum+=weight;if(std::abs(difference)>worst){worst=std::abs(difference);if(worst_triangle)*worst_triangle=triangle;if(worst_previous)*worst_previous=previous.cp[triangle];if(worst_current)*worst_current=current.cp[triangle];}
				}
				return weight_sum>0?std::sqrt(weighted_square/weight_sum):std::numeric_limits<double>::quiet_NaN();
			};
			double lag[4]={NAN,NAN,NAN,NAN};for(std::size_t q=0;q<history.size()&&q<4;++q)lag[q]=compare(history[history.size()-1-q]);
			std::size_t worst_triangle=0;double previous_cp=NAN,current_cp=NAN,worst_change=NAN;if(!history.empty())worst_change=compare(history.back(),&worst_triangle,&previous_cp,&current_cp);
			Vec3d centroid{};if(worst_triangle<wing.triangle_count()){for(int corner=0;corner<3;++corner){const std::uint32_t vertex=wing.indices[3*worst_triangle+corner];centroid.x+=wing.positions[3*vertex];centroid.y+=wing.positions[3*vertex+1];centroid.z+=wing.positions[3*vertex+2];}centroid=centroid/3.0;}
			std::printf("[paraglider-surface-time] step=%d t=%.9g range=%.6g normalized-rms-lag[1,2,3,4]=[%.6g %.6g %.6g %.6g] worst-triangle=%zu Cp=%.6g->%.6g centre=[%.6g %.6g %.6g]\n",step,time,current.range,lag[0],lag[1],lag[2],lag[3],worst_triangle,previous_cp,current_cp,centroid.x,centroid.y,centroid.z);
			history.push_back(std::move(current));while(history.size()>4)history.pop_front();
		}
	};

	double elapsed_ms(std::chrono::steady_clock::time_point begin,std::chrono::steady_clock::time_point end)
	{
		return std::chrono::duration<double,std::milli>(end-begin).count();
	}

	void report(int step,const ExternalAeroCore& core,const ExternalAeroStepStats& stats,double rho,AmrHostFields* diagnostic_fields,const TriangleBvh* bvh,const Aabb3d* wing_box,SurfaceTemporalDiagnostic* surface_temporal=nullptr,const TriMesh* wing=nullptr)
	{
		const AerodynamicLoads loads=core.aerodynamic_loads();const ExternalAeroConservationStats conservation=core.conservation_stats();
		if(step==0)std::printf("[paraglider-case] phase=initialization-impulse (pressure/load is not an evolved aerodynamic result)\n");
		if(surface_temporal&&wing)surface_temporal->sample(step,core.physical_time(),*wing,loads);
		std::printf("[paraglider-case] step=%d t=%.6f s dt=%.3e max|u_i|=%.5g [regular=%.5g special=%.5g CFL-source=%.5g] CFL=%.3f Fp=[%.8g %.8g %.8g] Fv=[%.8g %.8g %.8g] N div[max=%.3e rmsV=%.3e] flux-error[max=%.3e abs=%.3e net=%.3e] PCG=%d/%.3e step=%.2f ms projection=%.2f ms EB=%.2f ms\n",step,core.physical_time(),stats.dt,stats.max_abs_velocity,stats.max_abs_regular_velocity,stats.max_abs_special_velocity,stats.cfl_velocity,stats.effective_cfl,loads.pressure_force.x,loads.pressure_force.y,loads.pressure_force.z,loads.viscous_force.x,loads.viscous_force.y,loads.viscous_force.z,conservation.max_abs_divergence,conservation.volume_weighted_rms_divergence,conservation.max_integrated_flux_error,conservation.absolute_integrated_flux_error,conservation.net_integrated_flux_error,stats.pressure.iterations,stats.pressure.relative_residual,stats.gpu_step_ms,stats.projection_ms,stats.embedded_transport_ms);
		CompositeAmrFluxes fluxes;core.download_special_fluxes(fluxes);double coarse_fine_max=0,embedded_max=0;std::size_t embedded_edge=0;for(double value:fluxes.coarse_fine_velocity)coarse_fine_max=std::max(coarse_fine_max,std::abs(value));for(std::size_t edge=0;edge<fluxes.embedded_velocity.size();++edge)if(std::abs(fluxes.embedded_velocity[edge])>embedded_max){embedded_max=std::abs(fluxes.embedded_velocity[edge]);embedded_edge=edge;}const CompositeAmrPressureSystem& system=core.pressure_system();if(embedded_edge<system.embedded.size()){std::vector<double> pressure;core.download_pressure(pressure);const CoarseFinePressureConnection& connection=system.embedded[embedded_edge];const int lower=connection.direction>0?connection.coarse_dof:connection.fine_dof,upper=connection.direction>0?connection.fine_dof:connection.coarse_dof;const double delta_p=pressure[upper]-pressure[lower],correction=-(stats.dt/rho)*delta_p*pressure_gradient_factor(connection);std::printf("[paraglider-case]   special maxima: coarse/fine=%.5g EB=%.5g; EB edge=%zu axis=%d A=%.3e d=%.3e dn=%.3e Va=%.3e Vb=%.3e dp=%.5g Pa du_p=%.5g m/s centroid=[%.5g %.5g %.5g]\n",coarse_fine_max,embedded_max,embedded_edge,static_cast<int>(connection.axis),connection.open_area,connection.centre_distance,connection.normal_distance,system.volume[connection.coarse_dof],system.volume[connection.fine_dof],delta_p,correction,connection.face_centroid.x,connection.face_centroid.y,connection.face_centroid.z);
			const int bs=system.brick_size,last_level=static_cast<int>(system.hierarchy->levels().size())-1;const int structured_end=system.level_offset[last_level]+static_cast<int>(system.hierarchy->levels()[last_level].bricks.size())*bs*bs*bs;
			auto describe_node=[&](int node,const char* label)
			{
				int degree=0;double area_sum=0,absolute_flux=0,net_flux=0,wall_area=0;for(std::size_t edge=0;edge<system.embedded.size();++edge){const auto& incident=system.embedded[edge];const int lo=incident.direction>0?incident.coarse_dof:incident.fine_dof,hi=incident.direction>0?incident.fine_dof:incident.coarse_dof;if(node!=lo&&node!=hi)continue;const double flux=incident.open_area*fluxes.embedded_velocity[edge];++degree;area_sum+=incident.open_area;absolute_flux+=std::abs(flux);net_flux+=node==lo?flux:-flux;}for(const auto& patch:system.surface_patches)if(patch.plus_dof==node||patch.minus_dof==node)wall_area+=patch.area;const NearestSurfacePoint nearest=bvh?bvh->nearest(system.centroid[node]):NearestSurfacePoint{};std::printf("[paraglider-case]     %s node=%d %s degree=%d A-sum=%.3e |Q|=%.3e Qnet=%.3e wall-A=%.3e centre=[%.5g %.5g %.5g] fabric-distance=%.3e triangle=%u face=%u\n",label,node,node>=structured_end?"fragment":"regular",degree,area_sum,absolute_flux,net_flux,wall_area,system.centroid[node].x,system.centroid[node].y,system.centroid[node].z,nearest.found?nearest.distance:-1.0,nearest.triangle_id,nearest.source_face_id);
			};describe_node(lower,"lower");describe_node(upper,"upper");}
		if(diagnostic_fields&&bvh){const FieldMaximum maximum=locate_regular_maximum(core,*diagnostic_fields);const NearestSurfacePoint nearest=bvh->nearest(maximum.point);std::printf("[paraglider-case]   active regular maximum: %.6g m/s component=%c level=%d brick=%d ijk=[%d %d %d] point=[%.6g %.6g %.6g] fabric-distance=%.6g m triangle=%u\n",maximum.value,maximum.component>=0?"uvw"[maximum.component]:'?',maximum.level,maximum.brick,maximum.i,maximum.j,maximum.k,maximum.point.x,maximum.point.y,maximum.point.z,nearest.found?nearest.distance:-1.0,nearest.triangle_id);if(wing_box){const DirectionalStats direction=directional_stats(core,*diagnostic_fields,*wing_box);print_direction_region("active",direction.all);print_direction_region("background",direction.background);}}
	}
}

int main(int argc,char** argv)
{
	// Long imported-wing studies are commonly redirected to a background log. Keep
	// each diagnostic sample visible while the case is running instead of buffering
	// the entire multi-minute history until process exit.
	std::setvbuf(stdout,nullptr,_IONBF,0);
	std::string config_path="configs/planb_parakite.json",step_path_override;int steps=50,sample_every=5,max_levels=0;double projection_tolerance=0,min_volume_fraction=0,smagorinsky_cs=-1,target_physical_time=0,thin_y_fraction=-1,thin_y_width=0.125,aoa_degrees=0;bool locate_regular_max=false,surface_temporal=false,conservative_cell_momentum=true,no_smooth_wall=false,no_pressure_impulse=false,steps_explicit=false,thin_y_width_explicit=false,half_wing_override=false;
	for(int i=1;i<argc;++i)
	{
		const std::string argument=argv[i];
		if(argument=="--config"&&i+1<argc)config_path=argv[++i];
		else if(argument=="--step"&&i+1<argc)step_path_override=argv[++i];
		else if(argument=="--steps"&&i+1<argc){steps=std::atoi(argv[++i]);steps_explicit=true;}
		else if(argument=="--physical-time"&&i+1<argc)target_physical_time=std::atof(argv[++i]);
		else if(argument=="--sample-every"&&i+1<argc)sample_every=std::atoi(argv[++i]);
		else if(argument=="--projection-tolerance"&&i+1<argc)projection_tolerance=std::atof(argv[++i]);
		else if(argument=="--min-volume-fraction"&&i+1<argc)min_volume_fraction=std::atof(argv[++i]);
		else if(argument=="--smagorinsky-cs"&&i+1<argc)smagorinsky_cs=std::atof(argv[++i]);
		else if(argument=="--max-levels"&&i+1<argc)max_levels=std::atoi(argv[++i]);
		else if(argument=="--aoa"&&i+1<argc)aoa_degrees=std::atof(argv[++i]);
		else if(argument=="--thin-y-fraction"&&i+1<argc)thin_y_fraction=std::atof(argv[++i]);
		else if(argument=="--thin-y-width"&&i+1<argc){thin_y_width=std::atof(argv[++i]);thin_y_width_explicit=true;}
		else if(argument=="--half-wing")half_wing_override=true;
		else if(argument=="--locate-regular-max")locate_regular_max=true;
		else if(argument=="--surface-temporal")surface_temporal=true;
		else if(argument=="--conservative-cell-momentum")conservative_cell_momentum=true;
		else if(argument=="--staggered-momentum")conservative_cell_momentum=false;
		else if(argument=="--no-smooth-wall")no_smooth_wall=true;
		else if(argument=="--no-pressure-impulse")no_pressure_impulse=true;
		else{std::fprintf(stderr,"usage: paraglider_case_probe [--config file] [--step file] [--steps N] [--physical-time seconds] [--sample-every N] [--projection-tolerance value] [--min-volume-fraction value] [--smagorinsky-cs value] [--max-levels N] [--aoa degrees] [--half-wing] [--thin-y-fraction 0..1] [--thin-y-width metres] [--locate-regular-max] [--surface-temporal] [--conservative-cell-momentum|--staggered-momentum] [--no-smooth-wall] [--no-pressure-impulse]\n");return 2;}
	}
	if(thin_y_width_explicit&&thin_y_fraction<0)thin_y_fraction=0.5;if(steps<0||sample_every<1||target_physical_time<0||(thin_y_fraction>=0&&(thin_y_fraction<0.05||thin_y_fraction>0.95))||!(thin_y_width>0)){std::fprintf(stderr,"steps and physical-time must be nonnegative, sample-every positive, thin-Y fraction in [0.05,0.95], and thin-Y width positive\n");return 2;}if(target_physical_time>0&&!steps_explicit)steps=std::numeric_limits<int>::max();
	ParagliderConfig config;std::string error;if(!load_paraglider_config(config_path,config,&error)){std::fprintf(stderr,"[paraglider-case] %s\n",error.c_str());return 2;}if(!step_path_override.empty())config.step_path=step_path_override;if(projection_tolerance>0)config.solver.projection_tolerance=projection_tolerance;if(min_volume_fraction>0)config.amr.min_volume_fraction=min_volume_fraction;if(smagorinsky_cs>=0)config.solver.smagorinsky_cs=smagorinsky_cs;if(max_levels>0)config.amr.max_levels=max_levels;if(aoa_degrees!=0)config.placement=pitch_placement(config.placement,aoa_degrees);if(half_wing_override)config.domain.half_wing_symmetry=true;if(config.domain.half_wing_symmetry&&thin_y_fraction>=0){std::fprintf(stderr,"[paraglider-case] --half-wing and cropped-Y diagnosis are mutually exclusive\n");return 2;}
	const auto load_begin=std::chrono::steady_clock::now();TriMesh source=load_step_mesh(config.step_path,config.tessellation_deflection_mm,&error);if(source.empty()){std::fprintf(stderr,"[paraglider-case] STEP load failed: %s\n",error.c_str());return 2;}TriMesh wing;if(thin_y_fraction>=0){const int ratio=1<<std::max(0,config.amr.max_levels-1);const double h=config.amr.base_cell_size/ratio;const int layers=std::max(2,static_cast<int>(std::ceil(thin_y_width/h-1e-9)));if(layers>128){std::fprintf(stderr,"[paraglider-case] cropped-Y width requires %d cells; maximum is 128\n",layers);return 2;}const double width=layers*h;ModelPlacement orientation=config.placement;orientation.tx=orientation.ty=orientation.tz=0;const TriMesh oriented=placed_mesh(source,orientation);const double centre=oriented.bbox_min[1]+thin_y_fraction*(oriented.bbox_max[1]-oriented.bbox_min[1]);TriMesh clipped=clip_mesh_to_axis_slab(oriented,1,centre-width*0.5,centre+width*0.5);if(clipped.empty()){std::fprintf(stderr,"[paraglider-case] cropped-Y slab contains no triangles\n");return 2;}config.amr.base_cell_size=h;config.amr.max_levels=1;config.amr.brick_size=layers;config.domain.lateral_margin=0;config.reference.area=0;config.reference.length=0;config.placement=frame_wing_for_external_domain(clipped,ModelPlacement{},config.domain.upstream_margin,0,config.domain.vertical_margin,width);wing=placed_mesh(clipped,config.placement);std::fprintf(stderr,"[paraglider-thin-y] station=%.6g source-y=%.6g width=%.6g m h=%.6g layers=%d triangles=%zu\n",thin_y_fraction,centre,width,h,layers,wing.triangle_count());}else if(config.domain.half_wing_symmetry){const ModelPlacement full_frame=frame_wing_for_external_domain(source,config.placement,config.domain.upstream_margin,config.domain.lateral_margin,config.domain.vertical_margin,config.amr.base_cell_size*config.amr.brick_size);ModelPlacement orientation=config.placement;orientation.tx=orientation.ty=orientation.tz=0;const TriMesh oriented=placed_mesh(source,orientation);const double centre=0.5*(oriented.bbox_min[1]+oriented.bbox_max[1]);TriMesh clipped=clip_mesh_to_axis_slab(oriented,1,centre,oriented.bbox_max[1]);if(clipped.empty()){std::fprintf(stderr,"[paraglider-case] half-wing crop contains no triangles\n");return 2;}config.placement=frame_positive_y_half_for_external_domain(clipped,ModelPlacement{},config.domain.upstream_margin,config.domain.vertical_margin,config.amr.base_cell_size*config.amr.brick_size);config.reference.moment_origin=config.reference.moment_origin+Vec3d{config.placement.tx-full_frame.tx,config.placement.ty-full_frame.ty,config.placement.tz-full_frame.tz};wing=placed_mesh(clipped,config.placement);std::fprintf(stderr,"[paraglider-half-wing] centre=%.9g retained=%zu/%zu plane-y=%.9g tip-y=%.9g moment-origin=[%.9g %.9g %.9g]\n",centre,wing.triangle_count(),source.triangle_count(),wing.bbox_min[1],wing.bbox_max[1],config.reference.moment_origin.x,config.reference.moment_origin.y,config.reference.moment_origin.z);}else{config.placement=frame_wing_for_external_domain(source,config.placement,config.domain.upstream_margin,config.domain.lateral_margin,config.domain.vertical_margin,config.amr.base_cell_size*config.amr.brick_size);wing=placed_mesh(source,config.placement);}std::fprintf(stderr,"[paraglider-placement] t=[%.17g %.17g %.17g] M=[%.17g %.17g %.17g; %.17g %.17g %.17g; %.17g %.17g %.17g]\n",config.placement.tx,config.placement.ty,config.placement.tz,config.placement.m[0],config.placement.m[1],config.placement.m[2],config.placement.m[3],config.placement.m[4],config.placement.m[5],config.placement.m[6],config.placement.m[7],config.placement.m[8]);TriangleBvh bvh(wing);const auto load_end=std::chrono::steady_clock::now();
	try
	{
		const auto build_begin=std::chrono::steady_clock::now();ExternalAeroExecutionOptions execution;execution.conservative_cell_momentum=conservative_cell_momentum;execution.smooth_fabric_wall=!no_smooth_wall;execution.pressure_impulse=!no_pressure_impulse;ExternalAeroCore core(wing,bvh,config,execution);const auto build_end=std::chrono::steady_clock::now();ExternalAeroStepStats stats=core.initialize();
		std::size_t discarded_apertures=0;double discarded_area=0;for(const AmrEbLevelAtlas& level:core.embedded_boundary().levels){discarded_apertures+=level.topology.discarded_subgrid_apertures;discarded_area+=level.topology.discarded_subgrid_aperture_area;}
		std::printf("[paraglider-case] momentum=%s triangles=%zu bricks=%zu DOFs=%d EB-edges=%zu MUSCL=%d LS-full=%d wall-nodes=%d patches=%zu gauges=%zu clamped-face-interpolation=%d raw-traction-closure=%d/max-%.6g-1m discarded-apertures=%zu/%.6g-m2 load=%.2f ms preprocess=%.2f ms GPU=%.2f MiB tolerance=%.3e\n",core.uses_conservative_cell_momentum()?"conservative-cell":"staggered-reference",wing.triangle_count(),core.hierarchy().active_brick_count(),core.pressure_system().storage_size,core.pressure_system().embedded.size(),core.embedded_high_order_stencil_count(),core.embedded_least_squares_full_rank_count(),core.fabric_wall_node_count(),core.pressure_system().surface_patches.size(),core.pressure_system().gauges.size(),core.clamped_cell_flux_interpolation_count(),core.pressure_closure_correction_count(),core.max_pressure_closure_acceleration(),discarded_apertures,discarded_area,elapsed_ms(load_begin,load_end),elapsed_ms(build_begin,build_end),core.gpu_bytes()/(1024.0*1024.0),config.solver.projection_tolerance);
		std::unique_ptr<AmrHostFields> diagnostic_fields=locate_regular_max?std::make_unique<AmrHostFields>(core.hierarchy()):nullptr;
		SurfaceTemporalDiagnostic surface_history;
		const Aabb3d wing_box{{wing.bbox_min[0],wing.bbox_min[1],wing.bbox_min[2]},{wing.bbox_max[0],wing.bbox_max[1],wing.bbox_max[2]}};
		if(!stats.pressure.converged){std::fprintf(stderr,"[paraglider-case] initialization did not converge: iterations=%d residual=%.3e\n",stats.pressure.iterations,stats.pressure.relative_residual);return 3;}report(0,core,stats,config.freestream.rho,diagnostic_fields.get(),&bvh,&wing_box,surface_temporal?&surface_history:nullptr,&wing);
		int completed_steps=0;for(int step=1;step<=steps&&(!(target_physical_time>0)||core.physical_time()<target_physical_time);++step)
		{
			stats=core.step();completed_steps=step;if(!stats.pressure.converged){std::fprintf(stderr,"[paraglider-case] step %d did not converge: iterations=%d residual=%.3e\n",step,stats.pressure.iterations,stats.pressure.relative_residual);return 4;}const bool reached_time=target_physical_time>0&&core.physical_time()>=target_physical_time;
			if(step==1||step==steps||reached_time||step%sample_every==0)report(step,core,stats,config.freestream.rho,diagnostic_fields.get(),&bvh,&wing_box,surface_temporal?&surface_history:nullptr,&wing);
		}
		if(target_physical_time>0&&core.physical_time()<target_physical_time){std::fprintf(stderr,"[paraglider-case] physical-time target %.9g s not reached within %d steps (t=%.9g s)\n",target_physical_time,steps,core.physical_time());return 6;}const AerodynamicLoads summary_loads=core.aerodynamic_loads();const Vec3d summary_force=summary_loads.viscous_loads_valid?summary_loads.total_force:summary_loads.pressure_force;const ExternalAeroConservationStats summary_conservation=core.conservation_stats();const double summary_cd=summary_loads.force_coefficients_valid?summary_loads.cd:std::numeric_limits<double>::quiet_NaN(),summary_cl=summary_loads.force_coefficients_valid?summary_loads.cl:std::numeric_limits<double>::quiet_NaN();std::printf("[paraglider-case-summary] levels=%zu finest_h=%.9g steps=%d t=%.9g prestep_regular_max=%.9g prestep_special_max=%.9g Fx=%.9g Fy=%.9g Fz=%.9g Fpx=%.9g Fpy=%.9g Fpz=%.9g Fvx=%.9g Fvy=%.9g Fvz=%.9g Cd=%.9g Cl=%.9g div_max=%.9g div_rms=%.9g flux_net=%.9g gpu_mib=%.9g step_ms=%.9g projection_ms=%.9g pressure_iterations=%d pressure_residual=%.9g\n",core.hierarchy().levels().size(),core.hierarchy().finest_cell_size(),completed_steps,core.physical_time(),stats.max_abs_regular_velocity,stats.max_abs_special_velocity,summary_force.x,summary_force.y,summary_force.z,summary_loads.pressure_force.x,summary_loads.pressure_force.y,summary_loads.pressure_force.z,summary_loads.viscous_force.x,summary_loads.viscous_force.y,summary_loads.viscous_force.z,summary_cd,summary_cl,summary_conservation.max_abs_divergence,summary_conservation.volume_weighted_rms_divergence,summary_conservation.net_integrated_flux_error,core.gpu_bytes()/(1024.0*1024.0),stats.gpu_step_ms,stats.projection_ms,stats.pressure.iterations,stats.pressure.relative_residual);
	}
	catch(const std::exception& exception){std::fprintf(stderr,"[paraglider-case] failed: %s\n",exception.what());return 5;}
	return 0;
}
