#include "core/fluid/external_aero_core.h"

#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <vector>

namespace paracfd::core
{
	namespace
	{
		double elapsed_ms(std::chrono::steady_clock::time_point begin,std::chrono::steady_clock::time_point end){return std::chrono::duration<double,std::milli>(end-begin).count();}
	}

	ExternalAeroCore::ExternalAeroCore(const TriMesh& wing,const TriangleBvh& bvh,const ParagliderConfig& config,ExternalAeroExecutionOptions options):config_(config),source_triangle_count_(wing.triangle_count()),use_conservative_cell_momentum_(options.conservative_cell_momentum),enable_smooth_fabric_wall_(options.smooth_fabric_wall),enable_pressure_impulse_(options.pressure_impulse)
	{
		if(wing.empty()||bvh.empty())throw std::invalid_argument("external aerodynamic core requires a placed fabric mesh");hierarchy_=AmrHierarchy::build_static(automatic_flow_domain(wing,config_.domain),wing,bvh,config_.amr);EmbeddedBoundaryBuildOptions eb_options;eb_options.min_volume_fraction=config_.amr.min_volume_fraction;eb_options.min_aperture_area_fraction=config_.amr.min_aperture_area_fraction;eb_options.complex_subdivisions=config_.amr.complex_subdivisions;embedded_boundary_=build_amr_embedded_boundary_atlas(hierarchy_,wing,bvh,eb_options);if(!embedded_boundary_.ready_for_flow())throw std::runtime_error("external aerodynamic core has unresolved finest-level embedded-boundary topology");pressure_system_=build_composite_amr_pressure_system(hierarchy_,embedded_boundary_,true);fields_=std::make_unique<DeviceAmrFields>(hierarchy_);projection_=std::make_unique<DeviceCompositeAmrProjection>(pressure_system_,*fields_);if(use_conservative_cell_momentum_)cell_momentum_=std::make_unique<DeviceCompositeCellMomentumTransport>(pressure_system_,*fields_);else advection_=std::make_unique<DeviceAmrAdvection>(hierarchy_,bvh,2.5,&pressure_system_);
	}

	ExternalAeroCore::~ExternalAeroCore()=default;

	ExternalAeroStepStats ExternalAeroCore::project(bool warm_start,double dt,bool sync_coarse_fine)
	{
		if(!(dt>0))throw std::invalid_argument("external aerodynamic projection timestep");ExternalAeroStepStats stats;stats.dt=dt;stats.conservative_cell_momentum=use_conservative_cell_momentum_;const auto begin=std::chrono::steady_clock::now();if(sync_coarse_fine)projection_->sync_coarse_fine_from_fields();stats.pressure=projection_->project(static_cast<Real>(config_.freestream.rho),static_cast<Real>(stats.dt),config_.solver.projection_tolerance,config_.solver.projection_max_iterations,warm_start);const auto end=std::chrono::steady_clock::now();stats.projection_ms=elapsed_ms(begin,end);stats.gpu_step_ms=stats.projection_ms;stats.physical_time=physical_time_;return stats;
	}

	ExternalAeroStepStats ExternalAeroCore::initialize()
	{
		const Real speed=static_cast<Real>(config_.freestream.speed);fields_->initialize_freestream(speed);fields_->apply_external_aero_boundaries(speed);projection_->initialize_special_freestream(speed);if(cell_momentum_)cell_momentum_->initialize_freestream(speed);const double maximum=std::max(1e-9,std::abs(config_.freestream.speed)),dt=config_.solver.cfl*hierarchy_.finest_cell_size()/maximum;ExternalAeroStepStats stats=project(false,dt);if(cell_momentum_&&enable_pressure_impulse_&&stats.pressure.converged)cell_momentum_->apply_projected_pressure_gradient(projection_->pressure(),static_cast<Real>(dt),static_cast<Real>(config_.freestream.rho));stats.max_abs_regular_velocity=cell_momentum_?cell_momentum_->max_abs_velocity():fields_->max_abs_velocity();stats.max_abs_special_velocity=projection_->max_abs_special_velocity();stats.max_embedded_cfl_rate=projection_->max_embedded_cfl_rate();stats.max_abs_velocity=std::max(stats.max_abs_regular_velocity,stats.max_abs_special_velocity);stats.cfl_velocity=maximum;stats.effective_cfl=maximum*dt/hierarchy_.finest_cell_size();initialized_=stats.pressure.converged;stats.physical_time=physical_time_;return stats;
	}

	ExternalAeroStepStats ExternalAeroCore::step()
	{
		if(!initialized_)throw std::logic_error("external aerodynamic core must have a converged initialization before stepping");
		const double mass_flux_maximum=fields_->max_abs_velocity(),regular_maximum=cell_momentum_?cell_momentum_->max_abs_velocity():mass_flux_maximum,special_maximum=projection_->max_abs_special_velocity(),embedded_cfl_rate=projection_->max_embedded_cfl_rate(),h=hierarchy_.finest_cell_size(),cell_outflow_rate=cell_momentum_?cell_momentum_->max_outflow_rate(projection_->coarse_fine_velocity_device(),projection_->embedded_velocity_device(),true):0;
		if(!std::isfinite(regular_maximum)||!std::isfinite(special_maximum)||!std::isfinite(embedded_cfl_rate))throw std::runtime_error("non-finite regular or embedded-boundary velocity before timestep");
		const double cfl_rate=cell_momentum_?std::max({1e-9,std::abs(config_.freestream.speed)/h,cell_outflow_rate}):std::max({1e-9,std::abs(config_.freestream.speed)/h,mass_flux_maximum/h,embedded_cfl_rate}),dt_value=config_.solver.cfl/cfl_rate,cfl_velocity=cfl_rate*h;const Real dt=static_cast<Real>(dt_value),nu=static_cast<Real>(config_.freestream.nu),cs=static_cast<Real>(config_.solver.smagorinsky_cs),speed=static_cast<Real>(config_.freestream.speed),rho=static_cast<Real>(config_.freestream.rho);
		const auto step_begin=std::chrono::steady_clock::now(),advection_begin=step_begin;ExternalAeroStepStats stats;
		auto advection_end=step_begin,turbulence_end=step_begin,embedded_transport_end=step_begin;
		if(cell_momentum_)
		{
			cell_momentum_->step(projection_->coarse_fine_velocity_device(),projection_->embedded_velocity_device(),dt,true,speed);advection_end=std::chrono::steady_clock::now();cell_momentum_->diffuse_smagorinsky(nu,cs,dt);turbulence_end=std::chrono::steady_clock::now();if(enable_smooth_fabric_wall_){cell_momentum_->apply_smooth_fabric_wall_model(dt,nu);smooth_fabric_wall_applied_=true;}cell_momentum_->reconstruct_flux_increments(projection_->coarse_fine_velocity_device(),projection_->embedded_velocity_device());fields_->apply_external_aero_boundaries(speed);embedded_transport_end=std::chrono::steady_clock::now();stats=project(true,dt_value,false);if(enable_pressure_impulse_&&stats.pressure.converged)cell_momentum_->apply_projected_pressure_gradient(projection_->pressure(),dt,rho);
		}
		else
		{
			advection_->advect(*fields_,dt);advection_end=std::chrono::steady_clock::now();advection_->diffuse_smagorinsky(*fields_,nu,cs,dt);turbulence_end=std::chrono::steady_clock::now();fields_->apply_external_aero_boundaries(speed);projection_->transport_embedded_apertures(dt,nu,cs,false);projection_->apply_smooth_fabric_wall_model(dt,nu);smooth_fabric_wall_applied_=true;embedded_transport_end=std::chrono::steady_clock::now();stats=project(true,dt_value);
		}
		const auto step_end=std::chrono::steady_clock::now();
		stats.max_abs_regular_velocity=regular_maximum;stats.max_abs_special_velocity=special_maximum;stats.max_embedded_cfl_rate=embedded_cfl_rate;stats.max_abs_velocity=std::max(regular_maximum,special_maximum);stats.cfl_velocity=cfl_velocity;stats.effective_cfl=cfl_rate*dt_value;stats.advection_ms=elapsed_ms(advection_begin,advection_end);stats.turbulence_ms=elapsed_ms(advection_end,turbulence_end);stats.embedded_transport_ms=elapsed_ms(turbulence_end,embedded_transport_end);stats.les_applied=true;stats.embedded_transport_applied=true;stats.smooth_fabric_wall_applied=smooth_fabric_wall_applied_;stats.gpu_step_ms=elapsed_ms(step_begin,step_end);if(stats.pressure.converged)physical_time_+=stats.dt;stats.physical_time=physical_time_;return stats;
	}

	AerodynamicLoads ExternalAeroCore::pressure_loads(double reference_pressure) const
	{
		std::vector<Real> device_pressure;projection_->download_pressure(device_pressure);std::vector<double> pressure(device_pressure.size());for(std::size_t q=0;q<pressure.size();++q)pressure[q]=static_cast<double>(device_pressure[q]);return compute_pressure_loads(pressure_system_,pressure,source_triangle_count_,config_.freestream,config_.reference,reference_pressure);
	}
	AerodynamicLoads ExternalAeroCore::aerodynamic_loads(double reference_pressure) const
	{
		AerodynamicLoads loads=pressure_loads(reference_pressure);if(smooth_fabric_wall_applied_){std::vector<SmoothFabricWallPatchLoad> viscous;if(cell_momentum_)cell_momentum_->download_smooth_fabric_wall_loads(static_cast<Real>(config_.freestream.rho),viscous);else projection_->download_smooth_fabric_wall_loads(static_cast<Real>(config_.freestream.rho),viscous);accumulate_viscous_loads(loads,viscous,config_.freestream,config_.reference);}return loads;
	}

	void ExternalAeroCore::download_fields(AmrHostFields& host) const{fields_->download(host);}
	void ExternalAeroCore::download_special_fluxes(CompositeAmrFluxes& host) const{projection_->download_special_fluxes(host);}
	void ExternalAeroCore::download_pressure(std::vector<double>& host) const{std::vector<Real> values;projection_->download_pressure(values);host.assign(values.begin(),values.end());}
	ExternalAeroConservationStats ExternalAeroCore::conservation_stats() const
	{
		std::vector<Real> divergence;projection_->download_divergence(divergence);ExternalAeroConservationStats out;double weighted_square=0,total_volume=0;
		for(std::size_t q=0;q<divergence.size()&&q<pressure_system_.volume.size();++q)if(pressure_system_.active[q]&&pressure_system_.volume[q]>0)
		{
			const double value=static_cast<double>(divergence[q]),volume=pressure_system_.volume[q],absolute=std::abs(value),integrated=value*volume;
			if(absolute>out.max_abs_divergence){out.max_abs_divergence=absolute;out.worst_control_volume=volume;}out.max_integrated_flux_error=std::max(out.max_integrated_flux_error,std::abs(integrated));out.absolute_integrated_flux_error+=std::abs(integrated);out.net_integrated_flux_error+=integrated;weighted_square+=value*value*volume;total_volume+=volume;
		}
		out.volume_weighted_rms_divergence=total_volume>0?std::sqrt(weighted_square/total_volume):0;return out;
	}
	double ExternalAeroCore::max_abs_divergence() const{return conservation_stats().max_abs_divergence;}
	std::size_t ExternalAeroCore::gpu_bytes() const{return fields_->bytes()+projection_->bytes()+(advection_?advection_->bytes():0)+(cell_momentum_?cell_momentum_->bytes():0);}
	std::size_t ExternalAeroCore::protected_face_count() const{return advection_?advection_->protected_face_count():0;}
	std::size_t ExternalAeroCore::active_face_count() const{return advection_?advection_->active_face_count():0;}
	int ExternalAeroCore::embedded_high_order_stencil_count() const{return projection_->embedded_high_order_stencil_count();}
	int ExternalAeroCore::embedded_least_squares_full_rank_count() const{return projection_->embedded_least_squares_full_rank_count();}
	int ExternalAeroCore::fabric_wall_node_count() const{return projection_->fabric_wall_node_count();}
	int ExternalAeroCore::clamped_cell_flux_interpolation_count() const{return cell_momentum_?cell_momentum_->clamped_flux_interpolation_count():0;}
	int ExternalAeroCore::pressure_closure_correction_count() const{return cell_momentum_?cell_momentum_->pressure_closure_correction_count():0;}
	double ExternalAeroCore::max_pressure_closure_acceleration() const{return cell_momentum_?cell_momentum_->max_pressure_closure_acceleration():0;}
}
