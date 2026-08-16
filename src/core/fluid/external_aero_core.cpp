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

	ExternalAeroCore::ExternalAeroCore(const TriMesh& wing,const TriangleBvh& bvh,const ParagliderConfig& config):config_(config),source_triangle_count_(wing.triangle_count())
	{
		if(wing.empty()||bvh.empty())throw std::invalid_argument("external aerodynamic core requires a placed fabric mesh");hierarchy_=AmrHierarchy::build_static(automatic_flow_domain(wing,config_.domain),wing,bvh,config_.amr);EmbeddedBoundaryBuildOptions options;options.min_volume_fraction=config_.amr.min_volume_fraction;options.complex_subdivisions=config_.amr.complex_subdivisions;embedded_boundary_=build_amr_embedded_boundary_atlas(hierarchy_,wing,bvh,options);if(!embedded_boundary_.ready_for_flow())throw std::runtime_error("external aerodynamic core has unresolved finest-level embedded-boundary topology");pressure_system_=build_composite_amr_pressure_system(hierarchy_,embedded_boundary_,true);fields_=std::make_unique<DeviceAmrFields>(hierarchy_);advection_=std::make_unique<DeviceAmrAdvection>(hierarchy_,bvh);projection_=std::make_unique<DeviceCompositeAmrProjection>(pressure_system_,*fields_);
	}

	ExternalAeroCore::~ExternalAeroCore()=default;

	ExternalAeroStepStats ExternalAeroCore::project(bool warm_start)
	{
		ExternalAeroStepStats stats;stats.dt=config_.solver.cfl*hierarchy_.finest_cell_size()/std::max(1e-9,config_.freestream.speed);const auto begin=std::chrono::steady_clock::now();stats.pressure=projection_->project(static_cast<Real>(config_.freestream.rho),static_cast<Real>(stats.dt),config_.solver.projection_tolerance,config_.solver.projection_max_iterations,warm_start);const auto end=std::chrono::steady_clock::now();stats.projection_ms=elapsed_ms(begin,end);stats.gpu_step_ms=stats.projection_ms;stats.physical_time=physical_time_;return stats;
	}

	ExternalAeroStepStats ExternalAeroCore::initialize()
	{
		fields_->initialize_freestream(static_cast<Real>(config_.freestream.speed));fields_->apply_external_aero_boundaries(static_cast<Real>(config_.freestream.speed));projection_->initialize_special_freestream(static_cast<Real>(config_.freestream.speed));ExternalAeroStepStats stats=project(false);initialized_=stats.pressure.converged;stats.physical_time=physical_time_;return stats;
	}

	ExternalAeroStepStats ExternalAeroCore::step()
	{
		if(!initialized_)throw std::logic_error("external aerodynamic core must have a converged initialization before stepping");const Real dt=static_cast<Real>(config_.solver.cfl*hierarchy_.finest_cell_size()/std::max(1e-9,config_.freestream.speed));const auto step_begin=std::chrono::steady_clock::now(),advection_begin=step_begin;advection_->advect(*fields_,dt);const auto advection_end=std::chrono::steady_clock::now();advection_->diffuse_smagorinsky(*fields_,static_cast<Real>(config_.freestream.nu),static_cast<Real>(config_.solver.smagorinsky_cs),dt);const auto turbulence_end=std::chrono::steady_clock::now();fields_->apply_external_aero_boundaries(static_cast<Real>(config_.freestream.speed));ExternalAeroStepStats stats=project(true);const auto step_end=std::chrono::steady_clock::now();stats.advection_ms=elapsed_ms(advection_begin,advection_end);stats.turbulence_ms=elapsed_ms(advection_end,turbulence_end);stats.les_applied=true;stats.gpu_step_ms=elapsed_ms(step_begin,step_end);if(stats.pressure.converged)physical_time_+=stats.dt;stats.physical_time=physical_time_;return stats;
	}

	AerodynamicLoads ExternalAeroCore::pressure_loads(double reference_pressure) const
	{
		std::vector<Real> device_pressure;projection_->download_pressure(device_pressure);std::vector<double> pressure(device_pressure.size());for(std::size_t q=0;q<pressure.size();++q)pressure[q]=static_cast<double>(device_pressure[q]);return compute_pressure_loads(pressure_system_,pressure,source_triangle_count_,config_.freestream,config_.reference,reference_pressure);
	}

	void ExternalAeroCore::download_fields(AmrHostFields& host) const{fields_->download(host);}
	std::size_t ExternalAeroCore::gpu_bytes() const{return fields_->bytes()+advection_->bytes()+projection_->bytes();}
	std::size_t ExternalAeroCore::protected_face_count() const{return advection_->protected_face_count();}
	std::size_t ExternalAeroCore::active_face_count() const{return advection_->active_face_count();}
}
