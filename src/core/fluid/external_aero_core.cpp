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

	ExternalAeroCore::ExternalAeroCore(const TriMesh& wing,const TriangleBvh& bvh,
		const ParagliderConfig& config,ExternalAeroExecutionOptions options)
		:config_(config),source_triangle_count_(wing.triangle_count()),
		 use_conservative_cell_momentum_(options.conservative_cell_momentum),
		 enable_smooth_fabric_wall_(options.smooth_fabric_wall),
		 enable_pressure_impulse_(options.pressure_impulse)
	{
		if(wing.empty()||bvh.empty())
			throw std::invalid_argument("external aerodynamic core requires a placed fabric mesh");
		symmetry_plane_y_=wing.bbox_min[1];
		hierarchy_=AmrHierarchy::build_static(automatic_flow_domain(wing,config_.domain),
			wing,bvh,config_.amr,config_.domain.half_wing_symmetry);

		EmbeddedBoundaryBuildOptions eb_options;
		eb_options.min_volume_fraction=config_.amr.min_volume_fraction;
		eb_options.retain_signed_bracketed_small_roots_for_face_state=
			!use_conservative_cell_momentum_;
		eb_options.min_aperture_area_fraction=config_.amr.min_aperture_area_fraction;
		eb_options.complex_subdivisions=config_.amr.complex_subdivisions;
		eb_options.allow_unverified_same_fragment_patches_for_diagnostics=
			options.allow_unsafe_same_fragment_patches;
		eb_options.exact_cell_decomposer=options.exact_cell_decomposer;
		eb_options.use_exact_cell_decomposer_as_development_oracle=
			options.use_exact_cell_decomposer_as_development_oracle;
		embedded_boundary_=build_amr_embedded_boundary_atlas(hierarchy_,wing,bvh,eb_options);

		if(!embedded_boundary_.ready_for_flow())
		{
			std::size_t owned_count=0;
			std::string first_owned,first_any;
			std::vector<ExternalAeroPreprocessingProblem> diagnostics;
			auto describe=[&](const AmrEbLevelAtlas& level,const UnresolvedEbCell& problem,
				bool owned)
			{
				const auto coordinate=level.topology.grid.cell_coord(problem.parent_cell);
				const Vec3d centroid=level.topology.grid.cell_centroid(problem.parent_cell);
				std::string reason=(owned?"owned ":"topology-halo ")+std::string("level ")+
					std::to_string(level.level)+" cell ["+std::to_string(coordinate[0])+","+
					std::to_string(coordinate[1])+","+std::to_string(coordinate[2])+"] centroid ["+
					std::to_string(centroid.x)+","+std::to_string(centroid.y)+","+
					std::to_string(centroid.z)+"]: "+problem.reason;
				if(!problem.source_triangles.empty())
				{
					const std::uint32_t triangle=problem.source_triangles.front();
					reason+=problem.source_triangles_are_candidates?
						" (candidate triangle ":" (source triangle ";
					reason+=std::to_string(triangle);
					if(triangle<wing.source_face_ids.size())
						reason+=", CAD face "+std::to_string(wing.source_face_ids[triangle]);
					if(triangle<bvh.triangle_count())
					{
						const BvhTriangle& geometry=bvh.triangle(triangle);
						reason+=", vertices ["+std::to_string(geometry.a.x)+","+
							std::to_string(geometry.a.y)+","+std::to_string(geometry.a.z)+"] ["+
							std::to_string(geometry.b.x)+","+std::to_string(geometry.b.y)+","+
							std::to_string(geometry.b.z)+"] ["+std::to_string(geometry.c.x)+","+
							std::to_string(geometry.c.y)+","+std::to_string(geometry.c.z)+"]";
					}
					reason+=")";
				}
				return reason;
			};
			for(const AmrEbLevelAtlas& level:embedded_boundary_.levels)
				for(const UnresolvedEbCell& problem:level.topology.unresolved)
				{
					const bool owned=problem.parent_cell>=0&&
						problem.parent_cell<static_cast<int>(level.owned_cell.size())&&
						level.owned_cell[problem.parent_cell];
					const auto coordinate=level.topology.grid.cell_coord(problem.parent_cell);
					const Aabb3d cell_box=level.topology.grid.cell_box(
						coordinate[0],coordinate[1],coordinate[2]);
					ExternalAeroPreprocessingProblem diagnostic;
					diagnostic.level=level.level;
					diagnostic.cell_coordinate={coordinate[0],coordinate[1],coordinate[2]};
					diagnostic.cell_lo=cell_box.lo;
					diagnostic.cell_hi=cell_box.hi;
					diagnostic.centroid=level.topology.grid.cell_centroid(problem.parent_cell);
					diagnostic.owned=owned;
					diagnostic.reason=problem.reason;
					diagnostic.source_triangles=problem.source_triangles;
					diagnostic.source_triangles_are_candidates=
						problem.source_triangles_are_candidates;
					// Some stabilization failures are aggregate-level and therefore have no
					// direct source list. The exact BVH query recovers every fabric triangle
					// crossing the failed control volume without guessing from its centroid.
					if(diagnostic.source_triangles.empty())
					{
						diagnostic.source_triangles=bvh.query_aabb(cell_box);
						diagnostic.source_triangles_are_candidates=true;
					}
					diagnostic.source_triangles.erase(std::remove_if(
						diagnostic.source_triangles.begin(),diagnostic.source_triangles.end(),
						[&](std::uint32_t triangle){return triangle>=wing.triangle_count();}),
						diagnostic.source_triangles.end());
					std::sort(diagnostic.source_triangles.begin(),diagnostic.source_triangles.end());
					diagnostic.source_triangles.erase(std::unique(
						diagnostic.source_triangles.begin(),diagnostic.source_triangles.end()),
						diagnostic.source_triangles.end());
					for(const std::uint32_t triangle:diagnostic.source_triangles)
						if(triangle<wing.source_face_ids.size())
							diagnostic.source_face_ids.push_back(wing.source_face_ids[triangle]);
					std::sort(diagnostic.source_face_ids.begin(),diagnostic.source_face_ids.end());
					diagnostic.source_face_ids.erase(std::unique(
						diagnostic.source_face_ids.begin(),diagnostic.source_face_ids.end()),
						diagnostic.source_face_ids.end());
					diagnostics.push_back(std::move(diagnostic));
					owned_count+=owned;
					if(first_any.empty())first_any=describe(level,problem,owned);
					if(owned&&first_owned.empty())first_owned=describe(level,problem,true);
				}
			const std::size_t total=embedded_boundary_.unresolved_count();
			throw ExternalAeroPreprocessingError("external aerodynamic core has "+std::to_string(total)+
				" unresolved EB cell(s) ("+std::to_string(owned_count)+" owned, "+
				std::to_string(total-owned_count)+" topology halo); first: "+
				(first_owned.empty()?first_any:first_owned),std::move(diagnostics),owned_count);
		}

		CompositeAmrPressureBuildOptions pressure_options;
		pressure_options.pressure_outlet_xmax=true;
		if(options.use_qualitative_first_order_orthogonal_pressure)
			pressure_options.unsupported_nonorthogonal_correction=
				UnsupportedNonorthogonalCorrectionPolicy::qualitative_preview_first_order_orthogonal;
		pressure_system_=build_composite_amr_pressure_system(
			hierarchy_,embedded_boundary_,pressure_options);
		fields_=std::make_unique<DeviceAmrFields>(hierarchy_);
		projection_=std::make_unique<DeviceCompositeAmrProjection>(pressure_system_,*fields_);
		if(use_conservative_cell_momentum_)
			cell_momentum_=std::make_unique<DeviceCompositeCellMomentumTransport>(pressure_system_,*fields_);
		else
			advection_=std::make_unique<DeviceAmrAdvection>(hierarchy_,bvh,2.5,&pressure_system_);
	}

	ExternalAeroCore::~ExternalAeroCore()=default;

	ExternalAeroStepStats ExternalAeroCore::project(bool warm_start,double dt,bool sync_coarse_fine)
	{
		if(!(dt>0))throw std::invalid_argument("external aerodynamic projection timestep");ExternalAeroStepStats stats;stats.dt=dt;stats.conservative_cell_momentum=use_conservative_cell_momentum_;const auto begin=std::chrono::steady_clock::now();if(sync_coarse_fine)projection_->sync_coarse_fine_from_fields();stats.pressure=projection_->project(static_cast<Real>(config_.freestream.rho),static_cast<Real>(stats.dt),config_.solver.projection_tolerance,config_.solver.projection_max_iterations,warm_start);const auto end=std::chrono::steady_clock::now();stats.projection_ms=elapsed_ms(begin,end);stats.gpu_step_ms=stats.projection_ms;stats.physical_time=physical_time_;return stats;
	}

	ExternalAeroStepStats ExternalAeroCore::initialize()
	{
		const Real speed=static_cast<Real>(config_.freestream.speed);
		if(pressure_system_.freestream_connected.size()!=static_cast<std::size_t>(pressure_system_.storage_size))throw std::logic_error("external aerodynamic pressure components have no initialization classification");

		// A zero-thickness closed surface creates a valid, disconnected fluid component.
		// Initializing every MAC face/aperture to the external freestream gives that sealed
		// component an unphysical impulse and leaves a large, persistent cavity pressure.
		// Seed one velocity per actual control volume, classify it by connectivity to the
		// external boundary, and reconstruct both ordinary MAC and compact special fluxes.
		CompositeCellMomentumState state;
		state.x.assign(pressure_system_.storage_size,Real(0));
		state.y.assign(pressure_system_.storage_size,Real(0));
		state.z.assign(pressure_system_.storage_size,Real(0));
		for(int q=0;q<pressure_system_.storage_size;++q)
			if(pressure_system_.active[q]&&pressure_system_.freestream_connected[q])state.x[q]=speed;

		AmrHostFields initial_fields(hierarchy_);
		CompositeAmrFluxes initial_fluxes;
		reconstruct_composite_cell_fluxes_cpu(pressure_system_,state,initial_fields,initial_fluxes);

		// The cell-to-face reconstruction only visits connections between two pressure
		// DOFs. Supply the two physical X boundary faces explicitly; the boundary kernel
		// subsequently fills their ghosts and all free-slip boundary values.
		const int bs=pressure_system_.brick_size;
		for(int level_index=0;level_index<static_cast<int>(hierarchy_.levels().size());++level_index)
		{
			const AmrLevel& level=hierarchy_.levels()[level_index];
			AmrHostLevelFields& values=initial_fields.levels()[level_index];
			for(int brick=0;brick<static_cast<int>(level.bricks.size());++brick)
			{
				const BrickMetadata& meta=level.bricks[brick];
				if(!meta.active())continue;
				for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)
				{
					if(meta.flags&BRICK_XMIN)
					{
						const int dof=pressure_system_.dof(level_index,brick,0,j,k);
						values.u[values.layout.u_index(brick,0,j,k)]=pressure_system_.active[dof]&&pressure_system_.freestream_connected[dof]?speed:Real(0);
					}
					if(meta.flags&BRICK_XMAX)
					{
						const int dof=pressure_system_.dof(level_index,brick,bs-1,j,k);
						values.u[values.layout.u_index(brick,bs,j,k)]=pressure_system_.active[dof]&&pressure_system_.freestream_connected[dof]?speed:Real(0);
					}
				}
			}
		}

		fields_->upload(initial_fields);
		projection_->upload_special_fluxes(initial_fluxes);
		if(cell_momentum_)
		{
			cell_momentum_->upload_state(state);
			cell_momentum_->reconstruct_fluxes(projection_->coarse_fine_velocity_device(),projection_->embedded_velocity_device());
		}
		fields_->apply_external_aero_boundaries(speed);const double maximum=std::max(1e-9,std::abs(config_.freestream.speed)),dt=config_.solver.cfl*hierarchy_.finest_cell_size()/maximum;ExternalAeroStepStats stats=project(false,dt);if(cell_momentum_&&enable_pressure_impulse_&&stats.pressure.converged)cell_momentum_->apply_projected_pressure_gradient(projection_->pressure(),static_cast<Real>(dt),static_cast<Real>(config_.freestream.rho));stats.max_abs_regular_velocity=cell_momentum_?cell_momentum_->max_abs_velocity():fields_->max_abs_velocity();stats.max_abs_special_velocity=projection_->max_abs_special_velocity();stats.max_embedded_cfl_rate=projection_->max_embedded_cfl_rate();stats.max_abs_velocity=std::max(stats.max_abs_regular_velocity,stats.max_abs_special_velocity);stats.cfl_velocity=maximum;stats.effective_cfl=maximum*dt/hierarchy_.finest_cell_size();initialized_=stats.pressure.converged;stats.physical_time=physical_time_;return stats;
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
			cell_momentum_->step(projection_->coarse_fine_velocity_device(),projection_->embedded_velocity_device(),dt,true,speed);advection_end=std::chrono::steady_clock::now();cell_momentum_->diffuse_smagorinsky(nu,cs,dt);const double diffusion_rate=cell_momentum_->last_diffusion_rate();const int diffusion_substeps=cell_momentum_->last_diffusion_substeps();turbulence_end=std::chrono::steady_clock::now();if(enable_smooth_fabric_wall_){cell_momentum_->apply_smooth_fabric_wall_model(dt,nu);smooth_fabric_wall_applied_=true;}cell_momentum_->reconstruct_pressure_consistent_fluxes(projection_->pressure(),projection_->coarse_fine_velocity_device(),projection_->embedded_velocity_device());fields_->apply_external_aero_boundaries(speed);embedded_transport_end=std::chrono::steady_clock::now();stats=project(true,dt_value,false);stats.max_diffusion_rate=diffusion_rate;stats.diffusion_substeps=diffusion_substeps;if(enable_pressure_impulse_&&stats.pressure.converged)cell_momentum_->apply_projected_pressure_gradient(projection_->pressure(),dt,rho);
		}
		else
		{
			advection_->advect(*fields_,dt);advection_end=std::chrono::steady_clock::now();advection_->diffuse_smagorinsky(*fields_,nu,cs,dt);turbulence_end=std::chrono::steady_clock::now();fields_->apply_external_aero_boundaries(speed);projection_->transport_embedded_apertures(dt,nu,cs,false);if(enable_smooth_fabric_wall_){projection_->apply_smooth_fabric_wall_model(dt,nu);smooth_fabric_wall_applied_=true;}embedded_transport_end=std::chrono::steady_clock::now();stats=project(true,dt_value);
		}
		const auto step_end=std::chrono::steady_clock::now();
		stats.max_abs_regular_velocity=regular_maximum;stats.max_abs_special_velocity=special_maximum;stats.max_embedded_cfl_rate=embedded_cfl_rate;stats.max_abs_velocity=std::max(regular_maximum,special_maximum);stats.cfl_velocity=cfl_velocity;stats.effective_cfl=cfl_rate*dt_value;stats.advection_ms=elapsed_ms(advection_begin,advection_end);stats.turbulence_ms=elapsed_ms(advection_end,turbulence_end);stats.embedded_transport_ms=elapsed_ms(turbulence_end,embedded_transport_end);stats.les_applied=true;stats.embedded_transport_applied=true;stats.smooth_fabric_wall_applied=smooth_fabric_wall_applied_;stats.gpu_step_ms=elapsed_ms(step_begin,step_end);if(stats.pressure.converged)physical_time_+=stats.dt;stats.physical_time=physical_time_;return stats;
	}

	AerodynamicLoads ExternalAeroCore::unscaled_pressure_loads(double reference_pressure) const
	{
		std::vector<Real> device_pressure;projection_->download_pressure(device_pressure);std::vector<double> pressure(device_pressure.size());for(std::size_t q=0;q<pressure.size();++q)pressure[q]=static_cast<double>(device_pressure[q]);return compute_pressure_loads(pressure_system_,pressure,source_triangle_count_,config_.freestream,config_.reference,reference_pressure);
	}
	AerodynamicLoads ExternalAeroCore::pressure_loads(double reference_pressure) const
	{
		AerodynamicLoads loads=unscaled_pressure_loads(reference_pressure);if(config_.domain.half_wing_symmetry)reconstruct_y_symmetric_integrated_loads(loads,symmetry_plane_y_,config_.reference.moment_origin.y);return loads;
	}
	AerodynamicLoads ExternalAeroCore::aerodynamic_loads(double reference_pressure) const
	{
		AerodynamicLoads loads=unscaled_pressure_loads(reference_pressure);if(smooth_fabric_wall_applied_){std::vector<SmoothFabricWallPatchLoad> viscous;if(cell_momentum_)cell_momentum_->download_smooth_fabric_wall_loads(static_cast<Real>(config_.freestream.rho),viscous);else projection_->download_smooth_fabric_wall_loads(static_cast<Real>(config_.freestream.rho),viscous);accumulate_viscous_loads(loads,viscous,config_.freestream,config_.reference);}if(config_.domain.half_wing_symmetry)reconstruct_y_symmetric_integrated_loads(loads,symmetry_plane_y_,config_.reference.moment_origin.y);return loads;
	}

	void ExternalAeroCore::download_fields(AmrHostFields& host) const{fields_->download(host);}
	void ExternalAeroCore::download_special_fluxes(CompositeAmrFluxes& host) const{projection_->download_special_fluxes(host);}
	void ExternalAeroCore::download_pressure(std::vector<double>& host) const{std::vector<Real> values;projection_->download_pressure(values);host.assign(values.begin(),values.end());}
	bool ExternalAeroCore::download_cell_momentum_state(CompositeCellMomentumState& host) const{if(!cell_momentum_)return false;cell_momentum_->download_state(host);return true;}
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
