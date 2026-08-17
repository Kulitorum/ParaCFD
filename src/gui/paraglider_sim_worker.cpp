#include "gui/paraglider_sim_worker.h"

#include <QThread>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <stdexcept>

namespace paracfd::gui
{
	ParagliderSimWorker::ParagliderSimWorker(
		std::unique_ptr<paracfd::core::ExternalAeroCore> core)
		: core_(std::move(core))
	{
		display_amr_ = std::make_unique<paracfd::core::AmrHostFields>(core_->hierarchy());
		latest_flow_change_=std::numeric_limits<double>::infinity();
	}

	void ParagliderSimWorker::setPlaying(bool playing)
	{
		playing_.store(playing);if(playing){auto_paused_.store(false);settling_reset_.store(true);}
		else if(!auto_paused_.load())settling_reset_.store(true);
		std::lock_guard lock(snapshot_mutex_);snapshot_.playing=playing;snapshot_.auto_paused=auto_paused_.load();
		if(playing){snapshot_.settling_ready=false;snapshot_.mean_force_ready=false;snapshot_.pause_reason=SimulationPauseReason::None;snapshot_.settling_score=0;snapshot_.settling_force_drift=0;snapshot_.settling_force_rms=0;snapshot_.mean_force_drift=0;snapshot_.mean_force_rms=0;snapshot_.flow_throughs=0;}
		++snapshot_.generation;
	}

	void ParagliderSimWorker::stepOnce()
	{
		auto_paused_.store(false);settling_reset_.store(true);step_requests_.fetch_add(1);
		std::lock_guard lock(snapshot_mutex_);snapshot_.auto_paused=false;snapshot_.pause_reason=SimulationPauseReason::None;++snapshot_.generation;
	}

	void ParagliderSimWorker::configureAutoPause(bool enabled,double sensitivity)
	{
		auto_pause_enabled_.store(enabled);auto_pause_sensitivity_.store(std::clamp(sensitivity,0.0,1.0));settling_reset_.store(true);
		std::lock_guard lock(snapshot_mutex_);snapshot_.auto_pause_enabled=enabled;++snapshot_.generation;
	}

	void ParagliderSimWorker::configureSweepExit(bool enabled,double relative_mean_tolerance,double maximum_flow_throughs)
	{
		sweep_exit_enabled_.store(enabled);
		sweep_mean_tolerance_.store(std::clamp(relative_mean_tolerance,0.001,0.2));
		sweep_max_flow_throughs_.store(std::clamp(maximum_flow_throughs,1.0,20.0));
		settling_reset_.store(true);
	}

	bool ParagliderSimWorker::withFlowField(
		const std::function<void(const FlowField&)>& fn) const
	{
		std::lock_guard lock(flow_mutex_);
		if (!flow_ready_ || flow_u_.empty()) return false;
		FlowField field;
		field.u = flow_u_.data(); field.v = flow_v_.data(); field.w = flow_w_.data();
		field.p = flow_p_.data(); field.grid = flow_grid_; field.generation = flow_generation_;
		fn(field);
		return true;
	}

	bool ParagliderSimWorker::sampleAmrSlice(
		const SliceParams& params, std::vector<float>& values, FieldRange& range) const
	{
		using namespace paracfd::core;
		std::lock_guard lock(display_amr_mutex_);
		if (!display_amr_ready_ || params.nu <= 0 || params.nv <= 0) return false;
		const AmrHierarchy& hierarchy = core_->hierarchy();
		const Aabb3d& domain = hierarchy.domain();
		const Vec3d extent = domain.hi - domain.lo;
		values.assign(static_cast<std::size_t>(params.nu) * params.nv, 0.0f);
		range = {};
		range.field_min = std::numeric_limits<float>::max();
		range.field_max = std::numeric_limits<float>::lowest();

		auto inside = [](double value, double lo, double hi)
		{
			return std::clamp(value, lo, std::nextafter(hi, lo));
		};
		for (int b = 0; b < params.nv; ++b)
			for (int a = 0; a < params.nu; ++a)
			{
				float x, y, z;
				slice_vertex_world(params, a, b, x, y, z);
				const Vec3d world{
					inside(domain.lo.x + x, domain.lo.x, domain.lo.x + extent.x),
					inside(domain.lo.y + y, domain.lo.y, domain.lo.y + extent.y),
					inside(domain.lo.z + z, domain.lo.z, domain.lo.z + extent.z)};
				const BrickLocation location = hierarchy.locate_finest(world);
				if (!location.found()) continue;
				const AmrHostLevelFields& level = display_amr_->levels()[location.level];
				const BrickFieldLayout& layout = level.layout;
				const int brick = location.brick, i = location.cell.x, j = location.cell.y, k = location.cell.z;
				const double u = 0.5 * (static_cast<double>(level.u[layout.u_index(brick, i, j, k)])
					+ static_cast<double>(level.u[layout.u_index(brick, i + 1, j, k)]));
				const double v = 0.5 * (static_cast<double>(level.v[layout.v_index(brick, i, j, k)])
					+ static_cast<double>(level.v[layout.v_index(brick, i, j + 1, k)]));
				const double w = 0.5 * (static_cast<double>(level.w[layout.w_index(brick, i, j, k)])
					+ static_cast<double>(level.w[layout.w_index(brick, i, j, k + 1)]));
				const float speed = static_cast<float>(std::sqrt(u * u + v * v + w * w));
				float scalar = speed;
				switch (params.field)
				{
				case Field::VelU: scalar = static_cast<float>(u); break;
				case Field::VelV: scalar = static_cast<float>(v); break;
				case Field::VelW: scalar = static_cast<float>(w); break;
				case Field::Pressure:
				{
					const int pressure_dof=core_->pressure_system().pressure_dof_at_point(core_->embedded_boundary(),world);
					scalar=pressure_dof>=0&&pressure_dof<static_cast<int>(display_pressure_.size())
						?static_cast<float>(display_pressure_[pressure_dof]):0.0f;
					break;
				}
				case Field::SpeedMag: break;
				}
				if (!std::isfinite(scalar) || !std::isfinite(speed)) continue;
				values[static_cast<std::size_t>(b) * params.nu + a] = scalar;
				range.field_min = std::min(range.field_min, scalar);
				range.field_max = std::max(range.field_max, scalar);
				range.speed_max = std::max(range.speed_max, speed);
			}
		range.valid = range.field_min <= range.field_max;
		return true;
	}

	void ParagliderSimWorker::publishFlowField()
	{
		using namespace paracfd::core;
		const AmrHierarchy& hierarchy = core_->hierarchy();
		if (hierarchy.levels().empty()) return;
		std::unique_lock display_lock(display_amr_mutex_);
		core_->download_fields(*display_amr_);
		core_->download_pressure(display_pressure_);

		const Aabb3d& domain = hierarchy.domain();
		const double h = hierarchy.levels().front().h;
		const Vec3d extent = domain.hi - domain.lo;
		MacGrid grid;
		grid.nx = std::max(1, static_cast<int>(std::llround(extent.x / h)));
		grid.ny = std::max(1, static_cast<int>(std::llround(extent.y / h)));
		grid.nz = std::max(1, static_cast<int>(std::llround(extent.z / h)));
		grid.h = h;

		// Resample cell-centred values from the finest active brick. This copy is solely a
		// visualization product; the CFD fields remain staggered, FP32, and GPU resident.
		const std::size_t cells = static_cast<std::size_t>(grid.p_count());
		std::vector<double> uc(cells), vc(cells), wc(cells), pressure(cells);
		for (int k = 0; k < grid.nz; ++k)
			for (int j = 0; j < grid.ny; ++j)
				for (int i = 0; i < grid.nx; ++i)
				{
					const Vec3d world{domain.lo.x + (i + 0.5) * h,
						domain.lo.y + (j + 0.5) * h, domain.lo.z + (k + 0.5) * h};
					const BrickLocation location = hierarchy.locate_finest(world);
					if (!location.found()) continue;
					const AmrHostLevelFields& level = display_amr_->levels()[location.level];
					const BrickFieldLayout& layout = level.layout;
					const int bi = location.brick, li = location.cell.x, lj = location.cell.y, lk = location.cell.z;
					const std::size_t out = static_cast<std::size_t>(grid.pidx(i, j, k));
					uc[out] = 0.5 * (static_cast<double>(level.u[layout.u_index(bi, li, lj, lk)])
						+ static_cast<double>(level.u[layout.u_index(bi, li + 1, lj, lk)]));
					vc[out] = 0.5 * (static_cast<double>(level.v[layout.v_index(bi, li, lj, lk)])
						+ static_cast<double>(level.v[layout.v_index(bi, li, lj + 1, lk)]));
					wc[out] = 0.5 * (static_cast<double>(level.w[layout.w_index(bi, li, lj, lk)])
						+ static_cast<double>(level.w[layout.w_index(bi, li, lj, lk + 1)]));
					pressure[out] = static_cast<double>(level.p[layout.cell_index(bi, li, lj, lk)]);
				}

		std::vector<double> u(static_cast<std::size_t>(grid.u_count()));
		std::vector<double> v(static_cast<std::size_t>(grid.v_count()));
		std::vector<double> w(static_cast<std::size_t>(grid.w_count()));
		auto finite_or_zero = [](double value) { return std::isfinite(value) ? value : 0.0; };
		for (int k = 0; k < grid.nz; ++k) for (int j = 0; j < grid.ny; ++j) for (int i = 0; i <= grid.nx; ++i)
		{
			const int left = std::max(0, i - 1), right = std::min(grid.nx - 1, i);
			u[grid.uidx(i,j,k)] = finite_or_zero(0.5 * (uc[grid.pidx(left,j,k)] + uc[grid.pidx(right,j,k)]));
		}
		for (int k = 0; k < grid.nz; ++k) for (int j = 0; j <= grid.ny; ++j) for (int i = 0; i < grid.nx; ++i)
		{
			const int below = std::max(0, j - 1), above = std::min(grid.ny - 1, j);
			v[grid.vidx(i,j,k)] = finite_or_zero(0.5 * (vc[grid.pidx(i,below,k)] + vc[grid.pidx(i,above,k)]));
		}
		for (int k = 0; k <= grid.nz; ++k) for (int j = 0; j < grid.ny; ++j) for (int i = 0; i < grid.nx; ++i)
		{
			const int back = std::max(0, k - 1), front = std::min(grid.nz - 1, k);
			w[grid.widx(i,j,k)] = finite_or_zero(0.5 * (wc[grid.pidx(i,j,back)] + wc[grid.pidx(i,j,front)]));
		}
		for (double& value : pressure) value = finite_or_zero(value);
		display_amr_ready_ = true;
		display_lock.unlock();

		std::lock_guard lock(flow_mutex_);
		if(flow_u_.size()==u.size()&&flow_v_.size()==v.size()&&flow_w_.size()==w.size())
		{
			double difference2=0,scale2=0;auto accumulate=[&](const std::vector<double>& previous,const std::vector<double>& current){for(std::size_t q=0;q<current.size();++q){const double d=current[q]-previous[q];difference2+=d*d;scale2+=current[q]*current[q];}};accumulate(flow_u_,u);accumulate(flow_v_,v);accumulate(flow_w_,w);latest_flow_change_=std::sqrt(difference2/std::max(scale2,1e-30));
		}
		else latest_flow_change_=std::numeric_limits<double>::infinity();
		flow_grid_ = grid;
		flow_u_.swap(u); flow_v_.swap(v); flow_w_.swap(w); flow_p_.swap(pressure);
		++flow_generation_;
		flow_ready_ = true;
	}

	void ParagliderSimWorker::updateSettling(double physical_time,const paracfd::core::Vec3d& force,const paracfd::core::ExternalAeroConservationStats& conservation)
	{
		if(settling_reset_.exchange(false))
		{
			settling_history_.clear();mean_force_history_.clear();
			settling_consecutive_=mean_consecutive_=0;
			settling_score_=settling_force_drift_=settling_force_rms_=flow_throughs_=0;
			settling_epoch_time_=physical_time;settling_ready_=false;
			mean_convergence_={};pause_reason_=SimulationPauseReason::None;
		}
		const bool steady_exit_enabled=auto_pause_enabled_.load();
		const bool sweep_exit_enabled=sweep_exit_enabled_.load();
		if(!steady_exit_enabled&&!sweep_exit_enabled)return;
		const double speed=std::max(1e-9,std::abs(core_->config().freestream.speed));
		const double flow_time=(core_->hierarchy().domain().hi.x-core_->hierarchy().domain().lo.x)/speed;
		if(!(flow_time>0))return;
		flow_throughs_=std::max(0.0,physical_time-settling_epoch_time_)/flow_time;
		mean_force_history_.push_back({physical_time,force});
		const double mean_window=std::max(0.5,flow_time),mean_oldest=physical_time-mean_window;
		while(mean_force_history_.size()>2&&mean_force_history_[1].time<mean_oldest)mean_force_history_.pop_front();
		mean_convergence_=paracfd::core::assess_aerodynamic_mean_convergence(mean_force_history_,physical_time,mean_window);

		if(!std::isfinite(latest_flow_change_))return;
		settling_history_.push_back({physical_time,force,latest_flow_change_});
		const double window=std::max(0.25,0.5*flow_time),oldest=physical_time-window;
		while(!settling_history_.empty()&&settling_history_.front().time<oldest)settling_history_.pop_front();
		if(settling_history_.size()<10)
		{
			if(sweep_exit_enabled&&flow_throughs_>=sweep_max_flow_throughs_.load())
			{
				playing_.store(false);auto_paused_.store(true);pause_reason_=SimulationPauseReason::MaximumFlowThroughs;
			}
			return;
		}
		const double history_span=physical_time-settling_history_.front().time;
		if(!(history_span>0.1*window))return;
		const double analysis_window=std::min(window,history_span);
		const double split=physical_time-0.5*analysis_window;paracfd::core::Vec3d old_mean{},new_mean{};int old_count=0,new_count=0;
		for(const SettlingSample& sample:settling_history_){paracfd::core::Vec3d& mean=sample.time<split?old_mean:new_mean;mean.x+=sample.force.x;mean.y+=sample.force.y;mean.z+=sample.force.z;if(sample.time<split)++old_count;else ++new_count;}
		if(old_count<4||new_count<4)return;old_mean.x/=old_count;old_mean.y/=old_count;old_mean.z/=old_count;new_mean.x/=new_count;new_mean.y/=new_count;new_mean.z/=new_count;
		auto magnitude=[](const paracfd::core::Vec3d& value){return std::sqrt(value.x*value.x+value.y*value.y+value.z*value.z);};
		const double scale=std::max(1.0,0.5*(magnitude(old_mean)+magnitude(new_mean)));settling_force_drift_=magnitude(new_mean-old_mean)/scale;
		double square=0,flow_change=0;for(const SettlingSample& sample:settling_history_)if(sample.time>=split){const paracfd::core::Vec3d delta=sample.force-new_mean;square+=delta.x*delta.x+delta.y*delta.y+delta.z*delta.z;flow_change+=sample.flow_change;}
		settling_force_rms_=std::sqrt(square/new_count)/scale;flow_change/=new_count;
		const double sensitivity=auto_pause_sensitivity_.load(),range=std::pow(20.0,sensitivity);
		const double drift_tolerance=0.001*range,noise_tolerance=0.003*range,flow_tolerance=0.0005*range;
		settling_score_=std::max({settling_force_drift_/drift_tolerance,settling_force_rms_/noise_tolerance,flow_change/flow_tolerance});settling_ready_=true;
		const bool conservative=conservation.volume_weighted_rms_divergence<1e-3;
		const bool observation_complete=flow_throughs_>=kAutoPauseMinimumFlowThroughs&&history_span>=0.9*window;
		if(steady_exit_enabled&&observation_complete&&settling_score_<1.0&&conservative)++settling_consecutive_;else settling_consecutive_=0;
		if(sweep_exit_enabled&&flow_throughs_>=1.0&&mean_convergence_.ready&&conservative&&mean_convergence_.relative_drift<sweep_mean_tolerance_.load())++mean_consecutive_;else mean_consecutive_=0;

		SimulationPauseReason reason=SimulationPauseReason::None;
		// Sweeps always retain at least one complete force window, even if the
		// instantaneous-field criterion happens to trigger earlier.
		const bool mean_force_converged=mean_convergence_.ready&&mean_convergence_.relative_drift<sweep_mean_tolerance_.load();
		if(settling_consecutive_>=3&&(!sweep_exit_enabled||mean_force_converged))reason=SimulationPauseReason::Steady;
		else if(mean_consecutive_>=3)reason=SimulationPauseReason::MeanConverged;
		else if(sweep_exit_enabled&&flow_throughs_>=sweep_max_flow_throughs_.load())reason=SimulationPauseReason::MaximumFlowThroughs;
		if(reason!=SimulationPauseReason::None)
		{
			playing_.store(false);auto_paused_.store(true);pause_reason_=reason;
			const char* label=reason==SimulationPauseReason::Steady?"steady":reason==SimulationPauseReason::MeanConverged?"mean-converged":"maximum-flow-throughs";
			std::fprintf(stderr,"[paraglider] auto-pause: reason=%s, t=%.6g s, flow-throughs=%.3f, settle-score=%.4g, mean-drift=%.4g, mean-rms=%.4g\n",label,physical_time,flow_throughs_,settling_score_,mean_convergence_.relative_drift,mean_convergence_.current_rms_fraction);
		}
	}

	bool ParagliderSimWorker::latestSnapshot(std::uint64_t& generation,
		std::uint64_t& surface_generation, ParagliderDisplaySnapshot& out) const
	{
		std::lock_guard lock(snapshot_mutex_);
		if (snapshot_.generation == generation) return false;
		out.generation = snapshot_.generation;
		out.surface_generation = snapshot_.surface_generation;
		out.steps = snapshot_.steps;
		out.physical_time = snapshot_.physical_time;
		out.dt = snapshot_.dt;
		out.step_ms = snapshot_.step_ms;
		out.projection_ms = snapshot_.projection_ms;
		out.residual = snapshot_.residual;
		out.max_abs_regular_velocity = snapshot_.max_abs_regular_velocity;
		out.max_abs_special_velocity = snapshot_.max_abs_special_velocity;
		out.max_embedded_cfl_rate = snapshot_.max_embedded_cfl_rate;
		out.effective_cfl = snapshot_.effective_cfl;
		out.pressure_iterations = snapshot_.pressure_iterations;
		out.initialized = snapshot_.initialized;
		out.converged = snapshot_.converged;
		out.conservative_cell_momentum = snapshot_.conservative_cell_momentum;
		out.conservation_valid = snapshot_.conservation_valid;
		out.max_abs_divergence = snapshot_.max_abs_divergence;
		out.volume_weighted_rms_divergence = snapshot_.volume_weighted_rms_divergence;
		out.max_integrated_flux_error = snapshot_.max_integrated_flux_error;
		out.absolute_integrated_flux_error = snapshot_.absolute_integrated_flux_error;
		out.net_integrated_flux_error = snapshot_.net_integrated_flux_error;
		out.flow_change=snapshot_.flow_change;out.settling_score=snapshot_.settling_score;out.settling_force_drift=snapshot_.settling_force_drift;out.settling_force_rms=snapshot_.settling_force_rms;out.flow_throughs=snapshot_.flow_throughs;
		out.mean_force=snapshot_.mean_force;out.mean_force_drift=snapshot_.mean_force_drift;out.mean_force_rms=snapshot_.mean_force_rms;out.mean_force_ready=snapshot_.mean_force_ready;out.pause_reason=snapshot_.pause_reason;
		out.gpu_bytes = snapshot_.gpu_bytes;
		out.playing=snapshot_.playing;out.auto_pause_enabled=snapshot_.auto_pause_enabled;out.auto_paused=snapshot_.auto_paused;out.settling_ready=snapshot_.settling_ready;
		out.error = snapshot_.error;
		out.cp_min = snapshot_.cp_min;
		out.cp_max = snapshot_.cp_max;
		out.pressure_force = snapshot_.pressure_force;
		out.viscous_force = snapshot_.viscous_force;
		out.total_force = snapshot_.total_force;
		out.viscous_loads_valid = snapshot_.viscous_loads_valid;
		out.coefficients_valid = snapshot_.coefficients_valid;
		out.cd_pressure = snapshot_.cd_pressure;
		out.cs_pressure = snapshot_.cs_pressure;
		out.cl_pressure = snapshot_.cl_pressure;
		out.cd_viscous = snapshot_.cd_viscous;
		out.cs_viscous = snapshot_.cs_viscous;
		out.cl_viscous = snapshot_.cl_viscous;
		out.cd = snapshot_.cd;
		out.cs = snapshot_.cs;
		out.cl = snapshot_.cl;
		if (snapshot_.surface_generation != surface_generation)
		{
			out.cp_plus = snapshot_.cp_plus;
			out.cp_minus = snapshot_.cp_minus;
			out.delta_cp = snapshot_.delta_cp;
			out.side_cp_min = snapshot_.side_cp_min;
			out.side_cp_max = snapshot_.side_cp_max;
			surface_generation = snapshot_.surface_generation;
		}
		generation = snapshot_.generation;
		return true;
	}

	void ParagliderSimWorker::publish(const paracfd::core::ExternalAeroStepStats& stats,bool include_surface)
	{
		std::vector<float> cp_plus,cp_minus,delta_cp;
		float cp_min = 0.0f, cp_max = 0.0f,side_cp_min=0.0f,side_cp_max=0.0f;
		paracfd::core::AerodynamicLoads loads;
		paracfd::core::ExternalAeroConservationStats conservation;
		const bool have_surface = include_surface && stats.pressure.converged;
		if (have_surface)
		{
			loads = core_->aerodynamic_loads();
			conservation = core_->conservation_stats();
			const float missing=std::numeric_limits<float>::quiet_NaN();
			cp_plus.assign(loads.triangles.size(),missing);
			cp_minus.assign(loads.triangles.size(),missing);
			delta_cp.assign(loads.triangles.size(),missing);
			float maximum = 0.0f,side_maximum=0.0f;
			for (std::size_t triangle = 0; triangle < loads.triangles.size(); ++triangle)
			{
				const auto& source=loads.triangles[triangle];
				if(std::isfinite(source.cp_plus)){cp_plus[triangle]=static_cast<float>(source.cp_plus);side_maximum=std::max(side_maximum,std::abs(cp_plus[triangle]));}
				if(std::isfinite(source.cp_minus)){cp_minus[triangle]=static_cast<float>(source.cp_minus);side_maximum=std::max(side_maximum,std::abs(cp_minus[triangle]));}
				if(std::isfinite(source.delta_cp)){delta_cp[triangle]=static_cast<float>(source.delta_cp);maximum=std::max(maximum,std::abs(delta_cp[triangle]));}
			}
			maximum = std::max(maximum, 1e-5f);
			side_maximum = std::max(side_maximum, 1e-5f);
			cp_min = -maximum;
			cp_max = maximum;
			side_cp_min=-side_maximum;
			side_cp_max=side_maximum;
		}
		if(have_surface)updateSettling(stats.physical_time,loads.viscous_loads_valid?loads.total_force:loads.pressure_force,conservation);

		std::lock_guard lock(snapshot_mutex_);
		snapshot_.steps = steps_;
		snapshot_.physical_time = stats.physical_time;
		snapshot_.dt = stats.dt;
		snapshot_.step_ms = stats.gpu_step_ms;
		snapshot_.projection_ms = stats.projection_ms;
		snapshot_.pressure_iterations = stats.pressure.iterations;
		snapshot_.residual = stats.pressure.relative_residual;
		snapshot_.max_abs_regular_velocity = stats.max_abs_regular_velocity;
		snapshot_.max_abs_special_velocity = stats.max_abs_special_velocity;
		snapshot_.max_embedded_cfl_rate = stats.max_embedded_cfl_rate;
		snapshot_.effective_cfl = stats.effective_cfl;
		snapshot_.gpu_bytes = core_->gpu_bytes();
		snapshot_.initialized = core_->initialized();
		snapshot_.converged = stats.pressure.converged;
		snapshot_.conservative_cell_momentum = core_->uses_conservative_cell_momentum();
		snapshot_.playing=playing_.load();snapshot_.auto_pause_enabled=auto_pause_enabled_.load();snapshot_.auto_paused=auto_paused_.load();snapshot_.settling_ready=settling_ready_;
		snapshot_.flow_change=latest_flow_change_;snapshot_.settling_score=settling_score_;snapshot_.settling_force_drift=settling_force_drift_;snapshot_.settling_force_rms=settling_force_rms_;snapshot_.flow_throughs=flow_throughs_;
		snapshot_.mean_force=(mean_convergence_.previous_mean+mean_convergence_.current_mean)*0.5;snapshot_.mean_force_drift=mean_convergence_.relative_drift;snapshot_.mean_force_rms=mean_convergence_.current_rms_fraction;snapshot_.mean_force_ready=mean_convergence_.ready;snapshot_.pause_reason=pause_reason_;
		snapshot_.error.clear();
		if (have_surface)
		{
			snapshot_.cp_plus = std::move(cp_plus);
			snapshot_.cp_minus = std::move(cp_minus);
			snapshot_.delta_cp = std::move(delta_cp);
			snapshot_.cp_min = cp_min;
			snapshot_.cp_max = cp_max;
			snapshot_.side_cp_min=side_cp_min;
			snapshot_.side_cp_max=side_cp_max;
			snapshot_.pressure_force = loads.pressure_force;
			snapshot_.viscous_force = loads.viscous_force;
			snapshot_.total_force = loads.total_force;
			snapshot_.viscous_loads_valid = loads.viscous_loads_valid;
			snapshot_.coefficients_valid = loads.force_coefficients_valid;
			snapshot_.cd_pressure = loads.cd_pressure;
			snapshot_.cs_pressure = loads.cs_pressure;
			snapshot_.cl_pressure = loads.cl_pressure;
			snapshot_.cd_viscous = loads.cd_viscous;
			snapshot_.cs_viscous = loads.cs_viscous;
			snapshot_.cl_viscous = loads.cl_viscous;
			snapshot_.cd = loads.cd;
			snapshot_.cs = loads.cs;
			snapshot_.cl = loads.cl;
			snapshot_.conservation_valid = true;
			snapshot_.max_abs_divergence = conservation.max_abs_divergence;
			snapshot_.volume_weighted_rms_divergence = conservation.volume_weighted_rms_divergence;
			snapshot_.max_integrated_flux_error = conservation.max_integrated_flux_error;
			snapshot_.absolute_integrated_flux_error = conservation.absolute_integrated_flux_error;
			snapshot_.net_integrated_flux_error = conservation.net_integrated_flux_error;
			++snapshot_.surface_generation;
		}
		++snapshot_.generation;
	}

	void ParagliderSimWorker::run()
	{
		try
		{
			paracfd::core::ExternalAeroStepStats stats = core_->initialize();
			publishFlowField();
			publish(stats, true);
			if (!stats.pressure.converged)
			{
				char message[192];
				std::snprintf(message, sizeof(message),
					"initial pressure projection did not converge after %d iterations "
					"(relative residual %.6g; requested tolerance %.6g)",
					stats.pressure.iterations, stats.pressure.relative_residual,
					core_->config().solver.projection_tolerance);
				throw std::runtime_error(message);
			}
			while (!stop_.load())
			{
				bool advance = playing_.load();
				const int requests = step_requests_.load();
				if (!advance && requests > 0)
				{
					advance = true;
					step_requests_.fetch_sub(1);
				}
				if (!advance)
				{
					QThread::msleep(10);
					continue;
				}
				stats = core_->step();
				++steps_;
				// Surface loads require a deliberately throttled pressure download; scalar
				// solver telemetry remains available after every GPU step.
				const bool publish_fields = (steps_ % 10) == 0 || !stats.pressure.converged;
				if (publish_fields) publishFlowField();
				if (!stats.pressure.converged) playing_.store(false);
				publish(stats, publish_fields);
			}
		}
		catch (const std::exception& exception)
		{
			playing_.store(false);
			std::lock_guard lock(snapshot_mutex_);
			snapshot_.error = exception.what();
			snapshot_.converged = false;
			snapshot_.playing = false;
			++snapshot_.generation;
		}
		emit finished();
	}
}
