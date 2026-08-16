#include "gui/paraglider_sim_worker.h"

#include <QThread>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>

namespace paracfd::gui
{
	ParagliderSimWorker::ParagliderSimWorker(
		std::unique_ptr<paracfd::core::ExternalAeroCore> core)
		: core_(std::move(core))
	{
		display_amr_ = std::make_unique<paracfd::core::AmrHostFields>(core_->hierarchy());
	}

	bool ParagliderSimWorker::withFlowField(
		const std::function<void(const FlowField&)>& fn) const
	{
		std::lock_guard lock(flow_mutex_);
		if (!flow_ready_ || flow_u_.empty()) return false;
		FlowField field;
		field.u = flow_u_.data(); field.v = flow_v_.data(); field.w = flow_w_.data();
		field.p = flow_p_.data(); field.grid = flow_grid_;
		fn(field);
		return true;
	}

	void ParagliderSimWorker::publishFlowField()
	{
		using namespace paracfd::core;
		const AmrHierarchy& hierarchy = core_->hierarchy();
		if (hierarchy.levels().empty()) return;
		core_->download_fields(*display_amr_);

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

		std::lock_guard lock(flow_mutex_);
		flow_grid_ = grid;
		flow_u_.swap(u); flow_v_.swap(v); flow_w_.swap(w); flow_p_.swap(pressure);
		flow_ready_ = true;
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
		out.conservation_valid = snapshot_.conservation_valid;
		out.max_abs_divergence = snapshot_.max_abs_divergence;
		out.volume_weighted_rms_divergence = snapshot_.volume_weighted_rms_divergence;
		out.max_integrated_flux_error = snapshot_.max_integrated_flux_error;
		out.absolute_integrated_flux_error = snapshot_.absolute_integrated_flux_error;
		out.net_integrated_flux_error = snapshot_.net_integrated_flux_error;
		out.gpu_bytes = snapshot_.gpu_bytes;
		out.error = snapshot_.error;
		out.cp_min = snapshot_.cp_min;
		out.cp_max = snapshot_.cp_max;
		out.pressure_force = snapshot_.pressure_force;
		out.coefficients_valid = snapshot_.coefficients_valid;
		out.cd_pressure = snapshot_.cd_pressure;
		out.cs_pressure = snapshot_.cs_pressure;
		out.cl_pressure = snapshot_.cl_pressure;
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
			loads = core_->pressure_loads();
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
			snapshot_.coefficients_valid = loads.force_coefficients_valid;
			snapshot_.cd_pressure = loads.cd_pressure;
			snapshot_.cs_pressure = loads.cs_pressure;
			snapshot_.cl_pressure = loads.cl_pressure;
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
				publish(stats, publish_fields);
				if (!stats.pressure.converged) playing_.store(false);
			}
		}
		catch (const std::exception& exception)
		{
			std::lock_guard lock(snapshot_mutex_);
			snapshot_.error = exception.what();
			snapshot_.converged = false;
			++snapshot_.generation;
		}
		emit finished();
	}
}
