#include "gui/paraglider_sim_worker.h"

#include <QThread>

#include <algorithm>
#include <cmath>
#include <exception>

namespace paracfd::gui
{
	ParagliderSimWorker::ParagliderSimWorker(
		std::unique_ptr<paracfd::core::ExternalAeroCore> core)
		: core_(std::move(core))
	{
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
		out.pressure_iterations = snapshot_.pressure_iterations;
		out.initialized = snapshot_.initialized;
		out.converged = snapshot_.converged;
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
			out.delta_cp = snapshot_.delta_cp;
			surface_generation = snapshot_.surface_generation;
		}
		generation = snapshot_.generation;
		return true;
	}

	void ParagliderSimWorker::publish(const paracfd::core::ExternalAeroStepStats& stats,bool include_surface)
	{
		std::vector<float> delta_cp;
		float cp_min = 0.0f, cp_max = 0.0f;
		paracfd::core::AerodynamicLoads loads;
		const bool have_surface = include_surface && stats.pressure.converged;
		if (have_surface)
		{
			loads = core_->pressure_loads();
			delta_cp.resize(loads.triangles.size());
			float maximum = 0.0f;
			for (std::size_t triangle = 0; triangle < loads.triangles.size(); ++triangle)
			{
				const double value = loads.triangles[triangle].delta_cp;
				delta_cp[triangle] = std::isfinite(value) ? static_cast<float>(value) : 0.0f;
				maximum = std::max(maximum, std::abs(delta_cp[triangle]));
			}
			maximum = std::max(maximum, 1e-5f);
			cp_min = -maximum;
			cp_max = maximum;
		}

		std::lock_guard lock(snapshot_mutex_);
		snapshot_.steps = steps_;
		snapshot_.physical_time = stats.physical_time;
		snapshot_.dt = stats.dt;
		snapshot_.step_ms = stats.gpu_step_ms;
		snapshot_.projection_ms = stats.projection_ms;
		snapshot_.pressure_iterations = stats.pressure.iterations;
		snapshot_.residual = stats.pressure.relative_residual;
		snapshot_.initialized = core_->initialized();
		snapshot_.converged = stats.pressure.converged;
		snapshot_.error.clear();
		if (have_surface)
		{
			snapshot_.delta_cp = std::move(delta_cp);
			snapshot_.cp_min = cp_min;
			snapshot_.cp_max = cp_max;
			snapshot_.pressure_force = loads.pressure_force;
			snapshot_.coefficients_valid = loads.force_coefficients_valid;
			snapshot_.cd_pressure = loads.cd_pressure;
			snapshot_.cs_pressure = loads.cs_pressure;
			snapshot_.cl_pressure = loads.cl_pressure;
			++snapshot_.surface_generation;
		}
		++snapshot_.generation;
	}

	void ParagliderSimWorker::run()
	{
		try
		{
			paracfd::core::ExternalAeroStepStats stats = core_->initialize();
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
				publish(stats, (steps_ % 10) == 0 || !stats.pressure.converged);
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
