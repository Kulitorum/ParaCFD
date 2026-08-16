#pragma once

#include "core/fluid/external_aero_core.h"

#include <QObject>

#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace paracfd::gui
{
	struct ParagliderDisplaySnapshot
	{
		std::uint64_t generation = 0;
		std::uint64_t surface_generation = 0;
		long long steps = 0;
		double physical_time = 0,dt = 0,step_ms = 0,projection_ms = 0,residual = 0;
		double max_abs_regular_velocity = 0,max_abs_special_velocity = 0,effective_cfl = 0;
		int pressure_iterations = 0;
		bool initialized = false,converged = false;
		bool conservation_valid = false;
		double max_abs_divergence = 0,volume_weighted_rms_divergence = 0;
		double max_integrated_flux_error = 0,absolute_integrated_flux_error = 0;
		double net_integrated_flux_error = 0;
		std::size_t gpu_bytes = 0;
		std::string error;
		std::vector<float> delta_cp;
		float cp_min = -1,cp_max = 1;
		paracfd::core::Vec3d pressure_force{};
		bool coefficients_valid = false;
		double cd_pressure = 0,cs_pressure = 0,cl_pressure = 0;
	};

	class ParagliderSimWorker final : public QObject
	{
		Q_OBJECT
	public:
		explicit ParagliderSimWorker(std::unique_ptr<paracfd::core::ExternalAeroCore> core);
		~ParagliderSimWorker() override = default;
		void stop() { stop_.store(true); }
		void setPlaying(bool playing) { playing_.store(playing); }
		void stepOnce() { step_requests_.fetch_add(1); }
		bool latestSnapshot(std::uint64_t& generation, std::uint64_t& surface_generation,
			ParagliderDisplaySnapshot& out) const;

	public slots:
		void run();

	signals:
		void finished();

	private:
		void publish(const paracfd::core::ExternalAeroStepStats& stats,bool include_surface);
		std::unique_ptr<paracfd::core::ExternalAeroCore> core_;
		std::atomic<bool> stop_{false},playing_{false};
		std::atomic<int> step_requests_{0};
		mutable std::mutex snapshot_mutex_;
		ParagliderDisplaySnapshot snapshot_;
		long long steps_ = 0;
	};
}
