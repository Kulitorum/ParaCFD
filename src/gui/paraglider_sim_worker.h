#pragma once

#include "core/fluid/external_aero_core.h"
#include "gui/flow_particles.h"
#include "gui/slice_field.h"

#include <QObject>

#include <atomic>
#include <cstdint>
#include <deque>
#include <functional>
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
		double max_abs_regular_velocity = 0,max_abs_special_velocity = 0;
		double max_embedded_cfl_rate = 0,effective_cfl = 0;
		int pressure_iterations = 0;
		bool initialized = false,converged = false;
		bool conservation_valid = false;
		double max_abs_divergence = 0,volume_weighted_rms_divergence = 0;
		double max_integrated_flux_error = 0,absolute_integrated_flux_error = 0;
		double net_integrated_flux_error = 0;
		double flow_change = 0,settling_score = 0,settling_force_drift = 0,settling_force_rms = 0;
		double flow_throughs = 0;
		std::size_t gpu_bytes = 0;
		bool playing = false,auto_pause_enabled = false,auto_paused = false,settling_ready = false;
		std::string error;
		std::vector<float> cp_plus,cp_minus,delta_cp;
		float cp_min = -1,cp_max = 1;
		float side_cp_min = -1,side_cp_max = 1;
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
		void setPlaying(bool playing);
		void stepOnce();
		void configureAutoPause(bool enabled,double sensitivity);
		bool latestSnapshot(std::uint64_t& generation, std::uint64_t& surface_generation,
			ParagliderDisplaySnapshot& out) const;
		// Borrow the latest coarse, uniform display resampling of the live AMR fields.
		// The callback runs under a short mutex and must not retain the pointers.
		bool withFlowField(const std::function<void(const FlowField&)>& fn) const;
		// Sample the current plane directly from the finest active AMR bricks. Unlike withFlowField(),
		// this does not pass through the coarse 3-D particle/tracer snapshot.
		bool sampleAmrSlice(const SliceParams& params, std::vector<float>& values, FieldRange& range) const;

	public slots:
		void run();

	signals:
		void finished();

	private:
		void publish(const paracfd::core::ExternalAeroStepStats& stats,bool include_surface);
		void publishFlowField();
		void updateSettling(double physical_time,const paracfd::core::Vec3d& force,
			const paracfd::core::ExternalAeroConservationStats& conservation);
		struct SettlingSample{double time=0;paracfd::core::Vec3d force{};double flow_change=0;};
		std::unique_ptr<paracfd::core::ExternalAeroCore> core_;
		std::unique_ptr<paracfd::core::AmrHostFields> display_amr_;
		std::vector<double> display_pressure_; // includes appended two-sided EB fragment DOFs
		mutable std::mutex display_amr_mutex_;
		bool display_amr_ready_ = false;
		std::atomic<bool> stop_{false},playing_{false},auto_pause_enabled_{true},auto_paused_{false},settling_reset_{true};
		std::atomic<double> auto_pause_sensitivity_{0.65};
		std::atomic<int> step_requests_{0};
		mutable std::mutex snapshot_mutex_;
		ParagliderDisplaySnapshot snapshot_;
		mutable std::mutex flow_mutex_;
		paracfd::core::MacGrid flow_grid_{};
		std::vector<double> flow_u_, flow_v_, flow_w_, flow_p_;
		bool flow_ready_ = false;
		std::uint64_t flow_generation_ = 0;
		double latest_flow_change_ = 0;
		std::deque<SettlingSample> settling_history_;
		int settling_consecutive_=0;
		double settling_score_=0,settling_force_drift_=0,settling_force_rms_=0,flow_throughs_=0,settling_epoch_time_=0;
		bool settling_ready_=false;
		long long steps_ = 0;
	};
}
