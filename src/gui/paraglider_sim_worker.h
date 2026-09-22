#pragma once

#include "core/aero_convergence.h"
#include "core/fluid/external_aero_core.h"
#include "gui/flow_particles.h"
#include "gui/slice_field.h"

#include <QByteArray>
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
	// Playback records the independent source channels. Derived views are rebuilt on demand,
	// so selecting velocity also covers speed, vorticity, Q, arrows and tracers; selecting
	// pressure also covers Cp and the pressure gradient.
	struct RecordingOptions
	{
		bool velocity = true;
		bool pressure = true;
		bool surface_pressure = true;
		bool every_step = false;
		bool any() const { return velocity || pressure || surface_pressure; }
	};

	struct PlaybackDatasetInfo
	{
		RecordingOptions channels;
		std::size_t frame_count = 0;
		std::size_t compressed_bytes = 0;
		double first_time = 0, last_time = 0;
		bool recording = false;
	};

	struct PlaybackFrameSnapshot
	{
		std::size_t frame_index = 0, frame_count = 0;
		long long steps = 0;
		double physical_time = 0;
		bool has_surface_pressure = false;
		std::vector<float> cp_plus, cp_minus, delta_cp;
		std::vector<float> triangle_pressure_force_xyz;
		float cp_min = -1, cp_max = 1;
		float side_cp_min = -1, side_cp_max = 1;
	};

	struct FieldProbeSample
	{
		bool valid = false;
		double x = 0, y = 0, z = 0;
		double u = 0, v = 0, w = 0;
		double pressure_delta = 0, pressure_coefficient = 0;
		double pressure_gradient = 0, vorticity_magnitude = 0, q_criterion = 0;
	};

	struct ScalarVolume
	{
		int nx = 0, ny = 0, nz = 0;
		float h = 0.0f;
		std::uint64_t generation = 0;
		Field field = Field::SpeedMag;
		float minimum = 0.0f, maximum = 0.0f;
		float focus_abs = 1.0f; // robust 98th percentile |value| used by volume defaults
		float focus_positive = 1.0f; // robust positive range used by one-sided iso-surfaces
		std::vector<float> values;
		bool valid() const { return nx > 1 && ny > 1 && nz > 1 && values.size() == static_cast<std::size_t>(nx) * ny * nz; }
	};

	// Auto-pause is permitted only after the wake has travelled this fraction of
	// the streamwise domain. The score is computed earlier for live diagnostics.
	inline constexpr double kAutoPauseMinimumFlowThroughs = 0.5;

	using SimulationPauseReason=paracfd::core::AerodynamicRunExitReason;

	struct ParagliderDisplaySnapshot
	{
		std::uint64_t generation = 0;
		std::uint64_t surface_generation = 0;
		long long steps = 0;
		double physical_time = 0,dt = 0,step_ms = 0,projection_ms = 0,residual = 0;
		double max_diffusion_rate = 0;
		int diffusion_substeps = 1;
		double max_abs_regular_velocity = 0,max_abs_special_velocity = 0;
		double max_embedded_cfl_rate = 0,effective_cfl = 0;
		int pressure_iterations = 0;
		bool initialized = false,converged = false,conservative_cell_momentum = false;
		bool conservation_valid = false;
		double max_abs_divergence = 0,volume_weighted_rms_divergence = 0;
		double max_integrated_flux_error = 0,absolute_integrated_flux_error = 0;
		double net_integrated_flux_error = 0;
		double flow_change = 0,settling_score = 0,settling_force_drift = 0,settling_force_rms = 0;
		double flow_throughs = 0;
		paracfd::core::Vec3d mean_force{};
		paracfd::core::Vec3d mean_pressure_force{},mean_viscous_force{};
		double mean_force_drift = 0,mean_force_rms = 0;
		// Time-weighted force over the bounded final convergence window. This remains
		// useful at a maximum-flow exit even when the two half-window convergence
		// comparison is not yet complete.
		paracfd::core::Vec3d bounded_mean_force{};
		paracfd::core::Vec3d bounded_mean_pressure_force{},bounded_mean_viscous_force{};
		double bounded_mean_force_coverage = 0;
		std::size_t gpu_bytes = 0;
		bool playing = false,auto_pause_enabled = false,auto_paused = false,settling_ready = false,mean_force_ready = false;
		bool bounded_mean_force_ready = false,bounded_mean_force_complete = false;
		SimulationPauseReason pause_reason = SimulationPauseReason::None;
		std::string error;
		std::vector<float> cp_plus,cp_minus,delta_cp;
		std::vector<float> triangle_pressure_force_xyz;
		float cp_min = -1,cp_max = 1;
		float side_cp_min = -1,side_cp_max = 1;
		paracfd::core::Vec3d pressure_force{};
		paracfd::core::Vec3d viscous_force{},total_force{};
		bool viscous_loads_valid = false;
		bool coefficients_valid = false;
		double cd_pressure = 0,cs_pressure = 0,cl_pressure = 0;
		double cd_viscous = 0,cs_viscous = 0,cl_viscous = 0;
		double cd = 0,cs = 0,cl = 0;
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
		void configureConvergenceExit(double relative_mean_tolerance,double maximum_flow_throughs);
		void configureRecording(const RecordingOptions& options);
		PlaybackDatasetInfo playbackDatasetInfo() const;
		std::size_t playbackFrameAtTime(double physical_time,int direction) const;
		bool setPlaybackFrame(std::size_t index,PlaybackFrameSnapshot& snapshot);
		void endPlayback();
		bool playbackActive() const { return playback_active_.load(); }
		bool latestSnapshot(std::uint64_t& generation, std::uint64_t& surface_generation,
			ParagliderDisplaySnapshot& out) const;
		// Borrow the latest coarse, uniform display resampling of the live AMR fields.
		// The callback runs under a short mutex and must not retain the pointers.
		bool withFlowField(const std::function<void(const FlowField&)>& fn) const;
		// Sample the current plane directly from the finest active AMR bricks. Unlike withFlowField(),
		// this does not pass through the coarse 3-D particle/tracer snapshot.
		bool sampleAmrSlice(const SliceParams& params, std::vector<float>& values, FieldRange& range) const;
		// Finest-active AMR point probe. Derived gradients use neighbouring finest-visible samples.
		bool samplePoint(double x, double y, double z, FieldProbeSample& sample) const;
		// Coarse uniform visualization snapshot used for cached iso-surface extraction and ray casting.
		bool scalarVolume(Field field, ScalarVolume& volume) const;
		std::uint64_t visualizationGeneration() const;

	public slots:
		void run();

	signals:
		void finished();

	private:
		struct StateSample
		{
			bool valid = false;
			double u = 0, v = 0, w = 0, p = 0, h = 0;
		};
		bool sampleStateLocked(const paracfd::core::Vec3d& world, StateSample& sample,
			bool need_pressure) const;
		bool sampleDerivedLocked(const paracfd::core::Vec3d& world, Field field,
			FieldProbeSample& sample, bool all_derived) const;
		void publish(const paracfd::core::ExternalAeroStepStats& stats,bool include_surface,
			bool update_convergence=true);
		void publishFlowField();
		void captureRecordingFrame(double scheduled_time,long long steps);
		void captureScheduledRecordingFrames(double physical_time,long long steps);
		bool samplePlaybackSlice(const SliceParams& params,std::vector<float>& values,
			FieldRange& range) const;
		bool samplePlaybackPoint(double x,double y,double z,FieldProbeSample& sample) const;
		void updateSettling(double physical_time,const paracfd::core::AerodynamicLoads& loads,
			const paracfd::core::ExternalAeroConservationStats& conservation);
		struct SettlingSample{double time=0;paracfd::core::Vec3d force{};double flow_change=0;};
		std::unique_ptr<paracfd::core::ExternalAeroCore> core_;
		std::unique_ptr<paracfd::core::AmrHostFields> display_amr_;
		std::vector<double> display_pressure_; // includes appended two-sided EB fragment DOFs
		mutable std::mutex display_amr_mutex_;
		bool display_amr_ready_ = false;
		std::atomic<bool> stop_{false},playing_{false},auto_pause_enabled_{true},auto_paused_{false},settling_reset_{true};
		std::atomic<double> auto_pause_sensitivity_{0.65};
		std::atomic<double> mean_force_tolerance_{0.02},maximum_flow_throughs_{2.5};
		std::atomic<int> step_requests_{0};
		mutable std::mutex snapshot_mutex_;
		ParagliderDisplaySnapshot snapshot_;
		mutable std::mutex flow_mutex_;
		paracfd::core::MacGrid flow_grid_{};
		std::vector<double> flow_u_, flow_v_, flow_w_, flow_p_;
		bool flow_ready_ = false;
		std::uint64_t flow_generation_ = 0;
		struct RecordedFrame
		{
			double physical_time = 0;
			long long steps = 0;
			paracfd::core::MacGrid grid{};
			QByteArray velocity,pressure,surface_pressure;
			std::size_t surface_triangles = 0;
			std::size_t compressed_bytes = 0;
			float cp_min = -1,cp_max = 1,side_cp_min = -1,side_cp_max = 1;
		};
		mutable std::mutex recording_mutex_;
		RecordingOptions recording_options_{};
		std::vector<RecordedFrame> recorded_frames_;
		std::size_t recorded_compressed_bytes_ = 0;
		std::uint64_t recording_epoch_ = 0;
		double next_recording_time_ = 0;
		static constexpr double kRecordingInterval = 1.0 / 30.0;
		std::atomic<bool> recording_started_{false},playback_active_{false};
		paracfd::core::MacGrid playback_grid_{};
		std::vector<double> playback_u_,playback_v_,playback_w_,playback_p_;
		bool playback_has_velocity_ = false,playback_has_pressure_ = false;
		std::uint64_t playback_generation_ = (std::uint64_t{1} << 63);
		double latest_flow_change_ = 0;
		std::deque<SettlingSample> settling_history_;
		std::deque<paracfd::core::TimedAerodynamicForce> mean_force_history_,
			mean_pressure_force_history_,mean_viscous_force_history_;
		paracfd::core::AerodynamicMeanConvergence mean_convergence_,
			mean_pressure_convergence_,mean_viscous_convergence_;
		paracfd::core::Vec3d bounded_mean_force_{},bounded_mean_pressure_force_{},
			bounded_mean_viscous_force_{};
		double bounded_mean_force_coverage_=0;
		bool bounded_mean_force_ready_=false,bounded_mean_force_complete_=false;
		int settling_consecutive_=0,mean_consecutive_=0;
		double settling_score_=0,settling_force_drift_=0,settling_force_rms_=0,flow_throughs_=0,settling_epoch_time_=0;
		bool settling_ready_=false;
		std::atomic<SimulationPauseReason> pause_reason_{SimulationPauseReason::None};
		long long steps_ = 0;
	};
}
