// sim_worker.h — the simulation worker (G1). Owns the ChannelFluidCore and steps it on
// a dedicated QThread. It performs NO OpenGL.
//
// Rendering is decoupled from stepping: step() runs WITHOUT any lock (the worker owns
// the core exclusively), and after each step the worker copies the {u,v,w,p} fields into
// display snapshot buffers under display_mutex(). The main-thread slice kernel reads ONLY
// those snapshots (also under display_mutex()), on its own CUDA stream — so an expensive
// step never blocks a frame and the viewer runs at display rate. Control (play/pause/
// step/quit) is via thread-safe atomics; stats are pushed to the GUI via a queued signal
// (cobod-slicer mainwindow.cpp:1691 pattern).
#pragma once

#include "core/fluid/channel_core.h"
#include "gui/flow_particles.h" // FlowField
#include "gui/scene_io.h"       // CheckpointState / CheckpointStatePtr

#include <QObject>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <vector>

namespace windcfd::gui
{
	// Flow-diversion metrics sampled from a settled voxelized-obstacle run (headless assertions
	// for the --voxelize smoke path). Computed AFTER the worker thread has joined, so touching
	// the core from the caller's thread is race-free.
	struct DiversionReport
	{
		bool active = false;           // a voxel obstacle was present when sampled
		double inlet_U = 0.0;          // reference inlet speed [m/s]
		long long solid_cells = 0;     // obstacle cell count
		double mean_speed_solid = 0.0; // mean cell-centred |u| over solid cells (≈0: no penetration)
		double max_speed_fluid = 0.0;  // max cell-centred |u| over fluid cells (>U: diversion accel)
	};

	class SimWorker : public QObject
	{
		Q_OBJECT
	public:
		explicit SimWorker(std::unique_ptr<windcfd::core::ChannelFluidCore> core, QObject* parent = nullptr)
			: QObject(parent), core_(std::move(core)) {}
		~SimWorker() override;

		// --- Main-thread render access (always take display_mutex() first) --------
		std::mutex& display_mutex() { return disp_mtx_; }
		bool display_ready() const { return disp_ready_.load(); }
		const double* disp_u() const { return du_; }
		const double* disp_v() const { return dv_; }
		const double* disp_w() const { return dw_; }
		const double* disp_p() const { return dp_; }
		const unsigned char* disp_solid() const { return ds_; } // solid mask snapshot (auto-range excludes it)

		bool playing() const { return playing_.load(); }
		long long steps() const { return steps_.load(); }
		double sim_time() const { return sim_time_.load(); }

		// --- Live flow solid-mask snapshot (for the "Show voxels" overlay) -----------------------
		// The worker publishes the CURRENT flow solid mask (obstacle ∪ any runtime-marked solid,
		// g.p_count() bytes, 1=solid, g.pidx order) whenever it CHANGES: at load and on an obstacle
		// rebuild. maskGeneration() bumps on every publish so the viewer re-extracts the overlay
		// geometry only then — never a per-frame mask copy. copyMask() (main thread) returns the
		// latest mask + its grid, or false if none.
		std::uint64_t maskGeneration() const { return mask_gen_.load(); }
		bool copyMask(std::vector<unsigned char>& mask, windcfd::core::MacGrid& grid);

		// Factory that rebuilds the core with a given obstacle mask + surface mode. Set once
		// (captures the immutable SimRecipe). Invoked ONLY on the worker thread by run().
		void setRebuildFactory(std::function<std::unique_ptr<windcfd::core::ChannelFluidCore>(
			const std::vector<unsigned char>&, int)> f) { factory_ = std::move(f); }

		// Queue a core rebuild with a new obstacle mask (thread-safe from the main thread). The
		// worker picks it up at the top of its next loop iteration and rebuilds ON ITS THREAD —
		// the core is never constructed/replaced from the main/GL thread. Flow is re-initialised.
		void requestRebuild(std::vector<unsigned char> solid, int solid_mode);

		// Live-update the inlet (current) speed from the main thread WITHOUT a rebuild/reset: the worker
		// applies it to the running core at the top of its next loop iteration, and re-applies it after
		// any obstacle rebuild so a re-inject keeps the chosen current. The developing flow evolves
		// toward the new inlet speed. Used by the GUI "Input speed" control to test different currents.
		void setInletSpeed(double U) { inlet_speed_.store(U); inlet_override_.store(true); inlet_speed_dirty_.store(true); }

		// Live-switch the inlet profile between uniform (top-hat) and a log-law boundary-layer profile
		// (near-bed u→0 at the domain floor, a thinner wall-bounded inflow). Applied on the worker thread
		// and re-applied after a rebuild. `z0` = wall roughness, `bed_datum` = inlet floor elevation [m].
		// Used by the GUI "Inlet profile" toggle.
		void setInletProfile(bool loglaw, double z0, double bed_datum)
		{
			{ std::lock_guard<std::mutex> lk(inlet_prof_mtx_); inlet_prof_loglaw_ = loglaw; inlet_prof_z0_ = z0; inlet_prof_datum_ = bed_datum; }
			inlet_prof_override_.store(true);
			inlet_prof_dirty_.store(true);
		}

		// Live TIDAL REVERSAL driver (M9, RESEARCH §8). When `on`, the worker drives the inlet each step
		// with a reversing tide U_d(t) = U_max·s(t) (s: +1 plateau → cosine slack ramp → −1 plateau → …)
		// via ChannelFluidCore::set_inlet_speed(|U_d|) + set_flow_direction(sign), face-swapping the
		// inlet/outlet at slack where |U_d|→0. It takes over the manual "Input speed" while active.
		// `plateau_s`/`ramp_s` are the hold and compressed-slack durations [s]. Disabling restores steady
		// forward flow. Applied on the worker thread; no reset.
		void setTidalReversal(bool on, double Umax, double plateau_s, double ramp_s)
		{
			tidal_umax_.store(Umax); tidal_plateau_.store(plateau_s); tidal_ramp_.store(ramp_s);
			tidal_on_.store(on); tidal_dirty_.store(true);
		}
		bool tidalOn() const { return tidal_on_.load(); }
		double tidePhase() const { return tide_phase_.load(); } // current signed U_d(t) [m/s] (0 when off)

		// Sample flow-diversion metrics from the current core. MUST be called with the worker
		// thread stopped/joined (no concurrent step()). Reports whether the obstacle is present.
		void computeDiversion(DiversionReport& out) const;

		// --- Checkpoint / auto-save (scene_io) ---------------------------------------------------
		// Prime the step / sim-time counters before the thread starts (a scene restore continues from
		// the saved step, not t=0). Call BEFORE moveToThread/start, like setRebuildFactory.
		void primeCounters(long long steps, double sim_time) { steps_.store(steps); sim_time_.store(sim_time); }
		// Auto-save cadence: with n>0 the worker gathers a checkpoint and emits checkpointReady every n
		// steps (0 = off). Thread-safe; picked up at the next step boundary.
		void setAutosaveInterval(int n) { autosave_interval_.store(n < 0 ? 0 : n); }
		// Request a one-off checkpoint (e.g. a manual "Save Scene"): the worker gathers the full state on
		// its own thread (race-free) at the top of its next loop and emits checkpointReady(state, tag).
		// `tag` < 0 marks a manual save (auto-saves pass the step count). Serviced even while paused.
		void requestCheckpoint(long long tag) { snapshot_tag_.store(tag); snapshot_pending_.store(true); }

		// --- Flow-arrow host snapshot (main-thread particle advection) ---------------------------
		// When enabled, the worker copies {u,v,w,solid} to host front buffers at ~12 Hz (worker
		// thread, throttled; a no-op when disabled) so the viewer can advect the arrow tracers on
		// the CPU without touching the device or the render's display_mutex(). Off by default so a
		// viewer with arrows hidden pays nothing.
		void setWantHostFlow(bool on) { want_host_flow_.store(on); if (!on) flow_ready_.store(false); }

		// Throttle the display snapshot to at most once per `seconds` of wall time (the GUI "Fast sim"
		// control). publish_display() copies the fields + hard-syncs the device EVERY step, so when
		// nobody is watching, skipping it most steps hands those cycles back to the solver. 0 (default)
		// = publish every step (live 60 Hz view; G1 byte-identical). The view is still refreshed on
		// pause and on a manual single-step, so it is never stale. Thread-safe; read on the worker thread.
		void setDisplayInterval(double seconds) { display_interval_.store(seconds < 0.0 ? 0.0 : seconds); }

		// Run `fn` with a borrowed, read-only view of the latest host flow snapshot, under the flow
		// lock (so the worker cannot swap it mid-advection). Returns false (and does not call `fn`)
		// if no snapshot is ready. Cheap for the worker: it only contends during the O(1) buffer
		// swap. Main-thread caller only.
		bool withFlowField(const std::function<void(const FlowField&)>& fn);

	public slots:
		// Runs the step loop until stop(); connect to QThread::started.
		void run();
		void setPlaying(bool p) { playing_.store(p); }
		void stepOnce() { step_requests_.fetch_add(1); }
		void stop() { quit_.store(true); }

	signals:
		// Emitted periodically (queued) with running totals for the GUI status bar. simFps is the measured
		// simulation rate (wall-clock steps/second, lightly smoothed) — the sim throughput, distinct from
		// the render fps; it rises when "Fast sim" throttles the display.
		void stats(qint64 steps, double simTime, double dt, double simFps);

		// A checkpoint was gathered on the worker thread (queued → the GUI writes it to disk on the main
		// thread). `tag` < 0 = manual save; >= 0 = the step count of an auto-save. The shared_ptr keeps the
		// (large) host state alive until the writer is done; no big copy crosses the connection.
		void checkpointReady(windcfd::gui::CheckpointStatePtr state, qint64 tag);

	private:
		void alloc_display();     // lazy device snapshot buffers (needs the grid)
		void publish_display();   // copy core fields -> snapshots (under disp_mtx_)
		void free_display();
		void apply_pending_rebuild(); // worker-thread: swap in a new core if one was requested
		void maybe_publish_host_flow(); // worker-thread: throttled D2H of {u,v,w,solid} for arrows
		void emit_checkpoint(long long tag); // worker-thread: gather full state → emit checkpointReady
		void publish_mask(const std::vector<unsigned char>& mask); // worker-thread: snapshot the mask + bump gen

		std::unique_ptr<windcfd::core::ChannelFluidCore> core_;

		// Live flow solid-mask snapshot (host), republished only when the mask changes (see above).
		std::mutex mask_mtx_;
		std::vector<unsigned char> mask_snapshot_;
		windcfd::core::MacGrid mask_grid_{};
		std::atomic<std::uint64_t> mask_gen_{ 0 };

		// Pending core rebuild (obstacle re-injection). Guarded by rebuild_mtx_; flagged atomic.
		std::function<std::unique_ptr<windcfd::core::ChannelFluidCore>(const std::vector<unsigned char>&, int)> factory_;
		std::mutex rebuild_mtx_;
		std::vector<unsigned char> pending_mask_;
		int pending_mode_ = 0;
		std::atomic<bool> rebuild_pending_{ false };

		// Live inlet-speed override (GUI "Input speed"): applied on the worker thread at the top of the
		// loop and re-applied after a rebuild. inlet_override_ latches once the user has changed it, so
		// a subsequent obstacle re-inject keeps the chosen current instead of reverting to the recipe's.
		std::atomic<double> inlet_speed_{ 0.0 };
		std::atomic<bool> inlet_override_{ false };
		std::atomic<bool> inlet_speed_dirty_{ false };

		// Live inlet-profile override (GUI "Inlet profile" toggle), same apply-on-worker + re-apply-after-
		// rebuild pattern as the speed. Guarded by inlet_prof_mtx_ (three coupled values).
		std::mutex inlet_prof_mtx_;
		bool inlet_prof_loglaw_ = false;
		double inlet_prof_z0_ = 0.0, inlet_prof_datum_ = 0.0;
		std::atomic<bool> inlet_prof_override_{ false };
		std::atomic<bool> inlet_prof_dirty_{ false };

		// Live tidal-reversal driver (GUI "Tidal reversal"). Params are atomics; tidal_dirty_ flags an
		// enable/disable/param change (the worker resets its tide clock and, on disable, restores steady
		// forward flow). tide_clock_ + last tide-apply are WORKER-THREAD-ONLY (advanced by dt each step).
		// tide_phase_ = the current signed U_d(t) for the GUI status readout.
		std::atomic<bool> tidal_on_{ false };
		std::atomic<double> tidal_umax_{ 0.8 }, tidal_plateau_{ 30.0 }, tidal_ramp_{ 20.0 };
		std::atomic<bool> tidal_dirty_{ false };
		std::atomic<double> tide_phase_{ 0.0 };
		double tide_clock_ = 0.0; // worker-thread only: flow-time since the tide was (re)started [s]

		// Display snapshots (device), same layout/counts as the core MAC fields.
		std::mutex disp_mtx_;
		double* du_ = nullptr;
		double* dv_ = nullptr;
		double* dw_ = nullptr;
		double* dp_ = nullptr;
		unsigned char* ds_ = nullptr; // solid-cell snapshot (obstacle), for auto-range
		std::atomic<bool> disp_ready_{ false };

		// Display-snapshot throttle (GUI "Fast sim"): 0 ⇒ publish every step. When > 0 the worker skips
		// publish_display()/host-flow between wall-clock windows, always publishing on pause + on a manual
		// single-step. The timestamp + dirty flag below are touched ONLY on the worker thread.
		std::atomic<double> display_interval_{ 0.0 };
		std::chrono::steady_clock::time_point last_display_pub_{};
		bool display_dirty_ = false;

		std::atomic<bool> playing_{ true };
		std::atomic<bool> quit_{ false };
		std::atomic<int> step_requests_{ 0 };
		std::atomic<long long> steps_{ 0 };
		std::atomic<double> sim_time_{ 0.0 };

		// Checkpoint / auto-save (see requestCheckpoint / setAutosaveInterval). A manual request is a
		// one-shot flag + tag; auto-save fires whenever steps % interval == 0.
		std::atomic<bool> snapshot_pending_{ false };
		std::atomic<long long> snapshot_tag_{ -1 };
		std::atomic<int> autosave_interval_{ 0 };

		// Host flow snapshot for the CPU arrow advection (double-buffered: the worker D2H-copies into
		// the *back* buffers, then swaps under flow_mtx_ so the viewer always reads a complete front).
		std::atomic<bool> want_host_flow_{ false };
		std::atomic<bool> flow_ready_{ false };
		std::mutex flow_mtx_;
		std::vector<double> fu_front_, fv_front_, fw_front_, fu_back_, fv_back_, fw_back_;
		std::vector<unsigned char> fs_front_, fs_back_;
		windcfd::core::MacGrid flow_grid_front_{};
		int flow_sign_front_ = 1; // current inlet direction snapshotted with the flow (guarded by flow_mtx_)
		std::chrono::steady_clock::time_point last_flow_pub_{};
	};
}
