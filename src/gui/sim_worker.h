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
#include "core/sediment/seabed_engine.h"
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

namespace scour::gui
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
		explicit SimWorker(std::unique_ptr<scour::core::ChannelFluidCore> core, QObject* parent = nullptr)
			: QObject(parent), core_(std::move(core)) {}
		~SimWorker() override;

		// --- Main-thread render access (always take display_mutex() first) --------
		std::mutex& display_mutex() { return disp_mtx_; }
		bool display_ready() const { return disp_ready_.load(); }
		const double* disp_u() const { return du_; }
		const double* disp_v() const { return dv_; }
		const double* disp_w() const { return dw_; }
		const double* disp_p() const { return dp_; }
		const double* disp_c() const { return dc_; } // suspended concentration snapshot (0 when no sediment)
		const unsigned char* disp_solid() const { return ds_; } // solid mask snapshot (auto-range excludes it)

		bool playing() const { return playing_.load(); }
		long long steps() const { return steps_.load(); }
		double sim_time() const { return sim_time_.load(); }

		// --- Erodible seabed (the morphodynamic loop) ----------------------------------------
		// Attach the morphodynamic engine; the worker then runs the FULL loop (fluid → τ_b →
		// suspended → Exner → avalanche) after `spinup` flow-only steps, and re-masks the flow core
		// in place as the bed moves. Set ONCE before the thread starts. Owned by the worker.
		void setMorpho(std::unique_ptr<scour::core::SeabedMorpho> m, int spinup);
		bool hasBed() const { return morpho_ != nullptr; }

		// Convert a RUNNING fluid viewer into an erodible-seabed run IN PLACE, preserving the developed
		// flow (the GUI "Add sand"). Thread-safe from the main thread; applied on the worker thread at the
		// top of its next loop: the pre-built engine is attached and the initial sand+structure mask is
		// installed via ChannelFluidCore::update_solid (untouched cells keep their {u,v,w,p}; the now-solid
		// cells are zeroed) — NO factory rebuild, NO init_uniform, NO t=0 reset. The bed-inlet mask is
		// enabled and (loglaw) the inlet switches to the log-law BL profile (z0/bed_datum), as a seabed run
		// needs. Morphology begins `spinup` steps later so the re-masked flow settles first. Replaces any
		// current engine (re-adding sand). The engine must be pre-built on the caller's grid.
		void requestSeabedConversion(std::unique_ptr<scour::core::SeabedMorpho> engine,
			std::vector<unsigned char> initial_solid, int spinup, bool loglaw, double z0, double bed_datum);
		// Live structure update (PLAN G3.3 SINK coupling): after the main thread re-settles the dropped
		// pile against the eroded bed and re-voxelizes it, hand the new rigid structure mask here; the
		// worker installs it into the morpho engine (SeabedMorpho::set_structure) and re-masks the flow in
		// place (update_solid), preserving the developed flow. Thread-safe; applied at the loop top. A
		// no-op without a bed (the flag stays pending, harmless).
		void requestStructureUpdate(std::vector<unsigned char> structure)
		{
			{ std::lock_guard<std::mutex> lk(structure_mtx_); pending_structure_ = std::move(structure); }
			structure_pending_.store(true);
		}
		double bedZ0() const { return bed_z0_; }
		// Bumped every time the worker republishes the bed elevation (a morphology step). Lets the GUI
		// video recorder capture a frame only when the sand has actually moved (mirrors maskGeneration).
		std::uint64_t bedGeneration() const { return bed_gen_.load(); }
		// Copy the latest per-column bed elevation z_b [m] (main-thread render). Returns false if no
		// bed. Takes display_mutex() internally. `zb` is nx·ny, row-major (j·nx+i).
		bool copyBed(std::vector<float>& zb);
		// Copy the latest per-column NET bed-exchange rate [m/s] (delivery − pickup) for the bed-surface
		// "Exchange rate" colouring. Returns false if no bed / not yet published. Same nx·ny layout as zb.
		bool copyBedExchange(std::vector<float>& rate);

		// --- Live flow solid-mask snapshot (for the "Show voxels" overlay) -----------------------
		// The worker publishes the CURRENT flow solid mask (structure ∪ sand ∪ obstacle ∪ any
		// runtime-marked solid, g.p_count() bytes, 1=solid, g.pidx order) whenever it CHANGES: at
		// load, on an obstacle rebuild, and on each bed re-mask (update_solid). maskGeneration() bumps
		// on every publish so the viewer re-extracts the overlay geometry only then — never a per-frame
		// mask copy. copyMask() (main thread) returns the latest mask + its grid, or false if none.
		std::uint64_t maskGeneration() const { return mask_gen_.load(); }
		bool copyMask(std::vector<unsigned char>& mask, scour::core::MacGrid& grid);

		// Factory that rebuilds the core with a given obstacle mask + surface mode. Set once
		// (captures the immutable SimRecipe). Invoked ONLY on the worker thread by run().
		void setRebuildFactory(std::function<std::unique_ptr<scour::core::ChannelFluidCore>(
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
		void setSedBoundary(int mode) { sed_bc_.store(mode); sed_bc_dirty_.store(true); } // suspended x-BC (SedBoundary)

		// Live MORFAC (morphological acceleration) change (GUI "Bed speed-up"): applied on the worker thread
		// to the morpho engine. MORFAC multiplies morphological time only, so the bed evolves M× faster than
		// the flow clock — the "catch up faster" knob. A no-op without a bed, but the flag stays pending (the
		// apply short-circuits on a null engine) so a value dialed in before "Add sand" lands once the engine
		// attaches. ⚠ Physically M≤10 steady / ≤5 reversing (RESEARCH §7) — a value choice, not enforced here.
		void setMorfac(double m) { morfac_.store(m); morfac_dirty_.store(true); }

		// Live-switch the inlet profile between uniform (top-hat) and a log-law boundary-layer profile
		// (near-bed u→0 at the bed top, removing the top-hat leading-edge over-erosion). Applied on the
		// worker thread and re-applied after a rebuild. `z0` = grain roughness, `bed_datum` = inlet bed
		// elevation [m]. Used by the GUI "Inlet profile" toggle.
		void setInletProfile(bool loglaw, double z0, double bed_datum)
		{
			{ std::lock_guard<std::mutex> lk(inlet_prof_mtx_); inlet_prof_loglaw_ = loglaw; inlet_prof_z0_ = z0; inlet_prof_datum_ = bed_datum; }
			inlet_prof_override_.store(true);
			inlet_prof_dirty_.store(true);
		}

		// Live TIDAL REVERSAL driver (M9, RESEARCH §8). When `on`, the worker drives the inlet each step
		// with a reversing tide U_d(t) = U_max·s(t) (s: +1 plateau → cosine slack ramp → −1 plateau → …)
		// via ChannelFluidCore::set_inlet_speed(|U_d|) + set_flow_direction(sign) (+ engine inlet speed),
		// face-swapping the inlet/outlet at slack where |U_d|→0. It takes over the manual "Input speed"
		// while active. `plateau_s`/`ramp_s` are the hold and compressed-slack durations [s]. Disabling
		// restores steady forward flow. Applied on the worker thread; no reset. (⚠ MORFAC ≤ 5 in reversing
		// flow — RESEARCH §7; that is a run-setup choice, not enforced here.)
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
		// the saved step, not t=0). Call BEFORE moveToThread/start, like setMorpho/setRebuildFactory.
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
		void checkpointReady(scour::gui::CheckpointStatePtr state, qint64 tag);

	private:
		void alloc_display();     // lazy device snapshot buffers (needs the grid)
		void publish_display();   // copy core fields -> snapshots (under disp_mtx_)
		void free_display();
		void apply_pending_rebuild(); // worker-thread: swap in a new core if one was requested
		void apply_seabed_conversion(); // worker-thread: attach engine + update_solid, flow preserved
		void apply_structure_update();  // worker-thread: swap the rigid structure + re-mask (G3.3 sink)
		void maybe_publish_host_flow(); // worker-thread: throttled D2H of {u,v,w,solid} for arrows
		void emit_checkpoint(long long tag); // worker-thread: gather full state → emit checkpointReady
		void publish_mask(const std::vector<unsigned char>& mask); // worker-thread: snapshot the mask + bump gen
		// Classify a flow solid mask into a KIND mask (0=fluid, 1=erodible sediment/sand, 2=rigid solid:
		// structure/obstacle) and publish it, so the GUI voxel overlay can colour the two apart. With a
		// morpho engine, sand = solid AND NOT structure; without one (fluid viewer) every solid is rigid.
		void publish_solid_kind(const std::vector<unsigned char>& solid); // worker-thread

		std::unique_ptr<scour::core::ChannelFluidCore> core_;

		// Erodible seabed morphodynamic engine (null ⇒ fluid-only viewer, G1 unchanged). Stepped on
		// the worker thread after the flow spin-up; its bed elevation snapshot feeds the bed viz.
		std::unique_ptr<scour::core::SeabedMorpho> morpho_;
		int morpho_spinup_ = 0;
		double bed_z0_ = 0.0, morpho_V0_ = 0.0;
		std::vector<unsigned char> morpho_solid_; // remask scratch (worker thread)
		std::vector<float> disp_zb_;              // published bed elevation (under disp_mtx_)
		std::vector<float> disp_exch_;            // published net bed-exchange rate [m/s] (under disp_mtx_)

		// Live flow solid-mask snapshot (host), republished only when the mask changes (see above).
		std::mutex mask_mtx_;
		std::vector<unsigned char> mask_snapshot_;
		scour::core::MacGrid mask_grid_{};
		std::atomic<std::uint64_t> mask_gen_{ 0 };
		std::atomic<std::uint64_t> bed_gen_{ 0 }; // bumped on each bed-elevation republish (video capture trigger)

		// Pending core rebuild (obstacle re-injection). Guarded by rebuild_mtx_; flagged atomic.
		std::function<std::unique_ptr<scour::core::ChannelFluidCore>(const std::vector<unsigned char>&, int)> factory_;
		std::mutex rebuild_mtx_;
		std::vector<unsigned char> pending_mask_;
		int pending_mode_ = 0;
		std::atomic<bool> rebuild_pending_{ false };

		// Pending live seabed conversion (GUI "Add sand" — preserves the running flow). Guarded by
		// convert_mtx_; flagged atomic. Applied on the worker thread by apply_seabed_conversion().
		std::mutex convert_mtx_;
		std::unique_ptr<scour::core::SeabedMorpho> pending_engine_;
		std::vector<unsigned char> pending_convert_solid_;
		int pending_convert_spinup_ = 0;
		bool pending_convert_loglaw_ = false;
		double pending_convert_z0_ = 0.0, pending_convert_datum_ = 0.0;
		std::atomic<bool> convert_pending_{ false };

		// Pending live structure update (G3.3 sink coupling). Guarded by structure_mtx_; flagged atomic.
		std::mutex structure_mtx_;
		std::vector<unsigned char> pending_structure_;
		std::atomic<bool> structure_pending_{ false };

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

		// Live suspended-sediment boundary switch (GUI "Sediment BC"): applied on the worker thread to the
		// morpho engine. Default OPEN (open-sea equilibrium inflow); a no-op when no bed is active.
		std::atomic<int> sed_bc_{ scour::core::SED_BC_OPEN };
		std::atomic<bool> sed_bc_dirty_{ false };

		// Live MORFAC override (GUI "Bed speed-up"): applied on the worker thread to the morpho engine. Held
		// dirty until a bed is active (the apply short-circuits on a null engine), so a value chosen before
		// "Add sand" lands the step morphology attaches. morfac_ = 0 ⇒ never set (engine keeps its config M).
		std::atomic<double> morfac_{ 0.0 };
		std::atomic<bool> morfac_dirty_{ false };

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
		double* dc_ = nullptr; // suspended-sediment concentration snapshot (seabed runs; 0 otherwise)
		unsigned char* ds_ = nullptr; // solid-cell snapshot (bed + structure + obstacle), for auto-range
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
		scour::core::MacGrid flow_grid_front_{};
		int flow_sign_front_ = 1; // current inlet direction snapshotted with the flow (guarded by flow_mtx_)
		std::chrono::steady_clock::time_point last_flow_pub_{};
	};
}
