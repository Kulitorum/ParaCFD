// sim_worker.cpp — see sim_worker.h. The step loop publishes post-step device snapshots
// so the main-thread slice kernel can render at display rate without blocking on step().
#include "gui/sim_worker.h"

#include <QThread>

#include <cuda_runtime.h>

#include <cmath>
#include <cstdio>
#include <vector>

namespace scour::gui
{
	using scour::core::MacGrid;

	namespace
	{
		// Reversing tidal current U_d(t) = U_max·s(t): start at +U_max, hold a plateau, then cosine-ramp
		// to the opposite sign over `ramp` (the compressed slack, zero-crossing at its midpoint), repeat
		// (RESEARCH §8). Mirrors the gate's Tide schedule. ramp≤0 ⇒ an instantaneous flip at each plateau.
		double tide_u_signed(double t, double Umax, double plateau, double ramp)
		{
			double half = plateau + (ramp > 0.0 ? ramp : 0.0);
			if (half <= 0.0) return Umax;
			int n = (int)std::floor(t / half);
			double base = (n % 2 == 0) ? 1.0 : -1.0;
			double local = t - n * half;
			double s = (local < plateau || ramp <= 0.0) ? base : base * std::cos(3.14159265358979323846 * (local - plateau) / ramp);
			return Umax * s;
		}
	}

	SimWorker::~SimWorker() { free_display(); }

	void SimWorker::setMorpho(std::unique_ptr<scour::core::SeabedMorpho> m, int spinup)
	{
		morpho_ = std::move(m);
		morpho_spinup_ = spinup;
		bed_z0_ = morpho_ ? morpho_->sand_depth() : 0.0;
		if (morpho_) morpho_V0_ = morpho_->bed_volume() + morpho_->susp_volume();
	}

	bool SimWorker::copyBed(std::vector<float>& zb)
	{
		if (!morpho_) return false;
		std::lock_guard<std::mutex> lk(disp_mtx_);
		if (disp_zb_.empty()) return false;
		zb = disp_zb_;
		return true;
	}

	bool SimWorker::copyBedExchange(std::vector<float>& rate)
	{
		if (!morpho_) return false;
		std::lock_guard<std::mutex> lk(disp_mtx_);
		if (disp_exch_.empty()) return false;
		rate = disp_exch_;
		return true;
	}

	void SimWorker::publish_mask(const std::vector<unsigned char>& mask)
	{
		std::lock_guard<std::mutex> lk(mask_mtx_);
		mask_snapshot_ = mask;
		mask_grid_ = core_->grid();
		mask_gen_.fetch_add(1);
	}

	void SimWorker::publish_solid_kind(const std::vector<unsigned char>& solid)
	{
		// Kind codes: 0 = fluid, 1 = erodible sediment (sand), 2 = rigid solid (structure/obstacle).
		// The overlay's exposed-face test still treats any non-zero as solid, so the staircase is the
		// union; the per-cell code just picks the colour. Sand = solid AND NOT structure.
		const bool have_bed = (morpho_ != nullptr);
		const std::vector<unsigned char>* structure = have_bed ? &morpho_->structure_host() : nullptr;
		std::vector<unsigned char> kind(solid.size(), 0);
		for (std::size_t n = 0; n < solid.size(); ++n)
		{
			if (!solid[n]) continue;
			const bool is_struct = structure && n < structure->size() && (*structure)[n];
			kind[n] = is_struct ? 2 : (have_bed ? 1 : 2); // no bed ⇒ all solids are rigid (obstacle)
		}
		publish_mask(kind);
	}

	bool SimWorker::copyMask(std::vector<unsigned char>& mask, scour::core::MacGrid& grid)
	{
		std::lock_guard<std::mutex> lk(mask_mtx_);
		if (mask_snapshot_.empty()) return false;
		mask = mask_snapshot_;
		grid = mask_grid_;
		return true;
	}

	void SimWorker::requestRebuild(std::vector<unsigned char> solid, int solid_mode)
	{
		{
			std::lock_guard<std::mutex> lk(rebuild_mtx_);
			pending_mask_ = std::move(solid);
			pending_mode_ = solid_mode;
		}
		rebuild_pending_.store(true);
	}

	void SimWorker::requestSeabedConversion(std::unique_ptr<scour::core::SeabedMorpho> engine,
		std::vector<unsigned char> initial_solid, int spinup, bool loglaw, double z0, double bed_datum)
	{
		{
			std::lock_guard<std::mutex> lk(convert_mtx_);
			pending_engine_ = std::move(engine);
			pending_convert_solid_ = std::move(initial_solid);
			pending_convert_spinup_ = spinup;
			pending_convert_loglaw_ = loglaw;
			pending_convert_z0_ = z0;
			pending_convert_datum_ = bed_datum;
		}
		convert_pending_.store(true);
	}

	void SimWorker::apply_seabed_conversion()
	{
		if (!convert_pending_.exchange(false)) return;
		std::unique_ptr<scour::core::SeabedMorpho> engine;
		std::vector<unsigned char> solid;
		int spinup = 0; bool loglaw = false; double z0 = 0.0, datum = 0.0;
		{
			std::lock_guard<std::mutex> lk(convert_mtx_);
			engine = std::move(pending_engine_);
			solid.swap(pending_convert_solid_);
			spinup = pending_convert_spinup_;
			loglaw = pending_convert_loglaw_;
			z0 = pending_convert_z0_;
			datum = pending_convert_datum_;
		}
		if (!engine || !core_ || (int)solid.size() != core_->grid().p_count()) return;

		// Install the sand+structure mask IN PLACE — the developed flow is preserved in the untouched
		// cells; only the now-solid cells are zeroed. No factory rebuild, no init_uniform, no t=0 reset.
		core_->update_solid(solid);
		core_->set_bed_inlet_mask(true);              // the erodible bed now reaches the inlet plane
		if (loglaw) core_->set_inlet_profile(true, z0, datum); // seabed inlet BL (leading-edge fix)

		// Attach the (already-built) engine; morphology waits `spinup` MORE steps so the re-masked flow
		// settles first. Sim time/step count are NOT reset — the run continues from its developed state.
		morpho_ = std::move(engine);
		morpho_spinup_ = steps_.load() + spinup;
		bed_z0_ = morpho_->sand_depth();
		morpho_V0_ = morpho_->bed_volume() + morpho_->susp_volume();

		publish_display();               // reflect the zeroed sand cells immediately
		publish_solid_kind(solid);       // overlay: sand (orange) vs structure (dark-blue)
		std::fprintf(stderr, "[seabed] live conversion: sand+structure installed via update_solid "
			"(flow preserved), morphology begins at step %d\n", (int)morpho_spinup_);
	}

	void SimWorker::apply_structure_update()
	{
		if (!structure_pending_.exchange(false)) return;
		if (!morpho_ || !core_) return; // no bed yet: drop the update (the pile isn't erodible-coupled)
		std::vector<unsigned char> s;
		{ std::lock_guard<std::mutex> lk(structure_mtx_); s.swap(pending_structure_); }
		if ((int)s.size() != core_->grid().p_count()) return;

		// Swap the rigid structure in the engine → new flow solid mask (sand ∪ new structure); re-mask the
		// flow IN PLACE so the developed wake is preserved (untouched cells keep {u,v,w,p}).
		std::vector<unsigned char> new_solid = morpho_->set_structure(s);
		if ((int)new_solid.size() != core_->grid().p_count()) return;
		core_->update_solid(new_solid);
		publish_display();
		publish_solid_kind(new_solid); // overlay: re-extract the pile + bed at their new positions
		std::fprintf(stderr, "[seabed] structure updated (drop re-settle): flow re-masked in place\n");
	}

	void SimWorker::apply_pending_rebuild()
	{
		if (!rebuild_pending_.exchange(false)) return;
		std::vector<unsigned char> mask;
		int mode = 0;
		{
			std::lock_guard<std::mutex> lk(rebuild_mtx_);
			mask.swap(pending_mask_);
			mode = pending_mode_;
		}
		if (!factory_) return;
		auto nc = factory_(mask, mode); // constructed + flow-initialised ON THIS (worker) thread
		if (!nc) return;
		core_ = std::move(nc);          // old core (its CUDA buffers) freed here
		if (inlet_override_.load()) core_->set_inlet_speed(inlet_speed_.load()); // keep the chosen current
		if (inlet_prof_override_.load()) // keep the chosen inlet profile across a rebuild
		{
			bool loglaw; double z0, datum;
			{ std::lock_guard<std::mutex> lk(inlet_prof_mtx_); loglaw = inlet_prof_loglaw_; z0 = inlet_prof_z0_; datum = inlet_prof_datum_; }
			core_->set_inlet_profile(loglaw, z0, datum);
		}
		steps_.store(0);
		sim_time_.store(0.0);
		publish_display();              // show the fresh initial condition immediately
		publish_solid_kind(mask);       // overlay follows the new obstacle mask (rigid ⇒ dark-blue)
	}

	void SimWorker::computeDiversion(DiversionReport& out) const
	{
		out = DiversionReport{};
		if (!core_) return;
		MacGrid g = core_->grid();
		std::vector<double> u((size_t)g.u_count()), v((size_t)g.v_count()), w((size_t)g.w_count());
		std::vector<unsigned char> solid((size_t)g.p_count());
		cudaMemcpy(u.data(), core_->u_dev(), sizeof(double) * u.size(), cudaMemcpyDeviceToHost);
		cudaMemcpy(v.data(), core_->v_dev(), sizeof(double) * v.size(), cudaMemcpyDeviceToHost);
		cudaMemcpy(w.data(), core_->w_dev(), sizeof(double) * w.size(), cudaMemcpyDeviceToHost);
		cudaMemcpy(solid.data(), core_->solid_dev(), solid.size(), cudaMemcpyDeviceToHost);

		long long nsolid = 0;
		double sum_solid = 0.0, max_fluid = 0.0;
		for (int k = 0; k < g.nz; ++k)
			for (int j = 0; j < g.ny; ++j)
				for (int i = 0; i < g.nx; ++i)
				{
					// Cell-centred velocity = average of the two bracketing MAC faces.
					double uc = 0.5 * (u[g.uidx(i, j, k)] + u[g.uidx(i + 1, j, k)]);
					double vc = 0.5 * (v[g.vidx(i, j, k)] + v[g.vidx(i, j + 1, k)]);
					double wc = 0.5 * (w[g.widx(i, j, k)] + w[g.widx(i, j, k + 1)]);
					double sp = std::sqrt(uc * uc + vc * vc + wc * wc);
					if (solid[g.pidx(i, j, k)]) { ++nsolid; sum_solid += sp; }
					else if (sp > max_fluid) max_fluid = sp;
				}
		out.active = nsolid > 0;
		out.solid_cells = nsolid;
		out.mean_speed_solid = nsolid > 0 ? sum_solid / (double)nsolid : 0.0;
		out.max_speed_fluid = max_fluid;
	}

	void SimWorker::emit_checkpoint(long long tag)
	{
		if (!core_) return;
		auto s = std::make_shared<CheckpointState>();
		s->grid = core_->grid();
		s->steps = steps_.load();
		s->sim_time = sim_time_.load();
		s->bc = core_->bc();
		s->solid_mode = s->bc.solid_mode;
		s->bed_inlet_mask = core_->bed_inlet_mask_on();
		core_->copy_state_host(s->u, s->v, s->w, s->p); // D2H velocities + pressure (worker thread ⇒ race-free)
		s->solid.resize((size_t)s->grid.p_count());
		cudaMemcpy(s->solid.data(), core_->solid_dev(), s->solid.size(), cudaMemcpyDeviceToHost);
		if (morpho_)
		{
			s->has_bed = true;
			s->sp = morpho_->params();
			morpho_->save_state(s->bed_G, s->susp_c, s->bed_emax, s->bed_emay, s->morpho_steps);
			s->structure = morpho_->structure_host();
		}
		emit checkpointReady(s, (qint64)tag);
	}

	void SimWorker::maybe_publish_host_flow()
	{
		if (!want_host_flow_.load()) return;
		// Throttle to ~12 Hz: the flow evolves far slower than the 60 Hz repaint, so a coarse
		// refresh is plenty and keeps the D2H cost off the critical path.
		auto now = std::chrono::steady_clock::now();
		if (last_flow_pub_.time_since_epoch().count() != 0 &&
			std::chrono::duration<double>(now - last_flow_pub_).count() < 0.08)
			return;
		last_flow_pub_ = now;

		MacGrid g = core_->grid();
		fu_back_.resize((size_t)g.u_count());
		fv_back_.resize((size_t)g.v_count());
		fw_back_.resize((size_t)g.w_count());
		fs_back_.resize((size_t)g.p_count());
		// Read the LIVE core fields (this worker owns the core; no lock needed). The D2H copies are
		// default-stream and therefore ordered after the just-finished step()/morpho work.
		cudaMemcpy(fu_back_.data(), core_->u_dev(), sizeof(double) * fu_back_.size(), cudaMemcpyDeviceToHost);
		cudaMemcpy(fv_back_.data(), core_->v_dev(), sizeof(double) * fv_back_.size(), cudaMemcpyDeviceToHost);
		cudaMemcpy(fw_back_.data(), core_->w_dev(), sizeof(double) * fw_back_.size(), cudaMemcpyDeviceToHost);
		cudaMemcpy(fs_back_.data(), core_->solid_dev(), fs_back_.size(), cudaMemcpyDeviceToHost);
		const int fsgn = core_->flow_direction(); // inlet side, snapshotted with these fields (worker owns core)

		{
			std::lock_guard<std::mutex> lk(flow_mtx_);
			fu_front_.swap(fu_back_);
			fv_front_.swap(fv_back_);
			fw_front_.swap(fw_back_);
			fs_front_.swap(fs_back_);
			flow_grid_front_ = g;
			flow_sign_front_ = fsgn;
			flow_ready_.store(true);
		}
	}

	bool SimWorker::withFlowField(const std::function<void(const FlowField&)>& fn)
	{
		std::lock_guard<std::mutex> lk(flow_mtx_);
		if (!flow_ready_.load() || fu_front_.empty()) return false;
		FlowField f;
		f.u = fu_front_.data();
		f.v = fv_front_.data();
		f.w = fw_front_.data();
		f.solid = fs_front_.empty() ? nullptr : fs_front_.data();
		f.grid = flow_grid_front_;
		f.flow_sign = flow_sign_front_;
		fn(f);
		return true;
	}

	void SimWorker::alloc_display()
	{
		MacGrid g = core_->grid();
		cudaMalloc(&du_, sizeof(double) * (size_t)g.u_count());
		cudaMalloc(&dv_, sizeof(double) * (size_t)g.v_count());
		cudaMalloc(&dw_, sizeof(double) * (size_t)g.w_count());
		cudaMalloc(&dp_, sizeof(double) * (size_t)g.p_count());
		cudaMalloc(&dc_, sizeof(double) * (size_t)g.p_count());
		cudaMalloc(&ds_, (size_t)g.p_count());
		cudaMemset(du_, 0, sizeof(double) * (size_t)g.u_count());
		cudaMemset(dv_, 0, sizeof(double) * (size_t)g.v_count());
		cudaMemset(dw_, 0, sizeof(double) * (size_t)g.w_count());
		cudaMemset(dp_, 0, sizeof(double) * (size_t)g.p_count());
		cudaMemset(dc_, 0, sizeof(double) * (size_t)g.p_count());
		cudaMemset(ds_, 0, (size_t)g.p_count());
	}

	void SimWorker::free_display()
	{
		if (du_) cudaFree(du_);
		if (dv_) cudaFree(dv_);
		if (dw_) cudaFree(dw_);
		if (dp_) cudaFree(dp_);
		if (dc_) cudaFree(dc_);
		if (ds_) cudaFree(ds_);
		du_ = dv_ = dw_ = dp_ = dc_ = nullptr;
		ds_ = nullptr;
	}

	void SimWorker::publish_display()
	{
		MacGrid g = core_->grid();
		std::lock_guard<std::mutex> lk(disp_mtx_);
		cudaMemcpy(du_, core_->u_dev(), sizeof(double) * (size_t)g.u_count(), cudaMemcpyDeviceToDevice);
		cudaMemcpy(dv_, core_->v_dev(), sizeof(double) * (size_t)g.v_count(), cudaMemcpyDeviceToDevice);
		cudaMemcpy(dw_, core_->w_dev(), sizeof(double) * (size_t)g.w_count(), cudaMemcpyDeviceToDevice);
		cudaMemcpy(dp_, core_->p_dev(), sizeof(double) * (size_t)g.p_count(), cudaMemcpyDeviceToDevice);
		cudaMemcpy(ds_, core_->solid_dev(), (size_t)g.p_count(), cudaMemcpyDeviceToDevice); // solid mask for auto-range
		if (morpho_)
		{
			cudaMemcpy(dc_, morpho_->c_device(), sizeof(double) * (size_t)g.p_count(), cudaMemcpyDeviceToDevice); // suspended c
			morpho_->copy_zb_host(disp_zb_); // per-column bed elevation for the bed viz
			morpho_->copy_exchange_host(disp_exch_); // per-column net exchange rate for the bed-colour viz
			bed_gen_.fetch_add(1); // signal the GUI recorder that the bed snapshot advanced
		}
		cudaStreamSynchronize(0); // ensure snapshots complete before a render can read them
		disp_ready_.store(true);
	}

	void SimWorker::run()
	{
		alloc_display();
		publish_display(); // show the initial condition before the first step

		// Publish the initial flow solid mask once (a single D2H) so the voxel overlay shows the
		// starting bed/structure/obstacle; thereafter it is republished only on a mask change.
		{
			std::vector<unsigned char> m0((size_t)core_->grid().p_count());
			cudaMemcpy(m0.data(), core_->solid_dev(), m0.size(), cudaMemcpyDeviceToHost);
			publish_solid_kind(m0);
		}

		int since_stats = 0;
		double last_dt = core_->last_dt();
		// Sim-throughput (steps/second) measured on THIS thread across each stats window, so it reflects
		// the true step rate regardless of the display-publish cadence (it climbs under "Fast sim").
		std::chrono::steady_clock::time_point last_stats_time{};
		long long last_stats_steps = 0;
		double sim_fps = 0.0;
		while (!quit_.load())
		{
			// Obstacle re-injection: rebuild the core on THIS thread before anything else, so a
			// paused viewer still picks up the new mask (and never races the main/GL thread).
			if (rebuild_pending_.load()) apply_pending_rebuild();

			// Live "Add sand": convert a running fluid viewer to a seabed run in place (flow preserved).
			if (convert_pending_.load()) apply_seabed_conversion();

			// Live structure update (G3.3 sink coupling): install the re-settled pile + re-mask the flow.
			if (structure_pending_.load()) apply_structure_update();

			// Live current change (GUI "Input speed"): update the running core's inlet before stepping.
			// A paused viewer picks it up too (the wake evolves on the next play/step).
			if (inlet_speed_dirty_.exchange(false) && inlet_override_.load())
			{
				core_->set_inlet_speed(inlet_speed_.load());
				if (morpho_) morpho_->set_inlet_speed(inlet_speed_.load()); // re-derive the open-sea equilibrium inflow
			}

			// Live inlet-profile switch (GUI "Inlet profile": uniform ↔ boundary-layer).
			if (inlet_prof_dirty_.exchange(false))
			{
				bool loglaw; double z0, datum;
				{ std::lock_guard<std::mutex> lk(inlet_prof_mtx_); loglaw = inlet_prof_loglaw_; z0 = inlet_prof_z0_; datum = inlet_prof_datum_; }
				core_->set_inlet_profile(loglaw, z0, datum);
			}

			// Live suspended-sediment boundary switch (GUI "Sediment BC": open-sea / recycle / closed).
			if (sed_bc_dirty_.exchange(false) && morpho_) morpho_->set_boundary_mode(sed_bc_.load());

			// Live MORFAC change (GUI "Bed speed-up"): retune the morpho engine's acceleration so the bed
			// catches up faster. The `morpho_ &&` short-circuit keeps the flag pending until a bed exists,
			// so a value dialed in before "Add sand" still lands once the engine attaches this same iteration.
			if (morpho_ && morfac_dirty_.exchange(false)) morpho_->set_morfac(morfac_.load());

			// Live tidal reversal (GUI "Tidal reversal"): an enable/disable/param change restarts the tide
			// clock; on DISABLE, restore steady forward flow at the manual inlet speed (so the view returns
			// to a sensible +x current instead of freezing at whatever phase the tide was in).
			if (tidal_dirty_.exchange(false))
			{
				tide_clock_ = 0.0;
				if (!tidal_on_.load())
				{
					core_->set_flow_direction(1);
					double U = inlet_override_.load() ? inlet_speed_.load() : tidal_umax_.load();
					core_->set_inlet_speed(U);
					if (morpho_) morpho_->set_inlet_speed(U);
					tide_phase_.store(0.0);
				}
			}

			// Manual checkpoint (GUI "Save Scene"): gather the full state on THIS thread — race-free — and
			// emit it for the main thread to write. Serviced even while paused (it's before the do_step gate).
			if (snapshot_pending_.exchange(false)) emit_checkpoint(snapshot_tag_.load());

			bool do_step = playing_.load();
			int req = step_requests_.load();
			if (!do_step && req > 0) do_step = true;

			if (!do_step)
			{
				// Publish the final state once when we go idle (pause) even if the throttle window has
				// not elapsed, so a paused view is never stale ("Fast sim" mode).
				if (display_dirty_) { publish_display(); display_dirty_ = false; last_display_pub_ = std::chrono::steady_clock::now(); }
				last_stats_time = {}; since_stats = 0; // re-baseline sim-fps on resume (don't count paused wall time)
				QThread::msleep(3); // idle when paused (keep the event loop responsive)
				continue;
			}

			// Tidal reversal: drive the reversing current for THIS step (before step()), face-swapping the
			// inlet/outlet at slack. Overrides the manual inlet speed while active. Clock advances post-step.
			if (tidal_on_.load())
			{
				double Ud = tide_u_signed(tide_clock_, tidal_umax_.load(), tidal_plateau_.load(), tidal_ramp_.load());
				core_->set_inlet_speed(std::fabs(Ud));
				core_->set_flow_direction(Ud >= 0.0 ? 1 : -1);
				if (morpho_) morpho_->set_inlet_speed(std::fabs(Ud));
				tide_phase_.store(Ud);
			}

			last_dt = core_->step(); // heavy CUDA work — NOT under the display lock
			if (tidal_on_.load()) tide_clock_ += last_dt; // advance the tide by the flow-time actually stepped

			// Erodible seabed: after the flow spins up, advance the morphodynamic loop and re-mask the
			// flow core in place as the bed moves (flow preserved). NOT under the display lock.
			long long n = steps_.load();
			if (morpho_ && n >= morpho_spinup_)
			{
				if (morpho_->step(core_->u_dev(), core_->v_dev(), core_->w_dev(), core_->nut_dev(), last_dt, &morpho_solid_))
				{
					core_->update_solid(morpho_solid_);
					publish_solid_kind(morpho_solid_); // overlay follows the scouring bed's re-mask (sand/structure split)
				}
				if (n % 200 == 0)
				{
					double V = morpho_->bed_volume() + morpho_->susp_volume();
					double bin = morpho_->boundary_sand_in();  // net open-sea sand import [m³] (0 when CLOSED)
					double err = morpho_V0_ > 0 ? std::fabs(V - morpho_V0_ - bin) / morpho_V0_ : 0.0; // budget incl. boundary
					std::fprintf(stderr, "[seabed] step %lld: z_b[min=%.4f max=%.4f] (z_b0=%.3f) | max_ustar=%.4f | "
						"boundary sand=%.3e m3 | budget err=%.2e\n", n, morpho_->zb_min(), morpho_->zb_max(), bed_z0_, morpho_->max_ustar(), bin, err);
				}
			}
			// Display snapshot — throttled by the GUI "Fast sim" interval. publish_display() copies the
			// fields and hard-syncs the device; at a few steps/s that overhead is pure waste when the
			// interval says nobody is watching. Off (0) ⇒ every step (byte-identical to G1). A manual Step
			// always shows its result (req > 0); the idle branch above catches the last frame on pause.
			display_dirty_ = true;
			const double disp_iv = display_interval_.load();
			bool want_pub = (disp_iv <= 0.0) || (req > 0);
			if (!want_pub)
			{
				const auto now = std::chrono::steady_clock::now();
				want_pub = (last_display_pub_.time_since_epoch().count() == 0) ||
					(std::chrono::duration<double>(now - last_display_pub_).count() >= disp_iv);
			}
			if (want_pub)
			{
				publish_display();
				maybe_publish_host_flow(); // throttled host {u,v,w,solid} for the arrow tracers (no-op if off)
				last_display_pub_ = std::chrono::steady_clock::now();
				display_dirty_ = false;
			}

			steps_.fetch_add(1);
			sim_time_.store(sim_time_.load() + last_dt);
			if (req > 0) step_requests_.fetch_sub(1);

			// Auto-save a restart-point checkpoint every `autosave_interval_` steps (0 = off).
			if (const int iv = autosave_interval_.load(); iv > 0)
			{
				const long long n2 = steps_.load();
				if (n2 > 0 && (n2 % iv) == 0) emit_checkpoint(n2);
			}

			if (++since_stats >= 8)
			{
				const auto now = std::chrono::steady_clock::now();
				const long long nsteps = steps_.load();
				if (last_stats_time.time_since_epoch().count() != 0)
				{
					const double wall = std::chrono::duration<double>(now - last_stats_time).count();
					if (wall > 0.0)
					{
						const double inst = (double)(nsteps - last_stats_steps) / wall;
						sim_fps = (sim_fps > 0.0) ? 0.7 * sim_fps + 0.3 * inst : inst; // light EMA
					}
				}
				last_stats_time = now;
				last_stats_steps = nsteps;
				since_stats = 0;
				emit stats((qint64)nsteps, sim_time_.load(), last_dt, sim_fps);
			}
		}
		emit stats((qint64)steps_.load(), sim_time_.load(), last_dt, sim_fps);
	}
}
