// main_window.h — G1 top-level window. Hosts the SliceViewer (all GL on this, the main
// thread) and drives a SimWorker on a dedicated QThread (queued-connection hand-off,
// cobod-slicer mainwindow.cpp:1691 pattern). Controls live in the left settings dock;
// status bar: steps / sim-time / fps. Loading a STEP model auto-voxelizes it as the flow
// obstacle (no separate toggle).
#pragma once

#include "core/geometry/model_placement.h"
#include "core/geometry/tri_mesh.h"
#include "gui/scene_io.h"    // SceneFile / CheckpointStatePtr
#include "gui/sim_setup.h"
#include "gui/sim_worker.h" // DiversionReport
#ifdef SCOUR_HAVE_JOLT
#include "core/physics/settle.h" // SettleWorld — the drop/settle preparation phase (G3.2)
#endif

#include <QMainWindow>

#include <cstdint>
#include <memory>
#include <vector>

class QAction;
class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QGroupBox;
class QLabel;
class QMenu;
class QPushButton;
class QThread;
class QTimer;

namespace scour::core { class ChannelFluidCore; }

namespace scour::gui
{
	class SliceViewer;
	class VideoRecorder;
	class VideoSettingsDialog;

	class MainWindow : public QMainWindow
	{
		Q_OBJECT
	public:
		// Takes ownership of the initial core; `recipe` carries its domain (recipe.info) for the
		// viewer AND the ingredients to rebuild the core when a model is (de)injected as obstacle.
		// `scen` (optional) attaches the erodible-seabed morphodynamic engine + bed viz; when active,
		// the worker runs the full morphodynamic loop and the structure sits on the sand bed.
		MainWindow(std::unique_ptr<scour::core::ChannelFluidCore> core, const SimRecipe& recipe,
			SeabedScenario scen = {}, QWidget* parent = nullptr);
		~MainWindow() override;

		SliceViewer* viewer() { return viewer_; }

		// Load a STEP model (OpenCascade → metre-scale TriMesh), display it, and — in a fluid
		// viewer — auto-voxelize it into the flow as the solid obstacle, REPLACING the config
		// obstacle (the loaded model IS the obstacle; there is no separate toggle). In a seabed
		// scenario the bed owns the mask, so the model is only displayed. `noslip` selects
		// SOLID_NOSLIP; the shipping default is SOLID_FREESLIP (RESEARCH §3). Used by the File
		// menu and the --load-step CLI flag. Logs a summary; returns false on load failure.
		// Safe to call before the window is shown.
		bool loadStepFile(const QString& path, bool noslip = false);

		// Stop the worker (if running) and, if a voxel obstacle was injected, sample the flow-
		// diversion metrics. Idempotent. Call before reading diversion() from main().
		void finalize();
		const DiversionReport& diversion() const { return diversion_; }

		// Latest worker totals (for the scripted-gate summary in main.cpp).
		long long lastSteps() const { return last_steps_; }
		double lastSimTime() const { return last_sim_time_; }

	public:
		// Load a seabed scenario JSON at runtime (File menu): tear down the current sim and rebuild the
		// core + morphodynamic engine + bed viz on the (possibly new) grid. Returns false on error.
		bool startSeabedScenario(const QString& path);

		// Programmatically apply a new grid (domain + voxel/cell size h) and rebuild+reset the sim,
		// exactly as the dock's Apply button. Any non-positive argument keeps the current sim's value.
		// Used by the `--apply-grid` / `--apply-h` CLI smoke hook so the reset-at-new-grid feature is
		// testable headlessly (mirrors --load-step / --voxelize).
		void applyGridValues(double Lx, double Ly, double Lz, double h);

		// Programmatically set the inlet current speed U [m/s], exactly as the "Input speed" spin box:
		// applies live to the running flow (no reset). Used by the `--set-u` CLI smoke hook so the
		// different-currents feature is testable headlessly (mirrors --apply-h).
		void setInputSpeed(double U);

		// Programmatically add `depth` m of sand and rebuild as a sediment run, exactly as the "Add sand"
		// button. Used by the `--add-sand <m>` CLI smoke hook.
		void addSandMeters(double depth);

#ifdef SCOUR_HAVE_JOLT
		// Drop/settle PREPARATION PHASE (PLAN G3.2): load a STEP protection unit (e.g. XStone), scatter
		// `count` scaled copies above a `sandDepth`-m sand bed, ANIMATE them falling + settling under
		// gravity (Jolt SettleWorld, stepped on this main thread via a timer), then voxelize the settled
		// pile as the rigid structure and start the morphodynamic run. Drives the `--drop-step`/`--drop`
		// CLI smoke (and, later, a dock button). Returns false if the STEP fails to load.
		bool dropBlocks(const QString& stepPath, int count, unsigned seed,
			double scaleMin, double scaleMax, double sandDepth);
#endif

		// Programmatically set MORFAC (morphological acceleration) live, exactly as the "MORFAC (bed
		// speed-up)" control: the bed then evolves M× faster than the flow clock. Used by the `--morfac
		// <M>` CLI smoke hook (a no-op without an active seabed run).
		void setMorfacValue(double m);

		// Programmatically move/rotate/scale the loaded model's placement (the gizmo), then re-voxelize it
		// as the obstacle — exactly what a gizmo drag + Apply do. dx/dy/dz metres, rzDeg degrees about +Z,
		// scale a uniform factor (<= 0 ⇒ unchanged). Used by the `--model-place dx,dy,dz[,rz[,scale]]` hook.
		void nudgeModelPlacement(double dx, double dy, double dz, double rzDeg, double scale);

		// Programmatically enable/disable the live tidal-reversal driver, exactly as the "Tidal reversal"
		// controls: drives a reversing current U_d(t) (peak U_max, plateau/ramp seconds) live (no reset).
		// Used by the `--tidal Umax,plateau,ramp` CLI smoke hook.
		void setTidalReversal(bool on, double Umax, double plateau_s, double ramp_s);

		// Throttle on-screen graphics updates to at most once per `seconds` of wall time (the "Fast sim"
		// combo): the sim keeps stepping at full rate while the view — and the worker's per-step display
		// snapshot, which sits on the sim's critical path — refresh less often, trading smoothness for
		// sim speed when you're away from the screen. 0 = live 60 Hz. Re-applied across a rebuild/restore.
		// Used by the combo and the `--display-interval <s>` CLI smoke hook.
		void setDisplayThrottle(double seconds);

		// Start recording an MP4 of the run to `path`: a frame is captured whenever the sand bed changes
		// and streamed to ffmpeg (H.264). Used by the File menu ("Record Video…") and the `--record
		// <path>` CLI smoke hook. Returns false if ffmpeg is unavailable. Always finalized on close.
		bool startRecording(const QString& path);

		// --- Scene save / restore (scene_io) -----------------------------------------------------
		// Request a full-scene save to `path` (a .scn): the worker gathers the live state on its thread
		// and the write completes shortly after (onCheckpointReady). Returns false only on a bad setup.
		// Used by the "Save Scene As" menu and the `--save-scene <path>` CLI smoke hook.
		bool saveSceneNow(const QString& path);
		// Load + fully restore a .scn (setup + all field data), continuing from the saved step (not t=0).
		// Used by the "Load Scene" menu, Recent Files, and the `--load-scene <path>` CLI smoke hook.
		bool loadSceneFromPath(const QString& path);
		// Set the auto-save cadence in steps (0 = off). Used by the Auto-save menu + `--autosave <n>`.
		void setAutosave(int stepsInterval);

	protected:
		void closeEvent(QCloseEvent* e) override;

	private:
		// Voxelize the loaded model on the sim grid and inject it as the solid obstacle (replacing
		// the config obstacle), or — when off — restore the config obstacle. Rebuilds the core on
		// the worker thread. `noslip` selects SOLID_NOSLIP (default SOLID_FREESLIP, RESEARCH §3).
		// Internal: driven by loadStepFile (on), "Close model" (off) and the grid-Apply re-inject.
		bool setModelAsObstacle(bool on, bool noslip = false);

		void shutdownWorker();
		// `steps0`/`t0` prime the worker's counters (a scene restore resumes from the saved step).
		void spawnWorker(std::unique_ptr<scour::core::ChannelFluidCore> core, SeabedScenario scen,
			long long steps0 = 0, double t0 = 0.0);

		// --- Scene save / restore internals ------------------------------------------------------
		void saveSceneAs();     // File menu: pick a .scn path, then request a save
		void loadSceneDialog(); // File menu: pick a .scn (offering sibling restart points), then restore
		void restoreScene(SceneFile& sf); // tear down + rebuild core/engine + inject the saved state
		// Worker → main-thread hand-off: write the gathered state to disk. tag<0 = manual save (→ Recent
		// Files), tag>=0 = an auto-save at that step (→ <base>.<step>.scn, not Recent).
		void onCheckpointReady(scour::gui::CheckpointStatePtr state, qint64 tag);
		QString scene_base_path_;          // stem (no extension) of the current scene, for the auto-save series
		QString pending_manual_save_path_; // destination of an in-flight manual "Save Scene As"
		QString pending_recent_keep_;      // a just-saved scene optimistically in Recent Files before its
		                                   // async write lands — rebuildRecentMenu keeps it from being pruned
		int autosave_interval_ = 0;        // auto-save cadence in steps (0 = off); re-applied on each spawnWorker
		scour::core::TriMesh scene_mesh_;  // display mesh persisted into a saved scene (structure/obstacle model)
		scour::core::ModelPlacement scene_place_{}; // its display placement
		// Stop + join + delete the worker (and its thread) so a fresh core/engine can replace it. Shared
		// by the seabed-scenario load and the grid Apply (both fully rebuild the sim). No-op if none.
		void teardownWorkerForReload();

		// Build the left-side settings dock (grouped controls: Simulation / Visualization). Called
		// once from the constructor. The idea mirrors cobod-slicer's left panel (grouped, scrollable),
		// not its styling.
		void buildControlDock();

		// --- Editable domain + resolution (Simulation group) -------------------------------------
		// Apply rebuilds+resets the whole sim (fresh t=0) at the domain (Lx,Ly,Lz) and voxel/cell size
		// h in the spin boxes, preserving the current scenario + every physics parameter — ONLY the grid
		// (nx,ny,nz) changes. syncGridControls repopulates the boxes from the live sim; updateGridReadout
		// refreshes the live "→ nx×ny×nz = N cells" label (and guards against an over-fine grid).
		void applyGrid();
		void syncGridControls();
		void updateGridReadout();

		// --- Model-placement gizmo UI ------------------------------------------------------------
		// Enable the gizmo dock group + refresh its readout for the current model/structure (any gizmo-
		// editable model). Called on model load/close, seabed switch, grid Apply and scene restore.
		void updateGizmoUi();
		QGroupBox* gizmo_group_ = nullptr;
		QCheckBox* gizmo_enable_chk_ = nullptr; // enable the unified manipulator (all handles at once)
		QLabel* gizmo_info_ = nullptr;
		QDoubleSpinBox* lx_spin_ = nullptr;
		QDoubleSpinBox* ly_spin_ = nullptr;
		QDoubleSpinBox* lz_spin_ = nullptr;
		QDoubleSpinBox* h_spin_ = nullptr;
		QDoubleSpinBox* u_spin_ = nullptr; // inlet current speed U — live (no reset) + honoured by Apply
		QComboBox* inlet_profile_box_ = nullptr; // inlet profile: uniform (top-hat) ↔ boundary-layer (log-law)
		QComboBox* sed_bc_box_ = nullptr;        // suspended-sediment x-boundary: open-sea / recycle / closed
		QDoubleSpinBox* morfac_spin_ = nullptr;  // MORFAC (morphological acceleration) — live bed speed-up
		QComboBox* fast_sim_box_ = nullptr;      // "Fast sim": graphics-update throttle (Off / 0.5 / 1 / 3 s / paused)
		// Tidal reversal (live): enable + peak current / plateau-hold / slack-ramp durations. `syncTidal`
		// pushes the current control values to the worker; the enable box also greys the manual "Input
		// speed" (the tide drives the inlet while active).
		QCheckBox* tidal_chk_ = nullptr;
		QDoubleSpinBox* tidal_umax_spin_ = nullptr;
		QDoubleSpinBox* tidal_plateau_spin_ = nullptr;
		QDoubleSpinBox* tidal_ramp_spin_ = nullptr;
		void syncTidal();

		// "Add sand": fill the lower `sand_spin_` metres of every non-rigid column with sand (rigid
		// structure/obstacle stays solid) and reset the sim as an erodible-seabed run — turning even a
		// plain fluid viewer with a voxelized obstacle into a live sediment simulation. last_sp_ carries
		// the current scenario's sediment params (d50/α/…) across the rebuild; nominal defaults otherwise.
		void addSand();
		QDoubleSpinBox* sand_spin_ = nullptr;
		scour::core::SeabedParams last_sp_{};
		bool have_sp_ = false;
		// Sand depth [m] added via the "Add sand" button (0 ⇒ none). Non-zero marks a seabed run that was
		// CONVERTED from a fluid viewer at runtime (fluid config + runtime sand/model), as opposed to one
		// loaded from a real scenario FILE. A grid Apply uses it to rebuild the converted seabed correctly
		// (re-voxelize the model as the structure + re-apply the sand at the new grid) instead of falling
		// back to the fluid config's obstacle. Reset to 0 when a real scenario file is loaded.
		double added_sand_m_ = 0.0;

#ifdef SCOUR_HAVE_JOLT
		// Drop/settle preparation phase (G3.2). settle_world_ owns the live Jolt pile; drop_timer_ steps
		// it a few substeps per tick on the main thread and refreshes the block display; commitDrop()
		// voxelizes the settled pile + starts the seabed run. drop_mesh_ is the XStone display/voxelize
		// geometry; drop_placements_ the current per-block poses.
		void stepDropAnimation();
		void commitDrop();
		// G3.3 sink coupling: periodically re-settle the persistent pile onto the eroded bed and hand the
		// re-voxelized structure to the worker (flow preserved). Fires only when the bed has moved ≥ ½ cell.
		void maybeResettle();
		std::unique_ptr<scour::core::SettleWorld> settle_world_; // KEPT ALIVE after commit for re-settling
		QTimer* drop_timer_ = nullptr;      // initial fall animation
		QTimer* resettle_timer_ = nullptr;  // periodic sink re-settle (G3.3)
		scour::core::TriMesh drop_mesh_;
		std::vector<scour::core::ModelPlacement> drop_placements_;
		std::vector<scour::core::BlockInstance> drop_queue_; // units to release, sorted largest-first
		int drop_next_ = 0;         // index of the next unit to release
		int drop_release_ctr_ = 0;  // animation ticks until the next release
		std::vector<float> resettle_last_zb_; // bed baseline at the last re-settle (change trigger)
		double drop_sand_depth_ = 0.5;
		bool drop_active_ = false; // a drop-originated seabed run is live (gates maybeResettle)
		int drop_release_ticks_ = 12; // ticks between unit releases (~0.2 s at ~60 Hz); dock "Release cadence"
		void clearDropPile();        // stop the drop animation + sink re-settle and clear the pile display
#endif

		QLabel* grid_readout_ = nullptr;
		QPushButton* apply_btn_ = nullptr;
		QPushButton* play_btn_ = nullptr; // so a rebuild can honour the current play/pause state
		static constexpr long long kMaxCells = 40000000LL; // Apply guard: refuse an OOM-risking grid

		// Recent-files (feature 1): persisted via QSettings("COBOD","ScourProtection"). The submenu is
		// rebuilt on demand (pruning files that no longer exist), MRU-first, capped at kMaxRecent.
		void addRecentFile(const QString& path);
		void rebuildRecentMenu();
		static constexpr int kMaxRecent = 8;
		QMenu* recent_menu_ = nullptr;
		QCheckBox* slice_chk_ = nullptr; // "Show slice" (defaulted OFF in the seabed scenario)

		QTimer* repaint_timer_ = nullptr; // main-thread ~60 Hz repaint; interval widens under "Fast sim"
		double display_throttle_s_ = 0.0; // "Fast sim" graphics-update interval [s] (0 = live); re-applied on spawnWorker

		// Video recording (video_recorder.*): capture the viewer whenever the sand bed changes and stream
		// the frames to ffmpeg. The File-menu action opens a non-modal settings popup (openRecordDialog);
		// maybeCaptureFrame runs on the repaint tick and grabs a frame only when the bed generation advanced
		// with a real z_b change; updateRecordDialogStatus pushes the live frame count into the popup.
		void openRecordDialog();       // File menu: show the settings dialog (Start/Stop + params + live status)
		void updateRecordDialogStatus(); // refresh the dialog's "● REC — N frames" from the recorder
		void maybeAutoRecord();        // start recording when a seabed/drop run begins (if armed in the dialog)
		void maybeCaptureFrame();
		std::unique_ptr<VideoRecorder> recorder_;
		VideoSettingsDialog* video_dialog_ = nullptr; // non-modal; owned by `this` (Qt parent)
		QAction* record_action_ = nullptr;
		std::vector<float> last_capture_zb_;
		long long last_capture_step_ = -1000000000LL; // sim step of the last captured frame (first always captures)
		long long frame_step_interval_ = 20;          // capture at most one video frame per this many sim steps
		                                              // (bounds the clip length; lower ⇒ more/smoother frames)
		double frame_bed_eps_ = 0.0002;               // …and only if the bed moved > this [m] since — skips dead
		                                              // periods (spin-up / near-equilibrium) so a static bed writes nothing
		int rec_crf_ = 17;                            // H.264 quality passed to VideoRecorder::start (dialog-set)
		int rec_fps_ = 30;                            // playback fps passed to VideoRecorder::start (dialog-set)
		QString rec_preset_ = "fast";                 // x264 preset passed to VideoRecorder::start (dialog-set)
		bool auto_record_ = false;                    // auto-start recording when a seabed/drop run begins

		SliceViewer* viewer_ = nullptr;
		SimWorker* worker_ = nullptr;
		QThread* worker_thread_ = nullptr;
		QLabel* status_ = nullptr;
		QLabel* fps_label_ = nullptr;
		SimRecipe recipe_;
		scour::core::TriMesh model_mesh_; // CPU copy of the loaded model (for voxelization)
		bool model_injected_ = false;
		bool seabed_ = false; // erodible-seabed scenario active (disables the obstacle toggle)
		bool worker_down_ = false;
		DiversionReport diversion_;
		long long last_steps_ = 0;
		double last_sim_time_ = 0.0;
	};
}
