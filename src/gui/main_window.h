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

namespace windcfd::core { class ChannelFluidCore; }

namespace windcfd::gui
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
		MainWindow(std::unique_ptr<windcfd::core::ChannelFluidCore> core, const SimRecipe& recipe,
			QWidget* parent = nullptr);
		~MainWindow() override;

		SliceViewer* viewer() { return viewer_; }

		// Load a STEP model (OpenCascade → metre-scale TriMesh), display it, and auto-voxelize it into
		// the flow as the solid obstacle, REPLACING the config obstacle (the loaded model IS the
		// obstacle; there is no separate toggle). `noslip` selects SOLID_NOSLIP; the shipping default is
		// SOLID_FREESLIP (RESEARCH §3). Used by the File menu and the --load-step CLI flag. Logs a
		// summary; returns false on load failure. Safe to call before the window is shown.
		bool loadStepFile(const QString& path, bool noslip = false);

		// Load a 3D-printing CENTERLINE STEP (vertical wall ribbons) as a CENTERLINE — distinct from
		// loadStepFile's "STEP as mesh obstacle". The mesh is STORED (centerline_mesh_) so Build can
		// re-run without reloading; it is NOT voxelized as a mesh. On load, the domain is sized to the
		// footprint and a solid building (thickened walls + overhanging flat roof) is built + injected
		// once (buildBuilding). `noslip` selects SOLID_NOSLIP for the injected solid (default free-slip,
		// matching loadStepFile). Used by the File menu and the --load-centerline CLI flag. Logs the
		// tri count + bbox; returns false on load failure.
		bool loadCenterlineFile(const QString& path, bool noslip = false);

		// Stop the worker (if running) and, if a voxel obstacle was injected, sample the flow-
		// diversion metrics. Idempotent. Call before reading diversion() from main().
		void finalize();
		const DiversionReport& diversion() const { return diversion_; }

		// Latest wind loads sampled at shutdown (valid only if a building was present). For the headless
		// exit-log assertion in main().
		const windcfd::core::WindLoads& lastLoads() const { return last_loads_; }
		bool lastLoadsValid() const { return last_loads_valid_; }

		// Latest worker totals (for the scripted-gate summary in main.cpp).
		long long lastSteps() const { return last_steps_; }
		double lastSimTime() const { return last_sim_time_; }

	public:
		// Programmatically apply a new grid (domain + voxel/cell size h) and rebuild+reset the sim,
		// exactly as the dock's Apply button. Any non-positive argument keeps the current sim's value.
		// Used by the `--apply-grid` / `--apply-h` CLI smoke hook so the reset-at-new-grid feature is
		// testable headlessly (mirrors --load-step / --voxelize).
		void applyGridValues(double Lx, double Ly, double Lz, double h);

		// Programmatically press "Start Simulation": release the master run gate so the worker begins
		// advancing the flow (the sim is held/paused until this). Idempotent. Used by the headless/gate
		// path (--autoclose-ms) to auto-start after setup so loads develop, and by the dock's Start button.
		void startSimulation();

		// Programmatically set the inlet current speed U [m/s], exactly as the "Input speed" spin box:
		// applies live to the running flow (no reset). Used by the `--set-u` CLI smoke hook so the
		// different-currents feature is testable headlessly (mirrors --apply-h).
		void setInputSpeed(double U);

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

		// Start recording an MP4 of the run to `path`: one frame is captured every N simulation steps and
		// streamed to ffmpeg (H.264). Used by the File menu ("Record Video…") and the `--record <path>`
		// CLI smoke hook. Returns false if ffmpeg is unavailable. Always finalized on close.
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

		// (Re)generate the solid building from the stored centerline and inject it as the obstacle. Sizes
		// the domain around the sectioned footprint (with wind clearance) and rebuilds the core at that
		// grid through the SAME resize path as Apply (applyGrid / GridOverride), then voxelizes the
		// thickened walls + overhanging roof (core/geometry/building) and hands the mask to the worker
		// exactly as loadStepFile does. Driven by the "Build" button and once on centerline load. No-op
		// (with a status message) if no centerline is loaded.
		void buildBuilding();

		// Refresh the Building-group live wind-load coefficient readout from the worker's latest WindLoads
		// (Cd/Cl/Cs + Cp range). Called on the repaint tick; no-op until a building's loads are published.
		void updateWindLoadReadout();

		void shutdownWorker();
		// `steps0`/`t0` prime the worker's counters (a scene restore resumes from the saved step).
		void spawnWorker(std::unique_ptr<windcfd::core::ChannelFluidCore> core,
			long long steps0 = 0, double t0 = 0.0);

		// --- Scene save / restore internals ------------------------------------------------------
		void saveSceneAs();     // File menu: pick a .scn path, then request a save
		void loadSceneDialog(); // File menu: pick a .scn (offering sibling restart points), then restore
		void restoreScene(SceneFile& sf); // tear down + rebuild core + inject the saved state
		// Worker → main-thread hand-off: write the gathered state to disk. tag<0 = manual save (→ Recent
		// Files), tag>=0 = an auto-save at that step (→ <base>.<step>.scn, not Recent).
		void onCheckpointReady(windcfd::gui::CheckpointStatePtr state, qint64 tag);
		QString scene_base_path_;          // stem (no extension) of the current scene, for the auto-save series
		QString pending_manual_save_path_; // destination of an in-flight manual "Save Scene As"
		QString pending_recent_keep_;      // a just-saved scene optimistically in Recent Files before its
		                                   // async write lands — rebuildRecentMenu keeps it from being pruned
		int autosave_interval_ = 0;        // auto-save cadence in steps (0 = off); re-applied on each spawnWorker
		windcfd::core::TriMesh scene_mesh_;  // display mesh persisted into a saved scene (obstacle model)
		windcfd::core::ModelPlacement scene_place_{}; // its display placement
		// Stop + join + delete the worker (and its thread) so a fresh core can replace it. Shared by the
		// grid Apply and scene restore (both fully rebuild the sim). No-op if none.
		void teardownWorkerForReload();

		// Build the left-side settings dock (grouped controls: Simulation / Visualization). Called
		// once from the constructor. The idea mirrors cobod-slicer's left panel (grouped, scrollable),
		// not its styling.
		void buildControlDock();

		// --- Editable domain + resolution (Simulation group) -------------------------------------
		// Apply rebuilds+resets the whole sim (fresh t=0) at the domain (Lx,Ly,Lz) and voxel/cell size
		// h in the spin boxes, preserving every physics parameter — ONLY the grid (nx,ny,nz) changes.
		// syncGridControls repopulates the boxes from the live sim; updateGridReadout refreshes the live
		// "→ nx×ny×nz = N cells" label (and guards against an over-fine grid).
		void applyGrid();
		void syncGridControls();
		void updateGridReadout();

		// --- Model-placement gizmo UI ------------------------------------------------------------
		// Enable the gizmo dock group + refresh its readout for the current model (any gizmo-editable
		// model). Called on model load/close, grid Apply and scene restore.
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
		QComboBox* fast_sim_box_ = nullptr;      // "Fast sim": graphics-update throttle (Off / 0.5 / 1 / 3 s / paused)
		// Tidal reversal (live): enable + peak current / plateau-hold / slack-ramp durations. `syncTidal`
		// pushes the current control values to the worker; the enable box also greys the manual "Input
		// speed" (the tide drives the inlet while active).
		QCheckBox* tidal_chk_ = nullptr;
		QDoubleSpinBox* tidal_umax_spin_ = nullptr;
		QDoubleSpinBox* tidal_plateau_spin_ = nullptr;
		QDoubleSpinBox* tidal_ramp_spin_ = nullptr;
		void syncTidal();

		QLabel* grid_readout_ = nullptr;
		QPushButton* apply_btn_ = nullptr;
		QPushButton* play_btn_ = nullptr; // so a rebuild can honour the current play/pause state
		QPushButton* start_btn_ = nullptr; // "Start Simulation" — releases the master run gate (holds until pressed)
		QPushButton* step_btn_ = nullptr;  // single-step (enabled only once started)
		bool sim_started_ = false;         // master run gate mirror; re-applied to every (re)spawned worker
		static constexpr long long kMaxCells = 40000000LL; // Apply guard: refuse an OOM-risking grid

		// Recent-files (feature 1): persisted via QSettings("COBOD","WindCFD"). The submenu is
		// rebuilt on demand (pruning files that no longer exist), MRU-first, capped at kMaxRecent.
		void addRecentFile(const QString& path);
		void rebuildRecentMenu();
		static constexpr int kMaxRecent = 8;
		QMenu* recent_menu_ = nullptr;
		QCheckBox* slice_chk_ = nullptr; // "Show slice"

		QTimer* repaint_timer_ = nullptr; // main-thread ~60 Hz repaint; interval widens under "Fast sim"
		double display_throttle_s_ = 0.0; // "Fast sim" graphics-update interval [s] (0 = live); re-applied on spawnWorker

		// Video recording (video_recorder.*): capture the viewer once every N sim steps and stream the
		// frames to ffmpeg. The File-menu action opens a non-modal settings popup (openRecordDialog);
		// maybeCaptureFrame runs on the repaint tick and grabs a frame when the step cadence elapses;
		// updateRecordDialogStatus pushes the live frame count into the popup.
		void openRecordDialog();       // File menu: show the settings dialog (Start/Stop + params + live status)
		void updateRecordDialogStatus(); // refresh the dialog's "● REC — N frames" from the recorder
		void maybeCaptureFrame();
		std::unique_ptr<VideoRecorder> recorder_;
		VideoSettingsDialog* video_dialog_ = nullptr; // non-modal; owned by `this` (Qt parent)
		QAction* record_action_ = nullptr;
		long long last_capture_step_ = -1000000000LL; // sim step of the last captured frame (first always captures)
		long long frame_step_interval_ = 20;          // capture one video frame per this many sim steps
		                                              // (bounds the clip length; lower ⇒ more/smoother frames)
		int rec_crf_ = 17;                            // H.264 quality passed to VideoRecorder::start (dialog-set)
		int rec_fps_ = 30;                            // playback fps passed to VideoRecorder::start (dialog-set)
		QString rec_preset_ = "fast";                 // x264 preset passed to VideoRecorder::start (dialog-set)

		SliceViewer* viewer_ = nullptr;
		SimWorker* worker_ = nullptr;
		QThread* worker_thread_ = nullptr;
		QLabel* status_ = nullptr;
		QLabel* fps_label_ = nullptr;
		SimRecipe recipe_;
		windcfd::core::TriMesh model_mesh_; // CPU copy of the loaded model (for voxelization)

		// --- Building (centerline STEP → thickened walls + flat roof solid) ----------------------
		// Stored centerline surface + the dock's wall/roof params, so Build can re-voxelize + re-inject
		// the building without reloading the STEP. Distinct from model_mesh_ (the STEP-as-mesh obstacle).
		windcfd::core::TriMesh centerline_mesh_; // stored centerline (empty until a centerline is loaded)
		bool centerline_noslip_ = false;         // surface mode for the injected building solid
		QDoubleSpinBox* wall_thick_spin_ = nullptr;
		QDoubleSpinBox* wall_height_spin_ = nullptr;
		QDoubleSpinBox* corner_radius_spin_ = nullptr;
		QDoubleSpinBox* roof_overhang_spin_ = nullptr;
		QDoubleSpinBox* roof_thick_spin_ = nullptr;
		QPushButton* build_btn_ = nullptr;
		QLabel* load_readout_ = nullptr; // live wind-load coefficient readout (Cd/Cl/Cs + Cp range)

		bool model_injected_ = false;
		bool worker_down_ = false;
		DiversionReport diversion_;
		windcfd::core::WindLoads last_loads_{}; // wind loads captured at shutdown (headless exit-log)
		bool last_loads_valid_ = false;
		long long last_steps_ = 0;
		double last_sim_time_ = 0.0;
	};
}
