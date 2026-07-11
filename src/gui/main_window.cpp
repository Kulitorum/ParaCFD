// main_window.cpp — see main_window.h.
#include "gui/main_window.h"

#include "core/geometry/model_placement.h"
#include "core/geometry/step_import.h"
#include "core/geometry/voxelize.h"
#include "gui/slice_viewer.h"
#include "gui/sim_worker.h"
#include "gui/video_recorder.h"
#include "gui/video_settings_dialog.h"

#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QCheckBox>
#include <QCloseEvent>
#include <QComboBox>
#include <QDir>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QFrame>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QKeySequence>
#include <QLabel>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QMetaType>
#include <QPushButton>
#include <QQuaternion>
#include <QVector3D>
#include <QScrollArea>
#include <QSettings>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QStatusBar>
#include <QThread>
#include <QTimer>
#include <QVBoxLayout>
#include <QWidget>

#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

namespace scour::gui
{
	MainWindow::MainWindow(std::unique_ptr<scour::core::ChannelFluidCore> core, const SimRecipe& recipe,
		SeabedScenario scen, QWidget* parent)
		: QMainWindow(parent), recipe_(recipe)
	{
		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("ScourProtection — G1 slice viewer [%1: %2x%3x%4, h=%5 m]")
			.arg(QString::fromStdString(info.name))
			.arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h));

		viewer_ = new SliceViewer(this);
		viewer_->setInfo(info);
		setCentralWidget(viewer_);

		// --- File menu: load / clear a STEP model --------------------------------
		QMenu* fileMenu = menuBar()->addMenu("&File");
		QAction* loadScenario = fileMenu->addAction("Load seabed scenario…");
		connect(loadScenario, &QAction::triggered, this, [this] {
			QString fn = QFileDialog::getOpenFileName(this, "Load seabed scenario", QString(),
				"Scenario JSON (*.json);;All files (*)");
			if (!fn.isEmpty()) startSeabedScenario(fn);
		});
		fileMenu->addSeparator();
		QAction* openStep = fileMenu->addAction("Open STEP…");
		openStep->setShortcut(QKeySequence::Open);
		connect(openStep, &QAction::triggered, this, [this] {
			QString fn = QFileDialog::getOpenFileName(this, "Open STEP model", QString(),
				"STEP files (*.step *.stp);;All files (*)");
			if (!fn.isEmpty()) loadStepFile(fn);
		});
		QAction* closeStep = fileMenu->addAction("Close model");
		connect(closeStep, &QAction::triggered, this, [this] {
			setModelAsObstacle(false); // remove the model obstacle, restoring the config obstacle (if any)
			model_mesh_ = scour::core::TriMesh{};
			if (viewer_) viewer_->clearMesh();
			updateGizmoUi(); // no model ⇒ disable the placement gizmo
			statusBar()->showMessage("model closed", 3000);
		});

		// --- Scene save / restore (a full self-contained .scn: setup + model + all field data) ----
		fileMenu->addSeparator();
		QAction* saveScene = fileMenu->addAction("Save Scene As…");
		saveScene->setShortcut(QKeySequence::Save);
		connect(saveScene, &QAction::triggered, this, [this] { saveSceneAs(); });
		QAction* loadScene = fileMenu->addAction("Load Scene…");
		connect(loadScene, &QAction::triggered, this, [this] { loadSceneDialog(); });

		fileMenu->addSeparator();
		// Record an MP4 of the run: opens a settings popup (path / cadence / quality / fps / auto-record) with
		// Start/Stop + a live "● REC — N frames" status. A frame is captured on each sand-bed change (ffmpeg).
		record_action_ = fileMenu->addAction("Record Video…");
		connect(record_action_, &QAction::triggered, this, [this] { openRecordDialog(); });

		// Auto-save checkpoints: <scene>.<step>.scn every N steps (a restart point; not added to Recent).
		QMenu* autoMenu = fileMenu->addMenu("Auto-save checkpoints");
		QActionGroup* autoGroup = new QActionGroup(this);
		autoGroup->setExclusive(true);
		const int intervals[] = { 0, 500, 1000, 2000, 5000 };
		for (int iv : intervals)
		{
			QAction* a = autoMenu->addAction(iv == 0 ? QString("Off") : QString("Every %1 steps").arg(iv));
			a->setCheckable(true);
			a->setChecked(iv == autosave_interval_);
			autoGroup->addAction(a);
			connect(a, &QAction::triggered, this, [this, iv] { setAutosave(iv); });
		}

		// --- Recent files — STEP models AND saved scenes (persisted across sessions via QSettings) --
		recent_menu_ = fileMenu->addMenu("Recent Files");
		rebuildRecentMenu();

		// A loaded STEP model IS the obstacle: loadStepFile auto-voxelizes + injects it (replacing the
		// config obstacle), so there is no separate "Model as obstacle" toggle / Model menu any more.

		// --- Left control dock: grouped settings (Simulation / Visualization) --------------------
		// Replaces the old crowded top toolbar. The idea (not the styling) follows cobod-slicer's
		// left settings panel: a scrollable column of QGroupBox sections, one per settings area, so
		// simulation parameters (resolution, domain, dt, …) can slot into the Simulation group later.
		buildControlDock();

		// --- Status bar -----------------------------------------------------------
		status_ = new QLabel("step 0   t = 0.000 s");
		fps_label_ = new QLabel("-- fps");
		statusBar()->addWidget(status_, 1);
		statusBar()->addPermanentWidget(fps_label_);

		connect(viewer_, &SliceViewer::fpsUpdated, this,
			[this](double fps) { fps_label_->setText(QString("%1 fps").arg(fps, 0, 'f', 1)); });

		resize(1100, 720);

		recorder_ = std::make_unique<VideoRecorder>();

		// Continuous repaint at ~60 Hz on the main thread (decoupled from the sim step
		// rate; the viewer always shows the latest device state). The same tick also captures a video
		// frame when recording and the sand bed has changed (maybeCaptureFrame early-returns otherwise).
		repaint_timer_ = new QTimer(this);
		connect(repaint_timer_, &QTimer::timeout, this, [this] { viewer_->update(); maybeCaptureFrame(); updateRecordDialogStatus(); });
		repaint_timer_->start(16); // widened by setDisplayThrottle when "Fast sim" is engaged

		// Create + start the worker (and attach the seabed engine + bed viz if the scenario is active).
		spawnWorker(std::move(core), std::move(scen));
	}

	void MainWindow::buildControlDock()
	{
		QDockWidget* dock = new QDockWidget("Settings", this);
		dock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
		dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

		QWidget* panel = new QWidget;
		QVBoxLayout* col = new QVBoxLayout(panel);
		col->setContentsMargins(8, 8, 8, 8);
		col->setSpacing(10);

		// --- Simulation group ------------------------------------------------------------------
		QGroupBox* simGroup = new QGroupBox("Simulation (live — no reset)");
		QVBoxLayout* simCol = new QVBoxLayout(simGroup);
		QHBoxLayout* runRow = new QHBoxLayout;
		QPushButton* playBtn = new QPushButton("Pause");
		playBtn->setCheckable(true);
		playBtn->setChecked(true);
		playBtn->setToolTip("Play/pause the simulation stepping.");
		connect(playBtn, &QPushButton::toggled, this, [this, playBtn](bool on) {
			if (worker_) worker_->setPlaying(on);
			playBtn->setText(on ? "Pause" : "Play");
		});
		play_btn_ = playBtn; // a rebuild (scenario load / grid Apply) honours this play/pause state
		QPushButton* stepBtn = new QPushButton("Step");
		stepBtn->setToolTip("Advance a single step (while paused).");
		connect(stepBtn, &QPushButton::clicked, this, [this] { if (worker_) worker_->stepOnce(); });
		runRow->addWidget(playBtn);
		runRow->addWidget(stepBtn);
		simCol->addLayout(runRow);

		// Inlet current speed U — applies LIVE to the running flow (no reset): the wake/scour evolves
		// toward the new current so you can test different currents interactively. It is ALSO honoured
		// by Apply below (a grid rebuild keeps the chosen current). Distinct from the domain/h block,
		// which resets the sim at t=0.
		QFormLayout* speedForm = new QFormLayout;
		speedForm->setLabelAlignment(Qt::AlignLeft);
		u_spin_ = new QDoubleSpinBox;
		u_spin_->setRange(0.0, 10.0);
		u_spin_->setDecimals(3);
		u_spin_->setSingleStep(0.05);
		u_spin_->setSuffix(" m/s");
		u_spin_->setToolTip("Inlet current speed U. Applies live to the running flow (no reset) — the wake/scour evolves toward it; also used when you Apply a new grid.");
		connect(u_spin_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this](double v) {
			if (worker_) worker_->setInletSpeed(v);
			if (viewer_) viewer_->setReferenceU(v);
		});
		speedForm->addRow("Input speed U", u_spin_);

		// Inlet velocity profile: uniform (top-hat) vs boundary-layer (log-law). The BL profile tapers the
		// near-bed inflow to ~0 at the bed top, so the sediment's leading edge isn't over-eroded by the
		// top-hat's spurious full-U bed shear. Applies LIVE (the flow adjusts over a flow-through).
		inlet_profile_box_ = new QComboBox;
		inlet_profile_box_->addItems({ "Uniform (top-hat)", "Boundary layer (log-law)" });
		inlet_profile_box_->setToolTip("Inlet velocity profile. Boundary layer (log-law) → near-bed u≈0 at the bed: the physical fix for top-hat leading-edge scour. Applies live; sediment runs default to it.");
		connect(inlet_profile_box_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int idx) {
			if (!worker_) return;
			const bool loglaw = (idx == 1);
			const double z0 = (have_sp_ && last_sp_.d50 > 0.0) ? last_sp_.d50 / 12.0 : (0.2e-3 / 12.0);
			const double datum = have_sp_ ? last_sp_.sand_depth : 0.0;
			worker_->setInletProfile(loglaw, z0, datum);
		});
		speedForm->addRow("Inlet profile", inlet_profile_box_);

		// Suspended-sediment boundary on the streamwise (x) faces. Open sea: the inlet carries the Rouse
		// equilibrium load so sheltered pockets get sand to deposit and the outlet lets it leave; Recycle:
		// flux-matched recirculating flume (mass-conserving); Closed: the legacy zero-flux box. Live; the
		// switch only affects an active seabed run (no bed ⇒ no-op).
		sed_bc_box_ = new QComboBox;
		sed_bc_box_->addItems({ "Open sea (equilibrium)", "Recycle (flume)", "Closed (no exchange)" });
		sed_bc_box_->setToolTip("Suspended-sediment inlet/outlet. Open sea: incoming water carries the equilibrium suspended load (fills sheltered pockets), outlet lets it leave; Recycle: outlet load fed back, mass-conserving; Closed: no exchange (legacy). Applies live to an active seabed run.");
		connect(sed_bc_box_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int idx) {
			if (!worker_) return;
			int mode = (idx == 0) ? scour::core::SED_BC_OPEN : (idx == 1) ? scour::core::SED_BC_RECYCLE : scour::core::SED_BC_CLOSED;
			worker_->setSedBoundary(mode);
		});
		speedForm->addRow("Sediment BC", sed_bc_box_);
		simCol->addLayout(speedForm);

		// "Add sand" (LIVE — keeps the developed flow, NO t=0 reset): fill the lower N metres with sand
		// around any rigid solid and continue as an erodible-seabed sediment run, so even a plain fluid
		// viewer with a voxelized obstacle becomes a live scour sim. It stays in this live box precisely
		// because it does NOT reset — unlike the grid Apply, which is its own box below.
		auto* sandLine = new QFrame; sandLine->setFrameShape(QFrame::HLine); sandLine->setEnabled(false);
		simCol->addWidget(sandLine);
		QFormLayout* sandForm = new QFormLayout;
		sandForm->setLabelAlignment(Qt::AlignLeft);
		sand_spin_ = new QDoubleSpinBox;
		sand_spin_->setRange(0.0, 50.0);
		sand_spin_->setDecimals(3);
		sand_spin_->setSingleStep(0.10);
		sand_spin_->setValue(1.0);
		sand_spin_->setSuffix(" m");
		sand_spin_->setToolTip("Depth of sand to fill from the bottom up. Rigid structure/obstacle voxels stay solid; sand fills the rest.");
		sandForm->addRow("Sand depth", sand_spin_);
		simCol->addLayout(sandForm);
		QPushButton* addSandBtn = new QPushButton("Add sand (keep flow)");
		addSandBtn->setToolTip("Fill the lower N metres with sand around any solid and continue as an erodible-seabed sediment run. LIVE: the developed flow is PRESERVED — no t=0 reset (sim time/step continue). PRE-CALIBRATION / qualitative.");
		connect(addSandBtn, &QPushButton::clicked, this, [this] { addSand(); });
		simCol->addWidget(addSandBtn);

		// MORFAC (morphological acceleration): multiplies morphological time so the bed evolves M× faster
		// than the flow clock — the "catch up faster" knob for watching scour develop. Live: retunes the
		// running seabed engine (a no-op without a bed; a value chosen before "Add sand" is applied when the
		// engine attaches). ⚠ Physically valid to ~10 steady / ~5 reversing (RESEARCH §7); above that the
		// flow and bed decouple and the per-step |Δz_b| limiter clips, so higher isn't proportionally faster.
		QFormLayout* morfacForm = new QFormLayout;
		morfacForm->setLabelAlignment(Qt::AlignLeft);
		morfac_spin_ = new QDoubleSpinBox;
		morfac_spin_->setRange(1.0, 50.0);
		morfac_spin_->setDecimals(1);
		morfac_spin_->setSingleStep(1.0);
		morfac_spin_->setValue(5.0);
		morfac_spin_->setSuffix(QString::fromUtf8(" ×")); // "×"
		morfac_spin_->setToolTip("MORFAC — morphological acceleration. The bed evolves this many times faster than the flow clock, so scour 'catches up' sooner. Live: retunes the running seabed engine (no effect without a bed; a value set before 'Add sand' applies when the engine attaches). Physically valid to ~10 (steady) / ~5 (reversing tide); above that the flow and bed decouple and the per-step bed-change limiter clips, so higher isn't proportionally faster.");
		connect(morfac_spin_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
			[this](double v) { if (worker_) worker_->setMorfac(v); });
		morfacForm->addRow("MORFAC (bed speed-up)", morfac_spin_);
		simCol->addLayout(morfacForm);

		// Tidal reversal (LIVE, no reset): drive a reversing current U_d(t) that face-swaps the inlet/
		// outlet through a cosine slack ramp (RESEARCH §8). Watch a scour hole reorganize as the tide
		// flips. It takes over "Input speed U" while enabled (greyed below). ⚠ MORFAC ≤ 5 in reversing
		// flow (RESEARCH §7) — a run-setup choice; keep the scenario's morfac modest.
		auto* tideLine = new QFrame; tideLine->setFrameShape(QFrame::HLine); tideLine->setEnabled(false);
		simCol->addWidget(tideLine);
		tidal_chk_ = new QCheckBox("Tidal reversal (live)");
		tidal_chk_->setToolTip("Drive a reversing tidal current: +U_max plateau → cosine slack ramp → −U_max plateau → … , face-swapping the inlet/outlet at slack (RESEARCH §8). LIVE — no reset. Overrides 'Input speed U' while on. Use MORFAC ≤ 5 for reversing flow.");
		simCol->addWidget(tidal_chk_);
		QFormLayout* tideForm = new QFormLayout;
		tideForm->setLabelAlignment(Qt::AlignLeft);
		tidal_umax_spin_ = new QDoubleSpinBox; tidal_umax_spin_->setRange(0.0, 10.0); tidal_umax_spin_->setDecimals(3); tidal_umax_spin_->setSingleStep(0.05); tidal_umax_spin_->setSuffix(" m/s"); tidal_umax_spin_->setValue(0.8);
		tidal_umax_spin_->setToolTip("Peak (spring-tide) current U_max. The tide runs +U_max ↔ −U_max.");
		tidal_plateau_spin_ = new QDoubleSpinBox; tidal_plateau_spin_->setRange(1.0, 100000.0); tidal_plateau_spin_->setDecimals(1); tidal_plateau_spin_->setSingleStep(5.0); tidal_plateau_spin_->setSuffix(" s"); tidal_plateau_spin_->setValue(30.0);
		tidal_plateau_spin_->setToolTip("Plateau: how long the current holds at ±U_max before turning.");
		tidal_ramp_spin_ = new QDoubleSpinBox; tidal_ramp_spin_->setRange(1.0, 100000.0); tidal_ramp_spin_->setDecimals(1); tidal_ramp_spin_->setSingleStep(5.0); tidal_ramp_spin_->setSuffix(" s"); tidal_ramp_spin_->setValue(20.0);
		tidal_ramp_spin_->setToolTip("Slack ramp: the cosine turn duration between directions (the 'slack water'). RESEARCH §8 uses ≥ 3–6 flow-throughs.");
		tideForm->addRow("Peak U_max", tidal_umax_spin_);
		tideForm->addRow("Plateau", tidal_plateau_spin_);
		tideForm->addRow("Slack ramp", tidal_ramp_spin_);
		simCol->addLayout(tideForm);
		auto pushTide = [this] { syncTidal(); };
		connect(tidal_chk_, &QCheckBox::toggled, this, [this](bool) { syncTidal(); });
		connect(tidal_umax_spin_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [pushTide](double) { pushTide(); });
		connect(tidal_plateau_spin_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [pushTide](double) { pushTide(); });
		connect(tidal_ramp_spin_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [pushTide](double) { pushTide(); });

		col->addWidget(simGroup);

		// --- Domain & resolution group — the ONLY control that REBUILDS + RESETS the sim to t=0 --------
		// Editable domain + voxel/cell size. Apply rebuilds at the new grid (fresh t=0), preserving the
		// scenario + all physics — only nx,ny,nz change. Coarser (bigger h) = faster to iterate. h IS the
		// voxel size: one uniform grid = fluid cell = model/bed voxelization spacing. In its own box (vs
		// the live "Simulation" controls) so the reset is unmistakable.
		QGroupBox* gridGroup = new QGroupBox("Domain & resolution (Apply resets)");
		QVBoxLayout* gridCol = new QVBoxLayout(gridGroup);
		QFormLayout* gridForm = new QFormLayout;
		gridForm->setLabelAlignment(Qt::AlignLeft);
		auto makeDomainSpin = [this](const char* tip) {
			QDoubleSpinBox* s = new QDoubleSpinBox;
			s->setRange(0.10, 100.0);
			s->setDecimals(2);
			s->setSingleStep(0.5);
			s->setSuffix(" m");
			s->setToolTip(tip);
			connect(s, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this](double) { updateGridReadout(); });
			return s;
		};
		lx_spin_ = makeDomainSpin("Domain length in x (streamwise) [m].");
		ly_spin_ = makeDomainSpin("Domain length in y (span) [m].");
		lz_spin_ = makeDomainSpin("Domain height in z [m].");
		gridForm->addRow("Domain Lx", lx_spin_);
		gridForm->addRow("Domain Ly", ly_spin_);
		gridForm->addRow("Domain Lz", lz_spin_);

		h_spin_ = new QDoubleSpinBox;
		h_spin_->setRange(0.01, 1.0);
		h_spin_->setDecimals(3);
		h_spin_->setSingleStep(0.01);
		h_spin_->setSuffix(" m");
		h_spin_->setToolTip("Uniform grid spacing: the fluid cell size AND the resolution the model + bed are voxelized at. Smaller = finer & slower.");
		connect(h_spin_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this](double) { updateGridReadout(); });
		gridForm->addRow("Voxel / cell size h", h_spin_);
		gridCol->addLayout(gridForm);

		grid_readout_ = new QLabel;
		grid_readout_->setWordWrap(true);
		grid_readout_->setToolTip("Resolution + total cell count Apply would build (your gauge for speed/memory).");
		gridCol->addWidget(grid_readout_);

		apply_btn_ = new QPushButton("Apply (rebuild + reset)");
		apply_btn_->setToolTip("Rebuild the simulation at the domain + voxel size above and restart it at t=0 (keeps the scenario + all physics).");
		connect(apply_btn_, &QPushButton::clicked, this, [this] { applyGrid(); });
		gridCol->addWidget(apply_btn_);
		col->addWidget(gridGroup);

#ifdef SCOUR_HAVE_JOLT
		// --- Drop objects (preparation) — REBUILDS + RESETS to a fresh seabed run at t=0 ----------
		// Scatter scaled copies of a STEP protection unit above a sand bed, settle them under gravity
		// (Jolt), then voxelize the settled pile as the rigid structure and start the morphodynamic run
		// (largest-first timed release — an armour-layer pour). Like the grid Apply it RESETS to t=0, so it
		// sits in its own box. ⚠ Needs a domain sized for the units (e.g. g1_viewer_full's 10×10×5 m).
		// A grid Apply AFTER a drop PERSISTS the settled pile: applyGrid re-voxelizes the multi-block
		// placements at the new grid (the `have_pile` path) instead of reverting to the config obstacle.
		QGroupBox* dropGroup = new QGroupBox("Drop objects (preparation)");
		QVBoxLayout* dropCol = new QVBoxLayout(dropGroup);

		QHBoxLayout* dropPathRow = new QHBoxLayout;
		QLineEdit* dropPathEdit = new QLineEdit;
		{
			QSettings s("COBOD", "ScourProtection");
			dropPathEdit->setText(s.value("dropStepFile", "Experiments/V000 Code tests/XStone_Decomposed.stp").toString());
		}
		dropPathEdit->setToolTip("STEP protection unit to scatter + settle (default: XStone_Decomposed.stp, an 8-piece convex compound). Remembered across sessions.");
		QPushButton* dropBrowse = new QPushButton("Browse…");
		connect(dropBrowse, &QPushButton::clicked, this, [this, dropPathEdit] {
			const QString fn = QFileDialog::getOpenFileName(this, "Drop unit STEP", dropPathEdit->text(),
				"STEP files (*.step *.stp);;All files (*)");
			if (!fn.isEmpty()) dropPathEdit->setText(fn);
		});
		dropPathRow->addWidget(dropPathEdit, 1);
		dropPathRow->addWidget(dropBrowse);
		dropCol->addLayout(dropPathRow);

		QFormLayout* dropForm = new QFormLayout;
		dropForm->setLabelAlignment(Qt::AlignLeft);
		QSpinBox* dropCount = new QSpinBox; dropCount->setRange(1, 2000); dropCount->setValue(12);
		dropCount->setToolTip("Number of units to scatter + settle onto the bed.");
		dropForm->addRow("Count", dropCount);
		QDoubleSpinBox* dropSmin = new QDoubleSpinBox; dropSmin->setRange(0.05, 10.0); dropSmin->setDecimals(2); dropSmin->setSingleStep(0.1); dropSmin->setValue(0.6);
		dropSmin->setToolTip("Minimum random size scale applied to each unit.");
		dropForm->addRow("Size min", dropSmin);
		QDoubleSpinBox* dropSmax = new QDoubleSpinBox; dropSmax->setRange(0.05, 10.0); dropSmax->setDecimals(2); dropSmax->setSingleStep(0.1); dropSmax->setValue(1.0);
		dropSmax->setToolTip("Maximum random size scale applied to each unit.");
		dropForm->addRow("Size max", dropSmax);
		QDoubleSpinBox* dropSand = new QDoubleSpinBox; dropSand->setRange(0.0, 50.0); dropSand->setDecimals(3); dropSand->setSingleStep(0.1); dropSand->setSuffix(" m"); dropSand->setValue(0.5);
		dropSand->setToolTip("Depth of the sand bed the units drop onto — and the erodible bed for the ensuing scour run.");
		dropForm->addRow("Sand depth", dropSand);
		QSpinBox* dropSeed = new QSpinBox; dropSeed->setRange(0, 1000000000); dropSeed->setValue(1234);
		dropSeed->setToolTip("Random seed for the scatter (positions / orientations / sizes). Same seed ⇒ the same pile.");
		dropForm->addRow("Seed", dropSeed);
		QSpinBox* dropCadence = new QSpinBox; dropCadence->setRange(0, 600); dropCadence->setValue(drop_release_ticks_); dropCadence->setSuffix(" ticks");
		dropCadence->setToolTip("Animation ticks (~60 Hz) between releasing successive units — the staged pour. Larger = slower release; 0 = release all at once.");
		connect(dropCadence, QOverload<int>::of(&QSpinBox::valueChanged), this, [this](int v) { drop_release_ticks_ = v; });
		dropForm->addRow("Release cadence", dropCadence);
		dropCol->addLayout(dropForm);

		QPushButton* dropBtn = new QPushButton("Drop & settle");
		dropBtn->setToolTip("Scatter + settle the units under gravity, then start a fresh erodible-seabed run (RESETS to t=0) with the settled pile as the rigid structure. ⚠ Size the domain for the units (e.g. 10×10×5 m).");
		connect(dropBtn, &QPushButton::clicked, this, [this, dropPathEdit, dropCount, dropSeed, dropSmin, dropSmax, dropSand] {
			const QString path = dropPathEdit->text().trimmed();
			if (path.isEmpty()) { statusBar()->showMessage("drop: choose a STEP unit first", 4000); return; }
			QSettings("COBOD", "ScourProtection").setValue("dropStepFile", path); // remember last used
			dropBlocks(path, dropCount->value(), (unsigned)dropSeed->value(),
				dropSmin->value(), dropSmax->value(), dropSand->value());
		});
		dropCol->addWidget(dropBtn);
		QPushButton* clearDropBtn = new QPushButton("Clear pile / re-drop");
		clearDropBtn->setToolTip("Stop settling + remove the current pile so you can adjust the parameters and Drop & settle again. (Drop & settle also clears any prior pile automatically.)");
		connect(clearDropBtn, &QPushButton::clicked, this, [this] { clearDropPile(); });
		dropCol->addWidget(clearDropBtn);
		col->addWidget(dropGroup);
#endif

		// --- Visualization group ---------------------------------------------------------------
		QGroupBox* vizGroup = new QGroupBox("Visualization");
		QFormLayout* viz = new QFormLayout(vizGroup);
		viz->setLabelAlignment(Qt::AlignLeft);

		QComboBox* fieldBox = new QComboBox;
		fieldBox->addItems({ "Speed |u|", "u (x-vel)", "v (y-vel)", "w (z-vel)", "Pressure", "Concentration" });
		connect(fieldBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
			[this](int i) { if (viewer_) viewer_->setField(static_cast<Field>(i)); });
		viz->addRow("Field", fieldBox);

		QComboBox* axisBox = new QComboBox;
		axisBox->addItems({ "X-normal", "Y-normal", "Z-normal" });
		axisBox->setCurrentIndex(2); // default Z-normal (x-y wake plane)
		connect(axisBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
			[this](int i) { if (viewer_) viewer_->setAxis(static_cast<Axis>(i)); });
		viz->addRow("Slice plane", axisBox);

		QSlider* planeSlider = new QSlider(Qt::Horizontal);
		planeSlider->setRange(0, 100);
		planeSlider->setValue(50);
		connect(planeSlider, &QSlider::valueChanged, this,
			[this](int v) { if (viewer_) viewer_->setPlaneFraction(v / 100.0f); });
		viz->addRow("Plane pos", planeSlider);

		// Auto colour-range: rescale the colormap (slice + legend + arrows) to the live data each
		// update so an accelerating / blowing-up region stays visible instead of saturating to red.
		QCheckBox* autoRangeChk = new QCheckBox("Auto range");
		autoRangeChk->setChecked(viewer_ ? viewer_->autoRange() : true);
		autoRangeChk->setToolTip("Scale the colour range to the live data (slice, legend and arrows follow the min/max).");
		connect(autoRangeChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setAutoRange(on); });
		viz->addRow(autoRangeChk);

		// Fast sim: throttle graphics updates to reclaim time for the solver when you're not watching.
		// The worker's per-step display snapshot (a device copy + hard sync) and the arrow D2H are the
		// only "graphics" work on the sim's critical path; skipping them most steps is the win. The
		// repaint slows to match (nothing new to draw between snapshots). Off (default) = live 60 Hz,
		// so G1 and every gate are unchanged.
		fast_sim_box_ = new QComboBox;
		fast_sim_box_->addItems({ "Off (live)", "0.5 s", "1 s", "3 s", "Only when paused" });
		fast_sim_box_->setToolTip("Throttle on-screen updates to speed up the simulation when you're away. The sim keeps stepping at full rate; only the view refreshes less often. Off = live 60 Hz.");
		connect(fast_sim_box_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int i) {
			static const double kSecs[] = { 0.0, 0.5, 1.0, 3.0, 1e9 };
			setDisplayThrottle(kSecs[(i >= 0 && i < 5) ? i : 0]);
		});
		viz->addRow("Fast sim", fast_sim_box_);

		// Layer visibility toggles (each independent; all default ON).
		auto* line1 = new QFrame; line1->setFrameShape(QFrame::HLine); line1->setEnabled(false);
		viz->addRow(line1);
		QCheckBox* sliceChk = new QCheckBox("Show slice");
		sliceChk->setChecked(true);
		sliceChk->setToolTip("Show/hide the coloured slice plane. When shown, arrows draw black for contrast.");
		connect(sliceChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowSlice(on); });
		viz->addRow(sliceChk);
		slice_chk_ = sliceChk; // so the seabed scenario can default it OFF (the bed is the star there)
		QCheckBox* modelChk = new QCheckBox("Show model");
		modelChk->setChecked(true);
		modelChk->setToolTip("Show/hide the smooth STEP polygon mesh.");
		connect(modelChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowModel(on); });
		viz->addRow(modelChk);
		QCheckBox* solidVoxChk = new QCheckBox("Show solid voxels");
		solidVoxChk->setChecked(true);
		solidVoxChk->setToolTip("Show/hide the RIGID solid-voxel staircase (structure / obstacle), drawn dark-blue.");
		connect(solidVoxChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowSolidVoxels(on); });
		viz->addRow(solidVoxChk);
		QCheckBox* sedVoxChk = new QCheckBox("Show sediment voxels");
		sedVoxChk->setChecked(true);
		sedVoxChk->setToolTip("Show/hide the ERODIBLE sediment-voxel staircase (sand bed), drawn orange.");
		connect(sedVoxChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowSedimentVoxels(on); });
		viz->addRow(sedVoxChk);
		QCheckBox* bedChk = new QCheckBox("Show bed");
		bedChk->setChecked(true);
		bedChk->setToolTip("Show/hide the erodible seabed surface (the height-coloured z_b polygon + its Δz_b legend). Seabed scenario only.");
		connect(bedChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowBed(on); });
		viz->addRow(bedChk);
		QComboBox* bedColBox = new QComboBox;
		bedColBox->addItems({ "Elevation", "Exchange rate" });
		bedColBox->setToolTip("Colour the seabed surface by: Elevation (Δz_b — scour blue / deposit red), or "
			"Exchange rate = deposition − pickup (blue = sand being lifted / eroded, red = sand landing). "
			"Same surface, different colour source. Seabed scenario only.");
		connect(bedColBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
			[this](int i) { if (viewer_) viewer_->setBedColorMode(i); });
		viz->addRow("Bed colour", bedColBox);
		QCheckBox* axesChk = new QCheckBox("Show axes");
		axesChk->setChecked(true);
		axesChk->setToolTip("Show/hide the world-origin XYZ triad (X=red, Y=green, Z=blue) with metre ticks + corner gizmo.");
		connect(axesChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowAxes(on); });
		viz->addRow(axesChk);

		// Flow arrows.
		auto* line2 = new QFrame; line2->setFrameShape(QFrame::HLine); line2->setEnabled(false);
		viz->addRow(line2);
		QCheckBox* arrowsChk = new QCheckBox("Show arrows");
		arrowsChk->setChecked(true);
		arrowsChk->setToolTip("Animated particle-arrows streaming along the flow.");
		connect(arrowsChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowArrows(on); });
		viz->addRow(arrowsChk);

		QComboBox* arrowModeBox = new QComboBox;
		arrowModeBox->addItems({ "3D volume", "2D plane" });
		arrowModeBox->setToolTip("3D = advect through the whole volume; 2D = on the current slice plane.");
		connect(arrowModeBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
			[this](int i) { if (viewer_) viewer_->setArrowMode3D(i == 0); });
		viz->addRow("Arrow mode", arrowModeBox);

		QSlider* densitySlider = new QSlider(Qt::Horizontal);
		densitySlider->setRange(100, 8000);
		densitySlider->setValue(viewer_ ? viewer_->arrowDensity() : 1500);
		densitySlider->setToolTip("Number of flow-arrow particles.");
		connect(densitySlider, &QSlider::valueChanged, this,
			[this](int v) { if (viewer_) viewer_->setArrowDensity(v); });
		viz->addRow("Density", densitySlider);

		// Arrow speed multiplier [0,1] — scales ONLY the visual tracer motion so a fast current (8 m/s)
		// doesn't whip the arrows about (0 = frozen, 1 = default lively). Physics is untouched.
		QSlider* arrowSpeedSlider = new QSlider(Qt::Horizontal);
		arrowSpeedSlider->setRange(0, 100);
		arrowSpeedSlider->setValue(viewer_ ? (int)std::lround(viewer_->arrowSpeedMult() * 100.0f) : 100);
		arrowSpeedSlider->setToolTip("Visual speed of the flow arrows (0 = frozen, 1 = full). Slow it down for fast currents; does not change the simulation.");
		connect(arrowSpeedSlider, &QSlider::valueChanged, this,
			[this](int v) { if (viewer_) viewer_->setArrowSpeedMult(v / 100.0f); });
		viz->addRow("Arrow speed", arrowSpeedSlider);

		// Arrow width (thickness) — scales ONLY the lateral extent of the shaft + head, so lower = thinner /
		// less cluttered while the arrow LENGTH is kept. Slider [5,100] ⇒ width v/50 ∈ [0.1,2.0] (50 = the
		// glyph's own width). Defaults thinner than the glyph to declutter a dense field.
		QSlider* arrowWidthSlider = new QSlider(Qt::Horizontal);
		arrowWidthSlider->setRange(5, 100);
		arrowWidthSlider->setValue(viewer_ ? (int)std::lround(viewer_->arrowWidthMult() * 50.0f) : 25);
		arrowWidthSlider->setToolTip("Thickness of the flow arrows (shaft + head). Lower = thinner / less cluttered; the arrow length is unchanged.");
		connect(arrowWidthSlider, &QSlider::valueChanged, this,
			[this](int v) { if (viewer_) viewer_->setArrowWidthMult(v / 50.0f); });
		viz->addRow("Arrow width", arrowWidthSlider);

		// Flow tracers (inlet-seeded streamlines): long lines seeded on the inlet face that follow the
		// flow all the way across the domain, coloured by speed. Off by default (arrows are the default).
		auto* line3 = new QFrame; line3->setFrameShape(QFrame::HLine); line3->setEnabled(false);
		viz->addRow(line3);
		QCheckBox* tracersChk = new QCheckBox("Show tracers");
		tracersChk->setChecked(viewer_ ? viewer_->showTracers() : false);
		tracersChk->setToolTip("Inlet-seeded streamlines: long lines that start at the inlet and follow the flow across the domain, coloured by speed.");
		connect(tracersChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowTracers(on); });
		viz->addRow(tracersChk);

		QComboBox* tracerModeBox = new QComboBox;
		tracerModeBox->addItems({ "3D volume", "2D plane" });
		tracerModeBox->setToolTip("3D = seed the whole inlet plane, integrate the full 3D flow; 2D = seed the inlet edge of the current slice plane and stay on it.");
		connect(tracerModeBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
			[this](int i) { if (viewer_) viewer_->setTracerMode3D(i == 0); });
		viz->addRow("Tracer mode", tracerModeBox);

		QSlider* tracerGridSlider = new QSlider(Qt::Horizontal);
		tracerGridSlider->setRange(2, 120);
		tracerGridSlider->setValue(viewer_ ? viewer_->tracerGridDensity() : 12);
		tracerGridSlider->setToolTip("Number of inlet streamlines (seed count along the larger inlet dimension). Higher = more lines.");
		connect(tracerGridSlider, &QSlider::valueChanged, this,
			[this](int v) { if (viewer_) viewer_->setTracerGridDensity(v); });
		viz->addRow("Tracer count", tracerGridSlider);

		QSlider* tracerLenSlider = new QSlider(Qt::Horizontal);
		tracerLenSlider->setRange(50, 2000);
		tracerLenSlider->setValue(viewer_ ? viewer_->tracerTrail() : 600);
		tracerLenSlider->setToolTip("Maximum streamline length (integration steps). Longer = the line reaches further and winds up more inside eddies.");
		connect(tracerLenSlider, &QSlider::valueChanged, this,
			[this](int v) { if (viewer_) viewer_->setTracerTrail(v); });
		viz->addRow("Tracer length", tracerLenSlider);

		// Tracer width in PIXELS (screen-space ribbon). Slider [10,60] ⇒ width v/10 ∈ [1,6] px.
		QSlider* tracerWidthSlider = new QSlider(Qt::Horizontal);
		tracerWidthSlider->setRange(10, 60);
		tracerWidthSlider->setValue(viewer_ ? (int)std::lround(viewer_->tracerWidth() * 10.0f) : 20);
		tracerWidthSlider->setToolTip("Tracer line thickness in pixels (constant on screen regardless of zoom).");
		connect(tracerWidthSlider, &QSlider::valueChanged, this,
			[this](int v) { if (viewer_) viewer_->setTracerWidth(v / 10.0f); });
		viz->addRow("Tracer width", tracerWidthSlider);

		// "Boring" filter — hide short/straight tracers. Slider [990,1100] ⇒ boring v/1000 ∈ [0.99,1.1]
		// in fine 0.001 steps (QSlider is integer-only, so a coarse ×100 scale would step visibly).
		// 0.99 (min) shows all; higher only draws streamlines whose path length exceeds Lx·(this), so
		// straight crossers (≈ Lx) drop out and only the meandering (eddy/wake) lines remain.
		QSlider* tracerBoringSlider = new QSlider(Qt::Horizontal);
		tracerBoringSlider->setRange(990, 1100);
		tracerBoringSlider->setPageStep(10); // arrow = 0.001, page = 0.01
		tracerBoringSlider->setValue(viewer_ ? (int)std::lround(viewer_->tracerBoring() * 1000.0f) : 990);
		tracerBoringSlider->setToolTip("Hide 'boring' (short/straight) tracers. 0.99 (far left) shows all; higher only keeps streamlines whose path length exceeds this multiple of the domain length, so straight crossers drop out and only meandering lines (eddies/wake) remain.");
		connect(tracerBoringSlider, &QSlider::valueChanged, this,
			[this](int v) { if (viewer_) viewer_->setTracerBoring(v / 1000.0f); });
		// While the handle is held, bypass the 1 s anti-flicker retention so the filter tracks the drag
		// live; restore retention on release.
		connect(tracerBoringSlider, &QSlider::sliderPressed, this,
			[this] { if (viewer_) viewer_->setTracerBoringInstant(true); });
		connect(tracerBoringSlider, &QSlider::sliderReleased, this,
			[this] { if (viewer_) viewer_->setTracerBoringInstant(false); });
		viz->addRow("Hide boring", tracerBoringSlider);

		col->addWidget(vizGroup);

		// --- Clip plane group (see inside hollow structures) -----------------------------------
		// A movable plane that hides SOLIDS (STEP model + voxel solids + seabed) on the camera side so
		// the interior of a hollow structure is exposed; the flow slice + arrows are never clipped.
		QGroupBox* clipGroup = new QGroupBox("Clip plane (see inside)");
		QFormLayout* clip = new QFormLayout(clipGroup);
		clip->setLabelAlignment(Qt::AlignLeft);

		QCheckBox* clipChk = new QCheckBox("Enable clip plane");
		clipChk->setToolTip("Hide solids (model + voxels + seabed) between the camera and a movable plane, so you can see inside hollow structures. The flow slice and arrows stay visible.");
		connect(clipChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setClipEnabled(on); });
		clip->addRow(clipChk);

		QComboBox* clipModeBox = new QComboBox;
		clipModeBox->addItems({ "X-normal", "Y-normal", "Z-normal", "Face camera" });
		clipModeBox->setCurrentIndex(viewer_ ? viewer_->clipMode() : 2);
		clipModeBox->setToolTip("Plane orientation. X/Y/Z = axis-aligned (the slider shifts it along that axis). Face camera = perpendicular to the view (the slider pushes it into the scene along your line of sight).");
		connect(clipModeBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
			[this](int i) { if (viewer_) viewer_->setClipMode(i); });
		clip->addRow("Orientation", clipModeBox);

		QSlider* clipPosSlider = new QSlider(Qt::Horizontal);
		clipPosSlider->setRange(0, 100);
		clipPosSlider->setValue(50);
		clipPosSlider->setToolTip("Move the clip plane. Axis modes: position along the axis. Face camera: depth into the scene.");
		connect(clipPosSlider, &QSlider::valueChanged, this,
			[this](int v) { if (viewer_) viewer_->setClipFraction(v / 100.0f); });
		clip->addRow("Position", clipPosSlider);

		QCheckBox* clipFlipChk = new QCheckBox("Flip side");
		clipFlipChk->setToolTip("Swap which side of the plane is hidden.");
		connect(clipFlipChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setClipFlip(on); });
		clip->addRow(clipFlipChk);

		col->addWidget(clipGroup);

		// --- Model placement (gizmo) -----------------------------------------------------------
		// A single UNIFIED manipulator drawn on the model/structure: 3 translate arrows, 3 rotate rings and
		// 3 scale cubes (X red / Y green / Z blue) plus a grey uniform-scale centre — all shown at once, and
		// the handle you grab picks the operation (translate/rotate/scale). Works on a fluid-viewer model AND
		// a seabed structure; the next "Apply" (Domain & resolution) re-voxelizes it exactly where you placed it.
		QGroupBox* gizmoGroup = new QGroupBox("Model placement (gizmo)");
		QVBoxLayout* gizmoCol = new QVBoxLayout(gizmoGroup);
		gizmo_enable_chk_ = new QCheckBox("Enable manipulator (move + rotate + scale)");
		gizmo_enable_chk_->setToolTip("Show the transform gizmo on the model and drag its coloured handles: arrows = move along an axis, rings = rotate about an axis, cubes = scale (a grey centre cube scales uniformly) — X red, Y green, Z blue. A left-drag away from any handle still orbits the camera. 'Apply' (Domain & resolution) re-voxelizes the model where you placed it.");
		connect(gizmo_enable_chk_, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setGizmoEnabled(on); });
		gizmoCol->addWidget(gizmo_enable_chk_);
		QPushButton* resetPlaceBtn = new QPushButton("Reset placement");
		resetPlaceBtn->setToolTip("Return the model to the default centre-on-bed placement.");
		connect(resetPlaceBtn, &QPushButton::clicked, this, [this] { if (viewer_) viewer_->resetModelPlacement(); });
		gizmoCol->addWidget(resetPlaceBtn);
		gizmo_info_ = new QLabel("(load a STEP model)");
		gizmo_info_->setWordWrap(true);
		gizmo_info_->setStyleSheet("color:#9aa;");
		gizmoCol->addWidget(gizmo_info_);
		gizmo_group_ = gizmoGroup;
		gizmoGroup->setEnabled(false);
		if (viewer_) connect(viewer_, &SliceViewer::modelPlacementChanged, this, [this] { updateGizmoUi(); });
		col->addWidget(gizmoGroup);

		col->addStretch(1);

		// Scrollable so more groups/params can be added without overflowing the window.
		QScrollArea* scroll = new QScrollArea;
		scroll->setWidgetResizable(true);
		scroll->setWidget(panel);
		scroll->setMinimumWidth(250);
		dock->setWidget(scroll);
		addDockWidget(Qt::LeftDockWidgetArea, dock);

		syncGridControls(); // seed the domain/h boxes + readout from the current sim
	}

	void MainWindow::spawnWorker(std::unique_ptr<scour::core::ChannelFluidCore> core, SeabedScenario scen,
		long long steps0, double t0)
	{
		// The checkpoint hand-off carries a shared_ptr across a queued connection — register it once.
		static bool s_meta = false;
		if (!s_meta) { qRegisterMetaType<scour::gui::CheckpointStatePtr>("scour::gui::CheckpointStatePtr"); s_meta = true; }

		worker_ = new SimWorker(std::move(core)); // no parent: moved to worker_thread_
		if (play_btn_) worker_->setPlaying(play_btn_->isChecked()); // honour the current play/pause state
		worker_->primeCounters(steps0, t0);          // a scene restore resumes from the saved step (else 0)
		worker_->setAutosaveInterval(autosave_interval_); // keep auto-saving across a rebuild/restore
		worker_->setDisplayInterval(display_throttle_s_); // keep the "Fast sim" graphics throttle across a rebuild/restore
		// Rebuild factory: a re-inject rebuilds the core from an immutable recipe snapshot, entirely on
		// the worker thread (make_core is Qt-free and thread-safe). Unused in the seabed scenario.
		worker_->setRebuildFactory([recipe = recipe_](const std::vector<unsigned char>& solid, int mode)
			{ return make_core(recipe, solid, mode); });
		connect(worker_, &SimWorker::checkpointReady, this, &MainWindow::onCheckpointReady, Qt::QueuedConnection);
		worker_thread_ = new QThread(this);
		worker_->moveToThread(worker_thread_);
		connect(worker_thread_, &QThread::started, worker_, &SimWorker::run);
		connect(worker_, &SimWorker::stats, this, [this](qint64 steps, double t, double dt, double simFps) {
			last_steps_ = steps; last_sim_time_ = t;
			// Always show the live step/time/dt + the measured sim throughput (steps/s) — even in the seabed
			// scenario, where the sim is clearly progressing (the old !seabed_ guard froze this readout).
			// simFps is the SIM rate (climbs under "Fast sim"), distinct from the render fps on the right.
			// The pre-calibration disclaimer is appended as a marker so it is not lost.
			QString s = QString("step %1   t = %2 s   dt = %3 ms   %4 sim fps")
				.arg(steps).arg(t, 0, 'f', 3).arg(dt * 1e3, 0, 'f', 2).arg(simFps, 0, 'f', 1);
			if (worker_ && worker_->tidalOn())
			{
				const double u = worker_->tidePhase();
				const char* dir = qAbs(u) < 0.02 ? "slack" : (u > 0.0 ? "flood +x" : "ebb −x");
				s += QString("   |   tide %1 m/s (%2)").arg(u, 0, 'f', 2).arg(QString::fromUtf8(dir));
			}
			if (seabed_) s += "   |   seabed (pre-calibration, qualitative)";
			status_->setText(s);
		}, Qt::QueuedConnection);
		viewer_->setWorker(worker_);

		// --- Erodible-seabed scenario: attach the morphodynamic engine + bed viz -----------------
		if (scen.active)
		{
			seabed_ = true;
			last_sp_ = scen.params; have_sp_ = true; // remember the sediment params for a later "Add sand"
			// Reflect the scenario's MORFAC in the control (signals blocked so it shows, not re-pushes — the
			// engine already carries it). The user can then dial it live via setMorfac.
			if (morfac_spin_) { morfac_spin_->blockSignals(true); morfac_spin_->setValue(scen.params.morfac); morfac_spin_->blockSignals(false); }
			worker_->setMorpho(std::move(scen.engine), scen.spinup);
			viewer_->setBedSurface(recipe_.info.nx, recipe_.info.ny, scen.bed_z0);
			// The default horizontal slice (z≈mid-height) is an opaque sheet that would occlude the
			// bed from an overhead view. In the seabed scenario the bed IS the subject, so hide the
			// slice by default (unchecking the box drives viewer_->setShowSlice(false)); the user can
			// re-enable it for a flow cross-section. Fluid viewers keep the slice on.
			if (slice_chk_) slice_chk_->setChecked(false);
			if (scen.has_structure && !scen.structure_mesh.empty())
			{
				const auto& p = scen.structure_place;
				scene_mesh_ = scen.structure_mesh;           // persist for a scene save (before the move below)
				scene_place_ = p;
				// Seed the gizmo from the structure placement (NOT a fixed override matrix) so the seabed
				// structure stays MOVABLE: drag to reposition, then Apply re-voxelizes it into the sand at
				// the new placement (place-before-run). setMesh seeds the gizmo pivot from the bbox; then
				// setModelPlacement drives it to `p` (mesh_override_ off ⇒ hasModelPlacement() true ⇒ editable).
				viewer_->setMesh(scen.structure_mesh);
				viewer_->setModelPlacement(p);
			}
			if (status_) status_->setText("PRE-CALIBRATION DEMO — constants nominal, results qualitative");
		}
		else
			seabed_ = false;

		updateGizmoUi(); // gizmo is fluid-viewer-only (seabed structure is fixed)
		worker_thread_->start();
		maybeAutoRecord(); // if armed in the Record Video… dialog, start recording a fresh seabed/drop run
	}

	void MainWindow::teardownWorkerForReload()
	{
#ifdef SCOUR_HAVE_JOLT
		// A non-drop rebuild (scene load / Apply / Add sand / commit) ends the drop sink-coupling. commitDrop
		// re-arms it AFTER this teardown; the other callers leave it off (the pile is gone/replaced).
		drop_active_ = false;
		if (resettle_timer_) resettle_timer_->stop();
#endif
		// Stop + join + delete the worker AND its thread so a fresh core/engine can replace it with no
		// concurrency (the core is only ever (re)built while the worker thread is fully joined — never
		// racing a live worker, and construction touches no GL). Deleting the thread each time avoids a
		// QThread accumulating across repeated reloads.
		if (worker_) worker_->stop();
		if (worker_thread_)
		{
			worker_thread_->quit();
			worker_thread_->wait();
			delete worker_thread_;
			worker_thread_ = nullptr;
		}
		delete worker_;
		worker_ = nullptr;
		if (viewer_) viewer_->setWorker(nullptr);
		worker_down_ = false;
		model_injected_ = false;
	}

	bool MainWindow::startSeabedScenario(const QString& path)
	{
		teardownWorkerForReload();

		SimRecipe recipe; SeabedScenario scen; std::string warn;
		auto core = build_seabed_sim(path.toStdString(), recipe, scen, warn);
		if (!core)
		{
			statusBar()->showMessage(QString("seabed scenario load failed: %1").arg(QString::fromStdString(warn)), 6000);
			return false;
		}
		if (!warn.empty()) std::fprintf(stderr, "[seabed] %s\n", warn.c_str());
		recipe_ = recipe;

		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("ScourProtection — seabed scenario [%1: %2x%3x%4, h=%5 m]")
			.arg(QString::fromStdString(info.name)).arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h));
		if (viewer_)
		{
			viewer_->clearVoxelOverlay();
			viewer_->clearMesh();
			viewer_->setInfo(info); // re-frame the camera + rebuild slice geometry for the new grid
		}
		model_mesh_ = scour::core::TriMesh{};
		added_sand_m_ = 0.0; // a real scenario FILE: a grid Apply rebuilds it via build_seabed_sim, not the conversion path
		spawnWorker(std::move(core), std::move(scen));
		syncGridControls(); // the scenario's domain/h now drive the editable boxes
		statusBar()->showMessage("seabed scenario loaded", 4000);
		return true;
	}

	void MainWindow::applyGridValues(double Lx, double Ly, double Lz, double h)
	{
		if (!lx_spin_) return;
		const SimInfo& info = recipe_.info;
		{
			const QSignalBlocker b1(lx_spin_), b2(ly_spin_), b3(lz_spin_), b4(h_spin_), b5(u_spin_);
			lx_spin_->setValue(Lx > 0.0 ? Lx : info.Lx);
			ly_spin_->setValue(Ly > 0.0 ? Ly : info.Ly);
			lz_spin_->setValue(Lz > 0.0 ? Lz : info.Lz);
			h_spin_->setValue(h > 0.0 ? h : info.h);
			u_spin_->setValue(info.U);
		}
		updateGridReadout();
		applyGrid();
	}

	void MainWindow::setInputSpeed(double U)
	{
		if (!u_spin_ || U <= 0.0) return;
		u_spin_->setValue(U); // fires valueChanged → worker_->setInletSpeed + viewer_->setReferenceU (live, no reset)
	}

	void MainWindow::setMorfacValue(double m)
	{
		if (!morfac_spin_ || m <= 0.0) return;
		morfac_spin_->setValue(m);            // fires valueChanged → worker_->setMorfac when the value changes
		if (worker_) worker_->setMorfac(m);   // also apply directly, so an unchanged value still lands on the engine
	}

	// Push the current tidal-reversal control values to the worker (live, no reset). While the tide is on
	// it drives the inlet each step, so the manual "Input speed U" is greyed to avoid a fight; the viewer's
	// reference speed tracks U_max so the auto-range/arrow scale stays sensible through the reversal.
	void MainWindow::syncTidal()
	{
		if (!tidal_chk_) return;
		const bool on = tidal_chk_->isChecked();
		const double Umax = tidal_umax_spin_->value(), plateau = tidal_plateau_spin_->value(), ramp = tidal_ramp_spin_->value();
		if (worker_) worker_->setTidalReversal(on, Umax, plateau, ramp);
		if (on && viewer_) viewer_->setReferenceU(Umax);
		if (u_spin_) u_spin_->setEnabled(!on); // the tide owns the inlet speed while active
	}

	void MainWindow::setTidalReversal(bool on, double Umax, double plateau_s, double ramp_s)
	{
		if (!tidal_chk_) return;
		const QSignalBlocker b0(tidal_chk_), b1(tidal_umax_spin_), b2(tidal_plateau_spin_), b3(tidal_ramp_spin_);
		if (Umax > 0.0) tidal_umax_spin_->setValue(Umax);
		if (plateau_s > 0.0) tidal_plateau_spin_->setValue(plateau_s);
		if (ramp_s > 0.0) tidal_ramp_spin_->setValue(ramp_s);
		tidal_chk_->setChecked(on);
		syncTidal(); // blockers suppressed the per-control signals; push once here
	}

	void MainWindow::setDisplayThrottle(double seconds)
	{
		display_throttle_s_ = seconds < 0.0 ? 0.0 : seconds;
		if (worker_) worker_->setDisplayInterval(display_throttle_s_); // worker: throttle the per-step snapshot
		if (repaint_timer_)
		{
			// Repaint no faster than snapshots arrive, but clamp to ≤ 0.5 s so a paused / interactive frame
			// still appears promptly even at a long (or "only when paused") throttle. 0 ⇒ live 60 Hz.
			const double t = display_throttle_s_;
			int ms = (t <= 0.0) ? 16 : (int)(((t < 0.5) ? t : 0.5) * 1000.0 + 0.5);
			repaint_timer_->setInterval(ms < 16 ? 16 : ms);
		}
	}

	void MainWindow::addSandMeters(double depth)
	{
		if (sand_spin_ && depth > 0.0) sand_spin_->setValue(depth);
		addSand();
	}

	void MainWindow::syncGridControls()
	{
		if (!lx_spin_) return;
		const SimInfo& info = recipe_.info;
		const QSignalBlocker b1(lx_spin_), b2(ly_spin_), b3(lz_spin_), b4(h_spin_), b5(u_spin_);
		lx_spin_->setValue(info.Lx);
		ly_spin_->setValue(info.Ly);
		lz_spin_->setValue(info.Lz);
		h_spin_->setValue(info.h);
		u_spin_->setValue(info.U);
		if (inlet_profile_box_) // reflect the built inlet mode (seabed ⇒ log-law, fluid ⇒ uniform)
		{
			const QSignalBlocker bp(inlet_profile_box_);
			inlet_profile_box_->setCurrentIndex(recipe_.bc.inlet_mode == scour::core::INLET_LOGLAW ? 1 : 0);
		}
		if (sed_bc_box_) // reflect the scenario's suspended-sediment boundary (open-sea by default)
		{
			const QSignalBlocker bs(sed_bc_box_);
			int m = have_sp_ ? last_sp_.sed_bc : scour::core::SED_BC_OPEN;
			sed_bc_box_->setCurrentIndex(m == scour::core::SED_BC_OPEN ? 0 : m == scour::core::SED_BC_RECYCLE ? 1 : 2);
		}
		updateGridReadout();
	}

	void MainWindow::updateGridReadout()
	{
		if (!grid_readout_ || !lx_spin_) return;
		int nx = 0, ny = 0, nz = 0;
		grid_dims_for(lx_spin_->value(), ly_spin_->value(), lz_spin_->value(), h_spin_->value(), nx, ny, nz);
		const long long cells = (long long)nx * ny * nz;
		// Rough device-memory gauge (fluid + snapshot double fields, seabed adds more) — an order-of-
		// magnitude hint only, not an allocation contract.
		const double mb = cells * 240.0 / (1024.0 * 1024.0);
		const bool ok = cells > 0 && cells <= kMaxCells;
		QString txt = QString("→ %1 × %2 × %3 = %4 cells  (~%5 MB)")
			.arg(nx).arg(ny).arg(nz).arg(cells).arg(mb, 0, 'f', 0);
		if (!ok) txt += QString("\n⚠ too fine (> %1 M cells) — increase h").arg(kMaxCells / 1000000);
		grid_readout_->setText(txt);
		grid_readout_->setStyleSheet(ok ? QString() : QString("color:#c0392b;"));
		if (apply_btn_) apply_btn_->setEnabled(ok);
	}

	void MainWindow::updateGizmoUi()
	{
		// Enabled whenever a gizmo-editable model/structure exists — a fluid-viewer model, a seabed structure
		// (real scenario, Add-sand conversion, or restored scene), all of which are seeded into the gizmo. The
		// "Enable manipulator" checkbox (gizmo_on_) persists in the viewer across an Apply, so nothing to reset.
		const bool ok = viewer_ && viewer_->hasModelPlacement();
		if (gizmo_group_) gizmo_group_->setEnabled(ok);
		if (gizmo_enable_chk_ && viewer_ && gizmo_enable_chk_->isChecked() != viewer_->gizmoEnabled())
		{
			const QSignalBlocker b(gizmo_enable_chk_);
			gizmo_enable_chk_->setChecked(viewer_->gizmoEnabled()); // keep the box in sync with the viewer
		}
		if (!gizmo_info_) return;
		if (!ok) { gizmo_info_->setText("(load a STEP model)"); return; }
		const SliceViewer::ModelGizmoXform x = viewer_->modelXform();
		const QVector3D e = x.rot.toEulerAngles(); // (pitch, yaw, roll) = rotation about X, Y, Z [deg]
		gizmo_info_->setText(QString("pos (%1, %2, %3) m\nrot (%4, %5, %6)°   scale (%7, %8, %9)")
			.arg(x.t.x(), 0, 'f', 2).arg(x.t.y(), 0, 'f', 2).arg(x.t.z(), 0, 'f', 2)
			.arg(e.x(), 0, 'f', 0).arg(e.y(), 0, 'f', 0).arg(e.z(), 0, 'f', 0)
			.arg(x.scale.x(), 0, 'f', 2).arg(x.scale.y(), 0, 'f', 2).arg(x.scale.z(), 0, 'f', 2));
	}

	void MainWindow::nudgeModelPlacement(double dx, double dy, double dz, double rzDeg, double scale)
	{
		if (!viewer_ || !viewer_->hasModelPlacement())
		{
			statusBar()->showMessage("no gizmo-editable model loaded", 4000);
			return;
		}
		SliceViewer::ModelGizmoXform x = viewer_->modelXform();
		x.t += QVector3D((float)dx, (float)dy, (float)dz);
		if (rzDeg != 0.0) x.rot = QQuaternion::fromAxisAndAngle(0.0f, 0.0f, 1.0f, (float)rzDeg) * x.rot;
		if (scale > 0.0) x.scale *= (float)scale;
		viewer_->setModelXform(x);
		if (!seabed_) setModelAsObstacle(true); // re-voxelize where placed (fluid viewer)
		updateGizmoUi();
		std::fprintf(stderr, "[G1] model placement: t=(%.3f %.3f %.3f) m, rz=%.1f deg, scale=%.3f\n",
			x.t.x(), x.t.y(), x.t.z(), rzDeg, x.scale.x());
	}

	void MainWindow::applyGrid()
	{
		if (!lx_spin_) return;
		GridOverride ov;
		ov.active = true;
		ov.Lx = lx_spin_->value(); ov.Ly = ly_spin_->value(); ov.Lz = lz_spin_->value();
		ov.h = h_spin_->value();
		ov.U = u_spin_ ? u_spin_->value() : 0.0; // keep the "Input speed" the user dialed in across the rebuild
		if (ov.h <= 0.0 || ov.Lx <= 0.0 || ov.Ly <= 0.0 || ov.Lz <= 0.0)
		{
			statusBar()->showMessage("invalid domain/voxel size", 4000);
			return;
		}

		// Rebuild the SAME sim (scenario or fluid viewer) at the new grid: re-run the exact setup path
		// this recipe came from with only the domain + h overridden — every physics/scenario parameter
		// is preserved. A fluid viewer keeps its loaded STEP model (re-displayed, re-voxelized at the
		// new h if it was injected); a scenario re-voxelizes its own structure + re-inits the sand bed.
		// Three rebuild paths at the new grid:
		//  - real_scenario: loaded from a seabed scenario FILE → build_seabed_sim re-voxelizes its own
		//    structure + re-inits the sand from the config.
		//  - converted: a fluid viewer turned into a seabed run at runtime via "Add sand" (added_sand_m_>0)
		//    → the config has no structure/sand, so rebuild the seabed HERE with the loaded model (or the
		//    config obstacle) as the rigid structure + the stored sand depth. This is what keeps a resized
		//    model-plus-sand run from reverting to the config's default cylinder.
		//  - fluid: plain viewer → build_sim, then re-voxelize the loaded model as the obstacle (if any).
		// A committed drop/settle PILE (still active) must SURVIVE the resize: capture it BEFORE the
		// teardown below (which clears drop_active_) and route it through the converted-seabed path,
		// re-voxelizing the multi-block placements at the new grid instead of reverting to the config
		// obstacle. `have_pile` is a plain bool in both build configs (false when Jolt is compiled out).
#ifdef SCOUR_HAVE_JOLT
		const bool have_pile = drop_active_ && !drop_placements_.empty() && !drop_mesh_.empty();
#else
		const bool have_pile = false;
#endif
		const bool converted = added_sand_m_ > 0.0 || have_pile;
		const bool real_scenario = recipe_.is_scenario && !converted;
		const bool fluid = !real_scenario && !converted;
		const std::string src = recipe_.source_config;
		const bool had_model = !model_mesh_.empty();
		scour::core::TriMesh keep_mesh = had_model ? model_mesh_ : scour::core::TriMesh{};

		// Preserve the user's gizmo placement (move/rotate/scale) across the rebuild so the model/structure
		// is re-voxelized WHERE IT WAS PLACED, not re-centred. keep_place is the affine for voxelization;
		// keep_x restores the full gizmo state so the manipulator keeps editing from the same transform.
		// Captured from the VIEWER (not model_mesh_), so it works for a fluid obstacle AND a seabed structure
		// (which has no model_mesh_ but IS gizmo-editable). The "Enable manipulator" state persists in the viewer.
		const bool have_place = viewer_ && viewer_->hasModelPlacement();
		scour::core::ModelPlacement keep_place = have_place
			? viewer_->modelPlacement()
			: (had_model ? place_model_on_bed(keep_mesh, recipe_.info.Lx, recipe_.info.Ly) : scour::core::ModelPlacement{});
		SliceViewer::ModelGizmoXform keep_x = viewer_ ? viewer_->modelXform() : SliceViewer::ModelGizmoXform{};

		teardownWorkerForReload();

		SimRecipe recipe; SeabedScenario scen; std::string warn;
		std::unique_ptr<scour::core::ChannelFluidCore> core;
		if (real_scenario)
			// Pass the gizmo placement so a repositioned structure re-voxelizes where placed (place-before-run).
			core = build_seabed_sim(src, recipe, scen, warn, &ov, have_place ? &keep_place : nullptr);
		else if (converted)
		{
			using namespace scour::core;
			// Get the new grid + config obstacle at the new resolution (fluid recipe, not spawned).
			SimRecipe frec; std::string fw;
			build_sim(src, frec, fw, &ov);
			const MacGrid ng = frec.grid;
			const double U = (ov.U > 0.0) ? ov.U : frec.info.U;
			// Rigid structure at the new grid, in priority order:
			//  - a settled drop/settle PILE (have_pile): re-voxelize the SAME unit mesh at every resting
			//    block pose (world-fixed placements) so a resize KEEPS the pile — the persistence fix;
			//  - a single loaded model: voxelize it at the gizmo placement (moved/rotated/scaled);
			//  - neither: the config obstacle.
			ModelPlacement place = keep_place;
			std::vector<unsigned char> rigid;
			if (have_pile)
			{
#ifdef SCOUR_HAVE_JOLT
				int pile_cells = 0;
				rigid = voxelize_mesh_instances(drop_mesh_, ng, drop_placements_, &pile_cells);
				std::fprintf(stderr, "[drop] pile re-voxelized on resize: %d rigid cells across %zu blocks\n",
					pile_cells, drop_placements_.size());
#endif
			}
			else if (had_model)
				rigid = voxelize_mesh(keep_mesh, ng, place, nullptr, nullptr, nullptr, nullptr);
			else
				rigid = frec.base_solid;
			const SeabedParams* base = have_sp_ ? &last_sp_ : nullptr;
			core = build_seabed_from_structure(ng, U, added_sand_m_, rigid, base, src, /*spinup*/ 200,
				/*nu_fluid*/ 5.0e-3, recipe, scen, warn, /*build_core*/ true);
			// A single model is displayed as the structure mesh; a pile is displayed via its per-block
			// placements (re-attached after spawnWorker below), so leave scen.structure_mesh empty for it.
			if (core && had_model && !have_pile)
			{
				scen.structure_mesh = keep_mesh;   // display the model as the structure at its bed placement
				scen.structure_place = place;
				scen.has_structure = true;
			}
		}
		else
			core = build_sim(src, recipe, warn, &ov);
		if (!core)
		{
			statusBar()->showMessage(QString("apply failed: %1").arg(QString::fromStdString(warn)), 6000);
			return;
		}
		if (!warn.empty()) std::fprintf(stderr, "[G1] %s\n", warn.c_str());
		recipe_ = recipe;

		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("ScourProtection — %1 [%2: %3x%4x%5, h=%6 m]")
			.arg((real_scenario || converted) ? "seabed" : "G1 slice viewer")
			.arg(QString::fromStdString(info.name)).arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h));
		if (viewer_)
		{
			viewer_->clearVoxelOverlay();      // drop stale-size staircase geometry
			viewer_->clearMesh();              // spawnWorker re-attaches the scenario/structure mesh
			viewer_->setInfo(info);            // re-frame + re-size all grid-derived viewer geometry
		}
		model_mesh_ = keep_mesh;               // keep the loaded model (fluid re-inject / converted structure)
		spawnWorker(std::move(core), std::move(scen));

		// Fluid viewer with a loaded model: re-show it and re-voxelize it as the obstacle at the new h (it
		// IS the obstacle). The converted-seabed path already seated the model as the structure above.
		if (fluid && had_model && viewer_)
		{
			viewer_->setMesh(scour::core::TriMesh(model_mesh_)); // display copy (model_mesh_ retained)
			if (keep_x.valid) viewer_->setModelXform(keep_x);    // re-apply the user's placement (setMesh reset it)
			setModelAsObstacle(true);                            // re-voxelizes at viewer_->modelPlacement()
			updateGizmoUi();
		}

#ifdef SCOUR_HAVE_JOLT
		// A settled drop pile survives the resize: re-display it at its (world-fixed) block poses and
		// re-arm the G3.3 sink coupling so it keeps sinking into the developing scour hole. build_core
		// reset the flow/bed to t=0 (Apply always does), but the PILE persists as the rigid structure.
		if (have_pile && viewer_)
		{
			viewer_->setMesh(scour::core::TriMesh(drop_mesh_));
			viewer_->setBlockPlacements(drop_placements_);
			drop_active_ = true;              // teardownWorkerForReload cleared it; the pile is still here
			resettle_last_zb_.clear();
			if (resettle_timer_) resettle_timer_->start(2000);
		}
#endif

		syncGridControls(); // reflect the built grid (floor/clamp may differ from the typed value)
		statusBar()->showMessage(QString("grid applied: %1x%2x%3 @ h=%4 m — sim reset to t=0")
			.arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h), 5000);
	}

	void MainWindow::addSand()
	{
		if (!worker_ || !sand_spin_) return;
		const double X = sand_spin_->value();
		if (X <= 0.0) { statusBar()->showMessage("set a sand depth > 0 to add sand", 4000); return; }

		const scour::core::MacGrid g = recipe_.grid;

		// The RIGID structure = the kind-2 cells of the live flow mask (structure / obstacle), EXCLUDING
		// any existing sand (kind 1) so re-adding sand doesn't freeze the old bed. Works from a fluid
		// viewer (all solids rigid) or an existing seabed scenario (only the structure is rigid).
		std::vector<unsigned char> rigid((size_t)g.p_count(), 0);
		{
			std::vector<unsigned char> kind; scour::core::MacGrid mg;
			if (worker_->copyMask(kind, mg) && (int)kind.size() == g.p_count())
			{
				for (std::size_t n = 0; n < kind.size(); ++n) if (kind[n] == 2) rigid[n] = 1;
			}
			else if ((int)recipe_.base_solid.size() == g.p_count())
			{
				for (std::size_t n = 0; n < rigid.size(); ++n) if (recipe_.base_solid[n]) rigid[n] = 1;
			}
		}

		const scour::core::SeabedParams* base = have_sp_ ? &last_sp_ : nullptr;

		// Build ONLY the morphodynamic engine + initial sand+structure mask + recipe (build_core=false):
		// we do NOT construct a fresh flow core, and we do NOT tear down the worker. The developed flow is
		// PRESERVED — the worker installs the sand mask in place via ChannelFluidCore::update_solid, so you
		// can spin the flow up to steady state and then add sand without resetting it. The now-solid sand
		// cells are zeroed; every untouched fluid cell keeps its {u,v,w,p}. The loaded model (if any) stays
		// displayed where it is (sand fills around it); a procedural obstacle shows via the voxel overlay.
		SimRecipe recipe; SeabedScenario scen; std::string warn;
		build_seabed_from_structure(g, recipe_.info.U, X, rigid, base,
			recipe_.source_config, /*spinup*/ 200, /*nu_fluid*/ 5.0e-3, recipe, scen, warn, /*build_core*/ false);
		if (!scen.engine)
		{
			statusBar()->showMessage(QString("add sand failed: %1").arg(QString::fromStdString(warn)), 6000);
			return;
		}
		if (!warn.empty()) std::fprintf(stderr, "[seabed] %s\n", warn.c_str());

		// Adopt the seabed state on the GUI side (bed viz + params) WITHOUT resetting the view.
		seabed_ = true;
		added_sand_m_ = X; // mark this as a runtime conversion so a later grid Apply rebuilds it correctly
		last_sp_ = scen.params; have_sp_ = true;
		recipe_ = recipe; // provenance + inlet-profile combo (now log-law); grid is unchanged
		if (viewer_) viewer_->setBedSurface(recipe.info.nx, recipe.info.ny, scen.bed_z0);

		// Hand the engine + initial mask to the worker; it converts the running core in place (flow kept).
		worker_->requestSeabedConversion(std::move(scen.engine), recipe.base_solid, /*spinup*/ 200,
			/*loglaw*/ recipe.bc.inlet_mode == scour::core::INLET_LOGLAW, recipe.bc.z0, recipe.bc.bed_datum);

		setWindowTitle(QString("ScourProtection — seabed (added sand, flow kept) [%1x%2x%3, h=%4 m, sand=%5 m]")
			.arg(recipe.info.nx).arg(recipe.info.ny).arg(recipe.info.nz).arg(recipe.info.h).arg(X, 0, 'f', 2));
		syncGridControls();
		updateGizmoUi(); // now a seabed run ⇒ the model placement gizmo is disabled (structure fixed)
		maybeAutoRecord(); // "Add sand" begins a seabed run in place (no spawnWorker) — arm auto-record here too
		statusBar()->showMessage(QString("added %1 m sand — flow preserved (sediment run continues)").arg(X, 0, 'f', 2), 6000);
	}

	// ================= Drop/settle preparation phase (PLAN G3.2) =================================
#ifdef SCOUR_HAVE_JOLT
	bool MainWindow::dropBlocks(const QString& stepPath, int count, unsigned seed,
		double scaleMin, double scaleMax, double sandDepth)
	{
		using namespace scour::core;

		// Fresh drop: discard any prior pile + stop the sink re-settle before building the new one.
		drop_active_ = false;
		if (resettle_timer_) resettle_timer_->stop();
		settle_world_.reset();
		resettle_last_zb_.clear();

		// Load the XStone: the full (holed) mesh for DISPLAY + VOXELIZATION, and the pre-cut convex solids
		// for the Jolt collision compound (fallback: the whole mesh as a single convex hull).
		std::string err;
		TriMesh mesh = load_step_mesh(stepPath.toStdString(), 0.1, &err);
		if (mesh.empty())
		{
			statusBar()->showMessage(QString("drop: STEP load failed: %1").arg(QString::fromStdString(err)), 6000);
			std::fprintf(stderr, "[drop] STEP load failed: %s\n", err.c_str());
			return false;
		}
		std::vector<TriMesh> solids = load_step_solids(stepPath.toStdString(), 0.1, nullptr);
		ConvexShape shape;
		if (!solids.empty()) shape.pieces = std::move(solids);
		else shape.pieces.push_back(mesh); // bootstrap: whole-mesh hull
		std::fprintf(stderr, "[drop] XStone: %zu display tris, %zu convex pieces; dropping %d @ scale [%.2f,%.2f], sand=%.2f m\n",
			mesh.triangle_count(), shape.pieces.size(), count, scaleMin, scaleMax, sandDepth);

		// Scatter `count` units above the sand surface, centred on the domain, at random orientation/scale.
		const SimInfo& info = recipe_.info;
		const std::array<float, 3> bsz = mesh.bbox_size();
		const double unit_h = std::max({ (double)bsz[0], (double)bsz[1], (double)bsz[2] }) * scaleMax;
		const double footprint = std::min(info.Lx, info.Ly) * 0.6;
		// Drop the units from well ABOVE the sim domain (MH's ask): spawn 5 m higher than the domain
		// z-size, staggered over ~unit heights so they don't all release from the exact same z. The Jolt
		// settling world is decoupled from the CUDA fluid, so a spawn above Lz is fine — the units simply
		// fall further before landing on the sand at z = sandDepth. The ~10–13 m/s impact from that height
		// is caught by the CCD (LinearCast) + sub-stepping in settle.cpp, so nothing tunnels the floor.
		const double drop_lo = info.Lz + 5.0;              // 5 m above the top of the simulation space
		double drop_hi = drop_lo + 2.5 * unit_h;           // staggered release heights
		// Scatter is sorted LARGEST-FIRST; we release the units over time in that order (big ones drop and
		// settle first, smaller ones land among/on them later — an armour-layer pour).
		std::vector<BlockInstance> instances = scatter_blocks(count, 0.5 * info.Lx, 0.5 * info.Ly, footprint,
			drop_lo, drop_hi, scaleMin, scaleMax, bsz, seed);

		// Build an EMPTY settling world (just the ground); stepDropAnimation releases the units one at a
		// time. Land them on the sand surface at z = sandDepth.
		GroundSpec ground; ground.bed_z = sandDepth;
		SettleParams params; params.seed = seed; params.max_steps = 30000; // generous for a staged pour
		params.reserve_bodies = count; // size Jolt's pair/contact buffers for the FULL pile (blocks are
		                               // added incrementally) — else a big pour overflows + drops contacts.
		settle_world_ = std::make_unique<SettleWorld>(shape, std::vector<BlockInstance>{}, ground, params);
		drop_queue_ = std::move(instances);
		drop_next_ = 0;
		drop_release_ctr_ = 0; // release the first (largest) unit on the first tick
		drop_placements_.clear();
		drop_mesh_ = mesh;
		drop_sand_depth_ = sandDepth;

		// Pause the underlying flow and arm the display; the units appear as they are released.
		if (worker_) worker_->setPlaying(false);
		if (viewer_)
		{
			viewer_->clearVoxelOverlay();
			viewer_->setMesh(TriMesh(mesh)); // the drawn geometry (holed XStone)
			viewer_->clearBlockPlacements();
			viewer_->setShowModel(false);    // nothing shown until the first unit drops
		}
		statusBar()->showMessage(QString("dropping %1 XStone units (largest first) — settling…").arg(count), 0);

		// Animate the settle on the MAIN thread (cheap: tens of bodies); commitDrop() runs once it rests.
		if (!drop_timer_)
		{
			drop_timer_ = new QTimer(this);
			connect(drop_timer_, &QTimer::timeout, this, &MainWindow::stepDropAnimation);
		}
		drop_timer_->start(16); // ~60 Hz
		return true;
	}

	void MainWindow::stepDropAnimation()
	{
		if (!settle_world_) { if (drop_timer_) drop_timer_->stop(); return; }

		// Release the next (largest remaining) unit on a cadence — big units drop + settle first, then
		// the smaller ones land among/on them.
		if (drop_next_ < (int)drop_queue_.size() && --drop_release_ctr_ <= 0)
		{
			settle_world_->addBlock(drop_queue_[drop_next_++]);
			drop_release_ctr_ = drop_release_ticks_;
		}

		settle_world_->step(20); // ~20 physics substeps per tick
		drop_placements_ = settle_world_->placements();
		if (viewer_)
		{
			if (!drop_placements_.empty()) viewer_->setShowModel(true);
			viewer_->setBlockPlacements(drop_placements_);
		}

		// Done once every unit is released AND the whole pile has come to rest.
		if (drop_next_ >= (int)drop_queue_.size() && settle_world_->all_asleep())
		{
			if (drop_timer_) drop_timer_->stop();
			commitDrop();
		}
	}

	void MainWindow::commitDrop()
	{
		using namespace scour::core;
		if (!settle_world_) return;
		const SettleResult res = settle_world_->result();
		drop_placements_ = res.placements;
		// KEEP settle_world_ alive for the G3.3 sink coupling (re-settle as the bed erodes below the pile).
		std::fprintf(stderr, "[drop] settled: %d steps, %s, min_z=%.3f m, max_penetration=%.3f m\n",
			res.steps, res.all_asleep ? "at rest" : "budget hit", res.min_z, res.max_penetration);

		// Voxelize the settled pile into ONE rigid-structure mask (union over the placed blocks).
		const MacGrid g = recipe_.grid;
		int solid_cells = 0;
		std::vector<unsigned char> rigid = voxelize_mesh_instances(drop_mesh_, g, drop_placements_, &solid_cells);
		std::fprintf(stderr, "[drop] pile voxelized: %d rigid cells across %zu blocks\n", solid_cells, drop_placements_.size());

		// Build a FRESH seabed run (t=0): the pile is the rigid structure, `drop_sand_depth_` m of sand
		// fills below it — the same builder the "Add sand" path uses, with build_core=true for a full reset.
		const SeabedParams* base = have_sp_ ? &last_sp_ : nullptr;
		SimRecipe recipe; SeabedScenario scen; std::string warn;
		std::unique_ptr<ChannelFluidCore> core = build_seabed_from_structure(
			g, recipe_.info.U, drop_sand_depth_, rigid, base, recipe_.source_config,
			/*spinup*/ 200, /*nu_fluid*/ 5.0e-3, recipe, scen, warn, /*build_core*/ true);
		if (!core)
		{
			statusBar()->showMessage(QString("drop commit failed: %1").arg(QString::fromStdString(warn)), 6000);
			std::fprintf(stderr, "[drop] commit failed: %s\n", warn.c_str());
			return;
		}
		if (!warn.empty()) std::fprintf(stderr, "[drop] %s\n", warn.c_str());

		// Swap in the seabed sim. The pile is DISPLAYED via the block placements (not a single structure
		// mesh) + the worker's voxel overlay, so leave scen.structure_mesh empty.
		teardownWorkerForReload();
		recipe_ = recipe;
		added_sand_m_ = drop_sand_depth_;
		model_mesh_ = TriMesh{}; // the pile is not a single-placement obstacle
		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("ScourProtection — seabed (XStone drop) [%1x%2x%3, h=%4 m, %5 blocks]")
			.arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h).arg(drop_placements_.size()));
		if (viewer_)
		{
			viewer_->clearVoxelOverlay();
			viewer_->clearMesh();
			viewer_->setInfo(info);
		}
		spawnWorker(std::move(core), std::move(scen));
		// Re-attach the pile display AFTER spawnWorker (which set up the seabed bed viz).
		if (viewer_)
		{
			viewer_->setMesh(TriMesh(drop_mesh_));
			viewer_->setBlockPlacements(drop_placements_);
		}
		std::fprintf(stderr, "[drop] committed: seabed run started (%zu blocks, %d rigid cells, sand=%.2f m)\n",
			drop_placements_.size(), solid_cells, drop_sand_depth_);
		statusBar()->showMessage(QString("XStone pile settled (%1 units) — morphodynamic run started").arg(drop_placements_.size()), 6000);

		// G3.3 sink coupling: start the periodic re-settle so the pile sinks into the scour hole as the
		// bed erodes beneath it. It re-settles only when the bed has moved ≥ ½ cell since the last one.
		if (!resettle_timer_)
		{
			resettle_timer_ = new QTimer(this);
			connect(resettle_timer_, &QTimer::timeout, this, &MainWindow::maybeResettle);
		}
		resettle_last_zb_.clear();
		drop_active_ = true;
		resettle_timer_->start(2000);
	}

	void MainWindow::maybeResettle()
	{
		using namespace scour::core;
		if (!drop_active_ || !settle_world_ || !worker_ || drop_mesh_.empty()) return;

		std::vector<float> zb;
		if (!worker_->copyBed(zb)) return;
		const int nx = recipe_.info.nx, ny = recipe_.info.ny;
		if ((int)zb.size() != nx * ny) return;

		// First tick after commit just records the (flat) baseline; later ticks trigger on bed change.
		if ((int)resettle_last_zb_.size() != (int)zb.size()) { resettle_last_zb_ = zb; return; }
		const double h = recipe_.info.h;
		double dmax = 0.0;
		for (std::size_t n = 0; n < zb.size(); ++n) dmax = std::max(dmax, (double)std::fabs(zb[n] - resettle_last_zb_[n]));
		if (dmax < 0.5 * h) return; // not enough erosion yet
		resettle_last_zb_ = zb;

		// Downsample z_b to a coarse world-Z-up bed surface for the collision mesh (the scour hole is
		// smooth, so a ~64² patch captures it cheaply). Sample (si,sj) at column (si·stride, sj·stride).
		const int stride = std::max(1, std::max(nx, ny) / 64);
		BedField bed;
		bed.nx = (nx + stride - 1) / stride;
		bed.ny = (ny + stride - 1) / stride;
		bed.x0 = 0.5 * h; bed.y0 = 0.5 * h; bed.h = stride * h;
		bed.z.assign((size_t)bed.nx * bed.ny, 0.0f);
		for (int sj = 0; sj < bed.ny; ++sj)
			for (int si = 0; si < bed.nx; ++si)
			{
				const int i = std::min(nx - 1, si * stride), j = std::min(ny - 1, sj * stride);
				bed.z[(size_t)sj * bed.nx + si] = zb[(size_t)j * nx + i];
			}

		// Re-settle the pile onto the eroded bed (blocks that lost their sand support sink into the hole).
		const SettleResult res = settle_world_->resettleOnBed(bed);
		drop_placements_ = res.placements;
		if (viewer_) viewer_->setBlockPlacements(drop_placements_);

		// Re-voxelize the moved pile and hand the new structure to the worker (flow preserved in place).
		int solid_cells = 0;
		std::vector<unsigned char> rigid = voxelize_mesh_instances(drop_mesh_, recipe_.grid, drop_placements_, &solid_cells);
		worker_->requestStructureUpdate(std::move(rigid));
		std::fprintf(stderr, "[drop] re-settle: bed moved %.3f m (≥ %.3f) → pile sank (%d steps, min_z=%.3f), %d rigid cells\n",
			dmax, 0.5 * h, res.steps, res.min_z, solid_cells);
	}

	// Dock "Clear pile / re-drop": stop the fall animation + the sink re-settle and remove the pile display,
	// so the user can adjust parameters and Drop & settle again. Does NOT rebuild the sim — a fresh
	// Drop & settle (or an Apply) does that; this just tears down the transient drop/settle state.
	void MainWindow::clearDropPile()
	{
		drop_active_ = false;
		if (drop_timer_) drop_timer_->stop();
		if (resettle_timer_) resettle_timer_->stop();
		settle_world_.reset();
		drop_queue_.clear();
		drop_placements_.clear();
		drop_next_ = 0;
		drop_release_ctr_ = 0;
		resettle_last_zb_.clear();
		if (viewer_) viewer_->clearBlockPlacements();
		statusBar()->showMessage("pile cleared — adjust parameters and Drop & settle again", 4000);
	}
#endif

	// ================= Scene save / restore (scene_io) ==========================================

	void MainWindow::setAutosave(int stepsInterval)
	{
		autosave_interval_ = stepsInterval < 0 ? 0 : stepsInterval;
		if (worker_) worker_->setAutosaveInterval(autosave_interval_);
		statusBar()->showMessage(autosave_interval_ > 0
			? QString("auto-save every %1 steps → <scene>.<step>.scn").arg(autosave_interval_)
			: QString("auto-save off"), 4000);
	}

	void MainWindow::saveSceneAs()
	{
		QString start = scene_base_path_.isEmpty()
			? (QString::fromStdString(recipe_.info.name) + ".scn") : scene_base_path_;
		QString fn = QFileDialog::getSaveFileName(this, "Save Scene", start, "Scene files (*.scn);;All files (*)");
		if (fn.isEmpty()) return;
		if (!fn.endsWith(".scn", Qt::CaseInsensitive)) fn += ".scn";
		saveSceneNow(fn);
	}

	bool MainWindow::saveSceneNow(const QString& path)
	{
		if (!worker_) { statusBar()->showMessage("no active simulation to save", 4000); return false; }
		// The worker gathers the full device state on ITS thread (race-free) and emits checkpointReady;
		// onCheckpointReady then writes it. tag = -1 marks this as a manual save.
		pending_manual_save_path_ = path;
		worker_->requestCheckpoint(-1);

		// Register the scene in Recent Files + set the auto-save base IMMEDIATELY (don't wait for the async
		// write): the menu updates the instant you confirm, and the entry survives even if you close right
		// after saving. The .scn doesn't exist on disk yet, so `pending_recent_keep_` stops rebuildRecentMenu
		// from pruning it until the write lands (onCheckpointReady clears it). Numbered <name>.<step>.scn
		// checkpoints are auto-saves, not user scenes — never added.
		scene_base_path_ = path;
		if (!scour::gui::is_numbered_checkpoint(path.toStdString()))
		{
			pending_recent_keep_ = QFileInfo(path).absoluteFilePath();
			addRecentFile(path);
		}
		statusBar()->showMessage(QString("saving scene %1…").arg(QFileInfo(path).fileName()), 3000);
		return true;
	}

	void MainWindow::onCheckpointReady(scour::gui::CheckpointStatePtr state, qint64 tag)
	{
		if (!state) return;
		SceneFile sf;
		sf.def.recipe = recipe_;
		sf.def.has_bed = state->has_bed;
		sf.def.spinup = 0; // a restored run has already spun up
		sf.def.bed_z0 = state->has_bed ? state->sp.sand_depth : 0.0;
		sf.def.has_mesh = !scene_mesh_.empty();
		if (sf.def.has_mesh)
		{
			sf.def.mesh = scene_mesh_;
			// Persist the LIVE gizmo placement (move/rotate/scale) when the model is gizmo-editable, so a
			// restored scene resumes at exactly where the model was placed; scene_place_ (the default
			// centre-on-bed) is only the fallback for a seabed structure / no live placement.
			sf.def.place = (viewer_ && viewer_->hasModelPlacement()) ? viewer_->modelPlacement() : scene_place_;
		}
		sf.state = std::move(*state); // we hold the last reference — move the big host arrays into the writer

		const bool manual = (tag < 0);
		QString path;
		if (manual) { path = pending_manual_save_path_; pending_manual_save_path_.clear(); }
		else
		{
			QString base = scene_base_path_;
			if (base.isEmpty()) base = QDir::current().filePath(QString::fromStdString(recipe_.info.name) + ".scn");
			path = QString::fromStdString(checkpoint_path_for(base.toStdString(), (long long)tag));
		}
		if (path.isEmpty()) return;

		std::string warn;
		if (write_scene(path.toStdString(), sf, warn))
		{
			std::fprintf(stderr, "[scene] wrote %s (step %lld, %s)\n", path.toUtf8().constData(),
				(long long)sf.state.steps, sf.state.has_bed ? "seabed" : "fluid");
			statusBar()->showMessage(QString("%1 saved: %2")
				.arg(manual ? "scene" : "checkpoint").arg(QFileInfo(path).fileName()), 5000);
			if (manual && !pending_recent_keep_.isEmpty())
			{
				pending_recent_keep_.clear();  // the .scn is now on disk — normal pruning applies from here
				rebuildRecentMenu();           // refresh (the file exists now)
			}
		}
		else
		{
			std::fprintf(stderr, "[scene] save FAILED: %s\n", warn.c_str());
			statusBar()->showMessage(QString("scene save failed: %1").arg(QString::fromStdString(warn)), 6000);
			if (manual && !pending_recent_keep_.isEmpty())
			{
				pending_recent_keep_.clear();  // drop the optimistic entry: rebuild prunes the never-written file
				rebuildRecentMenu();
			}
		}
	}

	void MainWindow::loadSceneDialog()
	{
		QString fn = QFileDialog::getOpenFileName(this, "Load Scene", scene_base_path_,
			"Scene files (*.scn);;All files (*)");
		if (fn.isEmpty()) return;

		// Offer sibling restart points (the base scene + its auto-saved checkpoints) so the user can pick
		// a step to resume from, per "choose a restartpoint when loading".
		auto cps = scour::gui::list_checkpoints(fn.toStdString());
		if (cps.size() > 1)
		{
			QStringList items; int cur = 0;
			for (int i = 0; i < (int)cps.size(); ++i)
			{
				items << QString("step %1  —  %2").arg(cps[i].step)
					.arg(QFileInfo(QString::fromStdString(cps[i].path)).fileName());
				if (QFileInfo(QString::fromStdString(cps[i].path)) == QFileInfo(fn)) cur = i;
			}
			bool ok = false;
			QString pick = QInputDialog::getItem(this, "Restart point",
				"Choose a checkpoint to resume from:", items, cur, false, &ok);
			if (!ok) return;
			const int idx = items.indexOf(pick);
			if (idx >= 0) fn = QString::fromStdString(cps[idx].path);
		}
		loadSceneFromPath(fn);
	}

	bool MainWindow::loadSceneFromPath(const QString& path)
	{
		SceneFile sf; std::string warn;
		if (!scour::gui::read_scene(path.toStdString(), sf, warn))
		{
			std::fprintf(stderr, "[scene] load FAILED (%s): %s\n", path.toUtf8().constData(), warn.c_str());
			statusBar()->showMessage(QString("scene load failed: %1").arg(QString::fromStdString(warn)), 6000);
			return false;
		}
		restoreScene(sf);
		scene_base_path_ = path;
		if (!scour::gui::is_numbered_checkpoint(path.toStdString())) addRecentFile(path);
		statusBar()->showMessage(QString("scene restored: %1 @ step %2")
			.arg(QFileInfo(path).fileName()).arg(sf.state.steps), 6000);
		return true;
	}

	void MainWindow::restoreScene(SceneFile& sf)
	{
		SceneDefinition& d = sf.def;
		CheckpointState& st = sf.state;

		teardownWorkerForReload();
		recipe_ = d.recipe;

		// Build the core on the restored grid with the SAVED live solid mask + surface mode; then inject
		// the saved velocity/pressure fields, replacing make_core's init_uniform state. Set the live
		// inlet + bed-inlet state BEFORE load_state_host so its BC re-hold is consistent.
		std::vector<unsigned char> solid = std::move(st.solid);
		if ((int)solid.size() != recipe_.grid.p_count()) solid = recipe_.base_solid;
		auto core = make_core(recipe_, solid, st.solid_mode);
		core->set_bed_inlet_mask(st.bed_inlet_mask);
		if (st.bc.inlet_mode == scour::core::INLET_LOGLAW) core->set_inlet_profile(true, st.bc.z0, st.bc.bed_datum);
		core->set_inlet_speed(st.bc.U_inlet);
		core->load_state_host(st.u, st.v, st.w, st.p);

		// Seabed engine (if any): rebuild from params + saved structure, inject the bed state, resume
		// morphology immediately (spinup 0 — the flow is already developed).
		SeabedScenario scen;
		if (st.has_bed)
		{
			auto engine = std::make_unique<scour::core::SeabedMorpho>(recipe_.grid, st.sp, st.structure);
			engine->load_state(st.bed_G, st.susp_c, st.bed_emax, st.bed_emay, st.morpho_steps);
			scen.active = true;
			scen.engine = std::move(engine);
			scen.spinup = 0;
			scen.bed_z0 = st.sp.sand_depth;
			scen.params = st.sp;
			scen.has_structure = !st.structure.empty();
			if (d.has_mesh) { scen.structure_mesh = d.mesh; scen.structure_place = d.place; }
		}

		// GUI-side provenance so a later Apply / Add sand rebuilds correctly.
		seabed_ = st.has_bed;
		have_sp_ = st.has_bed;
		if (st.has_bed) last_sp_ = st.sp;
		added_sand_m_ = (st.has_bed && !recipe_.is_scenario) ? st.sp.sand_depth : 0.0;
		model_mesh_ = d.has_mesh ? d.mesh : scour::core::TriMesh{};
		scene_mesh_ = model_mesh_;
		scene_place_ = d.place;

		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("ScourProtection — %1 [%2: %3x%4x%5, h=%6 m] @ step %7")
			.arg(st.has_bed ? "seabed (restored)" : "restored")
			.arg(QString::fromStdString(info.name)).arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h).arg(st.steps));
		if (viewer_)
		{
			viewer_->clearVoxelOverlay();
			viewer_->clearMesh();
			viewer_->setInfo(info); // re-frame + re-size all grid-derived viewer geometry
		}

		spawnWorker(std::move(core), std::move(scen), st.steps, st.sim_time);

		// Fluid viewer with a display mesh: seabed meshes are attached inside spawnWorker; a fluid one
		// is attached here (the obstacle is already baked into the restored solid mask — no re-voxelize).
		if (!st.has_bed && d.has_mesh && viewer_)
		{
			viewer_->setMesh(scour::core::TriMesh(model_mesh_));
			viewer_->setModelPlacement(d.place); // restore the saved gizmo placement (move/rotate/scale)
			model_injected_ = true;
		}

		updateGizmoUi();
		syncGridControls();
		// Reflect the live inlet speed AFTER syncGridControls (which resets the box to the recipe U).
		if (u_spin_) { const QSignalBlocker b(u_spin_); u_spin_->setValue(st.bc.U_inlet); }
		if (viewer_) viewer_->setReferenceU(st.bc.U_inlet);
		if (worker_) worker_->setInletSpeed(st.bc.U_inlet); // latch so a later rebuild keeps the current
	}

	MainWindow::~MainWindow()
	{
		shutdownWorker();
	}

	bool MainWindow::loadStepFile(const QString& path, bool noslip)
	{
		QApplication::setOverrideCursor(Qt::WaitCursor);
		std::string err;
		scour::core::TriMesh mesh = scour::core::load_step_mesh(path.toStdString(), 0.1, &err);
		QApplication::restoreOverrideCursor();

		if (mesh.empty())
		{
			std::fprintf(stderr, "[G1] STEP load FAILED (%s): %s\n", path.toUtf8().constData(), err.c_str());
			statusBar()->showMessage(QString("STEP load failed: %1").arg(QString::fromStdString(err)), 6000);
			return false;
		}

		const auto& lo = mesh.bbox_min;
		const auto& hi = mesh.bbox_max;
		std::fprintf(stderr,
			"[G1] loaded STEP %s: %zu triangles, bbox=[%.4f %.4f %.4f]..[%.4f %.4f %.4f] m\n",
			path.toUtf8().constData(), mesh.triangle_count(),
			lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
		statusBar()->showMessage(QString("loaded %1 (%2 triangles)")
			.arg(QFileInfo(path).fileName()).arg(mesh.triangle_count()), 5000);

		model_mesh_ = mesh;            // keep a CPU copy for voxelization (viewer frees its own)
		scene_mesh_ = mesh;           // persist for a scene save (display mesh); placement recomputed below
		scene_place_ = scour::core::place_model_on_bed(mesh, recipe_.info.Lx, recipe_.info.Ly);
		if (viewer_) viewer_->setMesh(std::move(mesh));
		addRecentFile(path);          // remember it in the Recent Files menu (feature 1)

		// The loaded model IS the obstacle: voxelize it on the sim grid and inject it into the flow,
		// REPLACING the config obstacle (the default cylinder). In a seabed scenario the bed + structure
		// own the flow mask, so we only display the mesh there (setModelAsObstacle is a no-op then).
		if (!seabed_) setModelAsObstacle(true, noslip);
		updateGizmoUi(); // enable the placement gizmo for the freshly loaded (fluid) model
		return true;
	}

	void MainWindow::addRecentFile(const QString& path)
	{
		QSettings s("COBOD", "ScourProtection");
		QStringList files = s.value("recentStepFiles").toStringList();
		QString abs = QFileInfo(path).absoluteFilePath();
		files.removeAll(abs);
		files.prepend(abs);
		while (files.size() > kMaxRecent) files.removeLast();
		s.setValue("recentStepFiles", files);
		rebuildRecentMenu();
	}

	void MainWindow::rebuildRecentMenu()
	{
		if (!recent_menu_) return;
		recent_menu_->clear();
		QSettings s("COBOD", "ScourProtection");
		QStringList files = s.value("recentStepFiles").toStringList();

		// Prune entries whose file no longer exists (write the pruned list back). A scene that was just
		// saved but whose async write hasn't landed yet (pending_recent_keep_) is kept so it shows instantly.
		QStringList kept;
		for (const QString& f : files)
			if (QFileInfo::exists(f) || (!pending_recent_keep_.isEmpty() && f == pending_recent_keep_)) kept << f;
		if (kept != files) s.setValue("recentStepFiles", kept);

		if (kept.isEmpty())
		{
			QAction* none = recent_menu_->addAction("(no recent files)");
			none->setEnabled(false);
			return;
		}
		int n = 1;
		for (const QString& f : kept)
		{
			QAction* a = recent_menu_->addAction(QString("&%1  %2").arg(n++).arg(QFileInfo(f).fileName()));
			a->setToolTip(f);
			// A recent entry is either a saved scene (.scn → full restore) or a STEP model (→ obstacle).
			connect(a, &QAction::triggered, this, [this, f] {
				if (f.endsWith(".scn", Qt::CaseInsensitive)) loadSceneFromPath(f);
				else loadStepFile(f);
			});
		}
		recent_menu_->addSeparator();
		QAction* clear = recent_menu_->addAction("Clear Recent");
		connect(clear, &QAction::triggered, this, [this] {
			QSettings s2("COBOD", "ScourProtection");
			s2.remove("recentStepFiles");
			rebuildRecentMenu();
		});
	}

	bool MainWindow::setModelAsObstacle(bool on, bool noslip)
	{
		if (!worker_) return false;
		if (seabed_) return false; // the seabed scenario owns the flow mask (bed + structure)
		if (on && model_mesh_.empty())
		{
			statusBar()->showMessage("no model loaded — open a STEP first", 4000);
			return false;
		}

		if (on)
		{
			using namespace scour::core;
			// Voxelize where the user placed it: use the gizmo's live placement (move/rotate/scale) when a
			// model is loaded; otherwise the default centre-on-bed. The viewer's placement and its drawn mesh
			// share one transform, so the solid mask lands exactly under the mesh.
			ModelPlacement place = (viewer_ && viewer_->hasModelPlacement())
				? viewer_->modelPlacement()
				: place_model_on_bed(model_mesh_, recipe_.info.Lx, recipe_.info.Ly);
			double mesh_vol = 0.0, voxel_vol = 0.0;
			int thin = 0;
			std::vector<float> frac;
			std::vector<unsigned char> model_solid = voxelize_mesh(model_mesh_, recipe_.grid, place,
				&mesh_vol, &voxel_vol, &frac, &thin);

			// The loaded model REPLACES the config obstacle (the default cylinder): the injected mask is
			// the model voxels ALONE, not merged with recipe_.base_solid. "Close model" restores the
			// config obstacle via the `else` branch below.
			long long nsolid = 0;
			std::vector<unsigned char> obstacle(model_solid.size(), 0);
			for (std::size_t n = 0; n < obstacle.size(); ++n)
			{
				if (model_solid[n]) { obstacle[n] = 1; ++nsolid; }
			}

			const int mode = noslip ? SOLID_NOSLIP : SOLID_FREESLIP;
			worker_->requestRebuild(std::move(obstacle), mode);
			// The voxel overlay follows the LIVE flow mask the worker publishes after the rebuild.
			model_injected_ = true;

			QString vol = mesh_vol > 1e-12
				? QString("vol err %1%").arg(100.0 * std::abs(mesh_vol - voxel_vol) / mesh_vol, 0, 'f', 2)
				: QString("open shell (no vol gate)");
			statusBar()->showMessage(QString("model injected: %1 solid cells, %2%3 [%4]")
				.arg(nsolid).arg(vol)
				.arg(thin > 0 ? QString(", %1 thin cells kept").arg(thin) : QString())
				.arg(noslip ? "no-slip" : "free-slip"), 8000);
		}
		else
		{
			worker_->requestRebuild(recipe_.base_solid, recipe_.base_solid_mode);
			// Overlay follows the restored config mask once the worker republishes it (auto-sync).
			model_injected_ = false;
			statusBar()->showMessage("model obstacle removed (config obstacle restored)", 4000);
		}
		return true;
	}

	// ================= Video recording (video_recorder) =========================================
	bool MainWindow::startRecording(const QString& path)
	{
		if (!recorder_) recorder_ = std::make_unique<VideoRecorder>();
		if (recorder_->recording()) return true;
		const bool ok = recorder_->start(path, rec_crf_, rec_preset_, rec_fps_);
		last_capture_step_ = -1000000000LL; // capture the first frame immediately after recording starts
		last_capture_zb_.clear();
		// The menu item always opens the settings dialog (where Start/Stop live); a trailing ● just marks
		// that a recording is in progress — it is NOT a stop toggle.
		if (record_action_) record_action_->setText(ok ? QString::fromUtf8("Record Video…  \xE2\x97\x8F REC") : "Record Video…");
		if (ok) statusBar()->showMessage(
			QString("recording → %1 (h264 crf%2, one frame / %3 steps)").arg(path).arg(rec_crf_).arg(frame_step_interval_), 5000);
		updateRecordDialogStatus();
		return ok;
	}

	// File menu → "Record Video…": show the non-modal settings popup. Built lazily on first use; its Start
	// pulls the chosen params into the members and calls startRecording(), Stop finalizes the MP4, and the
	// repaint tick keeps its "● REC — N frames" status live via updateRecordDialogStatus().
	void MainWindow::openRecordDialog()
	{
		if (!video_dialog_)
		{
			video_dialog_ = new VideoSettingsDialog(this);
			connect(video_dialog_, &VideoSettingsDialog::startRequested, this, [this] {
				frame_step_interval_ = std::max<long long>(1, video_dialog_->cadenceSteps());
				frame_bed_eps_ = video_dialog_->bedEps();
				rec_crf_ = video_dialog_->crf();
				rec_preset_ = video_dialog_->preset();
				rec_fps_ = video_dialog_->fps();
				QString path = video_dialog_->path().trimmed();
				if (path.isEmpty()) path = "scour_run.mp4";
				if (!path.endsWith(".mp4", Qt::CaseInsensitive)) path += ".mp4";
				startRecording(path);
			});
			connect(video_dialog_, &VideoSettingsDialog::stopRequested, this, [this] {
				if (recorder_ && recorder_->recording())
				{
					recorder_->finish();
					if (record_action_) record_action_->setText("Record Video…");
					statusBar()->showMessage("recording stopped — MP4 finalized", 5000);
				}
				updateRecordDialogStatus();
			});
			connect(video_dialog_, &VideoSettingsDialog::autoRecordToggled, this, [this](bool on) { auto_record_ = on; });
		}
		// Seed with the live values (keep any path already resolved by the recorder).
		const QString shownPath = (recorder_ && !recorder_->path().isEmpty()) ? recorder_->path() : QString();
		video_dialog_->setValues(shownPath, (int)frame_step_interval_, rec_crf_, rec_preset_, rec_fps_, frame_bed_eps_, auto_record_);
		updateRecordDialogStatus();
		video_dialog_->show();
		video_dialog_->raise();
		video_dialog_->activateWindow();
	}

	void MainWindow::updateRecordDialogStatus()
	{
		if (!video_dialog_ || !video_dialog_->isVisible()) return;
		const bool rec = recorder_ && recorder_->recording();
		video_dialog_->setRecordingStatus(rec, rec ? recorder_->framesWritten() : 0);
	}

	// Auto-record hook: when a seabed/drop morphodynamic run begins (spawnWorker with an active bed), start
	// recording to the dialog's path if the user armed "Auto-record". A no-op without a bed, without the arm,
	// or (with a warning) without a path. Pulls the latest cadence/quality from the dialog first.
	void MainWindow::maybeAutoRecord()
	{
		if (!auto_record_ || !seabed_) return;
		if (recorder_ && recorder_->recording()) return;
		QString path = video_dialog_ ? video_dialog_->path().trimmed() : QString();
		if (path.isEmpty())
		{
			std::fprintf(stderr, "[video] auto-record armed but no output path set; skipping\n");
			statusBar()->showMessage("auto-record armed but no file path set (open Record Video… to set one)", 5000);
			return;
		}
		if (!path.endsWith(".mp4", Qt::CaseInsensitive)) path += ".mp4";
		if (video_dialog_)
		{
			frame_step_interval_ = std::max<long long>(1, video_dialog_->cadenceSteps());
			frame_bed_eps_ = video_dialog_->bedEps();
			rec_crf_ = video_dialog_->crf();
			rec_preset_ = video_dialog_->preset();
			rec_fps_ = video_dialog_->fps();
		}
		std::fprintf(stderr, "[video] auto-record: seabed run started -> %s\n", path.toUtf8().constData());
		startRecording(path);
	}

	void MainWindow::maybeCaptureFrame()
	{
		if (!recorder_ || !recorder_->recording() || !viewer_ || !worker_) return;
		// Capture at most ONE frame per frame_step_interval_ sim steps of evolution — a controllable cadence,
		// NOT one per step. Keying off the bed's max per-cell change alone failed: near active scour a cell
		// hits the |Δz_b| limiter (~2.5 mm/step), tripping any small threshold EVERY step, so ~one frame per
		// step (thousands ⇒ slow, duplicate-looking). The step cadence bounds the frame count regardless of
		// how fast the bed moves; the tiny bed-change check then skips dead periods (spin-up / equilibrium /
		// paused) so a static bed writes nothing.
		const long long step = worker_->steps();
		if (step - last_capture_step_ < frame_step_interval_) return;
		std::vector<float> zb;
		if (!worker_->copyBed(zb)) return; // no bed yet ⇒ nothing to record
		if (last_capture_zb_.size() == zb.size())
		{
			double dmax = 0.0;
			for (std::size_t n = 0; n < zb.size(); ++n)
				dmax = std::max(dmax, (double)std::fabs(zb[n] - last_capture_zb_[n]));
			if (dmax < frame_bed_eps_) return; // bed hasn't visibly moved since the last frame
		}
		last_capture_step_ = step;
		last_capture_zb_ = zb;
		recorder_->writeFrame(viewer_->grabFramebuffer());
	}

	void MainWindow::finalize()
	{
		if (recorder_) recorder_->finish(); // close the MP4 (idempotent) however the app exits
		shutdownWorker();
	}

	void MainWindow::closeEvent(QCloseEvent* e)
	{
		if (recorder_) recorder_->finish();
		shutdownWorker();
		QMainWindow::closeEvent(e);
	}

	void MainWindow::shutdownWorker()
	{
		if (worker_down_) return;
		worker_down_ = true;
		if (worker_) worker_->stop();
		if (worker_thread_)
		{
			worker_thread_->quit();
			worker_thread_->wait();
		}
		// Thread joined — sampling the core from this thread is now race-free.
		if (worker_ && model_injected_)
		{
			worker_->computeDiversion(diversion_);
			diversion_.inlet_U = recipe_.info.U;
		}
		delete worker_;
		worker_ = nullptr;
		// worker_thread_ is parented to `this`; Qt deletes it. Null the viewer's ref.
		if (viewer_) viewer_->setWorker(nullptr);
	}
}
