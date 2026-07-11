// main_window.cpp — see main_window.h.
#include "gui/main_window.h"

#include "core/geometry/building.h"
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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

namespace windcfd::gui
{
	MainWindow::MainWindow(std::unique_ptr<windcfd::core::ChannelFluidCore> core, const SimRecipe& recipe,
		QWidget* parent)
		: QMainWindow(parent), recipe_(recipe)
	{
		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("WindCFD — G1 slice viewer [%1: %2x%3x%4, h=%5 m]")
			.arg(QString::fromStdString(info.name))
			.arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h));

		viewer_ = new SliceViewer(this);
		viewer_->setInfo(info);
		setCentralWidget(viewer_);

		// --- File menu: load / clear a STEP model --------------------------------
		QMenu* fileMenu = menuBar()->addMenu("&File");
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
			model_mesh_ = windcfd::core::TriMesh{};
			if (viewer_) viewer_->clearMesh();
			initPlacementHistory(); // no model ⇒ clear the undo/redo history + disable the placement gizmo
			statusBar()->showMessage("model closed", 3000);
		});

		// Load a 3D-printing CENTERLINE STEP as a CENTERLINE (walls thickened + flat roof, built into a
		// solid obstacle) — distinct from "Open STEP…" which voxelizes the mesh directly. Set the
		// wall/roof params in the "Building" dock group, then press "Build" (also auto-built on load).
		QAction* openCenterline = fileMenu->addAction("Open centerline STEP…");
		connect(openCenterline, &QAction::triggered, this, [this] {
			QString fn = QFileDialog::getOpenFileName(this, "Open centerline STEP", QString(),
				"STEP files (*.step *.stp);;All files (*)");
			if (!fn.isEmpty()) loadCenterlineFile(fn);
		});

		// --- Scene save / restore (a full self-contained .scn: setup + model + all field data) ----
		fileMenu->addSeparator();
		QAction* saveScene = fileMenu->addAction("Save Scene As…");
		saveScene->setShortcut(QKeySequence::Save);
		connect(saveScene, &QAction::triggered, this, [this] { saveSceneAs(); });
		QAction* loadScene = fileMenu->addAction("Load Scene…");
		connect(loadScene, &QAction::triggered, this, [this] { loadSceneDialog(); });

		fileMenu->addSeparator();
		// Record an MP4 of the run: opens a settings popup (path / cadence / quality / fps) with Start/Stop
		// + a live "● REC — N frames" status. One frame is captured every N sim steps (ffmpeg).
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

		// --- Edit menu: model-placement Undo / Redo (30 levels) ----------------------------------
		// The gizmo can move/rotate/scale the model by accident (and previously that couldn't be undone —
		// only restarting fixed it). These give a visible, reversible history: Ctrl+Z / Ctrl+Shift+Z (and
		// Ctrl+Y), mirrored by the small buttons in the "Model placement" dock group. They only MOVE the
		// model — the user presses Build/Apply to re-voxelize the new pose (see the tooltips).
		QMenu* editMenu = menuBar()->addMenu("&Edit");
		undo_action_ = editMenu->addAction("Undo placement");
		undo_action_->setShortcut(QKeySequence::Undo); // Ctrl+Z
		undo_action_->setEnabled(false);
		connect(undo_action_, &QAction::triggered, this, [this] { undoPlacement(); });
		redo_action_ = editMenu->addAction("Redo placement");
		redo_action_->setShortcuts({ QKeySequence(QStringLiteral("Ctrl+Shift+Z")), QKeySequence(QStringLiteral("Ctrl+Y")) });
		redo_action_->setEnabled(false);
		connect(redo_action_, &QAction::triggered, this, [this] { redoPlacement(); });

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
		// frame when recording, on the fixed step cadence (maybeCaptureFrame early-returns otherwise).
		repaint_timer_ = new QTimer(this);
		connect(repaint_timer_, &QTimer::timeout, this, [this] { viewer_->update(); maybeCaptureFrame(); updateRecordDialogStatus(); updateWindLoadReadout(); });
		repaint_timer_->start(16); // widened by setDisplayThrottle when "Fast sim" is engaged

		// Create + start the worker.
		spawnWorker(std::move(core));
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

		// "Start Simulation" master gate: the worker is HELD (steps nothing) until this is pressed, so the
		// user can prepare the site first — load a centerline, place/rotate it with the gizmo, Build the
		// solid, and set the domain — all while the flow stays frozen. Display/repaint + Build/voxelize keep
		// working while held. Once started, the live Play/Pause + Step below take over.
		QPushButton* startBtn = new QPushButton("Start Simulation");
		startBtn->setToolTip("Begin advancing the flow. Prepare the site first (load a centerline, place/rotate it, Build, set the domain), then press Start.");
		connect(startBtn, &QPushButton::clicked, this, [this] { startSimulation(); });
		start_btn_ = startBtn;
		simCol->addWidget(startBtn);

		QHBoxLayout* runRow = new QHBoxLayout;
		QPushButton* playBtn = new QPushButton("Pause");
		playBtn->setCheckable(true);
		playBtn->setChecked(true);
		playBtn->setEnabled(false); // enabled once the sim is started (holds until then)
		playBtn->setToolTip("Play/pause the simulation stepping (available once the simulation is started).");
		connect(playBtn, &QPushButton::toggled, this, [this, playBtn](bool on) {
			if (worker_) worker_->setPlaying(on);
			playBtn->setText(on ? "Pause" : "Play");
		});
		play_btn_ = playBtn; // a rebuild (scenario load / grid Apply) honours this play/pause state
		QPushButton* stepBtn = new QPushButton("Step");
		stepBtn->setEnabled(false); // enabled once the sim is started
		stepBtn->setToolTip("Advance a single step (while paused). Available once the simulation is started.");
		connect(stepBtn, &QPushButton::clicked, this, [this] { if (worker_) worker_->stepOnce(); });
		step_btn_ = stepBtn;
		runRow->addWidget(playBtn);
		runRow->addWidget(stepBtn);
		simCol->addLayout(runRow);

		// Inlet current speed U — applies LIVE to the running flow (no reset): the wake evolves
		// toward the new current so you can test different currents interactively. It is ALSO honoured
		// by Apply below (a grid rebuild keeps the chosen current). Distinct from the domain/h block,
		// which resets the sim at t=0.
		QFormLayout* speedForm = new QFormLayout;
		speedForm->setLabelAlignment(Qt::AlignLeft);
		u_spin_ = new QDoubleSpinBox;
		u_spin_->setRange(0.0, 40.0);
		u_spin_->setDecimals(3);
		u_spin_->setSingleStep(0.05);
		u_spin_->setSuffix(" m/s");
		u_spin_->setToolTip("Inlet current speed U. Applies live to the running flow (no reset) — the wake evolves toward it; also used when you Apply a new grid.");
		connect(u_spin_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this](double v) {
			if (worker_) worker_->setInletSpeed(v);
			if (viewer_) viewer_->setReferenceU(v);
		});
		speedForm->addRow("Input speed U", u_spin_);

		// Inlet velocity profile: uniform (top-hat) vs boundary-layer (log-law). The BL profile tapers the
		// near-wall inflow to ~0 at the floor, a thinner wall-bounded inflow. Applies LIVE (the flow adjusts
		// over a flow-through).
		inlet_profile_box_ = new QComboBox;
		inlet_profile_box_->addItems({ "Uniform (top-hat)", "Boundary layer (log-law)" });
		inlet_profile_box_->setToolTip("Inlet velocity profile. Boundary layer (log-law) → near-bed u≈0 at the wall: a thinner wall-bounded inflow. Applies live.");
		connect(inlet_profile_box_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int idx) {
			if (!worker_) return;
			const bool loglaw = (idx == 1);
			const double z0 = 0.2e-3 / 12.0; // nominal wall roughness
			const double datum = 0.0;
			worker_->setInletProfile(loglaw, z0, datum);
		});
		speedForm->addRow("Inlet profile", inlet_profile_box_);
		simCol->addLayout(speedForm);

		// Tidal reversal (LIVE, no reset): drive a reversing current U_d(t) that face-swaps the inlet/
		// outlet through a cosine slack ramp (RESEARCH §8). Watch the wake reorganize as the tide flips.
		// It takes over "Input speed U" while enabled (greyed below).
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

		// --- Building group (centerline STEP → thickened walls + overhanging flat roof solid) ---------
		// Load a 3D-printing CENTERLINE STEP (File ▸ Open centerline STEP…), dial in the wall/roof params
		// here, then Build: the domain is sized around the sectioned footprint (with wind clearance) and
		// rebuilt through the SAME path as Apply, and a solid building is voxelized + injected as the flow
		// obstacle. Re-runs on the STORED centerline — no reload. Spin boxes seed from BuildingParams.
		QGroupBox* buildingGroup = new QGroupBox("Building (centerline → solid)");
		QVBoxLayout* buildingCol = new QVBoxLayout(buildingGroup);
		QFormLayout* buildingForm = new QFormLayout;
		buildingForm->setLabelAlignment(Qt::AlignLeft);
		const windcfd::core::BuildingParams bdefs; // seed the controls from the API defaults

		wall_thick_spin_ = new QDoubleSpinBox;
		wall_thick_spin_->setRange(0.05, 2.0);
		wall_thick_spin_->setDecimals(2);
		wall_thick_spin_->setSingleStep(0.05);
		wall_thick_spin_->setSuffix(" m");
		wall_thick_spin_->setValue(bdefs.wall_thickness);
		wall_thick_spin_->setToolTip("Full wall thickness: the centerline is thickened ± half of this. Also sets the build voxel size (~3 cells across the wall).");
		buildingForm->addRow("Wall thickness", wall_thick_spin_);

		wall_height_spin_ = new QDoubleSpinBox;
		wall_height_spin_->setRange(1.0, 50.0);
		wall_height_spin_->setDecimals(2);
		wall_height_spin_->setSingleStep(0.5);
		wall_height_spin_->setSuffix(" m");
		wall_height_spin_->setValue(bdefs.wall_height);
		wall_height_spin_->setToolTip("Wall height above the base (z = 0). The domain height is sized to walls + roof + clearance.");
		buildingForm->addRow("Wall height", wall_height_spin_);

		corner_radius_spin_ = new QDoubleSpinBox;
		corner_radius_spin_->setRange(0.0, 5.0);
		corner_radius_spin_->setDecimals(2);
		corner_radius_spin_->setSingleStep(0.05);
		corner_radius_spin_->setSuffix(" m");
		corner_radius_spin_->setValue(bdefs.corner_radius);
		corner_radius_spin_->setToolTip("Outer corner radius (0 = sharp).");
		buildingForm->addRow("Corner radius", corner_radius_spin_);

		roof_overhang_spin_ = new QDoubleSpinBox;
		roof_overhang_spin_->setRange(0.0, 5.0);
		roof_overhang_spin_->setDecimals(2);
		roof_overhang_spin_->setSingleStep(0.05);
		roof_overhang_spin_->setSuffix(" m");
		roof_overhang_spin_->setValue(bdefs.roof_overhang);
		roof_overhang_spin_->setToolTip("How far the flat roof extends beyond the outer wall face (0 = flush).");
		buildingForm->addRow("Roof overhang", roof_overhang_spin_);

		roof_thick_spin_ = new QDoubleSpinBox;
		roof_thick_spin_->setRange(0.0, 2.0);
		roof_thick_spin_->setDecimals(2);
		roof_thick_spin_->setSingleStep(0.05);
		roof_thick_spin_->setSuffix(" m");
		roof_thick_spin_->setValue(bdefs.roof_thickness);
		roof_thick_spin_->setToolTip("Flat roof slab thickness (0 = no roof).");
		buildingForm->addRow("Roof thickness", roof_thick_spin_);

		buildingCol->addLayout(buildingForm);

		build_btn_ = new QPushButton("Build");
		build_btn_->setToolTip("Size the domain around the building (with wind clearance), rebuild the sim at that grid, then voxelize the thickened walls + overhanging flat roof and inject them as the flow obstacle. Load a centerline via File ▸ Open centerline STEP… first.");
		connect(build_btn_, &QPushButton::clicked, this, [this] { buildBuilding(); });
		buildingCol->addWidget(build_btn_);

		// Live wind-load readout: force coefficients integrated from the pressure field over the building
		// surface (worker-side, ~every 30 steps). Dimensionless; updated on the repaint tick. Blank until a
		// building exists. Monospaced so the columns line up.
		load_readout_ = new QLabel("Wind loads: (build a building)");
		load_readout_->setTextInteractionFlags(Qt::TextSelectableByMouse);
		load_readout_->setStyleSheet("font-family: Consolas, monospace; font-size: 11px; color:#bcd;");
		load_readout_->setToolTip("Wind-load coefficients from the live pressure field (dimensionless): "
			"Cd = drag (+x, along wind), Cl = lift/uplift (+z), Cs = side (+y); Cp = surface pressure-coefficient range.");
		buildingCol->addWidget(load_readout_);
		col->addWidget(buildingGroup);

		// --- Visualization group ---------------------------------------------------------------
		QGroupBox* vizGroup = new QGroupBox("Visualization");
		QFormLayout* viz = new QFormLayout(vizGroup);
		viz->setLabelAlignment(Qt::AlignLeft);

		QComboBox* fieldBox = new QComboBox;
		fieldBox->addItems({ "Speed |u|", "u (x-vel)", "v (y-vel)", "w (z-vel)", "Pressure" });
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
		slice_chk_ = sliceChk;
		QCheckBox* modelChk = new QCheckBox("Show model");
		modelChk->setChecked(true);
		modelChk->setToolTip("Show/hide the smooth STEP polygon mesh.");
		connect(modelChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowModel(on); });
		viz->addRow(modelChk);
		QCheckBox* solidVoxChk = new QCheckBox("Show solid voxels");
		solidVoxChk->setChecked(true);
		solidVoxChk->setToolTip("Show/hide the solid-voxel staircase (obstacle), drawn dark-blue.");
		connect(solidVoxChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowSolidVoxels(on); });
		viz->addRow(solidVoxChk);
		QCheckBox* cpVoxChk = new QCheckBox("Colour building by Cp");
		cpVoxChk->setChecked(viewer_ ? viewer_->colourByCp() : true);
		cpVoxChk->setToolTip("Tint the building's voxel staircase by surface pressure coefficient Cp "
			"(blue = suction, red = pressure). Off = uniform dark-blue solid.");
		connect(cpVoxChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setColourByCp(on); });
		viz->addRow(cpVoxChk);
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
		// A movable plane that hides SOLIDS (STEP model + voxel solids) on the camera side so the interior
		// of a hollow structure is exposed; the flow slice + arrows are never clipped.
		QGroupBox* clipGroup = new QGroupBox("Clip plane (see inside)");
		QFormLayout* clip = new QFormLayout(clipGroup);
		clip->setLabelAlignment(Qt::AlignLeft);

		QCheckBox* clipChk = new QCheckBox("Enable clip plane");
		clipChk->setToolTip("Hide solids (model + voxels) between the camera and a movable plane, so you can see inside hollow structures. The flow slice and arrows stay visible.");
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
		// A single UNIFIED manipulator drawn on the model: 3 translate arrows, 3 rotate rings and 3 scale
		// cubes (X red / Y green / Z blue) plus a grey uniform-scale centre — all shown at once, and the
		// handle you grab picks the operation (translate/rotate/scale). The next "Apply" (Domain &
		// resolution) re-voxelizes it exactly where you placed it.
		QGroupBox* gizmoGroup = new QGroupBox("Model placement (gizmo)");
		QVBoxLayout* gizmoCol = new QVBoxLayout(gizmoGroup);
		gizmo_enable_chk_ = new QCheckBox("Enable manipulator (move + rotate + scale)");
		gizmo_enable_chk_->setToolTip("Show the transform gizmo on the model and drag its coloured handles: arrows = move along an axis, rings = rotate about an axis, cubes = scale (a grey centre cube scales uniformly) — X red, Y green, Z blue. A left-drag away from any handle still orbits the camera. 'Apply' (Domain & resolution) re-voxelizes the model where you placed it.");
		connect(gizmo_enable_chk_, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setGizmoEnabled(on); });
		gizmoCol->addWidget(gizmo_enable_chk_);

		// Undo / Redo the model placement (30 levels). An accidental gizmo move/rotate/scale is reversible
		// here (or via Ctrl+Z / Ctrl+Shift+Z) without restarting. They ONLY move the model — press Build (or
		// Apply) afterwards to re-voxelize the new pose.
		QHBoxLayout* undoRow = new QHBoxLayout;
		undo_btn_ = new QPushButton("Undo");
		undo_btn_->setToolTip("Undo the last model move/rotate/scale (Ctrl+Z). Only moves the model — press Build/Apply to re-voxelize the new pose.");
		undo_btn_->setEnabled(false);
		connect(undo_btn_, &QPushButton::clicked, this, [this] { undoPlacement(); });
		redo_btn_ = new QPushButton("Redo");
		redo_btn_->setToolTip("Redo the last undone placement change (Ctrl+Shift+Z / Ctrl+Y). Press Build/Apply to re-voxelize.");
		redo_btn_->setEnabled(false);
		connect(redo_btn_, &QPushButton::clicked, this, [this] { redoPlacement(); });
		undoRow->addWidget(undo_btn_);
		undoRow->addWidget(redo_btn_);
		gizmoCol->addLayout(undoRow);

		QPushButton* resetPlaceBtn = new QPushButton("Reset placement");
		resetPlaceBtn->setToolTip("Return the model to the default centre-on-bed placement (undoable with Ctrl+Z).");
		connect(resetPlaceBtn, &QPushButton::clicked, this, [this] { if (viewer_) viewer_->resetModelPlacement(); });
		gizmoCol->addWidget(resetPlaceBtn);
		gizmo_info_ = new QLabel("(load a STEP model)");
		gizmo_info_->setWordWrap(true);
		gizmo_info_->setStyleSheet("color:#9aa;");
		gizmoCol->addWidget(gizmo_info_);

		// Always-visible transform readout (monospaced): translation in metres, rotation decomposed as the
		// yaw about the vertical Z axis PLUS the full axis-angle, and (near-)uniform scale — so an accidental
		// gizmo rotation is immediately visible. Refreshed on every placement change + undo/redo/reset.
		placement_readout_ = new QLabel("(no model)");
		placement_readout_->setTextInteractionFlags(Qt::TextSelectableByMouse);
		placement_readout_->setStyleSheet("font-family: Consolas, monospace; font-size: 11px; color:#bcd;");
		placement_readout_->setToolTip("Live model placement: translation (m), rotation (yaw about Z + total axis-angle) and scale.");
		gizmoCol->addWidget(placement_readout_);

		gizmo_group_ = gizmoGroup;
		gizmoGroup->setEnabled(false);
		if (viewer_) connect(viewer_, &SliceViewer::modelPlacementChanged, this, [this] { onModelPlacementChanged(); });
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

	void MainWindow::spawnWorker(std::unique_ptr<windcfd::core::ChannelFluidCore> core,
		long long steps0, double t0)
	{
		// The checkpoint hand-off carries a shared_ptr across a queued connection — register it once.
		static bool s_meta = false;
		if (!s_meta) { qRegisterMetaType<windcfd::gui::CheckpointStatePtr>("windcfd::gui::CheckpointStatePtr"); s_meta = true; }

		worker_ = new SimWorker(std::move(core)); // no parent: moved to worker_thread_
		worker_->setStarted(sim_started_); // hold the fresh worker until "Start Simulation" (persists across a rebuild)
		if (play_btn_) worker_->setPlaying(play_btn_->isChecked()); // honour the current play/pause state
		worker_->primeCounters(steps0, t0);          // a scene restore resumes from the saved step (else 0)
		worker_->setAutosaveInterval(autosave_interval_); // keep auto-saving across a rebuild/restore
		worker_->setDisplayInterval(display_throttle_s_); // keep the "Fast sim" graphics throttle across a rebuild/restore
		worker_->setWindLoadRho(recipe_.pr.rho);          // dynamic-pressure density = the solver's rho (Cp/Cd consistency)
		// Rebuild factory: a re-inject rebuilds the core from an immutable recipe snapshot, entirely on
		// the worker thread (make_core is Qt-free and thread-safe).
		worker_->setRebuildFactory([recipe = recipe_](const std::vector<unsigned char>& solid, int mode)
			{ return make_core(recipe, solid, mode); });
		connect(worker_, &SimWorker::checkpointReady, this, &MainWindow::onCheckpointReady, Qt::QueuedConnection);
		worker_thread_ = new QThread(this);
		worker_->moveToThread(worker_thread_);
		connect(worker_thread_, &QThread::started, worker_, &SimWorker::run);
		connect(worker_, &SimWorker::stats, this, [this](qint64 steps, double t, double dt, double simFps) {
			last_steps_ = steps; last_sim_time_ = t;
			// Live step/time/dt + the measured sim throughput (steps/s). simFps is the SIM rate (climbs
			// under "Fast sim"), distinct from the render fps on the right.
			QString s = QString("step %1   t = %2 s   dt = %3 ms   %4 sim fps")
				.arg(steps).arg(t, 0, 'f', 3).arg(dt * 1e3, 0, 'f', 2).arg(simFps, 0, 'f', 1);
			if (worker_ && worker_->tidalOn())
			{
				const double u = worker_->tidePhase();
				const char* dir = qAbs(u) < 0.02 ? "slack" : (u > 0.0 ? "flood +x" : "ebb −x");
				s += QString("   |   tide %1 m/s (%2)").arg(u, 0, 'f', 2).arg(QString::fromUtf8(dir));
			}
			status_->setText(s);
		}, Qt::QueuedConnection);
		viewer_->setWorker(worker_);

		updateGizmoUi();
		worker_thread_->start();
	}

	void MainWindow::teardownWorkerForReload()
	{
		// Stop + join + delete the worker AND its thread so a fresh core can replace it with no
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

	void MainWindow::startSimulation()
	{
		if (sim_started_) return; // idempotent
		sim_started_ = true;
		if (worker_) worker_->setStarted(true); // release the master run gate — the worker begins stepping
		// Hand control to the live Play/Pause + Step controls, and make sure we resume playing.
		if (play_btn_)
		{
			play_btn_->setEnabled(true);
			if (!play_btn_->isChecked()) play_btn_->setChecked(true); // fires setPlaying(true)
			else if (worker_) worker_->setPlaying(true);
		}
		if (step_btn_) step_btn_->setEnabled(true);
		if (start_btn_) { start_btn_->setEnabled(false); start_btn_->setText("Simulation running"); }
		statusBar()->showMessage("simulation started", 3000);
		std::fprintf(stderr, "[G1] simulation started (run gate released)\n");
	}

	void MainWindow::setInputSpeed(double U)
	{
		if (!u_spin_ || U <= 0.0) return;
		u_spin_->setValue(U); // fires valueChanged → worker_->setInletSpeed + viewer_->setReferenceU (live, no reset)
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
		if (inlet_profile_box_) // reflect the built inlet mode (log-law vs uniform)
		{
			const QSignalBlocker bp(inlet_profile_box_);
			inlet_profile_box_->setCurrentIndex(recipe_.bc.inlet_mode == windcfd::core::INLET_LOGLAW ? 1 : 0);
		}
		updateGridReadout();
	}

	void MainWindow::updateGridReadout()
	{
		if (!grid_readout_ || !lx_spin_) return;
		int nx = 0, ny = 0, nz = 0;
		grid_dims_for(lx_spin_->value(), ly_spin_->value(), lz_spin_->value(), h_spin_->value(), nx, ny, nz);
		const long long cells = (long long)nx * ny * nz;
		// Rough device-memory gauge (fluid + snapshot double fields) — an order-of-
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
		// Enabled whenever a gizmo-editable model exists (a loaded fluid-viewer model or a restored scene,
		// both seeded into the gizmo). The "Enable manipulator" checkbox (gizmo_on_) persists in the viewer
		// across an Apply, so nothing to reset.
		const bool ok = viewer_ && viewer_->hasModelPlacement();
		if (gizmo_group_) gizmo_group_->setEnabled(ok);
		if (gizmo_enable_chk_ && viewer_ && gizmo_enable_chk_->isChecked() != viewer_->gizmoEnabled())
		{
			const QSignalBlocker b(gizmo_enable_chk_);
			gizmo_enable_chk_->setChecked(viewer_->gizmoEnabled()); // keep the box in sync with the viewer
		}

		// Undo/redo availability: a model must be loaded AND the respective stack non-empty.
		if (undo_action_) undo_action_->setEnabled(ok && !undo_.empty());
		if (redo_action_) redo_action_->setEnabled(ok && !redo_.empty());
		if (undo_btn_) undo_btn_->setEnabled(ok && !undo_.empty());
		if (redo_btn_) redo_btn_->setEnabled(ok && !redo_.empty());

		// Always-visible transform readout ("(no model)" when none is loaded).
		if (placement_readout_)
			placement_readout_->setText(ok ? formatPlacement(viewer_->modelXform()) : QStringLiteral("(no model)"));

		if (!gizmo_info_) return;
		if (!ok) { gizmo_info_->setText("(load a STEP model)"); return; }
		const SliceViewer::ModelGizmoXform x = viewer_->modelXform();
		const QVector3D e = x.rot.toEulerAngles(); // (pitch, yaw, roll) = rotation about X, Y, Z [deg]
		gizmo_info_->setText(QString("pos (%1, %2, %3) m\nrot (%4, %5, %6)°   scale (%7, %8, %9)")
			.arg(x.t.x(), 0, 'f', 2).arg(x.t.y(), 0, 'f', 2).arg(x.t.z(), 0, 'f', 2)
			.arg(e.x(), 0, 'f', 0).arg(e.y(), 0, 'f', 0).arg(e.z(), 0, 'f', 0)
			.arg(x.scale.x(), 0, 'f', 2).arg(x.scale.y(), 0, 'f', 2).arg(x.scale.z(), 0, 'f', 2));
	}

	// Decompose a gizmo transform into a readable, always-visible string. Rotation is reported two
	// complementary ways: the YAW about the vertical Z axis (how far the building has spun on its base —
	// the accidental rotation the user cares about) computed directly from the model's rotated local +X
	// axis, AND the full single axis-angle (unambiguous for any tilt). Scale collapses to one number when
	// (near-)uniform. Uses modelXform()'s rot/scale/t directly (no dependency on the affine decomposition).
	QString MainWindow::formatPlacement(const SliceViewer::ModelGizmoXform& x) const
	{
		constexpr float kRad2Deg = 57.2957795131f;

		// Yaw about Z: heading of the model's local +X axis projected into the world XY plane.
		const QVector3D fx = x.rot.rotatedVector(QVector3D(1, 0, 0));
		float yaw = std::atan2(fx.y(), fx.x()) * kRad2Deg;
		if (std::abs(yaw) < 0.05f) yaw = 0.0f; // clean up -0.0 / float noise

		// Total rotation as a single axis + angle.
		QVector3D axis; float ang = 0.0f;
		x.rot.getAxisAndAngle(&axis, &ang);
		if (ang > 180.0f) ang -= 360.0f; // report in (-180, 180]

		QString rotLine;
		if (std::abs(ang) < 0.05f)
			rotLine = QStringLiteral("  R  none");
		else
			rotLine = QString("  R  yaw(Z) %1°   axis(%2, %3, %4) %5°")
				.arg(yaw, 0, 'f', 1)
				.arg(axis.x(), 0, 'f', 2).arg(axis.y(), 0, 'f', 2).arg(axis.z(), 0, 'f', 2)
				.arg(ang, 0, 'f', 1);

		// Scale: single value when near-uniform, else per-axis.
		const float sx = x.scale.x(), sy = x.scale.y(), sz = x.scale.z();
		const float smax = std::max(sx, std::max(sy, sz));
		const bool uniform = (std::abs(sx - sy) + std::abs(sy - sz)) <= 1e-3f * std::max(1.0f, smax);
		const QString scaleLine = uniform
			? QString("  S  %1  (uniform)").arg(sx, 0, 'f', 3)
			: QString("  S  x%1  y%2  z%3").arg(sx, 0, 'f', 3).arg(sy, 0, 'f', 3).arg(sz, 0, 'f', 3);

		return QString("  T  (%1, %2, %3) m\n%4\n%5")
			.arg(x.t.x(), 0, 'f', 3).arg(x.t.y(), 0, 'f', 3).arg(x.t.z(), 0, 'f', 3)
			.arg(rotLine).arg(scaleLine);
	}

	// (Re)seat the undo baseline from the viewer's current placement and clear both stacks. Called on any
	// programmatic (re)placement that is NOT a user edit: model load/close, grid Apply, scene restore.
	void MainWindow::initPlacementHistory()
	{
		undo_.clear();
		redo_.clear();
		last_xform_ = (viewer_ && viewer_->hasModelPlacement())
			? viewer_->modelXform() : SliceViewer::ModelGizmoXform{};
		updateGizmoUi(); // refresh the readout + undo/redo enable state
	}

	// Record the last committed state onto the undo stack (capped at 30), drop the redo future, and adopt
	// the viewer's current placement as the new baseline. Shared by the drag/reset slot and the CLI nudge.
	void MainWindow::commitPlacementEdit()
	{
		if (!viewer_ || !viewer_->hasModelPlacement()) return;
		if (last_xform_.valid)
		{
			undo_.push_back(last_xform_);
			if ((int)undo_.size() > kUndoMax) undo_.pop_front();
			redo_.clear();
		}
		last_xform_ = viewer_->modelXform();
		updateGizmoUi();
	}

	// SliceViewer emits modelPlacementChanged on a gizmo drag-release, a Reset, or setModelPlacement (scene
	// restore). While we push a state back during undo/redo, restoring_placement_ suppresses recording so
	// the stacks aren't polluted (setModelXform does not currently emit the signal, but the guard is correct
	// regardless — reset/setModelPlacement DO emit). Otherwise the edit is committed to the history.
	void MainWindow::onModelPlacementChanged()
	{
		if (restoring_placement_) { updateGizmoUi(); return; } // programmatic restore — just refresh the readout
		commitPlacementEdit();
	}

	void MainWindow::undoPlacement()
	{
		if (undo_.empty() || !viewer_ || !viewer_->hasModelPlacement()) return;
		redo_.push_back(viewer_->modelXform());         // the current pose becomes the redo target
		if ((int)redo_.size() > kUndoMax) redo_.pop_front();
		const SliceViewer::ModelGizmoXform x = undo_.back();
		undo_.pop_back();
		restoring_placement_ = true;                    // re-entrancy guard (setModelXform may update())
		viewer_->setModelXform(x);
		restoring_placement_ = false;
		last_xform_ = x;                                // the restored state is the new baseline
		viewer_->update();
		updateGizmoUi();
		statusBar()->showMessage("placement undone — press Build/Apply to re-voxelize", 3000);
	}

	void MainWindow::redoPlacement()
	{
		if (redo_.empty() || !viewer_ || !viewer_->hasModelPlacement()) return;
		undo_.push_back(viewer_->modelXform());
		if ((int)undo_.size() > kUndoMax) undo_.pop_front();
		const SliceViewer::ModelGizmoXform x = redo_.back();
		redo_.pop_back();
		restoring_placement_ = true;
		viewer_->setModelXform(x);
		restoring_placement_ = false;
		last_xform_ = x;
		viewer_->update();
		updateGizmoUi();
		statusBar()->showMessage("placement redone — press Build/Apply to re-voxelize", 3000);
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
		setModelAsObstacle(true); // re-voxelize where placed
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

		// Rebuild the fluid viewer at the new grid: re-run build_sim with only the domain + h overridden —
		// every other physics parameter is preserved. A loaded STEP model is kept (re-displayed and
		// re-voxelized as the obstacle at the new h if it was injected).
		const std::string src = recipe_.source_config;
		const bool had_model = !model_mesh_.empty();
		windcfd::core::TriMesh keep_mesh = had_model ? model_mesh_ : windcfd::core::TriMesh{};

		// Preserve the user's gizmo placement (move/rotate/scale) across the rebuild so the model is
		// re-voxelized WHERE IT WAS PLACED, not re-centred. The "Enable manipulator" state persists in the viewer.
		SliceViewer::ModelGizmoXform keep_x = viewer_ ? viewer_->modelXform() : SliceViewer::ModelGizmoXform{};

		teardownWorkerForReload();

		SimRecipe recipe; std::string warn;
		std::unique_ptr<windcfd::core::ChannelFluidCore> core = build_sim(src, recipe, warn, &ov);
		if (!core)
		{
			statusBar()->showMessage(QString("apply failed: %1").arg(QString::fromStdString(warn)), 6000);
			return;
		}
		if (!warn.empty()) std::fprintf(stderr, "[G1] %s\n", warn.c_str());
		recipe_ = recipe;

		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("WindCFD — G1 slice viewer [%1: %2x%3x%4, h=%5 m]")
			.arg(QString::fromStdString(info.name)).arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h));
		if (viewer_)
		{
			viewer_->clearVoxelOverlay();      // drop stale-size staircase geometry
			viewer_->clearMesh();
			viewer_->setInfo(info);            // re-frame + re-size all grid-derived viewer geometry
		}
		model_mesh_ = keep_mesh;               // keep the loaded model (fluid re-inject)
		spawnWorker(std::move(core));

		// Fluid viewer with a loaded model: re-show it and re-voxelize it as the obstacle at the new h (it
		// IS the obstacle).
		if (had_model && viewer_)
		{
			viewer_->setMesh(windcfd::core::TriMesh(model_mesh_)); // display copy (model_mesh_ retained)
			if (keep_x.valid) viewer_->setModelXform(keep_x);    // re-apply the user's placement (setMesh reset it)
			setModelAsObstacle(true);                            // re-voxelizes at viewer_->modelPlacement()
			updateGizmoUi();
		}

		// Loaded centerline (building): re-display the placed house and RE-VOXELIZE it into the NEW
		// domain/resolution. buildBuilding no longer resizes the grid (the dependency is now this way round),
		// so an Apply that changes the domain re-runs the placed house's section→voxelize on the fresh grid.
		if (!centerline_mesh_.empty() && viewer_)
		{
			viewer_->setMesh(windcfd::core::TriMesh(centerline_mesh_)); // re-show (setMesh reset the gizmo)
			if (keep_x.valid) viewer_->setModelXform(keep_x);           // re-apply the user's placement
			updateGizmoUi();
			buildBuilding();                                            // re-voxelize the placed house into the new grid
		}

		syncGridControls(); // reflect the built grid (floor/clamp may differ from the typed value)
		statusBar()->showMessage(QString("grid applied: %1x%2x%3 @ h=%4 m — sim reset to t=0")
			.arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h), 5000);
	}

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
		if (!windcfd::gui::is_numbered_checkpoint(path.toStdString()))
		{
			pending_recent_keep_ = QFileInfo(path).absoluteFilePath();
			addRecentFile(path);
		}
		statusBar()->showMessage(QString("saving scene %1…").arg(QFileInfo(path).fileName()), 3000);
		return true;
	}

	void MainWindow::onCheckpointReady(windcfd::gui::CheckpointStatePtr state, qint64 tag)
	{
		if (!state) return;
		SceneFile sf;
		sf.def.recipe = recipe_;
		sf.def.has_mesh = !scene_mesh_.empty();
		if (sf.def.has_mesh)
		{
			sf.def.mesh = scene_mesh_;
			// Persist the LIVE gizmo placement (move/rotate/scale) when the model is gizmo-editable, so a
			// restored scene resumes at exactly where the model was placed; scene_place_ (the default
			// centre-on-bed) is only the fallback when there is no live placement.
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
			std::fprintf(stderr, "[scene] wrote %s (step %lld)\n", path.toUtf8().constData(),
				(long long)sf.state.steps);
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
		auto cps = windcfd::gui::list_checkpoints(fn.toStdString());
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
		if (!windcfd::gui::read_scene(path.toStdString(), sf, warn))
		{
			std::fprintf(stderr, "[scene] load FAILED (%s): %s\n", path.toUtf8().constData(), warn.c_str());
			statusBar()->showMessage(QString("scene load failed: %1").arg(QString::fromStdString(warn)), 6000);
			return false;
		}
		restoreScene(sf);
		scene_base_path_ = path;
		if (!windcfd::gui::is_numbered_checkpoint(path.toStdString())) addRecentFile(path);
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
		if (st.bc.inlet_mode == windcfd::core::INLET_LOGLAW) core->set_inlet_profile(true, st.bc.z0, st.bc.bed_datum);
		core->set_inlet_speed(st.bc.U_inlet);
		core->load_state_host(st.u, st.v, st.w, st.p);

		// GUI-side provenance so a later Apply rebuilds correctly.
		model_mesh_ = d.has_mesh ? d.mesh : windcfd::core::TriMesh{};
		scene_mesh_ = model_mesh_;
		scene_place_ = d.place;

		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("WindCFD — restored [%1: %2x%3x%4, h=%5 m] @ step %6")
			.arg(QString::fromStdString(info.name)).arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h).arg(st.steps));
		if (viewer_)
		{
			viewer_->clearVoxelOverlay();
			viewer_->clearMesh();
			viewer_->setInfo(info); // re-frame + re-size all grid-derived viewer geometry
		}

		spawnWorker(std::move(core), st.steps, st.sim_time);

		// Fluid viewer with a display mesh: attach it here (the obstacle is already baked into the restored
		// solid mask — no re-voxelize).
		if (d.has_mesh && viewer_)
		{
			viewer_->setMesh(windcfd::core::TriMesh(model_mesh_));
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
		windcfd::core::TriMesh mesh = windcfd::core::load_step_mesh(path.toStdString(), 0.1, &err);
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
		scene_place_ = windcfd::core::place_model_on_bed(mesh, recipe_.info.Lx, recipe_.info.Ly);
		if (viewer_) viewer_->setMesh(std::move(mesh));
		addRecentFile(path);          // remember it in the Recent Files menu (feature 1)

		// The loaded model IS the obstacle: voxelize it on the sim grid and inject it into the flow,
		// REPLACING the config obstacle (the default cylinder).
		setModelAsObstacle(true, noslip);
		updateGizmoUi(); // enable the placement gizmo for the freshly loaded model
		return true;
	}

	bool MainWindow::loadCenterlineFile(const QString& path, bool noslip)
	{
		QApplication::setOverrideCursor(Qt::WaitCursor);
		std::string err;
		windcfd::core::TriMesh mesh = windcfd::core::load_step_mesh(path.toStdString(), 0.1, &err);
		QApplication::restoreOverrideCursor();

		if (mesh.empty())
		{
			std::fprintf(stderr, "[G1] centerline STEP load FAILED (%s): %s\n", path.toUtf8().constData(), err.c_str());
			statusBar()->showMessage(QString("centerline load failed: %1").arg(QString::fromStdString(err)), 6000);
			return false;
		}

		const auto& lo = mesh.bbox_min;
		const auto& hi = mesh.bbox_max;
		std::fprintf(stderr,
			"[G1] loaded centerline STEP %s: %zu triangles, bbox=[%.4f %.4f %.4f]..[%.4f %.4f %.4f] m\n",
			path.toUtf8().constData(), mesh.triangle_count(),
			lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
		statusBar()->showMessage(QString("loaded centerline %1 (%2 triangles)")
			.arg(QFileInfo(path).fileName()).arg(mesh.triangle_count()), 5000);

		centerline_mesh_ = mesh; // keep the centerline so Build can re-run without reloading
		centerline_noslip_ = noslip;

		// A centerline defines a BUILDING obstacle, not a STEP-as-mesh obstacle: keep model_mesh_ EMPTY so
		// the Apply/Build paths use the centerline→building voxelizer (not voxelize_mesh over a stale mesh).
		model_mesh_ = windcfd::core::TriMesh{};

		// Display the centerline mesh in the viewer AND enable the placement gizmo for it, exactly like a
		// loaded STEP model: setMesh seeds a centre-on-bed default placement (in the CURRENT domain) that the
		// user can translate/rotate/scale, and Build voxelizes the TRANSFORMED model (see buildBuilding).
		scene_mesh_ = centerline_mesh_; // persist for a scene save (display mesh)
		scene_place_ = windcfd::core::place_model_on_bed(centerline_mesh_, recipe_.info.Lx, recipe_.info.Ly);
		if (viewer_) viewer_->setMesh(std::move(mesh)); // default centre-on-bed gizmo, valid before first paint
		addRecentFile(path);
		updateGizmoUi(); // enable the placement gizmo for the centerline

		buildBuilding(); // build + inject the solid once (into the current domain), right after loading
		return true;
	}

	void MainWindow::buildBuilding()
	{
		using namespace windcfd::core;
		if (centerline_mesh_.empty())
		{
			statusBar()->showMessage("no centerline loaded — File ▸ Open centerline STEP…", 5000);
			return;
		}
		if (!worker_) { statusBar()->showMessage("no active simulation to build into", 4000); return; }

		// 1) Place the centerline WHERE THE GIZMO PUT IT: the user can translate/rotate/scale the house, so
		//    voxelize the TRANSFORMED model. Read the live placement from the viewer (default centre-on-bed);
		//    placed_mesh applies it to every vertex so the section/voxelization land under the drawn mesh.
		const ModelPlacement place = (viewer_ && viewer_->hasModelPlacement())
			? viewer_->modelPlacement()
			: place_model_on_bed(centerline_mesh_, recipe_.info.Lx, recipe_.info.Ly);
		const TriMesh placed = placed_mesh(centerline_mesh_, place);

		// 2) Section the PLACED centerline to a 2D footprint at mid-height (NaN = robust mid-section). The
		//    footprint is already in domain xy coordinates (the placement positions the building), so no
		//    center_footprint is needed.
		Footprint fp = mesh_horizontal_section(placed, std::nan(""));
		if (fp.empty())
		{
			std::fprintf(stderr, "[G1] building: empty footprint (no closed loops from the centerline section)\n");
			statusBar()->showMessage("build failed: no closed loops from the centerline section", 6000);
			return;
		}

		// Building parameters from the dock (seeded from the BuildingParams defaults). Wall base on the floor.
		BuildingParams prm;
		if (wall_thick_spin_)    prm.wall_thickness = wall_thick_spin_->value();
		if (wall_height_spin_)   prm.wall_height    = wall_height_spin_->value();
		if (corner_radius_spin_) prm.corner_radius  = corner_radius_spin_->value();
		if (roof_overhang_spin_) prm.roof_overhang  = roof_overhang_spin_->value();
		if (roof_thick_spin_)    prm.roof_thickness = roof_thick_spin_->value();
		prm.base_z = 0.0;

		// 3) Voxelize the building into the CURRENT domain + resolution (the sim's live grid) — NO auto-resize.
		//    The user controls the simulation volume via the Domain size + voxel-size controls and Apply. Inject
		//    the mask as the obstacle through the SAME worker hand-off loadStepFile uses (rebuild off-thread).
		const MacGrid g = recipe_.grid;
		int solid = 0;
		std::vector<unsigned char> mask = voxelize_building(fp, prm, g, &solid);
		if (solid <= 0 || (int)mask.size() != g.p_count())
		{
			std::fprintf(stderr, "[G1] building: voxelize produced %d solid cells (mask %zu, expected %d) — "
				"is the building inside the current domain? adjust Domain size + Apply\n",
				solid, mask.size(), g.p_count());
			statusBar()->showMessage("build produced no solid cells (check domain size)", 6000);
			return;
		}

		const int mode = centerline_noslip_ ? SOLID_NOSLIP : SOLID_FREESLIP;
		worker_->requestRebuild(std::move(mask), mode); // core rebuild off-thread; overlay follows the published mask
		model_injected_ = true; // finalize() samples flow-diversion; the voxel overlay refreshes on the new mask

		std::fprintf(stderr,
			"[G1] building: current domain %.1f x %.1f x %.1f m @ h=%.3f -> %dx%dx%d = %d cells; "
			"%d solid cells (%.2f%% of domain) [%s]\n",
			g.nx * g.h, g.ny * g.h, g.nz * g.h, g.h, g.nx, g.ny, g.nz, g.p_count(), solid,
			100.0 * solid / std::max(1, g.p_count()), centerline_noslip_ ? "no-slip" : "free-slip");
		statusBar()->showMessage(QString("building built: %1 solid cells, grid %2×%3×%4 @ h=%5 m")
			.arg(solid).arg(g.nx).arg(g.ny).arg(g.nz).arg(g.h), 8000);
	}

	// Refresh the Building-group wind-load readout from the worker's latest published WindLoads. Called on
	// the ~60 Hz repaint tick (cheap: copies just the struct under the worker's loads mutex). Shows a hint
	// until a building is present and its first loads have been integrated.
	void MainWindow::updateWindLoadReadout()
	{
		if (!load_readout_) return;
		windcfd::core::WindLoads L;
		if (worker_ && worker_->latestLoads(L) && L.exposed_faces > 0)
		{
			load_readout_->setText(QString(
				"Wind loads (dimensionless)\n"
				"  Cd drag  (+x): %1\n"
				"  Cl uplift(+z): %2\n"
				"  Cs side  (+y): %3\n"
				"  Cp range: %4 … %5\n"
				"  exposed faces: %6")
				.arg(L.Cd, 7, 'f', 3).arg(L.Cl, 7, 'f', 3).arg(L.Cs, 7, 'f', 3)
				.arg(L.cp_min, 0, 'f', 2).arg(L.cp_max, 0, 'f', 2).arg(L.exposed_faces));
		}
		else
			load_readout_->setText("Wind loads: (build a building)");
	}

	void MainWindow::addRecentFile(const QString& path)
	{
		QSettings s("COBOD", "WindCFD");
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
		QSettings s("COBOD", "WindCFD");
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
			QSettings s2("COBOD", "WindCFD");
			s2.remove("recentStepFiles");
			rebuildRecentMenu();
		});
	}

	bool MainWindow::setModelAsObstacle(bool on, bool noslip)
	{
		if (!worker_) return false;
		if (on && model_mesh_.empty())
		{
			statusBar()->showMessage("no model loaded — open a STEP first", 4000);
			return false;
		}

		if (on)
		{
			using namespace windcfd::core;
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
				rec_crf_ = video_dialog_->crf();
				rec_preset_ = video_dialog_->preset();
				rec_fps_ = video_dialog_->fps();
				QString path = video_dialog_->path().trimmed();
				if (path.isEmpty()) path = "windcfd_run.mp4";
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
		}
		// Seed with the live values (keep any path already resolved by the recorder).
		const QString shownPath = (recorder_ && !recorder_->path().isEmpty()) ? recorder_->path() : QString();
		video_dialog_->setValues(shownPath, (int)frame_step_interval_, rec_crf_, rec_preset_, rec_fps_);
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

	void MainWindow::maybeCaptureFrame()
	{
		if (!recorder_ || !recorder_->recording() || !viewer_ || !worker_) return;
		// Capture at most ONE frame per frame_step_interval_ sim steps — a fixed time-lapse cadence that
		// bounds the clip length regardless of how fast the sim runs. A paused sim stops advancing the step
		// count, so it naturally writes nothing.
		const long long step = worker_->steps();
		if (step - last_capture_step_ < frame_step_interval_) return;
		last_capture_step_ = step;
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
		// Snapshot the last integrated wind loads for the headless exit-log (worker still alive here).
		if (worker_) last_loads_valid_ = worker_->latestLoads(last_loads_);
		delete worker_;
		worker_ = nullptr;
		// worker_thread_ is parented to `this`; Qt deletes it. Null the viewer's ref.
		if (viewer_) viewer_->setWorker(nullptr);
	}
}
