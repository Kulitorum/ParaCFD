// main_window.cpp — see main_window.h.
#include "gui/main_window.h"

#include "core/geometry/building.h"
#include "core/fluid/amr_eb.h"
#include "core/fluid/amr_grid.h"
#include "core/fluid/amr_pressure.h"
#include "core/fluid/external_aero_core.h"
#include "core/geometry/embedded_boundary.h"
#include "core/geometry/model_placement.h"
#include "core/geometry/step_import.h"
#include "core/geometry/triangle_bvh.h"
#include "core/geometry/voxelize.h"
#include "core/paraglider_config.h"
#include "gui/slice_viewer.h"
#include "gui/sim_worker.h"
#include "gui/paraglider_sim_worker.h"
#include "gui/video_recorder.h"
#include "gui/video_settings_dialog.h"

#include <cuda_runtime.h> // cudaMemGetInfo — size the Apply grid cap to the actual device VRAM

#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QCheckBox>
#include <QCloseEvent>
#include <QComboBox>
#include <QDir>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFile>
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

namespace paracfd::gui
{
	namespace
	{
		// Slurp a file's raw bytes (the source .stp) so a saved scene can embed the STEP itself, not
		// the derived triangulation. Returns empty on any read error (the caller degrades gracefully).
		std::vector<unsigned char> read_file_bytes(const QString& path)
		{
			QFile f(path);
			if (!f.open(QIODevice::ReadOnly)) return {};
			const QByteArray a = f.readAll();
			return std::vector<unsigned char>(a.constBegin(), a.constEnd());
		}
	}

	MainWindow::MainWindow(std::unique_ptr<paracfd::core::ChannelFluidCore> core, const SimRecipe& recipe,
		QWidget* parent)
		: QMainWindow(parent), recipe_(recipe)
	{
		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("ParaCFD — G1 slice viewer [%1: %2x%3x%4, h=%5 m]")
			.arg(QString::fromStdString(info.name))
			.arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h));

		viewer_ = new SliceViewer(this);
		viewer_->setInfo(info);
		pushGridToViewer(); // seed the grid-overlay cell-face coordinates (graded metrics / i·h)
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
			shutdownParagliderWorker();
			setModelAsObstacle(false); // remove the model obstacle, restoring the config obstacle (if any)
			model_mesh_ = paracfd::core::TriMesh{};
			if (viewer_) { viewer_->clearMesh(); viewer_->clearParagliderDebugBoxes(); }
			if (building_group_) building_group_->setVisible(true);
			if (start_btn_) { start_btn_->setEnabled(true); start_btn_->setText("Start Simulation"); }
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
		connect(repaint_timer_, &QTimer::timeout, this, [this] { updateParagliderReadout(); viewer_->update(); maybeCaptureFrame(); updateRecordDialogStatus(); updateWindLoadReadout(); updateAvgReadout(); });
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
			if (paraglider_worker_) paraglider_worker_->setPlaying(on);
			playBtn->setText(on ? "Pause" : "Play");
		});
		play_btn_ = playBtn; // a rebuild (scenario load / grid Apply) honours this play/pause state
		QPushButton* stepBtn = new QPushButton("Step");
		stepBtn->setEnabled(false); // enabled once the sim is started
		stepBtn->setToolTip("Advance a single step (while paused). Available once the simulation is started.");
		connect(stepBtn, &QPushButton::clicked, this, [this] { if (paraglider_worker_) paraglider_worker_->stepOnce(); else if (worker_) worker_->stepOnce(); });
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

		// --- Purpose-built paraglider controls -------------------------------------------------
		// These values rebuild the static block-AMR/EB hierarchy. They are kept separate from the
		// inherited channel controls while the legacy reference path still exists in this window.
		const paracfd::core::ParagliderConfig pgDefaults;
		QGroupBox* pgGroup = new QGroupBox("Paraglider CFD grid (Build resets)");
		paraglider_group_ = pgGroup;
		pgGroup->setVisible(false);
		QVBoxLayout* pgCol = new QVBoxLayout(pgGroup);
		QLabel* pgOrientationHint=new QLabel(QString::fromUtf8("WIND: +X  →   Bbox auto-aligns span/chord.\nLeading/trailing is ambiguous: flip 180° if needed."));
		pgOrientationHint->setWordWrap(true);pgOrientationHint->setStyleSheet("font-weight:600;color:#ff8b69;background:#20242a;padding:6px;");pgCol->addWidget(pgOrientationHint);
		QHBoxLayout* pgRotateRow=new QHBoxLayout;QPushButton* pgFlip=new QPushButton("Flip leading/trailing 180°");QPushButton* pgRotate=new QPushButton("Rotate +90° Z");pgFlip->setToolTip("The bbox identifies span and chord but not their sign. Flip the canopy when its leading edge points downwind; this rebuilds the CFD grid.");pgRotate->setToolTip("Manual axis correction for unusual STEP coordinates. Rotates the wing +90° about global Z and rebuilds the CFD grid.");connect(pgFlip,&QPushButton::clicked,this,[this]{nudgeModelPlacement(0,0,0,180,1);});connect(pgRotate,&QPushButton::clicked,this,[this]{nudgeModelPlacement(0,0,0,90,1);});pgRotateRow->addWidget(pgFlip);pgRotateRow->addWidget(pgRotate);pgCol->addLayout(pgRotateRow);
		QPushButton* pgBuild=new QPushButton("Build CFD Grid");pgBuild->setToolTip("Rebuild static AMR and zero-thickness embedded-boundary topology from the current STEP placement, then initialize a fresh GPU solver.");connect(pgBuild,&QPushButton::clicked,this,[this]{sim_started_=false;if(play_btn_)play_btn_->setChecked(false);buildParagliderPreviewGrid();});pgCol->addWidget(pgBuild);
		QFormLayout* pgForm = new QFormLayout;
		pgForm->setLabelAlignment(Qt::AlignLeft);
		auto pgDistance = [this](double value,double maximum,const char* tip)
		{
			QDoubleSpinBox* spin=new QDoubleSpinBox;spin->setRange(0.0,maximum);spin->setDecimals(3);spin->setSingleStep(0.25);spin->setSuffix(" m");spin->setValue(value);spin->setToolTip(tip);connect(spin,QOverload<double>::of(&QDoubleSpinBox::valueChanged),this,[this](double){updateParagliderGridReadout();});return spin;
		};
		pg_rho_spin_=new QDoubleSpinBox;pg_rho_spin_->setRange(0.1,10.0);pg_rho_spin_->setDecimals(4);pg_rho_spin_->setSingleStep(0.01);pg_rho_spin_->setValue(pgDefaults.freestream.rho);pg_rho_spin_->setSuffix(" kg/m³");
		pg_nu_spin_=new QDoubleSpinBox;pg_nu_spin_->setRange(1e-8,1e-2);pg_nu_spin_->setDecimals(8);pg_nu_spin_->setSingleStep(1e-6);pg_nu_spin_->setValue(pgDefaults.freestream.nu);pg_nu_spin_->setSuffix(" m²/s");
		pg_upstream_spin_=pgDistance(pgDefaults.domain.upstream_margin,50,"Distance from wing bbox to X-min freestream boundary.");
		pg_downstream_spin_=pgDistance(pgDefaults.domain.downstream_margin,100,"Distance from wing bbox to X-max pressure outlet.");
		pg_lateral_spin_=pgDistance(pgDefaults.domain.lateral_margin,50,"Far-field margin on both spanwise sides.");
		pg_vertical_spin_=pgDistance(pgDefaults.domain.vertical_margin,50,"Far-field margin above and below; there is no ground.");
		pg_base_h_spin_=pgDistance(pgDefaults.amr.base_cell_size,5,"Coarsest Cartesian cell size; every finer level is exactly 2:1.");pg_base_h_spin_->setRange(0.03125,5.0);pg_base_h_spin_->setSingleStep(0.03125);
		pg_levels_spin_=new QSpinBox;pg_levels_spin_->setRange(1,5);pg_levels_spin_->setValue(pgDefaults.amr.max_levels);connect(pg_levels_spin_,QOverload<int>::of(&QSpinBox::valueChanged),this,[this](int){updateParagliderGridReadout();});
		pg_wing_refine_spin_=pgDistance(pgDefaults.amr.wing_refinement_distance,20,"Refined volume surrounding the whole canopy.");
		pg_surface_refine_spin_=pgDistance(pgDefaults.amr.surface_refinement_distance,10,"Finest-brick distance from actual fabric triangles.");
		pg_wake_length_spin_=pgDistance(pgDefaults.amr.wake_length,100,"Static +X near-wake refinement length.");
		pg_wake_radius_spin_=pgDistance(pgDefaults.amr.wake_radius,50,"Static wake refinement radius around the wing bbox centre.");
		pg_cfl_spin_=new QDoubleSpinBox;pg_cfl_spin_->setRange(0.05,0.95);pg_cfl_spin_->setDecimals(2);pg_cfl_spin_->setSingleStep(0.05);pg_cfl_spin_->setValue(pgDefaults.solver.cfl);
		pg_cs_spin_=new QDoubleSpinBox;pg_cs_spin_->setRange(0.0,0.3);pg_cs_spin_->setDecimals(3);pg_cs_spin_->setSingleStep(0.01);pg_cs_spin_->setValue(pgDefaults.solver.smagorinsky_cs);
		pg_pressure_tolerance_spin_=new QDoubleSpinBox;pg_pressure_tolerance_spin_->setRange(1e-8,1e-2);pg_pressure_tolerance_spin_->setDecimals(8);pg_pressure_tolerance_spin_->setSingleStep(1e-5);pg_pressure_tolerance_spin_->setValue(pgDefaults.solver.projection_tolerance);
		pg_pressure_iterations_spin_=new QSpinBox;pg_pressure_iterations_spin_->setRange(20,5000);pg_pressure_iterations_spin_->setSingleStep(50);pg_pressure_iterations_spin_->setValue(pgDefaults.solver.projection_max_iterations);
		pg_reference_area_spin_=new QDoubleSpinBox;pg_reference_area_spin_->setRange(0,10000);pg_reference_area_spin_->setDecimals(3);pg_reference_area_spin_->setSingleStep(0.5);pg_reference_area_spin_->setSuffix(" m²");pg_reference_area_spin_->setToolTip("Set zero to withhold CL/CD/CS; ParaCFD does not invent a paraglider reference area.");
		pg_reference_length_spin_=new QDoubleSpinBox;pg_reference_length_spin_->setRange(0,1000);pg_reference_length_spin_->setDecimals(3);pg_reference_length_spin_->setSingleStep(0.1);pg_reference_length_spin_->setSuffix(" m");
		pg_surface_colour_box_=new QComboBox;
		pg_surface_colour_box_->addItems({"Visible side (Cp+ / Cp-)",QString::fromUtf8("Pressure difference ΔCp"),"Plus side Cp+","Minus side Cp-"});
		pg_surface_colour_box_->setCurrentIndex(1);
		pg_surface_colour_box_->setToolTip("Triangle winding defines the plus normal. Visible-side mode colours front faces by Cp+ and back faces by Cp-; reversing winding swaps the labels but not the physical pressure force.");
		connect(pg_surface_colour_box_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,[this](int){applyParagliderSurfaceColour();});
		pgForm->addRow("Air density rho",pg_rho_spin_);pgForm->addRow("Kinematic nu",pg_nu_spin_);
		pgForm->addRow("Upstream margin",pg_upstream_spin_);pgForm->addRow("Downstream margin",pg_downstream_spin_);pgForm->addRow("Lateral margin",pg_lateral_spin_);pgForm->addRow("Vertical margin",pg_vertical_spin_);
		pgForm->addRow("Base cell size",pg_base_h_spin_);pgForm->addRow("AMR levels",pg_levels_spin_);pgForm->addRow("Wing refine distance",pg_wing_refine_spin_);pgForm->addRow("Surface finest distance",pg_surface_refine_spin_);pgForm->addRow("Wake length",pg_wake_length_spin_);pgForm->addRow("Wake radius",pg_wake_radius_spin_);
		pgForm->addRow("CFL",pg_cfl_spin_);pgForm->addRow("Smagorinsky Cs",pg_cs_spin_);pgForm->addRow("Projection tolerance",pg_pressure_tolerance_spin_);pgForm->addRow("Pressure max iterations",pg_pressure_iterations_spin_);pgForm->addRow("Reference area",pg_reference_area_spin_);pgForm->addRow("Reference length",pg_reference_length_spin_);pgForm->addRow("Canopy colour",pg_surface_colour_box_);
		pgCol->addLayout(pgForm);
		pg_grid_readout_=new QLabel;pg_grid_readout_->setWordWrap(true);pg_grid_readout_->setStyleSheet("font-family: Consolas, monospace; font-size: 11px; color:#bcd;");pgCol->addWidget(pg_grid_readout_);
		col->addWidget(pgGroup);
		updateParagliderGridReadout();

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

		// --- Fine core (graded grid) — a UNIFORM h_fine core around the building inside a geometrically
		// graded coarse far field, so a rounded corner is resolvable without a fine grid everywhere. Enable
		// it, set h_fine / growth / margin, then Apply: the fine-core box AUTO-TRACKS the placed building
		// bbox + margin and the grid is regenerated at t=0. Off ⇒ the uniform grid above (opt-in).
		QGroupBox* fcGroup = new QGroupBox("Fine core (graded grid)");
		QVBoxLayout* fcCol = new QVBoxLayout(fcGroup);
		fine_core_chk_ = new QCheckBox("Enable graded fine core");
		fine_core_chk_->setToolTip("Refine a uniform h_fine core around the building and grade coarser to the domain edges. Apply regenerates the grid. Off ⇒ the uniform grid above.");
		fcCol->addWidget(fine_core_chk_);
		QFormLayout* fcForm = new QFormLayout;
		fcForm->setLabelAlignment(Qt::AlignLeft);
		fc_hfine_spin_ = new QDoubleSpinBox;
		fc_hfine_spin_->setRange(0.005, 1.0); fc_hfine_spin_->setDecimals(3); fc_hfine_spin_->setSingleStep(0.005); fc_hfine_spin_->setSuffix(" m"); fc_hfine_spin_->setValue(0.03);
		fc_hfine_spin_->setToolTip("Fine-core cell size h_fine. Rule: ≈ r/10 for a corner radius r (≥10 cells across it). h_fine ≥ r/4 is below the resolution floor (rounded ≈ sharp).");
		fc_growth_spin_ = new QDoubleSpinBox;
		fc_growth_spin_->setRange(1.01, 1.30); fc_growth_spin_->setDecimals(2); fc_growth_spin_->setSingleStep(0.01); fc_growth_spin_->setValue(1.15);
		fc_growth_spin_->setToolTip("Max adjacent-cell growth ratio in the graded far field (≤ 1.15 keeps near-2nd-order accuracy).");
		fc_margin_spin_ = new QDoubleSpinBox;
		fc_margin_spin_->setRange(0.0, 20.0); fc_margin_spin_->setDecimals(2); fc_margin_spin_->setSingleStep(0.25); fc_margin_spin_->setSuffix(" m"); fc_margin_spin_->setValue(1.0);
		fc_margin_spin_->setToolTip("Margin the fine core extends beyond the placed building bbox (captures the near-wake before grading coarsens). MUST enclose the building — windloads asserts it.");
		fcForm->addRow("Fine cell h_fine", fc_hfine_spin_);
		fcForm->addRow("Growth ratio", fc_growth_spin_);
		fcForm->addRow("Core margin", fc_margin_spin_);
		fcCol->addLayout(fcForm);
		connect(fine_core_chk_, &QCheckBox::toggled, this, [this](bool) { updateGridReadout(); });
		connect(fc_hfine_spin_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, [this](double) { updateGridReadout(); });
		col->addWidget(fcGroup);

		// --- Building group (centerline STEP → thickened walls + overhanging flat roof solid) ---------
		// Load a 3D-printing CENTERLINE STEP (File ▸ Open centerline STEP…), dial in the wall/roof params
		// here, then Build: the domain is sized around the sectioned footprint (with wind clearance) and
		// rebuilt through the SAME path as Apply, and a solid building is voxelized + injected as the flow
		// obstacle. Re-runs on the STORED centerline — no reload. Spin boxes seed from BuildingParams.
		QGroupBox* buildingGroup = new QGroupBox("Building (centerline → solid)");
		building_group_ = buildingGroup;
		QVBoxLayout* buildingCol = new QVBoxLayout(buildingGroup);
		QFormLayout* buildingForm = new QFormLayout;
		buildingForm->setLabelAlignment(Qt::AlignLeft);
		const paracfd::core::BuildingParams bdefs; // seed the controls from the API defaults

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

		// Roof toggle: cap the walls with a flat slab (default), or voxelize ONLY the wall/surface band —
		// an OPEN extruded surface with no top (e.g. a wing/airfoil profile). Off greys out the roof params.
		roof_chk_ = new QCheckBox("Roof (cap the walls)");
		roof_chk_->setChecked(bdefs.roof);
		roof_chk_->setToolTip("Cap the building with a flat roof slab (overhang + thickness above the walls). "
			"Uncheck to voxelize ONLY the wall/surface band — an open extruded surface, e.g. a wing/airfoil "
			"profile. Takes effect on the next Build/Apply.");
		connect(roof_chk_, &QCheckBox::toggled, this, [this](bool on) {
			if (roof_overhang_spin_) roof_overhang_spin_->setEnabled(on);
			if (roof_thick_spin_) roof_thick_spin_->setEnabled(on);
		});
		buildingForm->addRow("", roof_chk_);

		buildingCol->addLayout(buildingForm);

		build_btn_ = new QPushButton("Build");
		build_btn_->setToolTip("Voxelize the placed centerline model (thickened walls + overhanging flat roof, rounded or sharp corners) into the CURRENT domain and inject it as the flow obstacle. On a graded grid this REBUILDS the grid so the fine core re-centres on the building's current placement (resets to t=0). Set the domain size + voxel size via the Domain controls first; load a centerline via File ▸ Open centerline STEP…");
		connect(build_btn_, &QPushButton::clicked, this, [this] { buildBuilding(); });
		buildingCol->addWidget(build_btn_);

		// Solidify the building's sealed interior (the hollow room inside the walls, below the roof): a
		// flood-fill from the open boundary marks any fluid the wind can't reach as solid. ON by default —
		// the trapped interior is otherwise simulated at full per-cell cost (an ill-conditioned enclosed
		// pressure cavity) and adds spurious inner-wall load faces. Only genuinely sealed voids are filled
		// (an open/leaky structure stays fluid). Takes effect on the next Build/Apply.
		fill_interior_chk_ = new QCheckBox("Fill sealed interior");
		fill_interior_chk_->setChecked(true);
		fill_interior_chk_->setToolTip("Solidify the enclosed interior of the built solid (flood-fill any fluid the wind can't reach). "
			"Removes the trapped-fluid interior — cheaper pressure solve + cleaner surface loads. Only fills GENUINELY sealed "
			"cavities; an open or leaky structure stays fluid. Applies on the next Build/Apply.");
		buildingCol->addWidget(fill_interior_chk_);

		// Live wind-load readout: force coefficients integrated from the pressure field over the building
		// surface (worker-side, ~every 30 steps). Dimensionless; updated on the repaint tick. Blank until a
		// building exists. Monospaced so the columns line up.
		load_readout_ = new QLabel("Aerodynamic pressure loads: (open a STEP wing)");
		load_readout_->setTextInteractionFlags(Qt::TextSelectableByMouse);
		load_readout_->setWordWrap(true); // wrap instead of forcing the dock wider than the viewport
		load_readout_->setStyleSheet("font-family: Consolas, monospace; font-size: 11px; color:#1a4d7a;"); // dark blue: readable on the light dock
		load_readout_->setToolTip("Instantaneous aerodynamic pressure force from the paraglider solver. "
			"Axes are drag +X, side +Y, and lift +Z. Coefficients are shown only when a reference area is supplied; "
			"skin-friction force is not implemented yet.");

		// --- Converged time-averaged loads --------------------------------------------------------------
		// A single snapshot off an unsettled, turbulent flow is one random draw. Watch the flow spin up (the
		// instantaneous readout above stops trending), then Start averaging to accumulate mean / RMS / peak
		// loads + a time-averaged surface Cp (which the building recolours to). "Start in N h" defers the
		// start by wall-clock time so an overnight run can settle first and hand you a multi-hour average.
		QHBoxLayout* avgBtnRow = new QHBoxLayout;
		avg_start_btn_ = new QPushButton("Start averaging");
		avg_start_btn_->setToolTip("Begin accumulating converged time-averaged loads (mean / RMS / peak) and a time-averaged surface Cp from NOW. Press once the flow looks settled (the instantaneous Cd above has stopped trending). The building recolours to the averaged Cp.");
		connect(avg_start_btn_, &QPushButton::clicked, this, [this] { if (worker_) worker_->startAveraging(0.0); });
		avg_stop_btn_ = new QPushButton("Stop");
		avg_stop_btn_->setToolTip("Stop accumulating (freezes the averaged result so you can read it) or cancel a scheduled start.");
		connect(avg_stop_btn_, &QPushButton::clicked, this, [this] { if (worker_) worker_->stopAveraging(); });
		avgBtnRow->addWidget(avg_start_btn_);
		avgBtnRow->addWidget(avg_stop_btn_);
		buildingCol->addLayout(avgBtnRow);

		QHBoxLayout* avgSchedRow = new QHBoxLayout;
		avg_schedule_btn_ = new QPushButton("Start in…");
		avg_schedule_btn_->setToolTip("Deferred start for overnight runs: begin averaging after this many hours of WALL-CLOCK time. Start the sim, schedule +6 h, and wake up to a settled multi-hour average.");
		avg_delay_spin_ = new QDoubleSpinBox;
		avg_delay_spin_->setRange(0.1, 72.0);
		avg_delay_spin_->setDecimals(1);
		avg_delay_spin_->setSingleStep(0.5);
		avg_delay_spin_->setValue(6.0);
		avg_delay_spin_->setSuffix(" h");
		avg_delay_spin_->setToolTip("Wall-clock hours to wait before averaging starts (for overnight spin-up).");
		connect(avg_schedule_btn_, &QPushButton::clicked, this,
			[this] { if (worker_) worker_->startAveraging(avg_delay_spin_->value() * 3600.0); });
		avgSchedRow->addWidget(avg_schedule_btn_);
		avgSchedRow->addWidget(avg_delay_spin_);
		buildingCol->addLayout(avgSchedRow);

		avg_readout_ = new QLabel("Time-average: idle\n(press Start averaging when spun up)");
		avg_readout_->setTextInteractionFlags(Qt::TextSelectableByMouse);
		avg_readout_->setWordWrap(true); // wrap instead of forcing the dock wider than the viewport
		avg_readout_->setStyleSheet("font-family: Consolas, monospace; font-size: 11px; color:#1e6b2e;"); // dark green: readable on the light dock
		avg_readout_->setToolTip("Converged time-averaged loads over the current window: mean (the trustworthy coefficient), "
			"RMS (fluctuation), peak (min…max, design-critical), the averaged Cp range, and a convergence-drift hint "
			"(how far the running mean still moves — small ⇒ converged).");
		buildingCol->addWidget(avg_readout_);
		QGroupBox* aeroGroup = new QGroupBox("Paraglider aerodynamics");
		QVBoxLayout* aeroCol = new QVBoxLayout(aeroGroup);
		aeroCol->addWidget(load_readout_);
		col->addWidget(aeroGroup);
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
		QCheckBox* gridChk = new QCheckBox("Show grid on slice");
		gridChk->setChecked(false);
		gridChk->setToolTip("Draw the cell-boundary lines where the grid meets the slice plane. On a graded grid the lines "
			"bunch up in the fine h_fine core and spread out in the coarse far field — the mesh made visible.");
		connect(gridChk, &QCheckBox::toggled, this, [this](bool on) { if (viewer_) viewer_->setShowGrid(on); });
		viz->addRow(gridChk);

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

	void MainWindow::spawnWorker(std::unique_ptr<paracfd::core::ChannelFluidCore> core,
		long long steps0, double t0)
	{
		// The checkpoint hand-off carries a shared_ptr across a queued connection — register it once.
		static bool s_meta = false;
		if (!s_meta) { qRegisterMetaType<paracfd::gui::CheckpointStatePtr>("paracfd::gui::CheckpointStatePtr"); s_meta = true; }

		worker_ = new SimWorker(std::move(core)); // no parent: moved to worker_thread_
		worker_->setStarted(sim_started_); // hold the fresh worker until "Start Simulation" (persists across a rebuild)
		if (play_btn_) worker_->setPlaying(play_btn_->isChecked()); // honour the current play/pause state
		worker_->primeCounters(steps0, t0);          // a scene restore resumes from the saved step (else 0)
		worker_->setAutosaveInterval(autosave_interval_); // keep auto-saving across a rebuild/restore
		worker_->setDisplayInterval(display_throttle_s_); // keep the "Fast sim" graphics throttle across a rebuild/restore
		worker_->setWindLoadRho(recipe_.pr.rho);          // dynamic-pressure density = the solver's rho (Cp/Cd consistency)
		worker_->setMetrics(recipe_.metrics);             // graded fine-core metrics (null ⇒ uniform); host consumers use host_view()
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
		if (!model_mesh_.empty() && centerline_mesh_.empty())
		{
			if (!paraglider_worker_)
			{
				statusBar()->showMessage("build the paraglider CFD grid before starting", 5000);
				return;
			}
			if (sim_started_) return;
			sim_started_ = true;
			paraglider_worker_->setPlaying(true);
			if (play_btn_)
			{
				play_btn_->setEnabled(true);
				if (!play_btn_->isChecked()) play_btn_->setChecked(true);
			}
			if (step_btn_) step_btn_->setEnabled(true);
			if (start_btn_)
			{
				start_btn_->setEnabled(false);
				start_btn_->setText("Simulation running");
			}
			statusBar()->showMessage("paraglider GPU simulation started", 3000);
			return;
		}
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

	void MainWindow::startLoadAveraging(double delaySeconds)
	{
		if (worker_) worker_->startAveraging(delaySeconds);
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
		h_spin_->setValue(info.coarse_h); // the coarse voxel size (NOT h_fine on a graded grid)
		u_spin_->setValue(info.U);
		syncFineCoreControls(); // reflect the loaded graded fine-core state in the dock
		pushGridToViewer();     // keep the grid overlay's cell-face coords in sync (covers scene restore)
		if (inlet_profile_box_) // reflect the built inlet mode (log-law vs uniform)
		{
			const QSignalBlocker bp(inlet_profile_box_);
			inlet_profile_box_->setCurrentIndex(recipe_.bc.inlet_mode == paracfd::core::INLET_LOGLAW ? 1 : 0);
		}
		updateGridReadout();
	}

	// Dock → GridOverride: on Apply the GUI is authoritative for the fine core. The box AUTO-TRACKS the
	// placed building/model bbox + the margin (design lean: auto-track, overridable); with no model yet a
	// centred default is used. Off ⇒ graded=false ⇒ the uniform grid. The box MUST enclose the building —
	// windloads asserts bbox ⊆ the uniform fine core.
	void MainWindow::fillFineCoreOverride(GridOverride& ov) const
	{
		ov.set_fine_core = true;
		ov.graded = fine_core_chk_ && fine_core_chk_->isChecked();
		if (!ov.graded) return;
		paracfd::core::FineCoreSpec& fc = ov.fine_core;
		fc.Lx = ov.Lx; fc.Ly = ov.Ly; fc.Lz = ov.Lz; // complete the spec so GridMetrics::generate can use it directly
		fc.h_fine = fc_hfine_spin_ ? fc_hfine_spin_->value() : 0.03;
		fc.growth = fc_growth_spin_ ? fc_growth_spin_->value() : 1.15;
		const double margin = fc_margin_spin_ ? fc_margin_spin_->value() : 1.0;
		const paracfd::core::TriMesh& src = !centerline_mesh_.empty() ? centerline_mesh_ : model_mesh_;
		if (!src.empty())
		{
			const paracfd::core::ModelPlacement place = (viewer_ && viewer_->hasModelPlacement())
				? viewer_->modelPlacement()
				: paracfd::core::place_model_on_bed(src, ov.Lx, ov.Ly);
			const paracfd::core::TriMesh placed = paracfd::core::placed_mesh(src, place);
			// Wrap the placed bbox + margin, CLAMPED to the domain so the core never requests beyond it (a
			// core ≈ the domain size can't tile in exact h_fine cells → windloads would reject the surface).
			fc.x0 = std::max(0.0, (double)placed.bbox_min[0] - margin); fc.x1 = std::min(ov.Lx, (double)placed.bbox_max[0] + margin);
			fc.y0 = std::max(0.0, (double)placed.bbox_min[1] - margin); fc.y1 = std::min(ov.Ly, (double)placed.bbox_max[1] + margin);
			fc.z0 = 0.0; fc.z1 = std::min(ov.Lz, (double)placed.bbox_max[2] + margin); // ground → building top + margin (covers the roof)
		}
		else if (recipe_.graded && recipe_.fine_core.h_fine > 0.0)
		{
			// No model to auto-track, but a config/scene already defines the box — preserve it across Apply.
			fc.x0 = recipe_.fine_core.x0; fc.x1 = recipe_.fine_core.x1;
			fc.y0 = recipe_.fine_core.y0; fc.y1 = recipe_.fine_core.y1;
			fc.z0 = recipe_.fine_core.z0; fc.z1 = recipe_.fine_core.z1;
		}
		else
		{
			fc.x0 = 0.30 * ov.Lx; fc.x1 = 0.70 * ov.Lx; // centred default until a building is placed
			fc.y0 = 0.30 * ov.Ly; fc.y1 = 0.70 * ov.Ly;
			fc.z0 = 0.0; fc.z1 = 0.50 * ov.Lz;
		}
	}

	// Push the per-axis cumulative cell-face coordinates to the viewer's grid overlay: the graded metric
	// arrays (GridMetrics::xf/yf/zf) on a graded grid, else uniform i·h. The overlay draws the cell
	// boundaries where the grid meets the slice plane (bunched in the fine core, spread in the far field).
	void MainWindow::pushGridToViewer()
	{
		if (!viewer_) return;
		std::vector<double> xf, yf, zf;
		const SimInfo& info = recipe_.info;
		if (recipe_.graded && recipe_.metrics)
		{
			xf = recipe_.metrics->xf();
			yf = recipe_.metrics->yf();
			zf = recipe_.metrics->zf();
		}
		else
		{
			for (int i = 0; i <= info.nx; ++i) xf.push_back((double)i * info.h);
			for (int j = 0; j <= info.ny; ++j) yf.push_back((double)j * info.h);
			for (int k = 0; k <= info.nz; ++k) zf.push_back((double)k * info.h);
		}
		viewer_->setGridLines(xf, yf, zf);
	}

	// recipe_ → dock: reflect the loaded config/scene fine-core state so a subsequent Apply preserves it
	// (else pressing Apply with the checkbox off would silently drop a config/scene graded grid to uniform).
	void MainWindow::syncFineCoreControls()
	{
		if (!fine_core_chk_) return;
		const QSignalBlocker b0(fine_core_chk_), b1(fc_hfine_spin_), b2(fc_growth_spin_);
		fine_core_chk_->setChecked(recipe_.graded);
		if (recipe_.graded && recipe_.fine_core.h_fine > 0.0)
		{
			fc_hfine_spin_->setValue(recipe_.fine_core.h_fine);
			fc_growth_spin_->setValue(recipe_.fine_core.growth);
		}
	}

	// Apply grid cap sized to the actual device: ~80% of TOTAL VRAM at the ~240 B/cell gauge below. Using
	// total (not free) is deliberate — Apply/Build free the current grid BEFORE building the new one, so the
	// new grid effectively has the whole card. Cached (total VRAM is constant); kMaxCells is the fallback if
	// no CUDA device is queryable. ~240 B/cell is conservative (real steady ≈ 185 B/cell: core + MGPCG +
	// snapshots), so 80% of VRAM at 240 B/cell keeps healthy headroom for the OS / display / render buffers.
	long long MainWindow::maxCells()
	{
		if (max_cells_cache_ > 0) return max_cells_cache_;
		std::size_t freeB = 0, totalB = 0;
		if (cudaMemGetInfo(&freeB, &totalB) == cudaSuccess && totalB > 0)
		{
			vram_total_gb_ = (double)totalB / (1024.0 * 1024.0 * 1024.0);
			const double budget = 0.80 * (double)totalB; // 80% of total VRAM for the solver
			max_cells_cache_ = std::max(1000000LL, (long long)(budget / 240.0));
		}
		else
			max_cells_cache_ = kMaxCells; // no CUDA device visible → conservative fixed fallback
		return max_cells_cache_;
	}

	void MainWindow::updateGridReadout()
	{
		if (!grid_readout_ || !lx_spin_) return;
		int nx = 0, ny = 0, nz = 0;
		grid_dims_for(lx_spin_->value(), ly_spin_->value(), lz_spin_->value(), h_spin_->value(), nx, ny, nz);
		// Graded fine core: the actual cell count is set by the generator (fine core + graded far field), not
		// the coarse-uniform floor(L/h) — compute the real dims so the readout + the VRAM guard are truthful.
		bool graded = false;
		if (fine_core_chk_ && fine_core_chk_->isChecked())
		{
			GridOverride ov; ov.active = true;
			ov.Lx = lx_spin_->value(); ov.Ly = ly_spin_->value(); ov.Lz = lz_spin_->value(); ov.h = h_spin_->value();
			fillFineCoreOverride(ov);
			if (ov.graded && ov.fine_core.h_fine > 0.0 && ov.fine_core.h_fine < ov.h)
			{
				const paracfd::core::GridMetrics gm = paracfd::core::GridMetrics::generate(ov.fine_core);
				nx = gm.nx(); ny = gm.ny(); nz = gm.nz(); graded = true;
			}
		}
		const long long cells = (long long)nx * ny * nz;
		const long long cap = maxCells();
		// Rough device-memory gauge (fluid + snapshot double fields) — an order-of-magnitude hint, not an
		// allocation contract. Shown against the card's total VRAM so it's clear how much room is left.
		const double gb = cells * 240.0 / (1024.0 * 1024.0 * 1024.0);
		const bool ok = cells > 0 && cells <= cap;
		QString txt = QString("→ %1 × %2 × %3 = %4 cells%5\n  (~%6 GB")
			.arg(nx).arg(ny).arg(nz).arg(cells)
			.arg(graded ? QString("  (graded, h_fine=%1 m)").arg(fc_hfine_spin_->value(), 0, 'f', 3) : QString())
			.arg(gb, 0, 'f', 1);
		if (vram_total_gb_ > 0.0) txt += QString(" of %1 GB").arg(vram_total_gb_, 0, 'f', 0);
		txt += " VRAM)";
		if (!ok) txt += QString("\n⚠ too large (> %1 M cells ≈ 80%% of VRAM) — increase h").arg(cap / 1000000);
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
		if (!model_mesh_.empty() && centerline_mesh_.empty()) buildParagliderPreviewGrid();
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
		viewer_->setModelXform(x); // placement is consumed by the new static BVH/AMR/EB preprocessing
		buildParagliderPreviewGrid();
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
		fillFineCoreOverride(ov); // GUI drives the fine core on Apply: enable + h_fine/growth + auto-tracked box

		// Rebuild the fluid viewer at the new grid: re-run build_sim with only the domain + h overridden —
		// every other physics parameter is preserved. A loaded STEP model is kept (re-displayed and
		// re-voxelized as the obstacle at the new h if it was injected).
		const std::string src = recipe_.source_config;
		const bool had_model = !model_mesh_.empty();
		paracfd::core::TriMesh keep_mesh = had_model ? model_mesh_ : paracfd::core::TriMesh{};

		// Preserve the user's gizmo placement (move/rotate/scale) across the rebuild so the model is
		// re-voxelized WHERE IT WAS PLACED, not re-centred. The "Enable manipulator" state persists in the viewer.
		SliceViewer::ModelGizmoXform keep_x = viewer_ ? viewer_->modelXform() : SliceViewer::ModelGizmoXform{};

		teardownWorkerForReload();

		SimRecipe recipe; std::string warn;
		std::unique_ptr<paracfd::core::ChannelFluidCore> core = build_sim(src, recipe, warn, &ov);
		if (!core)
		{
			statusBar()->showMessage(QString("apply failed: %1").arg(QString::fromStdString(warn)), 6000);
			return;
		}
		if (!warn.empty()) std::fprintf(stderr, "[G1] %s\n", warn.c_str());
		recipe_ = recipe;

		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("ParaCFD — G1 slice viewer [%1: %2x%3x%4, h=%5 m]")
			.arg(QString::fromStdString(info.name)).arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h));
		if (viewer_)
		{
			viewer_->clearVoxelOverlay();      // drop stale-size staircase geometry
			viewer_->clearMesh();
			viewer_->setInfo(info);            // re-frame + re-size all grid-derived viewer geometry
			pushGridToViewer();                // refresh the grid overlay for the new (possibly graded) grid
		}
		model_mesh_ = keep_mesh;               // keep the loaded model (fluid re-inject)
		spawnWorker(std::move(core));

		// Paraglider STEP preview: re-show the smooth surface only. Never route an ordinary STEP
		// through the legacy solid/parity voxelizer; the new EB grid is rebuilt separately.
		if (had_model && viewer_)
		{
			viewer_->setMesh(paracfd::core::TriMesh(model_mesh_)); // display copy (model_mesh_ retained)
			if (keep_x.valid) viewer_->setModelXform(keep_x);    // re-apply the user's placement (setMesh reset it)
			viewer_->clearVoxelOverlay();
			model_injected_ = false;
			buildParagliderPreviewGrid();
			updateGizmoUi();
		}

		// Loaded centerline (building): re-display the placed house and RE-VOXELIZE it into the NEW
		// domain/resolution. buildBuilding no longer resizes the grid (the dependency is now this way round),
		// so an Apply that changes the domain re-runs the placed house's section→voxelize on the fresh grid.
		if (!centerline_mesh_.empty() && viewer_)
		{
			viewer_->setMesh(paracfd::core::TriMesh(centerline_mesh_)); // re-show (setMesh reset the gizmo)
			if (keep_x.valid) viewer_->setModelXform(keep_x);           // re-apply the user's placement
			updateGizmoUi();
			buildBuilding(false);                                       // voxelize into the ALREADY-regenerated grid (no re-reconstruct → no recursion)
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
		if (!paracfd::gui::is_numbered_checkpoint(path.toStdString()))
		{
			pending_recent_keep_ = QFileInfo(path).absoluteFilePath();
			addRecentFile(path);
		}
		statusBar()->showMessage(QString("saving scene %1…").arg(QFileInfo(path).fileName()), 3000);
		return true;
	}

	void MainWindow::onCheckpointReady(paracfd::gui::CheckpointStatePtr state, qint64 tag)
	{
		if (!state) return;
		SceneFile sf;
		sf.def.recipe = recipe_;
		sf.def.has_mesh = !scene_mesh_.empty();
		// Record whether the display mesh is a building CENTERLINE (so a restored scene routes it back into
		// centerline_mesh_ and Build can re-voxelize it) or a plain STEP-as-mesh obstacle. When a centerline is
		// loaded, model_mesh_ is deliberately empty and scene_mesh_ holds the centerline (see loadCenterline).
		sf.def.mesh_is_centerline = !centerline_mesh_.empty();
		// Embed the SOURCE STEP (the source of truth); write_scene then drops the derived triangulation and
		// regenerates it on load. Empty ⇒ legacy path embeds the mesh blobs instead.
		sf.def.step_data = step_data_;
		sf.def.step_name = step_source_name_;
		sf.def.step_deflection = 0.1; // matches the deflection the load paths pass to load_step_mesh
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
		auto cps = paracfd::gui::list_checkpoints(fn.toStdString());
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
		if (!paracfd::gui::read_scene(path.toStdString(), sf, warn))
		{
			std::fprintf(stderr, "[scene] load FAILED (%s): %s\n", path.toUtf8().constData(), warn.c_str());
			statusBar()->showMessage(QString("scene load failed: %1").arg(QString::fromStdString(warn)), 6000);
			return false;
		}
		restoreScene(sf);
		scene_base_path_ = path;
		if (!paracfd::gui::is_numbered_checkpoint(path.toStdString())) addRecentFile(path);
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
		if (st.bc.inlet_mode == paracfd::core::INLET_LOGLAW) core->set_inlet_profile(true, st.bc.z0, st.bc.bed_datum);
		core->set_inlet_speed(st.bc.U_inlet);
		core->load_state_host(st.u, st.v, st.w, st.p);

		// SOURCE OF TRUTH: the scene embeds the STEP — regenerate the triangulation from it here. The
		// stored mesh (d.mesh) is only a LEGACY fallback for scenes saved before STEP embedding. Keep the
		// bytes so a later re-save re-embeds them.
		step_data_ = d.step_data;
		step_source_name_ = d.step_name;
		paracfd::core::TriMesh model = d.mesh; // legacy scenes carry the triangulation directly
		if (!d.step_data.empty())
		{
			std::string err;
			paracfd::core::TriMesh regen = paracfd::core::load_step_mesh_from_memory(d.step_data, d.step_deflection, &err);
			if (!regen.empty()) model = std::move(regen);
			else std::fprintf(stderr, "[scene] STEP regenerate failed (%s) — falling back to stored mesh\n", err.c_str());
		}
		const bool have_model = !model.empty();

		// Route the model to the member each rebuild path reads: a building CENTERLINE feeds buildBuilding
		// (section → solid), a plain STEP feeds applyGrid (mesh re-voxelize) — otherwise Build reports "no
		// centerline loaded". For a centerline, model_mesh_ MUST stay empty (that is how applyGrid tells the
		// two apart); the injected surface mode comes from the restored solid mode.
		if (have_model && d.mesh_is_centerline)
		{
			centerline_mesh_ = model;
			centerline_noslip_ = (st.solid_mode == paracfd::core::SOLID_NOSLIP);
			model_mesh_ = paracfd::core::TriMesh{};
		}
		else
		{
			model_mesh_ = have_model ? model : paracfd::core::TriMesh{};
			centerline_mesh_ = paracfd::core::TriMesh{};
		}
		scene_mesh_ = have_model ? model : paracfd::core::TriMesh{};
		scene_place_ = d.place;

		const SimInfo& info = recipe_.info;
		setWindowTitle(QString("ParaCFD — restored [%1: %2x%3x%4, h=%5 m] @ step %6")
			.arg(QString::fromStdString(info.name)).arg(info.nx).arg(info.ny).arg(info.nz).arg(info.h).arg(st.steps));
		if (viewer_)
		{
			viewer_->clearVoxelOverlay();
			viewer_->clearMesh();
			viewer_->setInfo(info); // re-frame + re-size all grid-derived viewer geometry
			pushGridToViewer();     // refresh the grid overlay for the (possibly graded) grid
		}

		spawnWorker(std::move(core), st.steps, st.sim_time);

		// Fluid viewer with a display mesh: attach the regenerated model (the obstacle is already baked into
		// the restored solid mask — no re-voxelize). scene_mesh_ holds it for both the centerline and plain
		// cases (model_mesh_ is empty for a centerline).
		if (have_model && viewer_)
		{
			viewer_->setMesh(paracfd::core::TriMesh(scene_mesh_));
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
		shutdownParagliderWorker();
		shutdownWorker();
	}

	paracfd::core::ParagliderConfig MainWindow::paragliderConfigFromUi() const
	{
		paracfd::core::ParagliderConfig config;
		if (u_spin_) config.freestream.speed = u_spin_->value();
		if (pg_rho_spin_) config.freestream.rho = pg_rho_spin_->value();
		if (pg_nu_spin_) config.freestream.nu = pg_nu_spin_->value();
		if (pg_upstream_spin_) config.domain.upstream_margin = pg_upstream_spin_->value();
		if (pg_downstream_spin_) config.domain.downstream_margin = pg_downstream_spin_->value();
		if (pg_lateral_spin_) config.domain.lateral_margin = pg_lateral_spin_->value();
		if (pg_vertical_spin_) config.domain.vertical_margin = pg_vertical_spin_->value();
		if (pg_base_h_spin_) config.amr.base_cell_size = pg_base_h_spin_->value();
		if (pg_levels_spin_) config.amr.max_levels = pg_levels_spin_->value();
		if (pg_wing_refine_spin_) config.amr.wing_refinement_distance = pg_wing_refine_spin_->value();
		if (pg_surface_refine_spin_) config.amr.surface_refinement_distance = pg_surface_refine_spin_->value();
		if (pg_wake_length_spin_) config.amr.wake_length = pg_wake_length_spin_->value();
		if (pg_wake_radius_spin_) config.amr.wake_radius = pg_wake_radius_spin_->value();
		if (pg_cfl_spin_) config.solver.cfl = pg_cfl_spin_->value();
		if (pg_cs_spin_) config.solver.smagorinsky_cs = pg_cs_spin_->value();
		if (pg_pressure_tolerance_spin_) config.solver.projection_tolerance = pg_pressure_tolerance_spin_->value();
		if (pg_pressure_iterations_spin_) config.solver.projection_max_iterations = pg_pressure_iterations_spin_->value();
		if (pg_reference_area_spin_) config.reference.area = pg_reference_area_spin_->value();
		if (pg_reference_length_spin_) config.reference.length = pg_reference_length_spin_->value();
		return config;
	}

	void MainWindow::updateParagliderGridReadout()
	{
		if (!pg_grid_readout_) return;
		const paracfd::core::ParagliderConfig config = paragliderConfigFromUi();
		const int refinement = 1 << std::max(0, config.amr.max_levels - 1);
		const double finest = config.amr.base_cell_size / refinement;
		const double brick_width = config.amr.base_cell_size * config.amr.brick_size;
		pg_grid_readout_->setText(QString(
			"finest h = %1 m  (2:1 x%2)\n"
			"brick = %3³ cells, coarse width %4 m\n"
			"internal freestream axis = +X")
			.arg(finest, 0, 'g', 5).arg(refinement).arg(config.amr.brick_size)
			.arg(brick_width, 0, 'g', 5));
	}

	bool MainWindow::loadStepFile(const QString& path, bool noslip)
	{
		shutdownParagliderWorker();
		(void)noslip; // legacy CLI compatibility; zero-thickness fabric has no solid-wall mode here
		QApplication::setOverrideCursor(Qt::WaitCursor);
		std::string err;
		// A 2 mm display/preprocessing tessellation is already much finer than the current
		// 62.5 mm finest CFD cell and avoids million-triangle previews from the old 0.1 mm default.
		paracfd::core::TriMesh mesh = paracfd::core::load_step_mesh(path.toStdString(), 2.0, &err);
		QApplication::restoreOverrideCursor();

		if (mesh.empty())
		{
			std::fprintf(stderr, "[G1] STEP load FAILED (%s): %s\n", path.toUtf8().constData(), err.c_str());
			statusBar()->showMessage(QString("STEP load failed: %1").arg(QString::fromStdString(err)), 6000);
			return false;
		}
		// The paraglider core owns its own pooled AMR fields. Release the legacy channel
		// core before allocating them so a model load does not retain two full GPU solvers.
		teardownWorkerForReload();

		const auto& lo = mesh.bbox_min;
		const auto& hi = mesh.bbox_max;
		std::fprintf(stderr,
			"[G1] loaded STEP %s: %zu triangles, bbox=[%.4f %.4f %.4f]..[%.4f %.4f %.4f] m\n",
			path.toUtf8().constData(), mesh.triangle_count(),
			lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
		statusBar()->showMessage(QString("loaded paraglider surface %1 (%2 triangles) — build AMR/EB before CFD")
			.arg(QFileInfo(path).fileName()).arg(mesh.triangle_count()), 8000);

		model_mesh_ = mesh;            // keep a CPU copy for voxelization (viewer frees its own)
		if (building_group_) building_group_->setVisible(false);
		if (paraglider_group_) paraglider_group_->setVisible(true);
		// A paraglider session starts at the external-aero default, not the inherited
		// low-speed channel recipe. The user can still change it before rebuilding.
		if (u_spin_)
		{
			const QSignalBlocker blocker(u_spin_);
			u_spin_->setValue(paracfd::core::ParagliderConfig{}.freestream.speed);
		}
		if (viewer_ && u_spin_) viewer_->setReferenceU(u_spin_->value());
		scene_mesh_ = mesh;           // persist for a scene save (display mesh); placement recomputed below
		centerline_mesh_ = paracfd::core::TriMesh{}; // a plain STEP is the mesh obstacle, not a centerline
		step_data_ = read_file_bytes(path); // embed the SOURCE STEP in a saved scene (regenerates the mesh on load)
		step_source_name_ = QFileInfo(path).fileName().toStdString();
		scene_place_ = paracfd::core::place_model_on_bed(mesh, recipe_.info.Lx, recipe_.info.Ly);
		if (viewer_) viewer_->setMesh(std::move(mesh));
		// Place the face-only wing bbox at the configured aerodynamic margins. Standalone
		// STEP curve/wire entities never entered model_mesh_, so pilot lines cannot enlarge
		// this domain. Keep the gizmo transform editable rather than installing an override.
		if (viewer_ && viewer_->hasModelPlacement())
		{
			const paracfd::core::ParagliderConfig initial_config = paragliderConfigFromUi();
			const paracfd::core::DomainConfig& margins = initial_config.domain;
			const paracfd::core::AmrConfig& amr_defaults = initial_config.amr;
			// A bbox reliably identifies span versus chord for a deployed paraglider, but not
			// leading versus trailing edge. Both currently known exporters encode forward on
			// the negative chord axis, so use that convention and leave a prominent 180-degree
			// flip in the UI for files with the opposite sign.
			SliceViewer::ModelGizmoXform x=viewer_->modelXform();
			const double bbox_dx=model_mesh_.bbox_max[0]-model_mesh_.bbox_min[0];
			const double bbox_dy=model_mesh_.bbox_max[1]-model_mesh_.bbox_min[1];
			const double span_chord_ratio=std::max(bbox_dx,bbox_dy)/std::max(1e-9,std::min(bbox_dx,bbox_dy));
			if(span_chord_ratio>=1.25)
			{
				const bool span_is_x=bbox_dx>bbox_dy;
				const float yaw=span_is_x?90.0f:180.0f; // -Y or -X model-forward -> internal +X
				x.rot=QQuaternion::fromAxisAndAngle(0,0,1,yaw);
				viewer_->setModelXform(x);
				std::fprintf(stderr,"[paraglider-orientation] bbox %.3f x %.3f m: span=%c, assumed forward=-%c, yaw=%.0f deg (use 180-deg flip if LE/TE is reversed)\n",bbox_dx,bbox_dy,span_is_x?'X':'Y',span_is_x?'Y':'X',yaw);
			}
			else std::fprintf(stderr,"[paraglider-orientation] bbox XY ratio %.3f is ambiguous; keeping imported yaw (use placement controls)\n",span_chord_ratio);
			const paracfd::core::TriMesh initially_placed=paracfd::core::placed_mesh(model_mesh_,viewer_->modelPlacement());
			const double base_brick_width=amr_defaults.base_cell_size*amr_defaults.brick_size;
			const double requested_y=(initially_placed.bbox_max[1]-initially_placed.bbox_min[1])+2.0*margins.lateral_margin;
			const double requested_z=(initially_placed.bbox_max[2]-initially_placed.bbox_min[2])+2.0*margins.vertical_margin;
			const double padding_y=std::ceil(requested_y/base_brick_width)*base_brick_width-requested_y;
			const double padding_z=std::ceil(requested_z/base_brick_width)*base_brick_width-requested_z;
			x=viewer_->modelXform();
			x.t+=QVector3D(static_cast<float>(margins.upstream_margin-initially_placed.bbox_min[0]),static_cast<float>(margins.lateral_margin+0.5*padding_y-initially_placed.bbox_min[1]),static_cast<float>(margins.vertical_margin+0.5*padding_z-initially_placed.bbox_min[2]));
			viewer_->setModelXform(x);scene_place_=viewer_->modelPlacement();
		}
		addRecentFile(path);          // remember it in the Recent Files menu (feature 1)

		// Fundamental paraglider invariant: this is fluid/fabric/fluid, not a solid volume.
		// Clear any prior legacy model injection and keep the STEP as the visual/source surface.
		if (model_injected_) setModelAsObstacle(false);
		if (viewer_) viewer_->clearVoxelOverlay();
		model_injected_ = false;
		// Never let the legacy channel path appear to simulate this surface. The old core was
		// released above; Start is enabled only after the external-aero topology is ready.
		sim_started_ = false;
		if (play_btn_)
		{
			play_btn_->setChecked(false);
			play_btn_->setEnabled(false);
		}
		if (step_btn_) step_btn_->setEnabled(false);
		if (start_btn_)
		{
			start_btn_->setEnabled(false);
			start_btn_->setText("Building CFD grid...");
		}
		buildParagliderPreviewGrid();
		updateGizmoUi(); // enable the placement gizmo for the freshly loaded model
		return true;
	}

	void MainWindow::buildParagliderPreviewGrid()
	{
		using namespace paracfd::core;
		if(model_mesh_.empty()||!viewer_)return;
		shutdownParagliderWorker();
		sim_started_ = false;
		if (play_btn_) { play_btn_->setChecked(false); play_btn_->setEnabled(false); }
		if (step_btn_) step_btn_->setEnabled(false);
		if (start_btn_) { start_btn_->setEnabled(false); start_btn_->setText("Building CFD grid..."); }
		const TriMesh wing=placed_mesh(model_mesh_,viewer_->modelPlacement());
		TriangleBvh bvh(wing);
		const ParagliderConfig cfg=paragliderConfigFromUi();
		const Aabb3d requested=automatic_flow_domain(wing,cfg.domain);
		std::unique_ptr<ExternalAeroCore> external_core;
		try{external_core=std::make_unique<ExternalAeroCore>(wing,bvh,cfg);}catch(const std::exception& e){statusBar()->showMessage(QString("paraglider CFD grid failed: %1").arg(e.what()),12000);std::fprintf(stderr,"[paraglider-preview] external core failed: %s\n",e.what());if(start_btn_){start_btn_->setEnabled(false);start_btn_->setText("CFD grid failed");}return;}const AmrHierarchy& amr=external_core->hierarchy();
		std::vector<std::array<float,6>> brick_boxes,eb_boxes;std::size_t fragments=0,patches=0,apertures=0,unresolved=0,pressure_static=0;
		for(const AmrLevel& level:amr.levels())for(int brick_id=0;brick_id<(int)level.bricks.size();++brick_id)
		{
			const BrickMetadata& brick=level.bricks[brick_id];if(!brick.active())continue;const float width=amr.brick_size()*brick.h;brick_boxes.push_back({(float)brick.origin.x,(float)brick.origin.y,(float)brick.origin.z,(float)brick.origin.x+width,(float)brick.origin.y+width,(float)brick.origin.z+width});
		}
		const AmrEmbeddedBoundaryAtlas& atlas=external_core->embedded_boundary();
		for(const AmrEbLevelAtlas& level_atlas:atlas.levels)
		{
			const EmbeddedBoundary& eb=level_atlas.topology;const UniformEbGrid& grid=eb.grid;
			for(const FluidFragment& fragment:eb.fragments)if(level_atlas.owned_cell[fragment.parent_cell]){++fragments;if(fragment.pressure_static)++pressure_static;}
			for(const FaceAperture& aperture:eb.apertures)if(level_atlas.owned_cell[aperture.parent_face_cell])++apertures;
			for(const SurfacePatch& patch:eb.patches){const BrickLocation owner=amr.locate_finest(patch.centroid);if(owner.found()&&owner.level==level_atlas.level)++patches;}
			for(int cell:eb.irregular_cells)if(level_atlas.owned_cell[cell]){const auto q=grid.cell_coord(cell);const Aabb3d box=grid.cell_box(q[0],q[1],q[2]);eb_boxes.push_back({(float)box.lo.x,(float)box.lo.y,(float)box.lo.z,(float)box.hi.x,(float)box.hi.y,(float)box.hi.z});}
			for(const UnresolvedEbCell& problem:eb.unresolved)if(level_atlas.owned_cell[problem.parent_cell]){++unresolved;const auto q=grid.cell_coord(problem.parent_cell);const Aabb3d box=grid.cell_box(q[0],q[1],q[2]);eb_boxes.push_back({(float)box.lo.x,(float)box.lo.y,(float)box.lo.z,(float)box.hi.x,(float)box.hi.y,(float)box.hi.z});}
		}
		const CompositeAmrPressureSystem& pressure=external_core->pressure_system();const bool pressure_ready=true;const std::size_t pressure_dofs=pressure.storage_size,pressure_edges=pressure.coarse_fine.size()+pressure.embedded.size();
		const Vec3d padded_size=amr.domain().hi-amr.domain().lo;
		SimInfo display_info;
		display_info.h=amr.levels().front().h; display_info.coarse_h=display_info.h;
		display_info.nx=static_cast<int>(std::llround(padded_size.x/display_info.h));
		display_info.ny=static_cast<int>(std::llround(padded_size.y/display_info.h));
		display_info.nz=static_cast<int>(std::llround(padded_size.z/display_info.h));
		display_info.Lx=padded_size.x; display_info.Ly=padded_size.y; display_info.Lz=padded_size.z;
		display_info.U=cfg.freestream.speed; display_info.rho=cfg.freestream.rho; display_info.nu=cfg.freestream.nu;
		display_info.name="paraglider_amr_display";
		viewer_->setInfo(display_info); viewer_->setGridLines({}, {}, {});
		viewer_->setShowSlice(true);viewer_->setParagliderDebugBoxes(brick_boxes,eb_boxes);const Vec3d requested_size=requested.hi-requested.lo;setWindowTitle(QString("ParaCFD — paraglider geometry/AMR preview [%1 bricks, hmin=%2 m]").arg(amr.active_brick_count()).arg(amr.finest_cell_size(),0,'g',4));statusBar()->showMessage(QString("face-only domain %1 × %2 × %3 m; AMR %4 bricks; EB fragments %5, apertures %6, patches %7, unresolved %8, static pockets %9; pressure topology %10").arg(padded_size.x,0,'f',2).arg(padded_size.y,0,'f',2).arg(padded_size.z,0,'f',2).arg(amr.active_brick_count()).arg(fragments).arg(apertures).arg(patches).arg(unresolved).arg(pressure_static).arg(pressure_ready?QString("ready (%1 DOFs/%2 edges)").arg(pressure_dofs).arg(pressure_edges):QString("not ready")),15000);
		std::fprintf(stderr,"[paraglider-preview] face-only requested domain %.3f x %.3f x %.3f m; padded AMR %.3f x %.3f x %.3f m, %zu bricks; EB fragments=%zu apertures=%zu patches=%zu unresolved=%zu static=%zu pressure=%s (%zu DOFs, %zu special edges)\n",requested_size.x,requested_size.y,requested_size.z,padded_size.x,padded_size.y,padded_size.z,amr.active_brick_count(),fragments,apertures,patches,unresolved,pressure_static,pressure_ready?"ready":"not-ready",pressure_dofs,pressure_edges);
		spawnParagliderWorker(std::move(external_core));if(start_btn_){start_btn_->setEnabled(true);start_btn_->setText("Start Simulation");}
	}

	bool MainWindow::loadCenterlineFile(const QString& path, bool noslip)
	{
		shutdownParagliderWorker();
		if (building_group_) building_group_->setVisible(true);
		if (paraglider_group_) paraglider_group_->setVisible(false);
		QApplication::setOverrideCursor(Qt::WaitCursor);
		std::string err;
		paracfd::core::TriMesh mesh = paracfd::core::load_step_mesh(path.toStdString(), 0.1, &err);
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
		if(start_btn_){start_btn_->setEnabled(!sim_started_);start_btn_->setText(sim_started_?"Simulation running":"Start Simulation");}
		centerline_noslip_ = noslip;
		step_data_ = read_file_bytes(path); // embed the SOURCE STEP in a saved scene (regenerates the centerline on load)
		step_source_name_ = QFileInfo(path).fileName().toStdString();

		// A centerline defines a BUILDING obstacle, not a STEP-as-mesh obstacle: keep model_mesh_ EMPTY so
		// the Apply/Build paths use the centerline→building voxelizer (not voxelize_mesh over a stale mesh).
		model_mesh_ = paracfd::core::TriMesh{};

		// Display the centerline mesh in the viewer AND enable the placement gizmo for it, exactly like a
		// loaded STEP model: setMesh seeds a centre-on-bed default placement (in the CURRENT domain) that the
		// user can translate/rotate/scale, and Build voxelizes the TRANSFORMED model (see buildBuilding).
		scene_mesh_ = centerline_mesh_; // persist for a scene save (display mesh)
		scene_place_ = paracfd::core::place_model_on_bed(centerline_mesh_, recipe_.info.Lx, recipe_.info.Ly);
		if (viewer_) viewer_->setMesh(std::move(mesh)); // default centre-on-bed gizmo, valid before first paint
		addRecentFile(path);
		updateGizmoUi(); // enable the placement gizmo for the centerline

		buildBuilding(false); // build + inject the solid once into the CONFIGURED grid (no reconstruct on load;
		                      // the first interactive Build/Apply re-tracks a graded fine core to the placement)
		return true;
	}

	void MainWindow::buildBuilding(bool reconstruct_grid)
	{
		using namespace paracfd::core;
		// Build voxelizes the loaded model as a building. Prefer the centerline; fall back to a plain loaded
		// model (model_mesh_) so Build also works for a plain STEP and for LEGACY scenes restored before STEP
		// embedding (their model lands in model_mesh_, with no centerline flag to route it otherwise).
		const TriMesh& src_mesh = !centerline_mesh_.empty() ? centerline_mesh_ : model_mesh_;
		if (src_mesh.empty())
		{
			statusBar()->showMessage("no model loaded — File ▸ Open centerline STEP…", 5000);
			return;
		}
		if (!worker_) { statusBar()->showMessage("no active simulation to build into", 4000); return; }

		// GRADED grid: the fine core is anchored to the building, so an INTERACTIVE (re)build must RECONSTRUCT
		// the grid around the building's CURRENT placement — voxelizing into the existing core would leave the
		// fine zone at the OLD position after the model was moved (the old grid must not influence the new
		// build). Delegate to the full rebuild (applyGrid regenerates the metric grid with the fine core
		// re-tracked to the placement, resets to t=0, then calls buildBuilding(false) to voxelize into the
		// fresh grid). A uniform grid is position-independent, so it voxelizes into the current grid directly.
		// reconstruct_grid is true only for the interactive "Build" button (always post-startup, GL live); the
		// CLI/File-menu load and applyGrid's own re-voxelize pass false, so they voxelize into the built grid.
		if (reconstruct_grid && fine_core_chk_ && fine_core_chk_->isChecked())
		{
			applyGrid();
			return;
		}

		// 1) Place the model WHERE THE GIZMO PUT IT: the user can translate/rotate/scale the house, so
		//    voxelize the TRANSFORMED model. Read the live placement from the viewer (default centre-on-bed);
		//    placed_mesh applies it to every vertex so the section/voxelization land under the drawn mesh.
		const ModelPlacement place = (viewer_ && viewer_->hasModelPlacement())
			? viewer_->modelPlacement()
			: place_model_on_bed(src_mesh, recipe_.info.Lx, recipe_.info.Ly);
		const TriMesh placed = placed_mesh(src_mesh, place);

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
		if (roof_chk_)           prm.roof           = roof_chk_->isChecked(); // off => walls/surface only (wing profile)
		prm.base_z = 0.0;

		// 3) Voxelize the building into the CURRENT domain + resolution (the sim's live grid) — NO auto-resize.
		//    The user controls the simulation volume via the Domain size + voxel-size controls and Apply. Inject
		//    the mask as the obstacle through the SAME worker hand-off loadStepFile uses (rebuild off-thread).
		// HOST-view grid: on a graded grid recipe_.grid is the DEVICE view (device metric pointers); the
		// voxelizer runs on this (main) thread and must read host cell coordinates. (Graded-aware wall
		// voxelization itself is task 4.1; this only keeps the deref host-safe.)
		const MacGrid g = recipe_.metrics ? recipe_.metrics->host_view() : recipe_.grid;
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

		// Solidify the sealed interior (default on): the hollow room inside the walls, below the roof, is
		// trapped fluid — fill it so it isn't simulated as an ill-conditioned enclosed pressure cavity and
		// so the wind-load walk sees only the OUTER wall faces. Flood-fill only touches GENUINELY enclosed
		// voids; an open/leaky structure stays fluid (filled == 0). Must run BEFORE the mask is moved out.
		if (fill_interior_chk_ && fill_interior_chk_->isChecked())
		{
			const int filled = seal_enclosed_voids(mask, g);
			if (filled > 0)
			{
				solid += filled;
				std::fprintf(stderr, "[G1] building: sealed-interior fill solidified %d enclosed cells (%.2f%% of domain)\n",
					filled, 100.0 * filled / std::max(1, g.p_count()));
			}
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
		if (!load_readout_ || paraglider_worker_) return;
		paracfd::core::WindLoads L;
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

	// Refresh the time-average block: the run state (idle / overnight countdown / collecting / stopped) and,
	// once samples have been folded, the converged mean / RMS / peak coefficients, the elapsed window (in
	// seconds and in flow-through times L_x/U), the averaged Cp range and a convergence-drift hint. Called on
	// the repaint tick alongside updateWindLoadReadout; cheap (atomics + one small struct copy).
	void MainWindow::updateAvgReadout()
	{
		if (paraglider_worker_) return;
		if (!avg_readout_ || !worker_) return;
		using AvgPhase = SimWorker::AvgPhase;
		const AvgPhase ph = worker_->avgPhase();

		if (ph == AvgPhase::Pending)
		{
			const double rem = worker_->avgCountdownSeconds();
			const int h = (int)(rem / 3600.0);
			const int m = (int)((rem - 3600.0 * h) / 60.0);
			const int s = (int)(rem - 3600.0 * h - 60.0 * m);
			avg_readout_->setText(QString("Time-average: starts in %1h %2m %3s\n(overnight — flow keeps developing until then)")
				.arg(h).arg(m, 2, 10, QChar('0')).arg(s, 2, 10, QChar('0')));
			return;
		}

		paracfd::core::WindLoadStats S;
		if (worker_->latestAvgStats(S) && S.samples > 0)
		{
			const double ft = S.flow_through_time > 1e-9 ? S.duration / S.flow_through_time : 0.0;
			const QString head = (ph == AvgPhase::Collecting) ? "collecting" : "stopped";
			const QString drift = S.Cd_drift >= 0.0 ? QString("%1%").arg(S.Cd_drift * 100.0, 0, 'f', 2) : QString("—");
			// Kept to short lines (header on its own row, tight columns) so the block never forces the dock
			// wider than the viewport; word-wrap catches any overflow from unusually large transient values.
			avg_readout_->setText(QString(
				"Time-average — %1\n"
				"  %2 samples · %3 s · %4 flow-throughs\n"
				"        mean      rms    peak (min…max)\n"
				"  Cd %5 %6  %7…%8\n"
				"  Cl %9 %10  %11…%12\n"
				"  Cs %13 %14  %15…%16\n"
				"  Cp(avg): %17 … %18\n"
				"  drift (1st half vs all): %19")
				.arg(head).arg(S.samples).arg(S.duration, 0, 'f', 1).arg(ft, 0, 'f', 1)
				.arg(S.Cd_mean, 8, 'f', 3).arg(S.Cd_rms, 7, 'f', 3).arg(S.Cd_min, 0, 'f', 2).arg(S.Cd_max, 0, 'f', 2)
				.arg(S.Cl_mean, 8, 'f', 3).arg(S.Cl_rms, 7, 'f', 3).arg(S.Cl_min, 0, 'f', 2).arg(S.Cl_max, 0, 'f', 2)
				.arg(S.Cs_mean, 8, 'f', 3).arg(S.Cs_rms, 7, 'f', 3).arg(S.Cs_min, 0, 'f', 2).arg(S.Cs_max, 0, 'f', 2)
				.arg(S.cp_min, 0, 'f', 2).arg(S.cp_max, 0, 'f', 2)
				.arg(drift));
			return;
		}

		if (ph == AvgPhase::Collecting)
			avg_readout_->setText("Time-average: collecting… (waiting for the first load sample)");
		else
			avg_readout_->setText("Time-average: idle (press Start averaging when spun up)");
	}

	void MainWindow::addRecentFile(const QString& path)
	{
		QSettings s("COBOD", "ParaCFD");
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
		QSettings s("COBOD", "ParaCFD");
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
			QSettings s2("COBOD", "ParaCFD");
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
			using namespace paracfd::core;
			// Voxelize where the user placed it: use the gizmo's live placement (move/rotate/scale) when a
			// model is loaded; otherwise the default centre-on-bed. The viewer's placement and its drawn mesh
			// share one transform, so the solid mask lands exactly under the mesh.
			ModelPlacement place = (viewer_ && viewer_->hasModelPlacement())
				? viewer_->modelPlacement()
				: place_model_on_bed(model_mesh_, recipe_.info.Lx, recipe_.info.Ly);
			double mesh_vol = 0.0, voxel_vol = 0.0;
			int thin = 0;
			std::vector<float> frac;
			const MacGrid vg = recipe_.metrics ? recipe_.metrics->host_view() : recipe_.grid; // host view (device-view on graded)
			std::vector<unsigned char> model_solid = voxelize_mesh(model_mesh_, vg, place,
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

			// Solidify any sealed interior (same "Fill sealed interior" toggle as Build): a watertight solid
			// is already filled by voxelize_mesh's ray parity, so this is a no-op there; it seals a
			// non-watertight / open loaded shell that would otherwise trap fluid inside.
			if (fill_interior_chk_ && fill_interior_chk_->isChecked())
			{
				const int filled = seal_enclosed_voids(obstacle, vg);
				if (filled > 0)
				{
					nsolid += filled;
					std::fprintf(stderr, "[G1] model: sealed-interior fill solidified %d enclosed cells\n", filled);
				}
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
				if (path.isEmpty()) path = "paracfd_run.mp4";
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
		shutdownParagliderWorker();
		shutdownWorker();
	}

	void MainWindow::closeEvent(QCloseEvent* e)
	{
		if (recorder_) recorder_->finish();
		shutdownParagliderWorker();
		shutdownWorker();
		QMainWindow::closeEvent(e);
	}

	void MainWindow::spawnParagliderWorker(std::unique_ptr<paracfd::core::ExternalAeroCore> core)
	{
		shutdownParagliderWorker();
		paraglider_snapshot_generation_ = 0;
		paraglider_surface_generation_ = 0;
		paraglider_cp_plus_.clear();paraglider_cp_minus_.clear();paraglider_delta_cp_.clear();
		if(viewer_)viewer_->clearTriangleSurfaceColouring();
		paraglider_worker_ = new ParagliderSimWorker(std::move(core));
		if (viewer_) viewer_->setParagliderWorker(paraglider_worker_);
		paraglider_worker_->setPlaying(sim_started_ && play_btn_ && play_btn_->isChecked());
		paraglider_thread_ = new QThread(this);
		paraglider_worker_->moveToThread(paraglider_thread_);
		connect(paraglider_thread_, &QThread::started,
			paraglider_worker_, &ParagliderSimWorker::run);
		connect(paraglider_worker_, &ParagliderSimWorker::finished,
			paraglider_thread_, &QThread::quit);
		paraglider_thread_->start();
	}

	void MainWindow::shutdownParagliderWorker()
	{
		if (viewer_) viewer_->setParagliderWorker(nullptr);
		if (paraglider_worker_) paraglider_worker_->stop();
		if (paraglider_thread_)
		{
			paraglider_thread_->quit();
			paraglider_thread_->wait();
		}
		delete paraglider_worker_;
		paraglider_worker_ = nullptr;
		delete paraglider_thread_;
		paraglider_thread_ = nullptr;
		paraglider_snapshot_generation_ = 0;
		paraglider_surface_generation_ = 0;
	}

	void MainWindow::applyParagliderSurfaceColour()
	{
		if(!viewer_||paraglider_delta_cp_.empty())return;
		const int mode=pg_surface_colour_box_?pg_surface_colour_box_->currentIndex():1;
		if(mode==0)viewer_->setTriangleSurfaceCp(paraglider_cp_plus_,paraglider_cp_minus_,-paraglider_side_cp_range_,paraglider_side_cp_range_);
		else if(mode==2)viewer_->setTriangleSurfaceCp(paraglider_cp_plus_,paraglider_cp_plus_,-paraglider_side_cp_range_,paraglider_side_cp_range_);
		else if(mode==3)viewer_->setTriangleSurfaceCp(paraglider_cp_minus_,paraglider_cp_minus_,-paraglider_side_cp_range_,paraglider_side_cp_range_);
		else viewer_->setTriangleSurfaceCp(paraglider_delta_cp_,paraglider_delta_cp_,-paraglider_delta_cp_range_,paraglider_delta_cp_range_);
	}

	void MainWindow::updateParagliderReadout()
	{
		if (!paraglider_worker_) return;
		ParagliderDisplaySnapshot snapshot;
		if (!paraglider_worker_->latestSnapshot(paraglider_snapshot_generation_,
			paraglider_surface_generation_, snapshot)) return;
		if (!snapshot.error.empty())
		{
			if (status_) status_->setText(QString("paraglider CFD stopped: %1")
				.arg(QString::fromStdString(snapshot.error)));
			statusBar()->showMessage(QString("paraglider CFD error: %1")
				.arg(QString::fromStdString(snapshot.error)), 15000);
			if (play_btn_ && play_btn_->isChecked()) play_btn_->setChecked(false);
			return;
		}
		if (!snapshot.initialized)
		{
			if (status_) status_->setText("initializing paraglider pressure field on GPU...");
			return;
		}
		if (!snapshot.delta_cp.empty())
		{
			paraglider_cp_plus_=std::move(snapshot.cp_plus);
			paraglider_cp_minus_=std::move(snapshot.cp_minus);
			paraglider_delta_cp_=std::move(snapshot.delta_cp);
			paraglider_delta_cp_range_=std::max(std::abs(snapshot.cp_min),std::abs(snapshot.cp_max));
			paraglider_side_cp_range_=std::max(std::abs(snapshot.side_cp_min),std::abs(snapshot.side_cp_max));
			applyParagliderSurfaceColour();
		}
		if (status_)
		{
			status_->setText(QString("step %1   t = %2 s   dt = %3 ms   CFL %4   GPU %5 ms   pressure %6 ms / %7 it / r=%8%9")
				.arg(snapshot.steps).arg(snapshot.physical_time, 0, 'f', 3)
				.arg(snapshot.dt * 1e3, 0, 'f', 2).arg(snapshot.effective_cfl, 0, 'f', 2)
				.arg(snapshot.step_ms, 0, 'f', 1).arg(snapshot.projection_ms, 0, 'f', 1)
				.arg(snapshot.pressure_iterations).arg(snapshot.residual, 0, 'g', 3)
				.arg(snapshot.converged ? "" : "  NOT CONVERGED"));
		}
		last_steps_ = snapshot.steps;
		last_sim_time_ = snapshot.physical_time;
		if (load_readout_)
		{
			QString text = QString(
				"Paraglider pressure loads (instantaneous)\n"
				"  drag  Fx (+x): %1 N\n"
				"  side  Fy (+y): %2 N\n"
				"  lift  Fz (+z): %3 N\n"
				"  delta-Cp range: %4 ... %5")
				.arg(snapshot.pressure_force.x, 0, 'f', 3)
				.arg(snapshot.pressure_force.y, 0, 'f', 3)
				.arg(snapshot.pressure_force.z, 0, 'f', 3)
				.arg(snapshot.cp_min, 0, 'f', 3).arg(snapshot.cp_max, 0, 'f', 3);
			if (snapshot.coefficients_valid)
				text += QString("\n  Cd,p %1   Cs,p %2   Cl,p %3")
					.arg(snapshot.cd_pressure, 0, 'f', 4)
					.arg(snapshot.cs_pressure, 0, 'f', 4)
					.arg(snapshot.cl_pressure, 0, 'f', 4);
			else
				text += "\n  coefficients withheld: reference area not set";
			text += QString("\n  max velocity: regular %1 / compact EB %2 m/s"
				"\n  compact EB transport CFL rate: %3 1/s"
				"\n  persistent GPU estimate: %4 MiB")
				.arg(snapshot.max_abs_regular_velocity, 0, 'g', 5)
				.arg(snapshot.max_abs_special_velocity, 0, 'g', 5)
				.arg(snapshot.max_embedded_cfl_rate, 0, 'g', 5)
				.arg(snapshot.gpu_bytes / (1024.0 * 1024.0), 0, 'f', 1);
			if (snapshot.conservation_valid)
				text += QString("\n  divergence max/RMSV: %1 / %2 1/s"
					"\n  flux error max/sum/net: %3 / %4 / %5 m³/s")
					.arg(snapshot.max_abs_divergence, 0, 'g', 3)
					.arg(snapshot.volume_weighted_rms_divergence, 0, 'g', 3)
					.arg(snapshot.max_integrated_flux_error, 0, 'g', 3)
					.arg(snapshot.absolute_integrated_flux_error, 0, 'g', 3)
					.arg(snapshot.net_integrated_flux_error, 0, 'g', 3);
			if (snapshot.max_abs_special_velocity > std::max(50.0, 3.0 * snapshot.max_abs_regular_velocity))
				text += "\n  WARNING: compact EB velocity hotspot; loads are not trustworthy";
			text += "\n  pressure-only: skin friction is not implemented";
			load_readout_->setText(text);
		}
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
		if (worker_) last_avg_valid_ = worker_->latestAvgStats(last_avg_stats_); // converged loads (--average-now)
		delete worker_;
		worker_ = nullptr;
		// worker_thread_ is parented to `this`; Qt deletes it. Null the viewer's ref.
		if (viewer_) viewer_->setWorker(nullptr);
	}
}
