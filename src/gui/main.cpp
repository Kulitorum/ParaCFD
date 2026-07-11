// main.cpp — scour-gui entry point (milestone G1).
//
// Requests a GL 4.3 core-profile default surface format BEFORE QApplication (so the
// QOpenGLWidget's context and Qt's shared context match — cobod-slicer main.cpp:474
// pattern, bumped 3.3 -> 4.3 for clean CUDA-GL interop), builds a live M2 open-channel
// simulation from a JSON config, and shows the slice viewer.
//
// Flags (positional config path is also accepted):
//   --config <path>        viewer/sim JSON (default: configs/g1_viewer.json if present)
//   --autoclose-ms <N>     quit after N ms with exit 0 (scripted gate run)
//   --offscreen            request the Qt offscreen platform (best-effort; windowed is
//                          the reliable path for a real GL 4.3 context on Windows)
//   --load-step <path>     load a STEP model at startup (also available via File menu). The model
//                          is auto-voxelized into the flow as the solid obstacle (replacing the
//                          config obstacle); smoke-testable headlessly with --offscreen
//   --voxelize             DEPRECATED no-op: loading a model auto-injects it now (kept for scripts)
//   --noslip               use SOLID_NOSLIP for the loaded model (default: SOLID_FREESLIP,
//                          RESEARCH §3 production default)
//   --apply-h <h>          after startup, rebuild+reset the sim at voxel/cell size h (headless
//                          smoke of the dock's Apply; the domain is kept unless --apply-domain given)
//   --apply-domain Lx Ly Lz  pair with --apply-h to also change the domain on the Apply
//   --set-u <U>            after startup, set the inlet current speed U [m/s] LIVE (no reset) —
//                          headless smoke of the "Input speed" control
//   --add-sand <m>         after startup, fill the lower <m> metres with sand (around any solid) and
//                          continue as a sediment run PRESERVING the developed flow — headless smoke
//                          of the "Add sand" control
//   --model-place dx,dy,dz[,rz[,scale]]  after --load-step, move the model by (dx,dy,dz) m, rotate rz
//                          deg about +Z, scale by a uniform factor, then re-voxelize it as the obstacle —
//                          headless smoke of the placement gizmo + "voxelize where placed"
//   --morfac <M>           after startup, set MORFAC (morphological acceleration) live — the bed evolves
//                          M× faster than the flow clock. Headless smoke of the "MORFAC" control (a no-op
//                          without an active seabed run; pair with --scenario or --add-sand)
//   --display-interval <s> throttle graphics updates to every <s> wall-seconds (the "Fast sim" combo):
//                          the sim keeps stepping at full rate; the view + per-step snapshot refresh
//                          less often. 0 = live (default). Headless smoke of the Fast-sim control.
//   --tidal Umax,plateau,ramp  after startup, enable the live tidal-reversal driver (RESEARCH §8): a
//                          reversing current peaking at Umax [m/s], holding `plateau` s each direction
//                          with a `ramp` s cosine slack. Headless smoke of the "Tidal reversal" controls.
//   --clip <mode>,<frac>[,flip]  enable the clip plane at startup (headless smoke of the "Clip
//                          plane" controls): mode = x|y|z|camera (or 0..3), frac ∈ [0,1] the plane
//                          position; append ,flip to hide the other side. Hides solids only.
//   --tracers [density]    enable the grid-seeded streakline tracers at startup (headless smoke of the
//                          "Show tracers" controls): long lines that warp with the flow, coloured by
//                          speed. Optional density = seeds along the longest axis (default 10).
//   --bed-colour <src>     seabed surface colour source: elevation (default) | rate (net exchange =
//                          deposition − pickup); headless smoke of the "Bed colour" selector
//   --drop-step <path>     STEP protection unit (e.g. XStone_Decomposed.stp) for the drop/settle phase
//   --drop <count>,<seed>[,smin,smax,sand]  drop/settle PREPARATION PHASE (PLAN G3.2): scatter `count`
//                          scaled copies of --drop-step above a `sand`-m bed, settle them under gravity
//                          (Jolt), then start the morphodynamic run with the settled pile as structure.
//                          scale ∈ [smin,smax] (default 0.6..1.0), sand default 0.5 m. Give a generous
//                          --autoclose-ms so the settle animation + commit complete (needs a domain
//                          sized for the units, e.g. g1_viewer_full's 10x10x5 m)
//   --record <path.mp4>    record an MP4 of the run: a frame is captured whenever the sand bed changes
//                          and H.264-encoded via ffmpeg (QProcess pipe, no libav linkage). Finalized on
//                          close. Pair with --drop / a seabed run so there is an evolving bed to capture.
//                          Needs ffmpeg on PATH or via the SCOUR_FFMPEG env override.
//   --autosave <n>         auto-save a checkpoint (<scene>.<step>.scn) every n steps (0 = off)
//   --load-scene <path>    at startup, restore a saved .scn (setup + all data) — resumes from its step
//   --save-scene <path>    after startup, save the full scene to <path> (.scn) — headless save smoke
#include "core/cuda_probe.h"
#include "gui/main_window.h"
#include "gui/slice_viewer.h"
#include "gui/sim_setup.h"
#include "gui/video_recorder.h"

#include <QApplication>
#include <QColor>
#include <QElapsedTimer>
#include <QImage>
#include <QMessageBox>
#include <QSurfaceFormat>
#include <QTimer>

#include <nlohmann/json.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

using namespace scour::gui;

int main(int argc, char** argv)
{
	// --- Parse our flags (leave the rest for QApplication, e.g. -platform) ------
	std::string config;
	std::string step_path;
	std::string scenario_path;
	int autoclose_ms = 0;
	bool offscreen = false;
	bool voxelize = false;
	bool noslip = false;
	bool apply_grid = false;
	double apply_h = 0.0, apply_lx = 0.0, apply_ly = 0.0, apply_lz = 0.0;
	double set_u = 0.0;
	double add_sand = 0.0;
	double morfac = 0.0;
	double display_interval = 0.0;
	std::string tidal_arg;
	std::string clip_arg;
	bool tracers_on = false;
	int tracers_density = 0;
	std::string bed_colour_arg;
	std::string save_scene, load_scene;
	std::string model_place_arg;
	std::string drop_step_path, drop_arg;
	std::string record_path, record_selftest;
	int autosave = 0;
	for (int i = 1; i < argc; ++i)
	{
		std::string a = argv[i];
		if (a == "--config" && i + 1 < argc) config = argv[++i];
		else if (a == "--scenario" && i + 1 < argc) scenario_path = argv[++i];
		else if (a == "--autoclose-ms" && i + 1 < argc) autoclose_ms = std::atoi(argv[++i]);
		else if (a == "--load-step" && i + 1 < argc) step_path = argv[++i];
		else if (a == "--offscreen") offscreen = true;
		else if (a == "--voxelize") voxelize = true;
		else if (a == "--noslip") noslip = true;
		else if (a == "--apply-h" && i + 1 < argc) { apply_h = std::atof(argv[++i]); apply_grid = true; }
		else if (a == "--apply-domain" && i + 3 < argc) { apply_lx = std::atof(argv[++i]); apply_ly = std::atof(argv[++i]); apply_lz = std::atof(argv[++i]); apply_grid = true; }
		else if (a == "--set-u" && i + 1 < argc) set_u = std::atof(argv[++i]);
		else if (a == "--add-sand" && i + 1 < argc) add_sand = std::atof(argv[++i]);
		else if (a == "--model-place" && i + 1 < argc) model_place_arg = argv[++i];
		else if (a == "--morfac" && i + 1 < argc) morfac = std::atof(argv[++i]);
		else if (a == "--display-interval" && i + 1 < argc) display_interval = std::atof(argv[++i]);
		else if (a == "--tidal" && i + 1 < argc) tidal_arg = argv[++i];
		else if (a == "--clip" && i + 1 < argc) clip_arg = argv[++i];
		else if (a == "--tracers") { tracers_on = true; if (i + 1 < argc && argv[i + 1][0] != '-') tracers_density = std::atoi(argv[++i]); }
		else if (a == "--bed-colour" && i + 1 < argc) bed_colour_arg = argv[++i];
		else if (a == "--autosave" && i + 1 < argc) autosave = std::atoi(argv[++i]);
		else if (a == "--load-scene" && i + 1 < argc) load_scene = argv[++i];
		else if (a == "--save-scene" && i + 1 < argc) save_scene = argv[++i];
		else if (a == "--drop-step" && i + 1 < argc) drop_step_path = argv[++i];
		else if (a == "--drop" && i + 1 < argc) drop_arg = argv[++i];
		else if (a == "--record" && i + 1 < argc) record_path = argv[++i];
		else if (a == "--record-selftest" && i + 1 < argc) record_selftest = argv[++i];
		else if (!a.empty() && a[0] != '-') config = a; // positional config path
	}
	// A positional/--config file with a "sand_depth" key is a seabed scenario (unless --scenario set).
	if (scenario_path.empty() && !config.empty())
	{
		try
		{
			std::ifstream in(config); std::ostringstream ss; ss << in.rdbuf();
			auto j = nlohmann::json::parse(ss.str());
			if (j.contains("sand_depth") || j.contains("structure_step")) scenario_path = config;
		}
		catch (...) {}
	}
	if (config.empty()) config = "configs/g1_viewer.json";
	if (offscreen) qputenv("QT_QPA_PLATFORM", "offscreen");

	// --- GL 4.3 core default format BEFORE QApplication -------------------------
	{
		QSurfaceFormat fmt;
		fmt.setDepthBufferSize(24);
		fmt.setStencilBufferSize(0);
		fmt.setSamples(4);
		fmt.setSwapBehavior(QSurfaceFormat::DoubleBuffer);
		fmt.setProfile(QSurfaceFormat::CoreProfile);
		fmt.setVersion(4, 3);
		QSurfaceFormat::setDefaultFormat(fmt);
	}

	QApplication app(argc, argv);
	QApplication::setApplicationName("ScourProtection G1");

	// Headless smoke of the video pipe (offscreen has no real GL, so grabFramebuffer yields nothing —
	// this feeds SYNTHETIC frames straight through VideoRecorder to prove the ffmpeg pipe + MP4
	// finalization end-to-end, independent of the GL capture). Exits without launching the sim.
	if (!record_selftest.empty())
	{
		std::fprintf(stderr, "[G1] --record-selftest: piping 30 synthetic frames to ffmpeg -> '%s'\n", record_selftest.c_str());
		scour::gui::VideoRecorder rec;
		if (!rec.start(QString::fromStdString(record_selftest))) return 3;
		for (int f = 0; f < 30; ++f)
		{
			QImage img(640, 480, QImage::Format_RGBA8888);
			img.fill(QColor((f * 8) % 256, (f * 4) % 256, 128));
			rec.writeFrame(img);
			rec.flush(); // no event loop here to drain the pipe (the live GUI's event loop does that)
		}
		rec.finish();
		return 0;
	}

	// --- CUDA sanity (the whole sim runs on the GPU; it allocates device memory before the window
	// shows). Done AFTER QApplication so a fatal result can be shown as a GUI dialog: this app is
	// launched from a desktop/Start-menu icon with no console, so the old stderr-only message was
	// invisible (the console just flashed and closed). We do NOT bundle the NVIDIA driver — it is a
	// large, GPU-specific, separately-licensed package — so the actionable fix is to update it. -----
	{
		std::string err, gpu;
		int ccM = 0, ccm = 0;
		// 0 = usable, 1 = no device / driver too old, 2 = device present but too old to RUN our kernels.
		const int st = scour::core::cuda_probe_usable(&gpu, &ccM, &ccm, &err);
		if (st != 0)
		{
			std::fprintf(stderr, "[G1] GPU unusable (%d): %s (%s compute %d.%d)\n",
				st, err.c_str(), gpu.empty() ? "no device" : gpu.c_str(), ccM, ccm);
			QString title, msg;
			if (st == 2)
			{
				// A device is present but its compute capability predates this build's compiled archs
				// (sm_75/86/89). Updating the driver will NOT help — the GPU itself is too old.
				title = QStringLiteral("ScourProtection - GPU not supported");
				msg = QStringLiteral(
					"<b>ScourProtection cannot run on this GPU.</b><br><br>"
					"Your GPU (<b>%1</b>, CUDA compute capability %2.%3) is older than this build supports. "
					"It needs an NVIDIA GPU with <b>compute capability 7.5 (Turing) or newer</b> — a GeForce "
					"GTX&nbsp;16 / RTX&nbsp;20-series or later, or an equivalent Quadro / RTX&nbsp;A card. "
					"Updating the driver will not help; the GPU itself is too old.<br><br>"
					"<span style=\"color:gray;\">Technical detail: %4</span>")
					.arg(QString::fromStdString(gpu)).arg(ccM).arg(ccm)
					.arg(QString::fromStdString(err).toHtmlEscaped());
			}
			else
			{
				// No device, or the NVIDIA driver is too old for the CUDA 13.1 runtime this build links.
				const bool old_driver = err.find("insufficient") != std::string::npos
					|| err.find("driver version") != std::string::npos;
				// Rich text (HTML) so the NVIDIA driver page is a real clickable hyperlink; QMessageBox's
				// label has openExternalLinks set, so a click opens the default browser.
				const QString link = QStringLiteral(
					"<a href=\"https://www.nvidia.com/Download/index.aspx\">www.nvidia.com/Download</a>");
				title = QStringLiteral("ScourProtection - GPU required");
				msg = QStringLiteral(
					"<b>ScourProtection could not start.</b><br><br>"
					"It requires an NVIDIA (CUDA-capable) GPU with an up-to-date driver.<br><br>");
				if (old_driver)
					msg += QStringLiteral(
						"Your NVIDIA graphics driver is too old for this build (CUDA 13.1). Please update it to "
						"the latest Game Ready / Studio driver, then relaunch:<br><br>&nbsp;&nbsp;&nbsp;&nbsp;%1<br><br>").arg(link);
				else
					msg += QStringLiteral(
						"No usable NVIDIA GPU was detected. This application requires a CUDA-capable NVIDIA GPU "
						"and its driver; integrated-only graphics are not supported.<br><br>"
						"If you do have an NVIDIA GPU, install or update its driver from:<br><br>&nbsp;&nbsp;&nbsp;&nbsp;%1<br><br>").arg(link);
				msg += QStringLiteral("<span style=\"color:gray;\">Technical detail: %1</span>")
					.arg(QString::fromStdString(err).toHtmlEscaped());
			}
			if (!offscreen) // a headless/gate run has no display for a modal dialog
			{
				QMessageBox box(QMessageBox::Critical, title, msg, QMessageBox::Ok, nullptr);
				box.setTextFormat(Qt::RichText);
				box.setTextInteractionFlags(Qt::TextBrowserInteraction); // clickable link + selectable text
				box.exec();
			}
			return 2;
		}
		std::fprintf(stderr, "[G1] CUDA device 0: %s (compute %d.%d)\n", gpu.c_str(), ccM, ccm);
	}

	// --- Build the simulation ---------------------------------------------------
	SimRecipe recipe;
	std::string warn;
	SeabedScenario scen;
	std::unique_ptr<scour::core::ChannelFluidCore> core;
	if (!scenario_path.empty())
		core = build_seabed_sim(scenario_path, recipe, scen, warn);
	else
		core = build_sim(config, recipe, warn);
	const SimInfo& info = recipe.info;
	if (!warn.empty()) std::fprintf(stderr, "[G1] %s\n", warn.c_str());
	std::fprintf(stderr, "[G1] sim '%s': %dx%dx%d cells, h=%.3f m, U=%.3f m/s, nu=%.3e, %s\n",
		info.name.c_str(), info.nx, info.ny, info.nz, info.h, info.U, info.nu,
		scen.active ? "SEABED (morphodynamic)" : (info.cylinder ? "cylinder" : "empty channel"));

	MainWindow win(std::move(core), recipe, std::move(scen));
	win.show();

	// Optional auto-save cadence (applies to this and any rebuilt/restored worker).
	if (autosave > 0)
	{
		std::fprintf(stderr, "[G1] --autosave: checkpoint every %d steps\n", autosave);
		win.setAutosave(autosave);
	}

	// Optional restore of a saved scene at startup (tears down the config-built sim, resumes from the
	// saved step). Supersedes the config for the initial state.
	if (!load_scene.empty())
	{
		std::fprintf(stderr, "[G1] --load-scene: restoring '%s'\n", load_scene.c_str());
		if (!win.loadSceneFromPath(QString::fromStdString(load_scene)))
			std::fprintf(stderr, "[G1] scene restore failed\n");
	}

	// Optional drop/settle PREPARATION PHASE (PLAN G3.2): scatter + settle N STEP units above the sand
	// bed, then start the morphodynamic run with the settled pile as the rigid structure. The settle
	// animates during the event loop, so give a generous --autoclose-ms for the headless smoke.
#ifdef SCOUR_HAVE_JOLT
	if (!drop_arg.empty() && !drop_step_path.empty())
	{
		double c = 12, s = 1234, smin = 0.6, smax = 1.0, sand = 0.5;
		std::sscanf(drop_arg.c_str(), "%lf,%lf,%lf,%lf,%lf", &c, &s, &smin, &smax, &sand);
		const int count = (int)c;
		const unsigned seed = (unsigned)s;
		std::fprintf(stderr, "[G1] --drop: %d units of '%s' seed=%u scale=[%.2f,%.2f] sand=%.2f m\n",
			count, drop_step_path.c_str(), seed, smin, smax, sand);
		win.dropBlocks(QString::fromStdString(drop_step_path), count, seed, smin, smax, sand);
	}
	else if (!drop_arg.empty() || !drop_step_path.empty())
		std::fprintf(stderr, "[G1] --drop needs BOTH --drop-step <path> and --drop <count>,<seed>[,smin,smax,sand]\n");
#else
	if (!drop_arg.empty() || !drop_step_path.empty())
		std::fprintf(stderr, "[G1] --drop: built without Jolt (SCOUR_HAVE_JOLT off); drop/settle unavailable\n");
#endif

	// Optional MP4 recording of the run (frames captured on each sand-bed change; ffmpeg spawns on the
	// first captured frame). Armed before the event loop — pair with --drop / a seabed run to record.
	if (!record_path.empty())
	{
		std::fprintf(stderr, "[G1] --record: recording to '%s'\n", record_path.c_str());
		win.startRecording(QString::fromStdString(record_path));
	}

	// Optional startup STEP model (File menu does the same at runtime). The model is auto-voxelized
	// into the flow as the obstacle on load; --voxelize is retained as a no-op for script compat.
	if (!step_path.empty())
	{
		(void)voxelize; // auto-injected on load now; the flag no longer gates injection
		std::fprintf(stderr, "[G1] --load-step: loading + auto-injecting model as %s obstacle\n",
			noslip ? "no-slip" : "free-slip");
		if (!win.loadStepFile(QString::fromStdString(step_path), noslip))
			std::fprintf(stderr, "[G1] STEP load failed; no obstacle injected\n");
	}

	// Optional headless smoke of the placement gizmo: move/rotate/scale the loaded model, then re-voxelize
	// it as the obstacle (what a gizmo drag + Apply do). Applied after --load-step, before --add-sand and
	// --apply-grid so `--load-step … --model-place dx,dy,dz,rz,scale --apply-h …` exercises placed-voxelize.
	if (!model_place_arg.empty())
	{
		double dx = 0, dy = 0, dz = 0, rz = 0, sc = 0;
		std::sscanf(model_place_arg.c_str(), "%lf,%lf,%lf,%lf,%lf", &dx, &dy, &dz, &rz, &sc);
		std::fprintf(stderr, "[G1] --model-place: nudge t=(%.3f %.3f %.3f) rz=%.1f scale=%.3f, then re-voxelize\n",
			dx, dy, dz, rz, sc);
		win.nudgeModelPlacement(dx, dy, dz, rz, sc);
	}

	// Optional headless smoke of the "Add sand" control: convert to a sediment run with <m> m of sand.
	// Applied BEFORE --apply-grid so `--add-sand X --apply-domain …` exercises the converted-seabed resize
	// (a runtime seabed conversion re-voxelizing the model as the structure + re-applying the sand).
	if (add_sand > 0.0)
	{
		std::fprintf(stderr, "[G1] --add-sand: adding %.3f m of sand (flow preserved, sediment run continues)\n", add_sand);
		win.addSandMeters(add_sand);
	}

	// Optional headless smoke of the "MORFAC" control: accelerate the bed live (applied after --add-sand so
	// the engine exists; a no-op for a fluid-only viewer).
	if (morfac > 0.0)
	{
		std::fprintf(stderr, "[G1] --morfac: setting MORFAC to %.2fx (live; no-op without a bed)\n", morfac);
		win.setMorfacValue(morfac);
	}

	// Optional headless smoke of the dock's Apply: rebuild+reset at a new grid (domain + voxel size).
	if (apply_grid)
	{
		std::fprintf(stderr, "[G1] --apply-grid: rebuilding+resetting at Lx=%.3f Ly=%.3f Lz=%.3f h=%.4f (0 = keep)\n",
			apply_lx, apply_ly, apply_lz, apply_h);
		win.applyGridValues(apply_lx, apply_ly, apply_lz, apply_h);
	}

	// Optional headless smoke of the "Input speed" control: set the current live (no reset).
	if (set_u > 0.0)
	{
		std::fprintf(stderr, "[G1] --set-u: setting inlet current to %.3f m/s (live)\n", set_u);
		win.setInputSpeed(set_u);
	}

	// Optional headless smoke of the "Tidal reversal" controls: enable the reversing-current driver.
	if (!tidal_arg.empty())
	{
		double um = 0.8, pl = 30.0, rp = 20.0;
		std::sscanf(tidal_arg.c_str(), "%lf,%lf,%lf", &um, &pl, &rp);
		std::fprintf(stderr, "[G1] --tidal: reversing current U_max=%.3f plateau=%.1fs ramp=%.1fs (live)\n", um, pl, rp);
		win.setTidalReversal(true, um, pl, rp);
	}

	// Optional headless smoke of the "Clip plane" controls: enable a clip plane at the given mode/position.
	if (!clip_arg.empty() && win.viewer())
	{
		char mode[16] = { 0 };
		double frac = 0.5;
		char flip[16] = { 0 };
		std::sscanf(clip_arg.c_str(), "%15[^,],%lf,%15s", mode, &frac, flip);
		std::string ms = mode;
		int m = 2; // default Z
		if (ms == "x" || ms == "X" || ms == "0") m = 0;
		else if (ms == "y" || ms == "Y" || ms == "1") m = 1;
		else if (ms == "z" || ms == "Z" || ms == "2") m = 2;
		else if (ms == "camera" || ms == "cam" || ms == "3") m = 3;
		const bool do_flip = (std::string(flip).find("flip") != std::string::npos);
		std::fprintf(stderr, "[G1] --clip: enabling clip plane mode=%d frac=%.3f flip=%d (solids only)\n", m, frac, (int)do_flip);
		win.viewer()->setClipMode(m);
		win.viewer()->setClipFraction((float)frac);
		win.viewer()->setClipFlip(do_flip);
		win.viewer()->setClipEnabled(true);
	}

	// Optional headless smoke of the "Show tracers" controls: enable the grid-seeded streaklines.
	if (tracers_on && win.viewer())
	{
		const int dens = tracers_density > 0 ? tracers_density : 10;
		std::fprintf(stderr, "[G1] --tracers: enabling grid streaklines (density=%d, 3D)\n", dens);
		win.viewer()->setTracerGridDensity(dens);
		win.viewer()->setTracerMode3D(true);
		win.viewer()->setShowTracers(true);
	}

	// Optional headless smoke of the bed-surface colour source (elevation vs net exchange rate).
	if (!bed_colour_arg.empty() && win.viewer())
	{
		int mode = (bed_colour_arg == "rate" || bed_colour_arg == "exchange" || bed_colour_arg == "1") ? 1 : 0;
		std::fprintf(stderr, "[G1] --bed-colour: bed colour source = %s\n", mode == 1 ? "exchange rate" : "elevation");
		win.viewer()->setBedColorMode(mode);
	}

	// Optional headless smoke of the "Fast sim" graphics throttle (worker snapshot + repaint cadence).
	if (display_interval > 0.0)
	{
		std::fprintf(stderr, "[G1] --display-interval: throttling graphics to every %.2f s (fast sim)\n", display_interval);
		win.setDisplayThrottle(display_interval);
	}

	// Optional headless save smoke: after the sim has run a little, save the full scene to disk. The
	// worker gathers the state on its thread and the write completes on the main thread's event loop.
	if (!save_scene.empty())
	{
		const int delay = autoclose_ms > 0 ? std::max(800, autoclose_ms / 2) : 3000;
		std::fprintf(stderr, "[G1] --save-scene: saving '%s' after %d ms\n", save_scene.c_str(), delay);
		const QString sp = QString::fromStdString(save_scene);
		QTimer::singleShot(delay, &app, [&win, sp] { win.saveSceneNow(sp); });
	}

	QElapsedTimer wall;
	wall.start();

	if (autoclose_ms > 0)
	{
		std::fprintf(stderr, "[G1] scripted run: auto-close in %d ms\n", autoclose_ms);
		QTimer::singleShot(autoclose_ms, &app, &QApplication::quit);
	}

	int rc = app.exec();

	// Stop the worker + (if a model was injected) sample the flow-diversion metrics before the
	// window is destroyed, so the assertions below can be printed for the headless smoke path.
	win.finalize();

	double secs = wall.elapsed() / 1000.0;
	long long frames = win.viewer() ? win.viewer()->framesRendered() : 0;
	double fps = (secs > 0) ? frames / secs : 0.0;
	long long sim_steps = win.lastSteps();
	std::fprintf(stderr, "[G1] rendered %lld frames in %.2f s = %.1f fps (avg); sim advanced %lld steps to t=%.3f s\n",
		frames, secs, fps, sim_steps, win.lastSimTime());

	// --- Voxel-obstacle flow-diversion assertions (headless --voxelize check) ----
	const DiversionReport& dr = win.diversion();
	if (dr.active)
	{
		const double pen_tol = 1e-3 * dr.inlet_U;
		const bool pass_pen = dr.mean_speed_solid < pen_tol;      // no penetration into the solid
		const bool pass_wake = dr.max_speed_fluid > dr.inlet_U;   // diversion accelerates the flow
		std::fprintf(stderr,
			"[G1] voxel-flow: %lld solid cells; mean|u| in solid = %.3e m/s (< %.3e = 1e-3 U) [%s]; "
			"max fluid speed = %.3f m/s (> U = %.3f) [%s]\n",
			dr.solid_cells, dr.mean_speed_solid, pen_tol, pass_pen ? "PASS" : "FAIL",
			dr.max_speed_fluid, dr.inlet_U, pass_wake ? "PASS" : "FAIL");
		std::printf("G1-voxel: %s — %lld solid cells; mean|u|_solid=%.2e (<%.2e); max|u|_fluid=%.3f (>U=%.3f)\n",
			(pass_pen && pass_wake) ? "PASS" : "FAIL",
			dr.solid_cells, dr.mean_speed_solid, pen_tol, dr.max_speed_fluid, dr.inlet_U);
	}

	if (autoclose_ms > 0)
		std::printf("G1: %s — %lld frames, %.1f fps over %.2f s; %lld sim steps\n",
			rc == 0 ? "PASS" : "FAIL", frames, fps, secs, sim_steps);
	return rc;
}
