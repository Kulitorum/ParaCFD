// slice_viewer.h — G1 slice-plane viewer widget. A QOpenGLWidget (GL 4.3 core) that
// renders an axis-aligned slice of the live simulation: CUDA writes RGBA vertex colours
// straight into a GL-registered VBO (zero copy, slice_gl.cu) and the widget draws the
// coloured grid mesh with a slicer-style orbit/pan/zoom camera. ALL GL happens here on
// the main thread (PARACFD_ASSERT_GL_THREAD guards every entry point); the widget only
// READS the worker's post-step device snapshots, under the worker's display_mutex().
#pragma once

#include "core/geometry/model_placement.h"
#include "core/geometry/step_import.h"
#include "gui/camera.h"
#include "gui/flow_particles.h"
#include "gui/flow_tracers.h"
#include "gui/sim_setup.h"
#include "gui/slice_field.h"

#include <QElapsedTimer>
#include <QMatrix4x4>
#include <QOpenGLFunctions_4_3_Core>
#include <QOpenGLShaderProgram>
#include <QOpenGLWidget>
#include <QPoint>
#include <QPointF>
#include <QQuaternion>
#include <QVector3D>
#include <QVector4D>

#include <cstdint>

class QPainter;

namespace paracfd::gui
{
	class SimWorker;

	class SliceViewer : public QOpenGLWidget, protected QOpenGLFunctions_4_3_Core
	{
		Q_OBJECT
	public:
		explicit SliceViewer(QWidget* parent = nullptr);
		~SliceViewer() override;

		// Wire the data source + domain metadata (call before the first paint). Re-pushes the
		// arrow host-flow request to the (possibly new) worker.
		void setWorker(SimWorker* w);
		void setInfo(const SimInfo& info);

		// Live-update ONLY the reference speed used for the fixed per-field colour range + the arrow
		// speed-scale fallback (both matter when auto-range is off). Follows a live "Input speed"
		// change without re-framing the camera or rebuilding geometry (unlike setInfo). Main thread.
		void setReferenceU(double U);

		// Live view controls (main thread).
		void setField(Field f) { field_ = f; updateRange(); range_valid_ = false; update(); }
		void setAxis(Axis a) { axis_ = a; setDefaultPlane(); geometry_dirty_ = true; grid_dirty_ = true; arrows_.reset(); tracers_.reset(); update(); }
		void setPlaneFraction(float frac); // 0..1 along the current axis

		// Auto colour-range: when ON (default) the slice colormap, legend and arrow speed scale follow
		// a live GPU min/max reduction of the displayed field; when OFF, the fixed per-field range
		// (updateRange) is used. Main thread.
		void setAutoRange(bool on);
		bool autoRange() const { return auto_range_; }

		// --- Layer visibility toggles (feature 2) --------------------------------
		void setShowSlice(bool on) { show_slice_ = on; update(); }
		void setShowModel(bool on) { show_model_ = on; update(); }
		// The voxel-mask overlay is the exposed-face staircase of the rigid solid mask (obstacle).
		void setShowSolidVoxels(bool on) { show_solid_vox_ = on; update(); }
		void setShowAxes(bool on) { show_axes_ = on; update(); }
		// Grid overlay: draw the cell-boundary LINES where the grid intersects the current slice plane, so
		// the graded mesh is visible (lines cluster in the fine h_fine core, spread out in the coarse far
		// field). Uses the per-axis face coordinates from setGridLines (uniform i·h if none set). Main thread.
		void setShowGrid(bool on) { show_grid_ = on; update(); }
		bool showGrid() const { return show_grid_; }
		// Per-axis cumulative cell-face coordinates (metres) for the grid overlay — the graded metric arrays
		// (GridMetrics::xf/yf/zf) or i·h for a uniform grid. Empty ⇒ fall back to info_'s uniform spacing.
		void setGridLines(const std::vector<double>& xf, const std::vector<double>& yf, const std::vector<double>& zf);
		bool showSlice() const { return show_slice_; }
		bool showModel() const { return show_model_; }
		bool showSolidVoxels() const { return show_solid_vox_; }
		bool showAxes() const { return show_axes_; }

		// Colour the building's voxel-overlay faces by the owning cell's surface pressure coefficient Cp
		// (published by the worker; mapped through colormap.h over a symmetric range). ON (default) tints
		// the staircase blue=suction … red=pressure; OFF reverts to the uniform dark-blue solid colour.
		// Toggling re-uploads the overlay (recolour) on the next paint. Main thread.
		void setColourByCp(bool on);
		bool colourByCp() const { return colour_by_cp_; }

		// --- Clip plane (see inside hollow structures) ---------------------------
		// A single movable plane that hides SOLIDS on the camera side — the STEP mesh and the voxel-solid
		// staircase — so the interior of a hollow structure is exposed. The flow
		// slice and the arrows are NEVER clipped, so the flow field inside the revealed cavity stays on
		// screen. Modes: axis-aligned X/Y/Z (the position slider shifts it along that axis, auto-oriented
		// to hide the camera side) or "Face camera" (normal = view direction; the slider pushes it into
		// the scene along the line of sight). Flip swaps the hidden half. Implemented with gl_ClipDistance.
		void setClipEnabled(bool on) { clip_enabled_ = on; update(); }
		void setClipMode(int m) { clip_mode_ = (m < 0) ? 0 : (m > 3 ? 3 : m); update(); } // 0=X 1=Y 2=Z 3=Camera
		void setClipFraction(float f) { clip_frac_ = f < 0.0f ? 0.0f : (f > 1.0f ? 1.0f : f); update(); }
		void setClipFlip(bool on) { clip_flip_ = on; update(); }
		bool clipEnabled() const { return clip_enabled_; }
		int clipMode() const { return clip_mode_; }

		// --- Flow arrows (features 4 + 5) ----------------------------------------
		void setShowArrows(bool on);
		void setArrowMode3D(bool three_d); // false ⇒ arrows live on the current slice plane
		void setArrowDensity(int n);       // particle count (applied live)
		// Visual advection-speed multiplier in [0,1] (1 = the default lively motion, 0 = frozen). Scales
		// ONLY the tracer animation so fast currents (e.g. 8 m/s) don't whip the arrows about — physics
		// is untouched. Applied live.
		void setArrowSpeedMult(float m) { arrow_speed_mult_ = m < 0.0f ? 0.0f : (m > 1.0f ? 1.0f : m); update(); }
		// Lateral thickness scale for the arrow glyph (shaft + head): 1 = the glyph's own width, lower =
		// thinner / less cluttered. The arrow LENGTH is unaffected (declutter knob for a dense field).
		void setArrowWidthMult(float m) { arrow_width_ = m < 0.05f ? 0.05f : (m > 3.0f ? 3.0f : m); update(); }
		bool showArrows() const { return show_arrows_; }
		bool arrowMode3D() const { return arrow_3d_; }
		int arrowDensity() const { return arrow_density_; }
		float arrowSpeedMult() const { return arrow_speed_mult_; }
		float arrowWidthMult() const { return arrow_width_; }

		// --- Flow tracers (inlet-seeded streamlines: long lines coloured by speed) ----------------
		// A grid of seeds on the INLET face each integrate a streamline through the live velocity field
		// all the way across the domain (until it leaves an edge, meets a solid, stalls, or curls into
		// an eddy). Re-integrated every frame so they warp with the wake; a scrolling brightness pulse
		// reads as flow. Like the arrows they read the worker's host-flow snapshot and run on the main
		// thread; the DRAW is camera-facing RIBBONS (glMultiDrawArrays(GL_TRIANGLE_STRIP)) so the width
		// is real geometry (glLineWidth is capped by most core-profile drivers). Never clipped.
		void setShowTracers(bool on);
		void setTracerMode3D(bool three_d); // false ⇒ seed + integrate on the current slice plane
		void setTracerGridDensity(int d);   // inlet seed count along the larger inlet dimension
		void setTracerTrail(int n);         // max streamline length [integration steps]
		// Ribbon thickness in PIXELS (screen-space, constant regardless of zoom). A few px reads well.
		void setTracerWidth(float w) { tracer_width_ = w < 0.5f ? 0.5f : (w > 12.0f ? 12.0f : w); update(); }
		// "Boring" filter [0.99,1.1]: hide streamlines whose path length is below this multiple of the
		// domain length Lx. 0.99 (the minimum) disables the filter (show all); higher hides all but the
		// meandering (eddy/wake) lines. A straight crosser has path length ≈ Lx.
		void setTracerBoring(float m) { tracer_boring_ = m < 0.99f ? 0.99f : (m > 1.1f ? 1.1f : m); update(); }
		// While the boring slider is being dragged, bypass the anti-flicker retention so the filter
		// updates instantly; retention resumes on release.
		void setTracerBoringInstant(bool on) { tracer_boring_instant_ = on; update(); }
		bool showTracers() const { return show_tracers_; }
		bool tracerMode3D() const { return tracer_3d_; }
		int tracerGridDensity() const { return tracer_grid_density_; }
		int tracerTrail() const { return tracer_trail_; }
		float tracerWidth() const { return tracer_width_; }
		float tracerBoring() const { return tracer_boring_; }

		// --- Legend metadata (feature 3 queries by the overlay) ------------------
		Field field() const { return field_; }
		float fieldMin() const { return vmin_; }
		float fieldMax() const { return vmax_; }

		// Loaded STEP model (metres). setMesh takes ownership; the GL upload is deferred to
		// the next paint (main thread) so it is safe to call before the context exists (CLI
		// --load-step) or from a menu action. clearMesh removes it from the view.
		void setMesh(paracfd::core::TriMesh mesh);
		void clearMesh();
		bool hasMesh() const { return has_mesh_; }

		// Solid-cell overlay: the staircase voxelization of the model, drawn as the exposed
		// voxel-surface faces over/under the smooth mesh so the user can judge resolution. The
		// mask is the ChannelBC cell field (1=solid). Geometry build is deferred to the next
		// paint (main thread). clearVoxelOverlay hides it.
		void setVoxelOverlay(const std::vector<unsigned char>& solid, paracfd::core::MacGrid grid);
		void clearVoxelOverlay();

		// Draw the loaded model at an explicit translate (metres) instead of the auto bed placement.
		void setMeshTranslate(double tx, double ty, double tz);
		// Draw the loaded model under an explicit AFFINE placement (rotation·scale + translation),
		// overriding the gizmo transform.
		void setMeshPlacement(const paracfd::core::ModelPlacement& p);

		// --- Interactive model-placement gizmo -----------------------------------
		// A single UNIFIED 3D manipulator: all handles are drawn at once — 3 translate arrows, 3 rotate
		// rings and 3 scale cubes (one per axis) plus a uniform-scale centre cube — and whichever handle
		// you grab picks the operation dynamically (no mode switching). The transform is the model's world
		// placement — world(v) = t + rot·(scale ⊙ (v − pivot)) about its bbox centre — exposed to the
		// voxelizer as a ModelPlacement so the model is voxelized exactly where it is drawn (on the next
		// "Apply" / obstacle-inject).
		struct ModelGizmoXform
		{
			QVector3D pivot, t, scale{ 1, 1, 1 };
			QQuaternion rot;
			bool valid = false;
		};
		void setGizmoEnabled(bool on);
		bool gizmoEnabled() const { return gizmo_on_; }
		bool gizmoEditable() const { return has_mesh_ && gz_valid_ && !mesh_override_; }
		// True once a gizmo placement exists (a model is loaded and not seated by an override) —
		// so the voxelizer can use modelPlacement() instead of the default place_model_on_bed, even
		// before the first paint (the CLI --load-step path voxelizes immediately). Unlike gizmoEditable
		// it does not require the GL upload to have happened.
		bool hasModelPlacement() const { return gz_valid_ && !mesh_override_; }
		void resetModelPlacement();                                  // back to centre-on-bed
		paracfd::core::ModelPlacement modelPlacement() const;          // the affine the voxelizer uses
		void setModelPlacement(const paracfd::core::ModelPlacement& p);// restore an affine (scene load)
		ModelGizmoXform modelXform() const;                          // full gizmo state (Apply capture/restore)
		void setModelXform(const ModelGizmoXform& x);

		double averageFps() const { return avg_fps_; }
		long long framesRendered() const { return frames_total_; }

	signals:
		void fpsUpdated(double fps);
		// Emitted when the gizmo finishes editing the model placement (drag release / reset), so the
		// host can refresh its readout. The re-voxelization happens on the next "Apply", not per-drag.
		void modelPlacementChanged();

	protected:
		void initializeGL() override;
		void resizeGL(int w, int h) override;
		void paintGL() override;

		void mousePressEvent(QMouseEvent* e) override;
		void mouseMoveEvent(QMouseEvent* e) override;
		void mouseReleaseEvent(QMouseEvent* e) override;
		void wheelEvent(QWheelEvent* e) override;

	private:
		void buildSliceGeometry(); // (re)generate static positions for axis/plane
		void buildBoxGeometry();
		void uploadMesh();        // push pending_mesh_ into GL buffers (main thread)
		void uploadVoxelOverlay(); // build exposed voxel-surface faces into GL (main thread)
		void updateRange();
		void applyAutoRange(const FieldRange& fr); // EMA-fold a live reduction into [vmin_,vmax_]+speed scale
		void setDefaultPlane();
		SliceParams currentParams() const;

		void buildArrowGlyph();   // static unit-arrow instance-base mesh (main thread)
		void updateArrows();      // advect the tracers + upload instance data (main thread)
		void drawArrows(const QMatrix4x4& mvp); // instanced draw of the arrow field
		void buildTracerBuffers(); // wire the streakline VAO/VBOs into the slice shader (main thread)
		void updateTracers();      // advect the streaklines + upload line geometry (main thread)
		void drawTracers(const QMatrix4x4& mvp); // glMultiDrawArrays line-strip draw of the streaklines
		void updateWantHostFlow(); // ask the worker for the host-flow snapshot iff arrows OR tracers are on
		void drawLegendWith(QPainter& p); // QPainter colorbar overlay (feature 3)

		// Clip plane: world-space plane (n.x, n.y, n.z, d); the KEPT half-space is dot(pos,n)+d >= 0, so
		// solids on the camera side are cut. Returns an all-keep plane when disabled / no domain yet.
		QVector4D computeClipPlane() const;
		void drawClipPlaneViz(const QMatrix4x4& mvp, const QVector4D& plane); // faint quad + bright outline

		// --- Model-placement gizmo internals -------------------------------------
		QMatrix4x4 modelMatrix() const;   // effective model transform (override matrix or gizmo TRS)
		void computeDefaultXform();       // centre-on-bed gz_* from the stored bbox + info_ (place_on_bed)
		void drawGizmo(const QMatrix4x4& mvp);            // all manipulator glyphs at once (dynamic flat lines)
		bool pickGizmo(const QPointF& pos, int& op, int& axis) const; // nearest handle → (op,axis); false if none
		void beginGizmoDrag(int op, int axis, const QPointF& pos);    // capture the drag reference for (op,axis)
		void dragGizmo(const QPointF& cur, const QPointF& prev); // apply one drag step to the active handle
		void mouseRay(const QPointF& pos, QVector3D& orig, QVector3D& dir) const; // world-space pick ray
		QVector3D gizmoAxisDir(int a) const;              // world direction of local axis a (rot·e_a)
		float gizmoSize() const;                          // handle length [m], scaled to camera distance

		void buildGridGeometry(); // cell-boundary lines where the grid meets the current slice plane (main thread)
		void drawGrid(const QMatrix4x4& mvp); // draw the grid-overlay lines (flat colour) when show_grid_
		void buildAxesGeometry(); // world-origin XYZ triad + metre ticks (main thread; needs the domain)
		void buildCornerGizmo();  // static camera-aligned orientation triad (unit axes)
		void drawAxes(const QMatrix4x4& mvp);  // world triad (depth-tested) + corner gizmo (on top)
		void drawAxesLabels(QPainter& p);      // X/Y/Z + metre tick labels via the QPainter overlay
		bool projectPoint(const QMatrix4x4& mvp, const QVector3D& w, QPointF& px) const; // world→screen px

		SimWorker* worker_ = nullptr;
		SimInfo info_;
		bool have_info_ = false;

		Camera camera_;
		QPoint last_mouse_;
		Qt::MouseButton drag_btn_ = Qt::NoButton;

		// Slice view state. Default Z-normal = the x-y plane (best view of the cylinder
		// wake); matches the axis combo's default selection in MainWindow.
		Axis axis_ = Axis::Z;
		float plane_frac_ = 0.5f;
		Field field_ = Field::SpeedMag;
		float vmin_ = 0.0f, vmax_ = 1.0f;

		// Auto colour-range (default ON). A throttled GPU reduction over the device snapshot fills
		// FieldRange; applyAutoRange folds it into [vmin_,vmax_] (expand-fast / shrink-slow EMA) and
		// auto_speed_max_ (the arrow speed scale). range_valid_ resets on a field/domain/toggle change
		// so the first sample snaps rather than lagging a stale range.
		bool auto_range_ = true;
		bool range_valid_ = false;
		float auto_speed_max_ = 1.0f;   // smoothed max |u| over fluid cells [m/s] (arrow colour scale)
		int range_ctr_ = 0;             // paint counter → reduce every kRangeEvery frames (~12 Hz)
		int range_log_ctr_ = 0;         // throttles the [vmin,vmax] diagnostic line
		static constexpr int kRangeEvery = 5;

		// Fixed mesh resolution (nearest-cell sampling => oversampling is smooth & cheap).
		static constexpr int NRES = 200;

		// GL objects.
		QOpenGLShaderProgram prog_;
		unsigned int slice_vao_ = 0, pos_vbo_ = 0, color_vbo_ = 0, idx_ebo_ = 0;
		unsigned int box_vao_ = 0, box_vbo_ = 0;
		int index_count_ = 0, box_vertex_count_ = 0;
		void* cuda_res_ = nullptr; // registered color_vbo_
		bool gl_ready_ = false;
		bool geometry_dirty_ = true;

		// STEP model (lit triangle mesh). Own shader + VAO/VBOs; model matrix places it on
		// the domain floor centred in x/y. Upload deferred via mesh_upload_pending_.
		QOpenGLShaderProgram mesh_prog_;
		unsigned int mesh_vao_ = 0, mesh_pos_vbo_ = 0, mesh_norm_vbo_ = 0, mesh_idx_ebo_ = 0;
		int mesh_index_count_ = 0;
		bool has_mesh_ = false;
		bool mesh_upload_pending_ = false;
		paracfd::core::TriMesh pending_mesh_;

		// Model transform. The model matrix comes from the gizmo TRS (gz_* below); an optional explicit
		// override matrix can seat the model directly instead.
		QMatrix4x4 mesh_override_mat_; // explicit model matrix (override seating)
		bool mesh_override_ = false;   // true ⇒ use mesh_override_mat_ instead of the gizmo transform

		// Interactive gizmo state: world(v) = t + rot·(scale ⊙ (v − pivot)). gz_bbox_* is the model-local
		// bounding box kept so a "reset placement" can recompute the centre-on-bed default without the mesh
		// (the CPU copy is freed after upload). gizmo_on_ toggles the whole unified manipulator; the active/
		// hovered (op,axis) below selects which of the 3×3 (+uniform) functions is being driven. gz_drag_n_
		// pins the active world axis for the drag; gz_rot_last_ang_ is the previous plane angle of a rotate.
		QVector3D gz_pivot_{ 0, 0, 0 }, gz_t_{ 0, 0, 0 }, gz_scale_{ 1, 1, 1 };
		QQuaternion gz_rot_;
		bool gz_valid_ = false;
		QVector3D gz_bbox_min_{ 0, 0, 0 }, gz_bbox_max_{ 0, 0, 0 };
		bool gizmo_on_ = false;
		// Active/hovered handle = (operation, axis). op: 0 translate, 1 rotate, 2 scale (per-axis), 3 uniform
		// scale (axis ignored); -1 = none. All handles are drawn together and the grabbed one picks the op.
		int gz_hover_op_ = -1, gz_hover_axis_ = -1;
		int gz_drag_op_ = -1, gz_drag_axis_ = -1;
		QVector3D gz_drag_n_{ 0, 0, 0 };
		float gz_rot_last_ang_ = 0.0f;
		unsigned int gizmo_vao_ = 0, gizmo_vbo_ = 0; // dynamic manipulator line geometry

		// Voxel-overlay (staircase mask) GL objects. Exposed solid-cell faces, lit shader.
		unsigned int vox_vao_ = 0, vox_pos_vbo_ = 0, vox_norm_vbo_ = 0, vox_color_vbo_ = 0;
		int vox_vertex_count_ = 0;     // total exposed-face vertices
		bool has_vox_ = false;
		bool vox_upload_pending_ = false;
		std::vector<unsigned char> pending_vox_solid_; // retained after upload so a Cp update can re-colour
		paracfd::core::MacGrid pending_vox_grid_;
		// Live flow-mask overlay: the last mask generation pulled from the worker. The overlay is
		// re-extracted (exposed faces only) whenever the worker republishes a changed mask.
		std::uint64_t vox_mask_gen_ = 0;

		// Cp colouring of the voxel overlay (feature: colour building by surface pressure coefficient).
		// vox_cell_cp_ is the worker's per-solid-cell mean Cp (g.p_count(), NaN off the surface); the
		// symmetric range [vox_cp_lo_,vox_cp_hi_] is auto-tracked from the published cp_min/cp_max. The
		// overlay carries a per-vertex colour VBO built from these when colour_by_cp_ is on. vox_load_gen_
		// mirrors the worker's loadGeneration() so a fresh Cp field triggers a recolour (like the mask).
		bool colour_by_cp_ = true;
		std::vector<float> vox_cell_cp_;
		float vox_cp_lo_ = -1.5f, vox_cp_hi_ = 1.0f;
		bool has_vox_cp_ = false;      // per-vertex Cp colours are built + valid for the current overlay
		std::uint64_t vox_load_gen_ = 0;

		// Layer visibility (feature 2). All default ON.
		bool show_slice_ = true;
		bool show_model_ = true;
		bool show_solid_vox_ = true; // rigid solid voxels (obstacle, dark-blue)

		// Clip plane (see inside hollow structures). Cuts SOLIDS only (mesh + voxel solids); the
		// flow slice + arrows stay visible. clip_mode_: 0=X 1=Y 2=Z 3=Camera-facing; clip_frac_ ∈ [0,1]
		// sweeps the plane (axis position, or view-depth about the camera target); clip_flip_ swaps the
		// hidden side. clip_vao_/vbo_ hold the translucent plane-visualisation quad (4 corners, per-paint).
		bool clip_enabled_ = false;
		int clip_mode_ = 2;       // default Z-normal
		float clip_frac_ = 0.5f;
		bool clip_flip_ = false;
		unsigned int clip_vao_ = 0, clip_vbo_ = 0;

		// XYZ axis triad: coloured axes from the WORLD ORIGIN (0,0,0 = the x=0/y=0/z=0 domain corner,
		// the inlet-side bed corner; place_model_on_bed centres the model in x/y at z=0) with metre
		// ticks + labels, plus a camera-aligned orientation gizmo in the bottom-left. Default ON.
		bool show_axes_ = true;
		bool axes_dirty_ = true;                       // rebuild world-triad geometry on a domain change

		// Grid overlay (cell-boundary lines on the current slice plane; shows the graded mesh). Default OFF.
		// grid_*f_ are the per-axis cumulative face coordinates (metres); empty ⇒ uniform i·h from info_.
		bool show_grid_ = false;
		bool grid_dirty_ = true;                       // rebuild on a plane/axis/domain/metric change
		std::vector<float> grid_xf_, grid_yf_, grid_zf_;
		unsigned int grid_vao_ = 0, grid_vbo_ = 0;
		int grid_vertex_count_ = 0;
		unsigned int axis_vao_ = 0, axis_vbo_ = 0;     // world triad + tick marks (per-axis ranges below)
		int axis_vert_off_[3] = { 0, 0, 0 }, axis_vert_cnt_[3] = { 0, 0, 0 };
		float axis_tick_step_ = 1.0f;                  // metre spacing between ticks (shared by geom + labels)
		unsigned int corner_vao_ = 0, corner_vbo_ = 0; // static unit-axis gizmo
		QMatrix4x4 corner_mvp_;                         // gizmo transform (set in drawAxes, read by labels)

		// Flow-arrow field (features 4+5). Instanced unit-arrow glyph + per-particle instance VBO.
		FlowParticles arrows_;
		bool show_arrows_ = true;
		bool arrow_3d_ = true;         // 3D volume advection vs 2D on the slice plane
		int arrow_density_ = 1500;     // particle count
		float arrow_speed_mult_ = 1.0f; // visual advection-speed multiplier [0,1] (1 = default, 0 = frozen)
		float arrow_width_ = 0.5f;      // lateral thickness scale (1 = glyph width); default thinner to declutter
		QOpenGLShaderProgram arrow_prog_;
		unsigned int arrow_vao_ = 0, arrow_glyph_vbo_ = 0, arrow_inst_vbo_ = 0;
		int arrow_glyph_verts_ = 0;
		int arrow_inst_capacity_ = 0;  // #particles the instance VBO is sized for
		int arrow_draw_count_ = 0;     // particles to draw this frame
		float arrow_len_ = 0.3f;       // base glyph length [m]
		QElapsedTimer arrow_clock_;    // per-frame dt for advection
		bool arrow_clock_started_ = false;
		bool arrows_logged_ = false;   // one-shot "arrows live" diagnostic

		// Flow-tracer field (inlet-seeded streamlines). CPU integration + a camera-facing RIBBON
		// geometry buffer drawn through its own tracer shader (real thickness). Default OFF (arrows are
		// the default flow viz). Density = inlet seed count; the ribbons update every frame.
		FlowTracers tracers_;
		bool show_tracers_ = false;
		bool tracer_3d_ = true;          // 3D inlet-plane seeding vs 2D on the slice plane
		int tracer_grid_density_ = 12;   // inlet seed count along the larger inlet dimension
		int tracer_trail_ = 600;         // max streamline length [integration steps]
		float tracer_width_ = 2.0f;      // ribbon thickness in pixels (screen-space)
		float tracer_boring_ = 0.99f;    // hide streamlines shorter than Lx·this (0.99 = off … 1.1)
		bool tracer_boring_instant_ = false; // true while the boring slider is dragged (bypass the hold)
		QOpenGLShaderProgram tracer_prog_;
		unsigned int tracer_vao_ = 0, tracer_vbo_ = 0; // interleaved ribbon geometry (stride 11 floats)
		int tracer_vbo_capacity_ = 0;    // floats the VBO is sized for
		int tracer_strips_ = 0;          // triangle strips to draw this frame
		std::vector<int> tracer_firsts_, tracer_counts_; // glMultiDrawArrays offsets + lengths
		QElapsedTimer tracer_clock_;     // per-frame dt for the temporal boring-filter hold
		bool tracer_clock_started_ = false;
		bool tracers_logged_ = false;    // one-shot "tracers live" diagnostic

		// fps bookkeeping.
		QElapsedTimer fps_timer_;
		int fps_frames_ = 0;
		double avg_fps_ = 0.0;
		long long frames_total_ = 0;
	};
}
