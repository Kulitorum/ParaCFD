// slice_viewer.cpp — see slice_viewer.h. GL 4.3 core rendering of a CUDA-filled slice.
#include "gui/slice_viewer.h"

#include "core/geometry/model_placement.h"
#include "gui/colormap.h"
#include "gui/gl_thread_check.h"
#include "gui/paraglider_sim_worker.h"
#include "gui/slice_gl.h"
#include "gui/sim_worker.h"

#include <QMatrix3x3>
#include <QMouseEvent>
#include <QPainter>
#include <QVector2D>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <vector>

namespace paracfd::gui
{
	namespace
	{
		std::vector<float> debug_box_lines(const std::vector<std::array<float, 6>>& boxes)
		{
			std::vector<float> lines;lines.reserve(boxes.size()*12*2*3);
			static constexpr int edges[12][2]={{0,1},{0,2},{0,4},{1,3},{1,5},{2,3},{2,6},{3,7},{4,5},{4,6},{5,7},{6,7}};
			for(const auto& b:boxes)
			{
				const float p[8][3]={{b[0],b[1],b[2]},{b[3],b[1],b[2]},{b[0],b[4],b[2]},{b[3],b[4],b[2]},{b[0],b[1],b[5]},{b[3],b[1],b[5]},{b[0],b[4],b[5]},{b[3],b[4],b[5]}};
				for(const auto& edge:edges)for(int endpoint:edge)lines.insert(lines.end(),{p[endpoint][0],p[endpoint][1],p[endpoint][2]});
			}
			return lines;
		}

		const char* kVert = R"(#version 430 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec4 aColor;
uniform mat4 uMVP;
uniform int  uFlat;
uniform vec4 uColor;
uniform vec4 uClipPlane; // (n, d); keep dot(pos,n)+d >= 0. Only clips when GL_CLIP_DISTANCE0 is enabled.
out vec4 vColor;
void main()
{
	vColor = (uFlat != 0) ? uColor : aColor;
	gl_ClipDistance[0] = dot(vec4(aPos, 1.0), uClipPlane);
	gl_Position = uMVP * vec4(aPos, 1.0);
}
)";
		const char* kFrag = R"(#version 430 core
in vec4 vColor;
out vec4 fragColor;
void main() { fragColor = vColor; }
)";

		// Lit shader for the loaded STEP model: transforms by uMVP (= VP * model) for the
		// clip position, and lights with a camera headlight (two-sided so winding never
		// leaves a facet black). uModel is translation-only, so its 3x3 == identity and the
		// vertex normal doubles as the world normal.
		const char* kMeshVert = R"(#version 430 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNormal;
layout(location=2) in vec3 aColor; // per-vertex Cp colour (voxel overlay only; ignored unless uUseVertexColor)
uniform mat4 uMVP;
uniform mat4 uModel;
uniform vec4 uClipPlane; // (n, d); keep dot(world,n)+d >= 0. Only clips when GL_CLIP_DISTANCE0 is enabled.
out vec3 vN;
out vec3 vWorld;
out vec3 vColor;
void main()
{
	vec4 wp = uModel * vec4(aPos, 1.0);
	vWorld = wp.xyz;
	vN = mat3(uModel) * aNormal;
	vColor = aColor;
	gl_ClipDistance[0] = dot(vec4(wp.xyz, 1.0), uClipPlane);
	gl_Position = uMVP * vec4(aPos, 1.0);
}
)";
		const char* kMeshFrag = R"(#version 430 core
in vec3 vN;
in vec3 vWorld;
in vec3 vColor;
uniform vec3 uEye;
uniform vec4 uBaseColor;
uniform int  uUseVertexColor; // 0 = flat uBaseColor (STEP model / uniform voxels); !=0 = per-vertex Cp colour
uniform int  uUseTriangleColor; // STEP surface: gl_PrimitiveID -> per-triangle, per-side Cp SSBO
struct TriangleSideColour { vec4 plus; vec4 minus; };
layout(std430, binding=3) readonly buffer TriangleColours { TriangleSideColour uTriangleColor[]; };
out vec4 fragColor;
void main()
{
	vec3 N = normalize(vN);
	vec3 L = normalize(uEye - vWorld);
	float diff = max(abs(dot(N, L)), 0.0); // two-sided
	float ambient = 0.28;
	vec3 base;
	if (uUseTriangleColor != 0)
		base = gl_FrontFacing ? uTriangleColor[gl_PrimitiveID].plus.rgb :
			uTriangleColor[gl_PrimitiveID].minus.rgb;
	else
		base = (uUseVertexColor != 0) ? vColor : uBaseColor.rgb;
	vec3 c = base * (ambient + 0.72 * diff);
	fragColor = vec4(c, uBaseColor.a);
}
)";

		// Flow-arrow instanced shader (feature 4). Per-instance {pos,dir,speed,alpha}; a unit 2D
		// arrow glyph (aGlyph) is oriented so its shaft follows the velocity direction and its flat
		// plane faces the camera (billboard), length-scaled by speed. Coloured BY SPEED through the
		// SAME 5-stop ramp as the slice (GLSL mirror of gui/colormap.h) so arrows match the legend.
		const char* kArrowVert = R"(#version 430 core
layout(location=0) in vec2 aGlyph;   // base arrow (shaft + head), shaft along +x in [0,1]
layout(location=1) in vec3 aPos;     // instance world position
layout(location=2) in vec3 aDir;     // instance unit velocity direction
layout(location=3) in float aSpeed;  // normalised speed [0,1]
layout(location=4) in float aAlpha;  // fade
uniform mat4 uMVP;
uniform vec3 uEye;
uniform float uLen;                  // base glyph length [m]
uniform float uWidth;                // lateral thickness scale (1 = the glyph's own width; length unaffected)
out float vSpeed;
out float vAlpha;
void main()
{
	vSpeed = aSpeed;
	vAlpha = aAlpha;
	vec3 T = (dot(aDir, aDir) > 1e-8) ? normalize(aDir) : vec3(1.0, 0.0, 0.0);
	vec3 E = normalize(uEye - aPos);
	vec3 N = cross(T, E);                          // in-screen perpendicular to the shaft
	if (dot(N, N) < 1e-8) N = cross(T, vec3(0.0, 0.0, 1.0));
	if (dot(N, N) < 1e-8) N = cross(T, vec3(0.0, 1.0, 0.0));
	N = normalize(N);
	float L = uLen * (0.35 + 0.9 * aSpeed);        // length grows with speed
	// Shaft length scales with L; the lateral (across-shaft) extent scales with L*uWidth, so lowering
	// uWidth makes the shaft + head THINNER while the arrow length stays the same (declutter knob).
	vec3 world = aPos + T * (aGlyph.x * L) + N * (aGlyph.y * L * uWidth);
	gl_Position = uMVP * vec4(world, 1.0);
}
)";
		const char* kArrowFrag = R"(#version 430 core
in float vSpeed;
in float vAlpha;
uniform int uColorMode;    // 0 = speed colormap, 1 = solid uSolidColor
uniform vec3 uSolidColor;
out vec4 fragColor;
// MIRROR of gui/colormap.h scour_colormap() — keep the stops in sync with that header.
vec3 cmap(float t)
{
	t = clamp(t, 0.0, 1.0);
	vec3 c0 = vec3(0.23, 0.30, 0.75);
	vec3 c1 = vec3(0.10, 0.65, 0.90);
	vec3 c2 = vec3(0.20, 0.78, 0.30);
	vec3 c3 = vec3(0.98, 0.90, 0.20);
	vec3 c4 = vec3(0.85, 0.18, 0.15);
	float x = t * 4.0;
	int i = int(x);
	if (i > 3) i = 3;
	float f = x - float(i);
	vec3 lo = (i == 0) ? c0 : (i == 1) ? c1 : (i == 2) ? c2 : c3;
	vec3 hi = (i == 0) ? c1 : (i == 1) ? c2 : (i == 2) ? c3 : c4;
	return mix(lo, hi, f);
}
void main()
{
	if (vAlpha <= 0.01) discard;
	vec3 c = (uColorMode == 1) ? uSolidColor : cmap(vSpeed);
	fragColor = vec4(c, vAlpha);
}
)";

		// Flow-tracer ribbon shader (inlet streamlines). Each streamline point is emitted as two
		// vertices (aPos.w = side ±1) sharing a tangent aTan; the vertex is offset ±side along the
		// line's SCREEN-SPACE perpendicular by uHalfPx PIXELS, so the ribbon has a constant pixel
		// thickness regardless of zoom/depth (glLineWidth is capped on core profiles). Colour is baked
		// per-vertex on the CPU and passed straight through; depth (clip.z/w) is preserved so the line
		// is still depth-tested against solids.
		const char* kTracerVert = R"(#version 430 core
layout(location=0) in vec4 aPos;   // xyz world position, w = side (-1 / +1)
layout(location=1) in vec3 aTan;   // world line tangent
layout(location=2) in vec4 aColor; // rgba (speed colour, end-faded alpha)
uniform mat4 uMVP;
uniform vec2 uViewport;            // framebuffer size [px]
uniform float uHalfPx;             // half ribbon width [px]
out vec4 vColor;
void main()
{
	vColor = aColor;
	vec4 clip = uMVP * vec4(aPos.xyz, 1.0);
	vec4 clipT = uMVP * vec4(aPos.xyz + aTan, 1.0);
	// Line direction in pixel space; the ribbon expands along its perpendicular.
	vec2 s0 = clip.xy / clip.w;
	vec2 s1 = clipT.xy / clipT.w;
	vec2 dir = (s1 - s0) * uViewport;
	dir = (dot(dir, dir) > 1e-12) ? normalize(dir) : vec2(1.0, 0.0);
	vec2 nrm = vec2(-dir.y, dir.x);
	vec2 offNdc = nrm * (aPos.w * uHalfPx) * 2.0 / uViewport; // pixels → NDC
	clip.xy += offNdc * clip.w;                               // keep constant px under perspective divide
	gl_Position = clip;
}
)";
		const char* kTracerFrag = R"(#version 430 core
in vec4 vColor;
out vec4 fragColor;
void main()
{
	if (vColor.a <= 0.01) discard;
	fragColor = vColor;
}
)";
	}

	SliceViewer::SliceViewer(QWidget* parent) : QOpenGLWidget(parent)
	{
		setFocusPolicy(Qt::StrongFocus);
		setMouseTracking(true); // hover-highlight the gizmo handle under the cursor (no button held)
		arrows_.set_count(arrow_density_);
	}

	void SliceViewer::setWorker(SimWorker* w)
	{
		worker_ = w;
		vox_mask_gen_ = 0;      // re-pull the voxel overlay from the (new) worker's live mask
		vox_load_gen_ = 0;      // re-pull the Cp field from the (new) worker
		vox_cell_cp_.clear();
		has_vox_cp_ = false;
		clearVoxelOverlay();    // drop any stale mask geometry until the new one is published
		updateWantHostFlow();
	}

	void SliceViewer::setParagliderWorker(ParagliderSimWorker* w)
	{
		paraglider_worker_ = w;
		range_valid_ = false;
		arrows_.reset();
		tracers_.reset();
		update();
	}

	void SliceViewer::setColourByCp(bool on)
	{
		if (colour_by_cp_ == on) return;
		colour_by_cp_ = on;
		// Re-extract the overlay so it picks up (or drops) the per-vertex Cp colours. The mask geometry is
		// retained in pending_vox_solid_, so this is a pure recolour when a building is present.
		if (has_vox_ || !pending_vox_solid_.empty()) vox_upload_pending_ = true;
		update();
	}

	// The worker only D2H-publishes the host {u,v,w,solid} snapshot when asked; both the arrows AND
	// the tracers advect from it, so it is needed iff EITHER is on (a viewer with both hidden pays
	// nothing).
	void SliceViewer::updateWantHostFlow()
	{
		if (worker_) worker_->setWantHostFlow(show_arrows_ || show_tracers_);
	}

	void SliceViewer::setShowArrows(bool on)
	{
		show_arrows_ = on;
		updateWantHostFlow();
		if (on) arrows_.reset(); // re-seed so they don't reappear frozen
		update();
	}

	void SliceViewer::setShowTracers(bool on)
	{
		show_tracers_ = on;
		updateWantHostFlow();
		if (on) tracers_.reset(); // re-seed the grid so streaks don't reappear frozen
		update();
	}

	void SliceViewer::setTracerMode3D(bool three_d)
	{
		if (tracer_3d_ == three_d) return;
		tracer_3d_ = three_d;
		tracers_.reset(); // 2D<->3D changes the seed lattice, rebuild it
		update();
	}

	void SliceViewer::setTracerGridDensity(int d)
	{
		d = std::max(2, std::min(d, 120));
		if (d == tracer_grid_density_) return;
		tracer_grid_density_ = d;
		tracers_.reset(); // the seed count changed, rebuild the lattice
		update();
	}

	void SliceViewer::setTracerTrail(int n)
	{
		tracer_trail_ = std::max(20, std::min(n, 3000)); // max streamline length; applied next frame
		update();
	}

	void SliceViewer::setArrowMode3D(bool three_d)
	{
		if (arrow_3d_ == three_d) return;
		arrow_3d_ = three_d;
		arrows_.reset(); // 2D<->3D changes the seeding domain, respawn all
		update();
	}

	void SliceViewer::setArrowDensity(int n)
	{
		n = std::max(50, std::min(n, 16000));
		if (n == arrow_density_) return;
		arrow_density_ = n;
		arrows_.set_count(n); // reallocates + respawns
		update();
	}

	SliceViewer::~SliceViewer()
	{
		if (!gl_ready_) return;
		makeCurrent();
		if (cuda_res_) slice_gl_unregister(cuda_res_);
		if (pos_vbo_) glDeleteBuffers(1, &pos_vbo_);
		if (color_vbo_) glDeleteBuffers(1, &color_vbo_);
		if (idx_ebo_) glDeleteBuffers(1, &idx_ebo_);
		if (box_vbo_) glDeleteBuffers(1, &box_vbo_);
		if (grid_vbo_) glDeleteBuffers(1, &grid_vbo_);
		if (mesh_pos_vbo_) glDeleteBuffers(1, &mesh_pos_vbo_);
		if (mesh_norm_vbo_) glDeleteBuffers(1, &mesh_norm_vbo_);
		if (mesh_idx_ebo_) glDeleteBuffers(1, &mesh_idx_ebo_);
		if (mesh_triangle_colour_ssbo_) glDeleteBuffers(1, &mesh_triangle_colour_ssbo_);
		if (amr_debug_vbo_) glDeleteBuffers(1, &amr_debug_vbo_);
		if (eb_debug_vbo_) glDeleteBuffers(1, &eb_debug_vbo_);
		if (vox_pos_vbo_) glDeleteBuffers(1, &vox_pos_vbo_);
		if (vox_norm_vbo_) glDeleteBuffers(1, &vox_norm_vbo_);
		if (vox_color_vbo_) glDeleteBuffers(1, &vox_color_vbo_);
		if (arrow_glyph_vbo_) glDeleteBuffers(1, &arrow_glyph_vbo_);
		if (arrow_inst_vbo_) glDeleteBuffers(1, &arrow_inst_vbo_);
		if (tracer_vbo_) glDeleteBuffers(1, &tracer_vbo_);
		if (axis_vbo_) glDeleteBuffers(1, &axis_vbo_);
		if (corner_vbo_) glDeleteBuffers(1, &corner_vbo_);
		if (clip_vbo_) glDeleteBuffers(1, &clip_vbo_);
		if (gizmo_vbo_) glDeleteBuffers(1, &gizmo_vbo_);
		if (slice_vao_) glDeleteVertexArrays(1, &slice_vao_);
		if (box_vao_) glDeleteVertexArrays(1, &box_vao_);
		if (grid_vao_) glDeleteVertexArrays(1, &grid_vao_);
		if (mesh_vao_) glDeleteVertexArrays(1, &mesh_vao_);
		if (amr_debug_vao_) glDeleteVertexArrays(1, &amr_debug_vao_);
		if (eb_debug_vao_) glDeleteVertexArrays(1, &eb_debug_vao_);
		if (vox_vao_) glDeleteVertexArrays(1, &vox_vao_);
		if (arrow_vao_) glDeleteVertexArrays(1, &arrow_vao_);
		if (tracer_vao_) glDeleteVertexArrays(1, &tracer_vao_);
		if (axis_vao_) glDeleteVertexArrays(1, &axis_vao_);
		if (corner_vao_) glDeleteVertexArrays(1, &corner_vao_);
		if (clip_vao_) glDeleteVertexArrays(1, &clip_vao_);
		if (gizmo_vao_) glDeleteVertexArrays(1, &gizmo_vao_);
		doneCurrent();
	}

	void SliceViewer::setInfo(const SimInfo& info)
	{
		info_ = info;
		have_info_ = true;
		camera_.frameDomain((float)info.Lx, (float)info.Ly, (float)info.Lz);
		setDefaultPlane();
		updateRange();
		range_valid_ = false; // re-snap the auto-range on the new domain
		auto_speed_max_ = 1.8f * (float)info.U;
		geometry_dirty_ = true;
		axes_dirty_ = true; // rebuild the world triad + ticks for the new domain extents
		grid_dirty_ = true; // rebuild the grid overlay for the new domain (a graded grid re-pushes faces via setGridLines)
		// Base arrow length ~2% of the domain diagonal (a few cells) — visible over both the wake
		// slice and the model without cluttering. Re-seed the tracers for the new domain.
		float diag = std::sqrt((float)(info.Lx * info.Lx + info.Ly * info.Ly + info.Lz * info.Lz));
		arrow_len_ = std::max(4.0f * (float)info.h, 0.02f * diag);
		arrows_.reset();
		tracers_.reset(); // re-seed the streakline grid for the new domain
	}

	void SliceViewer::setReferenceU(double U)
	{
		if (!have_info_) return;
		info_.U = U;
		if (!auto_range_) updateRange(); // refresh the fixed per-field range to the new current
		update();
	}

	void SliceViewer::setDefaultPlane()
	{
		plane_frac_ = 0.5f;
	}

	void SliceViewer::setPlaneFraction(float frac)
	{
		plane_frac_ = std::min(1.0f, std::max(0.0f, frac));
		geometry_dirty_ = true;
		grid_dirty_ = true; // the grid overlay lives on the plane → moves with it
		update();
	}

	void SliceViewer::setGridLines(const std::vector<double>& xf, const std::vector<double>& yf, const std::vector<double>& zf)
	{
		grid_xf_.assign(xf.begin(), xf.end()); // double → float cell-face coordinates (metres)
		grid_yf_.assign(yf.begin(), yf.end());
		grid_zf_.assign(zf.begin(), zf.end());
		grid_dirty_ = true;
		update();
	}

	void SliceViewer::setMesh(paracfd::core::TriMesh mesh)
	{
		// Aerodynamic data belongs to a specific tessellation. Never let a newly loaded
		// STEP accidentally inherit colours merely because it has the same triangle count.
		clearTriangleSurfaceColouring();
		// Seed the gizmo transform (centre-on-bed default) from the model bbox NOW — host-side, so
		// modelPlacement() is valid even before the first paint (the CLI --load-step path voxelizes
		// immediately). The GL upload is still deferred to uploadMesh().
		gz_bbox_min_ = QVector3D(mesh.bbox_min[0], mesh.bbox_min[1], mesh.bbox_min[2]);
		gz_bbox_max_ = QVector3D(mesh.bbox_max[0], mesh.bbox_max[1], mesh.bbox_max[2]);
		computeDefaultXform();
		gz_drag_op_ = -1; gz_drag_axis_ = -1;
		fabric_bvh_dirty_ = true;
		gz_hover_op_ = -1; gz_hover_axis_ = -1;
		fabric_mesh_ = mesh;
		fabric_bvh_dirty_ = true;
		pending_mesh_ = std::move(mesh);
		mesh_upload_pending_ = true;
		update(); // uploaded on the next paint (main thread, context current)
	}

	void SliceViewer::clearMesh()
	{
		has_mesh_ = false;
		mesh_upload_pending_ = false;
		pending_mesh_ = paracfd::core::TriMesh{};
		fabric_mesh_ = paracfd::core::TriMesh{};
		clearTriangleSurfaceColouring();
		fabric_bvh_.reset();
		fabric_bvh_dirty_ = false;
		mesh_index_count_ = 0;
		mesh_override_ = false; // a fresh model (fluid viewer) is gizmo-editable again
		gz_valid_ = false;
		gz_drag_op_ = -1; gz_drag_axis_ = -1;
		gz_hover_op_ = -1; gz_hover_axis_ = -1;
		update();
	}

	void SliceViewer::setVoxelOverlay(const std::vector<unsigned char>& solid, paracfd::core::MacGrid grid)
	{
		pending_vox_solid_ = solid;
		pending_vox_grid_ = grid;
		vox_upload_pending_ = true;
		update(); // built on the next paint (main thread, context current)
	}

	void SliceViewer::clearVoxelOverlay()
	{
		has_vox_ = false;
		vox_upload_pending_ = false;
		pending_vox_solid_.clear();
		vox_vertex_count_ = 0;
		update();
	}

	void SliceViewer::setMeshTranslate(double tx, double ty, double tz)
	{
		mesh_override_mat_ = QMatrix4x4();
		mesh_override_mat_.translate((float)tx, (float)ty, (float)tz);
		mesh_override_ = true;
		fabric_bvh_dirty_ = true;
		update();
	}

	void SliceViewer::setMeshPlacement(const paracfd::core::ModelPlacement& p)
	{
		QMatrix4x4 m;
		m.setRow(0, QVector4D((float)p.m[0], (float)p.m[1], (float)p.m[2], (float)p.tx));
		m.setRow(1, QVector4D((float)p.m[3], (float)p.m[4], (float)p.m[5], (float)p.ty));
		m.setRow(2, QVector4D((float)p.m[6], (float)p.m[7], (float)p.m[8], (float)p.tz));
		m.setRow(3, QVector4D(0, 0, 0, 1));
		mesh_override_mat_ = m;
		mesh_override_ = true;
		fabric_bvh_dirty_ = true;
		update();
	}

	// --- Model-placement gizmo ---------------------------------------------------------------------

	QMatrix4x4 SliceViewer::modelMatrix() const
	{
		if (mesh_override_) return mesh_override_mat_;
		// world(v) = T(t) · R · S · T(-pivot) · v = t + rot·(scale ⊙ (v − pivot)).
		QMatrix4x4 M;
		M.translate(gz_t_);
		M.rotate(gz_rot_);
		M.scale(gz_scale_);
		M.translate(-gz_pivot_);
		return M;
	}

	paracfd::core::ModelPlacement SliceViewer::modelPlacement() const
	{
		paracfd::core::ModelPlacement p;
		const QMatrix4x4 M = modelMatrix();
		p.m[0] = M(0, 0); p.m[1] = M(0, 1); p.m[2] = M(0, 2);
		p.m[3] = M(1, 0); p.m[4] = M(1, 1); p.m[5] = M(1, 2);
		p.m[6] = M(2, 0); p.m[7] = M(2, 1); p.m[8] = M(2, 2);
		p.tx = M(0, 3); p.ty = M(1, 3); p.tz = M(2, 3);
		return p;
	}

	void SliceViewer::computeDefaultXform()
	{
		const QVector3D c = 0.5f * (gz_bbox_min_ + gz_bbox_max_);
		gz_pivot_ = c;
		gz_rot_ = QQuaternion();
		gz_scale_ = QVector3D(1, 1, 1);
		// Centre in x/y over the domain and drop the model's min-z vertex onto the bed (z=0). Falls back
		// to the model's own centre (no domain shift) when the domain is not yet known (Lx = min+max ⇒ the
		// x/y translate is a no-op), matching place_model_on_bed's fallback.
		const float Lx = have_info_ ? (float)info_.Lx : (gz_bbox_min_.x() + gz_bbox_max_.x());
		const float Ly = have_info_ ? (float)info_.Ly : (gz_bbox_min_.y() + gz_bbox_max_.y());
		gz_t_ = QVector3D(0.5f * Lx, 0.5f * Ly, c.z() - gz_bbox_min_.z());
		gz_valid_ = true;
		fabric_bvh_dirty_ = true;
	}

	void SliceViewer::setGizmoEnabled(bool on)
	{
		gizmo_on_ = on;
		gz_drag_op_ = -1; gz_drag_axis_ = -1;
		gz_hover_op_ = -1; gz_hover_axis_ = -1;
		update();
	}

	void SliceViewer::resetModelPlacement()
	{
		if (!gz_valid_) return;
		computeDefaultXform();
		gz_drag_op_ = -1; gz_drag_axis_ = -1;
		update();
		emit modelPlacementChanged();
	}

	SliceViewer::ModelGizmoXform SliceViewer::modelXform() const
	{
		ModelGizmoXform x;
		x.pivot = gz_pivot_; x.t = gz_t_; x.scale = gz_scale_; x.rot = gz_rot_; x.valid = gz_valid_;
		return x;
	}

	void SliceViewer::setModelXform(const ModelGizmoXform& x)
	{
		if (!x.valid) return;
		gz_pivot_ = x.pivot; gz_t_ = x.t; gz_scale_ = x.scale; gz_rot_ = x.rot; gz_valid_ = true;
		mesh_override_ = false; // an explicit gizmo state supersedes any prior override
		fabric_bvh_dirty_ = true;
		update();
	}

	void SliceViewer::setModelPlacement(const paracfd::core::ModelPlacement& p)
	{
		if (!gz_valid_) return; // need a model (its pivot) first
		// Decompose the affine linear part M = R·S: each column col_i = M·e_i = scale_i · (R's i-th column),
		// so scale_i = |col_i| and R = [col_i / scale_i]. Exact for the R·S placements the gizmo produces.
		const QVector3D c0((float)p.m[0], (float)p.m[3], (float)p.m[6]);
		const QVector3D c1((float)p.m[1], (float)p.m[4], (float)p.m[7]);
		const QVector3D c2((float)p.m[2], (float)p.m[5], (float)p.m[8]);
		const float s0 = c0.length(), s1 = c1.length(), s2 = c2.length();
		gz_scale_ = QVector3D(s0 > 1e-9f ? s0 : 1.0f, s1 > 1e-9f ? s1 : 1.0f, s2 > 1e-9f ? s2 : 1.0f);
		const QVector3D r0 = s0 > 1e-9f ? c0 / s0 : QVector3D(1, 0, 0);
		const QVector3D r1 = s1 > 1e-9f ? c1 / s1 : QVector3D(0, 1, 0);
		const QVector3D r2 = s2 > 1e-9f ? c2 / s2 : QVector3D(0, 0, 1);
		const float rd[9] = { r0.x(), r1.x(), r2.x(), r0.y(), r1.y(), r2.y(), r0.z(), r1.z(), r2.z() };
		gz_rot_ = QQuaternion::fromRotationMatrix(QMatrix3x3(rd));
		// world(v) = M·v + t; we want it as gz_t_ + R·S·(v − pivot) ⇒ gz_t_ = M·pivot + t.
		const QVector3D Mpiv(
			(float)(p.m[0] * gz_pivot_.x() + p.m[1] * gz_pivot_.y() + p.m[2] * gz_pivot_.z()),
			(float)(p.m[3] * gz_pivot_.x() + p.m[4] * gz_pivot_.y() + p.m[5] * gz_pivot_.z()),
			(float)(p.m[6] * gz_pivot_.x() + p.m[7] * gz_pivot_.y() + p.m[8] * gz_pivot_.z()));
		gz_t_ = QVector3D((float)p.tx, (float)p.ty, (float)p.tz) + Mpiv;
		mesh_override_ = false;
		fabric_bvh_dirty_ = true;
		update();
		emit modelPlacementChanged();
	}

	void SliceViewer::setTriangleSurfaceCp(const std::vector<float>& plus_values,
		const std::vector<float>& minus_values, float range_min, float range_max)
	{
		if (!(range_max > range_min)) { range_min = -1.0f; range_max = 1.0f; }
		if (plus_values.size() != minus_values.size())
		{
			clearTriangleSurfaceColouring();
			return;
		}
		mesh_triangle_colours_.resize(plus_values.size() * 8);
		const float span = range_max - range_min;
		for (std::size_t triangle = 0; triangle < plus_values.size(); ++triangle)
		{
			for (std::size_t side = 0; side < 2; ++side)
			{
				float r = 0.5f, g = 0.5f, b = 0.5f;
				const float cp = side == 0 ? plus_values[triangle] : minus_values[triangle];
				if (std::isfinite(cp)) scour_colormap((cp - range_min) / span, r, g, b);
				const std::size_t offset = 8 * triangle + 4 * side;
				mesh_triangle_colours_[offset + 0] = r;
				mesh_triangle_colours_[offset + 1] = g;
				mesh_triangle_colours_[offset + 2] = b;
				mesh_triangle_colours_[offset + 3] = 1.0f;
			}
		}
		mesh_triangle_colour_upload_pending_ = true;
		update();
	}

	void SliceViewer::clearTriangleSurfaceColouring()
	{
		mesh_triangle_colours_.clear();
		mesh_triangle_colour_upload_pending_ = true;
		has_mesh_triangle_colours_ = false;
		update();
	}

	void SliceViewer::setParagliderDebugBoxes(const std::vector<std::array<float,6>>& amr_bricks,
		const std::vector<std::array<float,6>>& eb_cells)
	{
		pending_amr_debug_lines_=debug_box_lines(amr_bricks);pending_eb_debug_lines_=debug_box_lines(eb_cells);
		paraglider_debug_upload_pending_=true;update();
	}

	void SliceViewer::clearParagliderDebugBoxes()
	{
		pending_amr_debug_lines_.clear();pending_eb_debug_lines_.clear();amr_debug_vertex_count_=0;eb_debug_vertex_count_=0;
		paraglider_debug_upload_pending_=true;update();
	}

	void SliceViewer::ensureFabricBvh()
	{
		if (!fabric_bvh_dirty_) return;
		fabric_bvh_dirty_ = false;
		if (fabric_mesh_.empty()) { fabric_bvh_.reset(); return; }
		paracfd::core::TriMesh placed = paracfd::core::placed_mesh(fabric_mesh_, modelPlacement());
		fabric_bvh_ = std::make_unique<paracfd::core::TriangleBvh>(placed);
	}

	QVector3D SliceViewer::gizmoAxisDir(int a) const
	{
		const QVector3D e = (a == 0) ? QVector3D(1, 0, 0) : (a == 1) ? QVector3D(0, 1, 0) : QVector3D(0, 0, 1);
		return gz_rot_.rotatedVector(e).normalized();
	}

	float SliceViewer::gizmoSize() const { return (0.16f / 3.0f) * camera_.distance(); } // 1/3 of the old size

	void SliceViewer::mouseRay(const QPointF& pos, QVector3D& orig, QVector3D& dir) const
	{
		const QMatrix4x4 inv = camera_.viewProjection().inverted();
		const float nx = (float)(pos.x() / std::max(1, width())) * 2.0f - 1.0f;
		const float ny = 1.0f - (float)(pos.y() / std::max(1, height())) * 2.0f;
		QVector4D pn = inv * QVector4D(nx, ny, -1.0f, 1.0f);
		QVector4D pf = inv * QVector4D(nx, ny, 1.0f, 1.0f);
		if (pn.w() != 0.0f) pn /= pn.w();
		if (pf.w() != 0.0f) pf /= pf.w();
		orig = pn.toVector3D();
		const QVector3D far = pf.toVector3D();
		QVector3D d = far - orig;
		const float L = d.length();
		dir = (L > 1e-9f) ? d / L : QVector3D(0, 0, 1);
	}

	namespace
	{
		// Screen-space distance from point p to segment [a,b] (all in pixels).
		float dist_to_seg(const QPointF& p, const QPointF& a, const QPointF& b)
		{
			const QPointF ab = b - a, ap = p - a;
			const double L2 = ab.x() * ab.x() + ab.y() * ab.y();
			double t = (L2 > 1e-6) ? (ap.x() * ab.x() + ap.y() * ab.y()) / L2 : 0.0;
			t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
			const QPointF q = a + t * ab, d = p - q;
			return (float)std::sqrt(d.x() * d.x() + d.y() * d.y());
		}
		// A stable in-plane orthonormal basis (u,v) for a plane with unit normal n, with u × v = n so the
		// signed angle atan2(rel·v, rel·u) increases right-handed about +n.
		void plane_basis(const QVector3D& n, QVector3D& u, QVector3D& v)
		{
			const QVector3D ref = (std::fabs(n.z()) < 0.9f) ? QVector3D(0, 0, 1) : QVector3D(1, 0, 0);
			u = QVector3D::crossProduct(n, ref);
			const float L = u.length();
			u = (L > 1e-6f) ? u / L : QVector3D(1, 0, 0);
			v = QVector3D::crossProduct(n, u);
		}
		constexpr double kTwoPi = 6.283185307179586;
	}

	// Unified gizmo handle layout (all in the model's local/rotated frame, scaled by g = gizmoSize()):
	//   translate arrow shaft [0, kTrans·g] + head; scale cube at kScale·g; rotate ring radius kRing·g;
	//   uniform-scale cube at the centre. The radii are chosen so the three per-axis handles don't overlap
	//   (arrow tip < scale cube < ring radius laterally) ⇒ nearest-in-screen picking disambiguates cleanly.
	namespace
	{
		constexpr float kTransLen = 0.70f; // translate arrow tip distance (·g)
		constexpr float kRingRad = 1.00f;  // rotate ring radius (·g)
		constexpr float kScalePos = 1.25f; // scale cube distance along the axis (·g) — OUTSIDE the ring so the
		                                   // enlarged, easy-to-grab cube is clear of the arrow + ring handles
		constexpr float kCubeHalf = 0.075f;// scale-cube half-size (·g) — 1.5× the rest, for an easier pick target
	}

	bool SliceViewer::pickGizmo(const QPointF& pos, int& op, int& axis) const
	{
		op = -1; axis = -1;
		if (!gizmoEditable() || !gizmo_on_) return false;
		const QMatrix4x4 mvp = camera_.viewProjection();
		QPointF O;
		if (!projectPoint(mvp, gz_t_, O)) return false;
		const float g = gizmoSize();

		// Tier 1 — compact POINT handles (uniform-scale centre + per-axis scale cubes): precise targets
		// that take priority so they win over the translate shaft / rotate rings passing nearby. The scale
		// handles use a 1.5× pick radius (they are drawn 1.5× larger and sit clear of the other handles).
		{
			float best = 12.0f; int bop = -1, bax = -1;
			const QPointF dc = pos - O;
			const float rc = (float)std::sqrt(dc.x() * dc.x() + dc.y() * dc.y());
			if (rc < 10.0f && rc < best) { best = rc; bop = 3; bax = -1; }
			for (int a = 0; a < 3; ++a)
			{
				QPointF sp;
				if (!projectPoint(mvp, gz_t_ + gizmoAxisDir(a) * (kScalePos * g), sp)) continue;
				const QPointF d = pos - sp;
				const float r = (float)std::sqrt(d.x() * d.x() + d.y() * d.y());
				if (r < 12.0f && r < best) { best = r; bop = 2; bax = a; }
			}
			if (bop >= 0) { op = bop; axis = bax; return true; }
		}

		// Tier 2 — EXTENDED handles (translate arrows + rotate rings): nearest wins within the pick radius.
		float best = 11.0f; int bop = -1, bax = -1;
		for (int a = 0; a < 3; ++a) // translate: segment [O, arrow tip]
		{
			QPointF tip;
			if (!projectPoint(mvp, gz_t_ + gizmoAxisDir(a) * (kTransLen * g), tip)) continue;
			const float d = dist_to_seg(pos, O, tip);
			if (d < best) { best = d; bop = 0; bax = a; }
		}
		const int N = 48;
		for (int a = 0; a < 3; ++a) // rotate: ring in the plane of the other two axes
		{
			const QVector3D u = gizmoAxisDir((a + 1) % 3), v = gizmoAxisDir((a + 2) % 3);
			QPointF prev; bool hp = false;
			for (int i = 0; i <= N; ++i)
			{
				const float th = (float)(kTwoPi * i / N);
				QPointF s;
				if (!projectPoint(mvp, gz_t_ + (std::cos(th) * u + std::sin(th) * v) * (kRingRad * g), s)) { hp = false; continue; }
				if (hp) { const float d = dist_to_seg(pos, prev, s); if (d < best) { best = d; bop = 1; bax = a; } }
				prev = s; hp = true;
			}
		}
		op = bop; axis = bax;
		return op >= 0;
	}

	void SliceViewer::beginGizmoDrag(int op, int axis, const QPointF& pos)
	{
		gz_drag_op_ = op;
		gz_drag_axis_ = axis;
		gz_drag_n_ = (axis >= 0 && axis < 3) ? gizmoAxisDir(axis) : QVector3D(0, 0, 1); // pinned for the drag
		if (op == 1 && axis >= 0 && axis < 3) // rotate: seed the plane angle
		{
			QVector3D o, d;
			mouseRay(pos, o, d);
			const float dn = QVector3D::dotProduct(d, gz_drag_n_);
			if (std::fabs(dn) > 1e-5f)
			{
				QVector3D u, v;
				plane_basis(gz_drag_n_, u, v);
				const float t = QVector3D::dotProduct(gz_t_ - o, gz_drag_n_) / dn;
				const QVector3D rel = (o + d * t) - gz_t_;
				gz_rot_last_ang_ = std::atan2(QVector3D::dotProduct(rel, v), QVector3D::dotProduct(rel, u));
			}
			else
				gz_rot_last_ang_ = 0.0f;
		}
	}

	void SliceViewer::dragGizmo(const QPointF& cur, const QPointF& prev)
	{
		if (gz_drag_op_ < 0 || !gz_valid_) return;
		const QMatrix4x4 mvp = camera_.viewProjection();
		QPointF O;
		const bool haveO = projectPoint(mvp, gz_t_, O);
		const int a = gz_drag_axis_;

		if (gz_drag_op_ == 0) // translate along the pinned axis
		{
			QPointF tipS;
			if (!haveO || !projectPoint(mvp, gz_t_ + gz_drag_n_, tipS)) return; // project 1 m along the axis
			const QPointF sd = tipS - O;
			const double len2 = sd.x() * sd.x() + sd.y() * sd.y();
			if (len2 < 1e-3) return; // axis nearly parallel to the view — ill-conditioned
			const QPointF dp = cur - prev;
			const float move = (float)((dp.x() * sd.x() + dp.y() * sd.y()) / len2);
			gz_t_ += gz_drag_n_ * move;
		}
		else if (gz_drag_op_ == 1) // rotate about the pinned axis (ray → rotation plane → signed angle)
		{
			QVector3D u, v;
			plane_basis(gz_drag_n_, u, v);
			QVector3D o, d;
			mouseRay(cur, o, d);
			const float dn = QVector3D::dotProduct(d, gz_drag_n_);
			if (std::fabs(dn) < 1e-5f) return;
			const float t = QVector3D::dotProduct(gz_t_ - o, gz_drag_n_) / dn;
			const QVector3D rel = (o + d * t) - gz_t_;
			float ang = std::atan2(QVector3D::dotProduct(rel, v), QVector3D::dotProduct(rel, u));
			float dAng = ang - gz_rot_last_ang_;
			while (dAng > 3.14159265f) dAng -= (float)kTwoPi;
			while (dAng < -3.14159265f) dAng += (float)kTwoPi;
			gz_rot_last_ang_ = ang;
			gz_rot_ = QQuaternion::fromAxisAndAngle(gz_drag_n_, dAng * float(180.0 / 3.141592653589793)) * gz_rot_;
			gz_rot_.normalize();
		}
		else if (gz_drag_op_ == 2) // scale along the pinned axis (motion projected onto the screen axis)
		{
			QPointF tipS;
			if (!haveO || !projectPoint(mvp, gz_t_ + gz_drag_n_, tipS)) return;
			const QPointF sd = tipS - O;
			const float len = (float)std::sqrt(sd.x() * sd.x() + sd.y() * sd.y());
			if (len < 1e-2f) return;
			const QPointF dp = cur - prev;
			const float along = (float)((dp.x() * sd.x() + dp.y() * sd.y()) / len);
			const float f = std::exp(along * 0.012f);
			QVector3D s = gz_scale_;
			if (a == 0) s.setX(s.x() * f); else if (a == 1) s.setY(s.y() * f); else s.setZ(s.z() * f);
			gz_scale_ = s;
			gz_scale_ = QVector3D(std::clamp(gz_scale_.x(), 0.02f, 50.0f),
				std::clamp(gz_scale_.y(), 0.02f, 50.0f), std::clamp(gz_scale_.z(), 0.02f, 50.0f));
		}
		else if (gz_drag_op_ == 3) // uniform scale (radial distance from the centre)
		{
			if (!haveO) return;
			const QPointF rp = prev - O, rc = cur - O;
			float lp = (float)std::sqrt(rp.x() * rp.x() + rp.y() * rp.y());
			const float lc = (float)std::sqrt(rc.x() * rc.x() + rc.y() * rc.y());
			if (lp < 3.0f) lp = 3.0f;
			gz_scale_ *= (lc / lp);
			gz_scale_ = QVector3D(std::clamp(gz_scale_.x(), 0.02f, 50.0f),
				std::clamp(gz_scale_.y(), 0.02f, 50.0f), std::clamp(gz_scale_.z(), 0.02f, 50.0f));
		}
	}

	void SliceViewer::drawGizmo(const QMatrix4x4& mvp)
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gizmoEditable() || !gizmo_on_) return;
		if (!gizmo_vao_) { glGenVertexArrays(1, &gizmo_vao_); glGenBuffers(1, &gizmo_vbo_); }

		const float g = gizmoSize();
		const QVector3D O = gz_t_;
		std::vector<float> all;
		struct Blk { int off, cnt; QVector4D col; bool active; };
		std::vector<Blk> blks;

		auto push = [&](const QVector3D& p0, const QVector3D& p1) {
			all.push_back(p0.x()); all.push_back(p0.y()); all.push_back(p0.z());
			all.push_back(p1.x()); all.push_back(p1.y()); all.push_back(p1.z());
		};
		auto beginBlk = [&]() { return (int)all.size() / 3; };
		auto emitBlk = [&](int off, const QVector4D& col, bool active) {
			const int cnt = (int)all.size() / 3 - off;
			if (cnt > 0) blks.push_back({ off, cnt, col, active });
		};
		static const int kCubeE[12][2] = { {0,1},{0,2},{0,4},{1,3},{1,5},{2,3},{2,6},{3,7},{4,5},{4,6},{5,7},{6,7} };
		auto pushCube = [&](const QVector3D& c, const QVector3D& ax, const QVector3D& ay, const QVector3D& az, float s) {
			QVector3D v[8];
			for (int i = 0; i < 8; ++i)
				v[i] = c + ax * ((i & 1) ? s : -s) + ay * ((i & 2) ? s : -s) + az * ((i & 4) ? s : -s);
			for (auto& e : kCubeE) push(v[e[0]], v[e[1]]);
		};

		const QVector4D axcol[3] = {
			QVector4D(0.90f, 0.28f, 0.24f, 1.0f), QVector4D(0.30f, 0.78f, 0.32f, 1.0f), QVector4D(0.32f, 0.55f, 0.95f, 1.0f) };
		const QVector4D grey(0.85f, 0.85f, 0.88f, 1.0f);
		// A handle is "active" when it is being dragged, or hovered while nothing is dragged.
		auto hot = [&](int op, int ax) { return (gz_drag_op_ == op && gz_drag_axis_ == ax) || (gz_drag_op_ < 0 && gz_hover_op_ == op && gz_hover_axis_ == ax); };

		for (int a = 0; a < 3; ++a)
		{
			const QVector3D dir = gizmoAxisDir(a), u = gizmoAxisDir((a + 1) % 3), v = gizmoAxisDir((a + 2) % 3);
			// translate arrow (shaft + arrowhead)
			{
				const int off = beginBlk();
				const QVector3D tip = O + dir * (kTransLen * g);
				push(O, tip);
				const QVector3D base = tip - dir * (0.14f * g), p1 = u * (0.05f * g), p2 = v * (0.05f * g);
				push(tip, base + p1); push(tip, base - p1); push(tip, base + p2); push(tip, base - p2);
				emitBlk(off, axcol[a], hot(0, a));
			}
			// rotate ring
			{
				const int off = beginBlk();
				const int N = 48;
				QVector3D prev; bool hp = false;
				for (int i = 0; i <= N; ++i)
				{
					const float th = (float)(kTwoPi * i / N);
					const QVector3D w = O + (std::cos(th) * u + std::sin(th) * v) * (kRingRad * g);
					if (hp) push(prev, w);
					prev = w; hp = true;
				}
				emitBlk(off, axcol[a], hot(1, a));
			}
			// scale cube (on the axis, past the arrowhead)
			{
				const int off = beginBlk();
				pushCube(O + dir * (kScalePos * g), dir, u, v, kCubeHalf * g);
				emitBlk(off, axcol[a], hot(2, a));
			}
		}
		// uniform-scale centre cube (world-aligned; 1.5× like the per-axis scale cubes)
		{
			const int off = beginBlk();
			pushCube(O, QVector3D(1, 0, 0), QVector3D(0, 1, 0), QVector3D(0, 0, 1), 0.083f * g);
			emitBlk(off, grey, (gz_drag_op_ == 3) || (gz_drag_op_ < 0 && gz_hover_op_ == 3));
		}
		if (all.empty()) return;

		glBindVertexArray(gizmo_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, gizmo_vbo_);
		glBufferData(GL_ARRAY_BUFFER, all.size() * sizeof(float), all.data(), GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);

		const QVector4D hotcol(1.0f, 0.85f, 0.20f, 1.0f);
		glDisable(GL_DEPTH_TEST); // manipulator always on top so it stays grabbable
		glDisable(GL_CLIP_DISTANCE0);
		prog_.bind();
		prog_.setUniformValue("uMVP", mvp);
		prog_.setUniformValue("uFlat", 1);
		prog_.setUniformValue("uClipPlane", QVector4D(0, 0, 0, 1));
		for (const Blk& b : blks)
		{
			glLineWidth(b.active ? 3.4f : 1.8f);
			prog_.setUniformValue("uColor", b.active ? hotcol : b.col);
			glDrawArrays(GL_LINES, b.off, b.cnt);
		}
		glBindVertexArray(0);
		prog_.release();
		glLineWidth(1.0f);
		glEnable(GL_DEPTH_TEST);
	}

	void SliceViewer::updateRange()
	{
		float U = have_info_ ? (float)info_.U : 1.0f;
		float rhoU2 = have_info_ ? (float)(info_.rho * info_.U * info_.U) : 1.0f;
		switch (field_)
		{
		case Field::SpeedMag: vmin_ = 0.0f; vmax_ = 1.8f * U; break;
		case Field::VelU: vmin_ = -0.5f * U; vmax_ = 1.7f * U; break;
		case Field::VelV: vmin_ = -0.9f * U; vmax_ = 0.9f * U; break;
		case Field::VelW: vmin_ = -0.9f * U; vmax_ = 0.9f * U; break;
		case Field::Pressure: vmin_ = -0.7f * rhoU2; vmax_ = 0.7f * rhoU2; break;
		}
		if (vmax_ <= vmin_) vmax_ = vmin_ + 1.0f;
	}

	void SliceViewer::setAutoRange(bool on)
	{
		auto_range_ = on;
		range_valid_ = false;     // re-snap when toggled back on
		if (!on) updateRange();   // restore the fixed per-field colour range
		update();
	}

	// Fold a live GPU reduction into the displayed range. The scale EXPANDS instantly (so an
	// anomalous fast region / blow-up lights up the very next update — never hidden or capped) and
	// CONTRACTS slowly (a light EMA) so a settled field's legend does not flicker. Magnitude fields
	// pin vmin at 0; signed fields (u,v,w,p) track both data ends. auto_speed_max_ (arrow colour
	// scale) tracks the domain's max |u| the same way, so slice + legend + arrows stay consistent.
	void SliceViewer::applyAutoRange(const FieldRange& fr)
	{
		if (!fr.valid || !std::isfinite(fr.field_min) || !std::isfinite(fr.field_max)) return;

		float tmin = (field_ == Field::SpeedMag) ? 0.0f : fr.field_min;
		float tmax = fr.field_max;
		if (tmax - tmin < 1e-6f) tmax = tmin + 1e-6f;
		float tspd = std::isfinite(fr.speed_max) ? std::max(fr.speed_max, 1e-3f) : auto_speed_max_;

		if (!range_valid_)
		{
			vmin_ = tmin; vmax_ = tmax; auto_speed_max_ = tspd; // first sample: snap directly
			range_valid_ = true;
		}
		else
		{
			const float shrink = 0.05f; // ~2 s contraction time constant at the ~12 Hz update rate
			vmax_ = (tmax >= vmax_) ? tmax : vmax_ + (tmax - vmax_) * shrink;
			vmin_ = (tmin <= vmin_) ? tmin : vmin_ + (tmin - vmin_) * shrink;
			auto_speed_max_ = (tspd >= auto_speed_max_) ? tspd : auto_speed_max_ + (tspd - auto_speed_max_) * shrink;
			if (vmax_ - vmin_ < 1e-6f) vmax_ = vmin_ + 1e-6f;
		}

		// Diagnostic: prove the range is tracking the live data (throttled to ~every 2 s).
		if (++range_log_ctr_ >= 24)
		{
			range_log_ctr_ = 0;
			const char* fn = (field_ == Field::SpeedMag) ? "|u|" : (field_ == Field::VelU) ? "u"
				: (field_ == Field::VelV) ? "v" : (field_ == Field::VelW) ? "w" : "p";
			std::fprintf(stderr, "[G1] auto-range %s: vmin=%.4g vmax=%.4g (raw[%.4g,%.4g]) speed_max=%.4g\n",
				fn, vmin_, vmax_, fr.field_min, fr.field_max, auto_speed_max_);
		}
	}

	SliceParams SliceViewer::currentParams() const
	{
		SliceParams sp;
		sp.grid.nx = info_.nx; sp.grid.ny = info_.ny; sp.grid.nz = info_.nz; sp.grid.h = info_.h;
		// TRUE domain extent (graded: xf[nx] via info_.Lx; uniform: nx·h). Positions the slice quad correctly
		// on a graded grid without dereferencing metric pointers here (the host geometry path has none). The
		// device sampler additionally gets the metric ARRAYS via worker_->displayGrid() (see paintGL).
		sp.Lx = (float)info_.Lx; sp.Ly = (float)info_.Ly; sp.Lz = (float)info_.Lz;
		sp.axis = axis_;
		float L = (axis_ == Axis::X) ? (float)info_.Lx : (axis_ == Axis::Y) ? (float)info_.Ly : (float)info_.Lz;
		sp.plane_pos = plane_frac_ * L;
		sp.nu = NRES; sp.nv = NRES;
		sp.field = field_;
		sp.vmin = vmin_; sp.vmax = vmax_;
		return sp;
	}

	void SliceViewer::initializeGL()
	{
		PARACFD_ASSERT_GL_THREAD();
		initializeOpenGLFunctions();

		const char* ver = reinterpret_cast<const char*>(glGetString(GL_VERSION));
		std::fprintf(stderr, "[G1] GL context: %s\n", ver ? ver : "(null)");

		glClearColor(0.09f, 0.10f, 0.12f, 1.0f);
		glEnable(GL_DEPTH_TEST);

		if (!prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kVert)
			|| !prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kFrag)
			|| !prog_.link())
		{
			std::fprintf(stderr, "[G1] shader error: %s\n", prog_.log().toUtf8().constData());
		}

		if (!mesh_prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kMeshVert)
			|| !mesh_prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kMeshFrag)
			|| !mesh_prog_.link())
		{
			std::fprintf(stderr, "[G1] mesh shader error: %s\n", mesh_prog_.log().toUtf8().constData());
		}

		if (!arrow_prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kArrowVert)
			|| !arrow_prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kArrowFrag)
			|| !arrow_prog_.link())
		{
			std::fprintf(stderr, "[G1] arrow shader error: %s\n", arrow_prog_.log().toUtf8().constData());
		}

		if (!tracer_prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kTracerVert)
			|| !tracer_prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kTracerFrag)
			|| !tracer_prog_.link())
		{
			std::fprintf(stderr, "[G1] tracer shader error: %s\n", tracer_prog_.log().toUtf8().constData());
		}

		glGenVertexArrays(1, &slice_vao_);
		glGenVertexArrays(1, &box_vao_);
		glGenVertexArrays(1, &grid_vao_);
		glGenVertexArrays(1, &mesh_vao_);
		glGenVertexArrays(1, &amr_debug_vao_);
		glGenVertexArrays(1, &eb_debug_vao_);
		glGenVertexArrays(1, &vox_vao_);
		glGenVertexArrays(1, &arrow_vao_);
		glGenVertexArrays(1, &tracer_vao_);
		glGenVertexArrays(1, &axis_vao_);
		glGenVertexArrays(1, &corner_vao_);
		glGenVertexArrays(1, &clip_vao_);
		glGenBuffers(1, &pos_vbo_);
		glGenBuffers(1, &color_vbo_);
		glGenBuffers(1, &idx_ebo_);
		glGenBuffers(1, &box_vbo_);
		glGenBuffers(1, &grid_vbo_);
		glGenBuffers(1, &mesh_pos_vbo_);
		glGenBuffers(1, &mesh_norm_vbo_);
		glGenBuffers(1, &mesh_idx_ebo_);
		glGenBuffers(1, &mesh_triangle_colour_ssbo_);
		glGenBuffers(1, &amr_debug_vbo_);
		glGenBuffers(1, &eb_debug_vbo_);
		glGenBuffers(1, &vox_pos_vbo_);
		glGenBuffers(1, &vox_norm_vbo_);
		glGenBuffers(1, &vox_color_vbo_);
		glGenBuffers(1, &arrow_glyph_vbo_);
		glGenBuffers(1, &arrow_inst_vbo_);
		glGenBuffers(1, &tracer_vbo_);
		glGenBuffers(1, &axis_vbo_);
		glGenBuffers(1, &corner_vbo_);
		glGenBuffers(1, &clip_vbo_);

		// Static index buffer for the (NRES x NRES) grid triangulation (axis-independent).
		std::vector<unsigned int> idx;
		idx.reserve((NRES - 1) * (NRES - 1) * 6);
		for (int b = 0; b < NRES - 1; ++b)
			for (int a = 0; a < NRES - 1; ++a)
			{
				unsigned int v00 = b * NRES + a, v10 = v00 + 1, v01 = v00 + NRES, v11 = v01 + 1;
				idx.push_back(v00); idx.push_back(v10); idx.push_back(v11);
				idx.push_back(v00); idx.push_back(v11); idx.push_back(v01);
			}
		index_count_ = (int)idx.size();
		glBindVertexArray(slice_vao_);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, idx_ebo_);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx.size() * sizeof(unsigned int), idx.data(), GL_STATIC_DRAW);

		// Colour VBO: NRES*NRES float4, written by CUDA every frame. Allocate then register.
		glBindBuffer(GL_ARRAY_BUFFER, color_vbo_);
		glBufferData(GL_ARRAY_BUFFER, (size_t)NRES * NRES * 4 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(1);

		// Position VBO wired into the same VAO (filled in buildSliceGeometry()).
		glBindBuffer(GL_ARRAY_BUFFER, pos_vbo_);
		glBufferData(GL_ARRAY_BUFFER, (size_t)NRES * NRES * 3 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);
		glBindVertexArray(0);

		cuda_res_ = slice_gl_register(color_vbo_);
		if (!cuda_res_)
			std::fprintf(stderr, "[G1] WARNING: CUDA-GL interop registration failed; slice will be blank.\n");

		gl_ready_ = true;
		buildBoxGeometry();
		buildArrowGlyph();
		buildTracerBuffers();
		buildCornerGizmo();
		if (have_info_) { buildSliceGeometry(); buildAxesGeometry(); }
		fps_timer_.start();
	}

	void SliceViewer::buildArrowGlyph()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_) return;
		// Unit arrow in the local (x=along-shaft, y=lateral) frame: a thin shaft rectangle plus a
		// triangular head, filled triangles so it reads without relying on GL line width.
		const float g[] = {
			// shaft (two triangles)
			0.00f, -0.05f, 0.70f, -0.05f, 0.70f, 0.05f,
			0.00f, -0.05f, 0.70f, 0.05f, 0.00f, 0.05f,
			// head (one triangle)
			0.62f, -0.16f, 1.00f, 0.00f, 0.62f, 0.16f
		};
		arrow_glyph_verts_ = (int)(sizeof(g) / sizeof(float) / 2);
		glBindVertexArray(arrow_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, arrow_glyph_vbo_);
		glBufferData(GL_ARRAY_BUFFER, sizeof(g), g, GL_STATIC_DRAW);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);
		// Per-instance attributes (pos3,dir3,speed,alpha) = 8 floats, divisor 1. Buffer sized on demand.
		glBindBuffer(GL_ARRAY_BUFFER, arrow_inst_vbo_);
		const int stride = 8 * sizeof(float);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)0);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
		glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, stride, (void*)(6 * sizeof(float)));
		glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, stride, (void*)(7 * sizeof(float)));
		glEnableVertexAttribArray(1);
		glEnableVertexAttribArray(2);
		glEnableVertexAttribArray(3);
		glEnableVertexAttribArray(4);
		glVertexAttribDivisor(1, 1);
		glVertexAttribDivisor(2, 1);
		glVertexAttribDivisor(3, 1);
		glVertexAttribDivisor(4, 1);
		glBindVertexArray(0);
	}

	void SliceViewer::updateArrows()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_ || !show_arrows_ || (!worker_ && !paraglider_worker_) || !have_info_)
		{
			arrow_draw_count_ = 0;
			return;
		}
		float dt = 0.0f;
		if (!arrow_clock_started_) { arrow_clock_.start(); arrow_clock_started_ = true; }
		else dt = (float)(arrow_clock_.restart() / 1000.0);

		ArrowView view;
		view.three_d = arrow_3d_;
		view.axis = (int)axis_;
		{
			float L = (axis_ == Axis::X) ? (float)info_.Lx : (axis_ == Axis::Y) ? (float)info_.Ly : (float)info_.Lz;
			view.plane_pos = plane_frac_ * L;
		}
		// Colour + normalise arrows by |vel| over the live speed max (auto-range) so an arrow's
		// colour reads on the same scale as the |u| slice/legend; falls back to a fixed 1.8U when
		// auto-range is off. Advection is a modest multiple of real time so the motion is legible on
		// the larger domains without desyncing from the physics.
		view.speed_scale = auto_range_ ? std::max(auto_speed_max_, 1e-3f)
			: 1.8f * (have_info_ ? (float)info_.U : 1.0f);
		view.gain = 1.5f * arrow_speed_mult_; // user speed multiplier [0,1] (visual only)

		ensureFabricBvh();
		auto consume = [&](const FlowField& f) { FlowField ff = f; ff.fabric = fabric_bvh_.get(); arrows_.advance(dt, view, ff); };
		const bool got = paraglider_worker_ ? paraglider_worker_->withFlowField(consume)
			: worker_->withFlowField(consume);
		if (!got) { arrow_draw_count_ = 0; return; }

		const std::vector<float>& inst = arrows_.instance_data();
		int n = arrows_.count();
		arrow_draw_count_ = n;
		if (n == 0) return;

		if (!arrows_logged_)
		{
			arrows_logged_ = true;
			std::fprintf(stderr, "[G1] flow arrows live: %d particles, %s mode\n", n, arrow_3d_ ? "3D" : "2D");
		}

		glBindVertexArray(arrow_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, arrow_inst_vbo_);
		if (n > arrow_inst_capacity_)
		{
			glBufferData(GL_ARRAY_BUFFER, (size_t)n * 8 * sizeof(float), inst.data(), GL_DYNAMIC_DRAW);
			arrow_inst_capacity_ = n;
		}
		else
		{
			glBufferSubData(GL_ARRAY_BUFFER, 0, (size_t)n * 8 * sizeof(float), inst.data());
		}
		glBindVertexArray(0);
	}

	void SliceViewer::drawArrows(const QMatrix4x4& mvp)
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!show_arrows_ || arrow_draw_count_ <= 0) return;
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDepthMask(GL_FALSE); // blend arrows without fighting each other's depth; still occluded by solids
		arrow_prog_.bind();
		arrow_prog_.setUniformValue("uMVP", mvp);
		arrow_prog_.setUniformValue("uEye", camera_.eye());
		arrow_prog_.setUniformValue("uLen", arrow_len_);
		arrow_prog_.setUniformValue("uWidth", arrow_width_);
		// Over a visible colour slice the speed-coloured arrows blend into the field, so draw them
		// solid WHITE for contrast (reads over any colormap stop); when the slice is hidden, colour
		// them by speed (matches the legend) so they stay informative against the dark background.
		arrow_prog_.setUniformValue("uColorMode", show_slice_ ? 1 : 0);
		arrow_prog_.setUniformValue("uSolidColor", QVector3D(1.0f, 1.0f, 1.0f));
		glBindVertexArray(arrow_vao_);
		glDrawArraysInstanced(GL_TRIANGLES, 0, arrow_glyph_verts_, arrow_draw_count_);
		glBindVertexArray(0);
		arrow_prog_.release();
		glDepthMask(GL_TRUE);
		glDisable(GL_BLEND);
	}

	// --- Flow tracers (inlet-seeded streamlines, drawn as camera-facing ribbons) -------------------
	// The streamlines are drawn through their own tracer shader as triangle-strip RIBBONS: an
	// interleaved VBO of [pos.xyz side | tan.xyz | rgba] (stride 11 floats). buildTracerBuffers()
	// wires that layout onto the tracer VAO; the VBO is (re)sized on demand in updateTracers().
	void SliceViewer::buildTracerBuffers()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_) return;
		const int stride = FlowTracers::kStride * sizeof(float);
		glBindVertexArray(tracer_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, tracer_vbo_);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, stride, (void*)0);                   // pos.xyz + side
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)(4 * sizeof(float))); // tangent
		glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, stride, (void*)(7 * sizeof(float))); // rgba
		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);
		glEnableVertexAttribArray(2);
		glBindVertexArray(0);
	}

	void SliceViewer::updateTracers()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_ || !show_tracers_ || (!worker_ && !paraglider_worker_) || !have_info_)
		{
			tracer_strips_ = 0;
			return;
		}
		float dt = 0.0f;
		if (!tracer_clock_started_) { tracer_clock_.start(); tracer_clock_started_ = true; }
		else dt = (float)(tracer_clock_.restart() / 1000.0);

		TracerView view;
		view.three_d = tracer_3d_;
		view.axis = (int)axis_;
		{
			float L = (axis_ == Axis::X) ? (float)info_.Lx : (axis_ == Axis::Y) ? (float)info_.Ly : (float)info_.Lz;
			view.plane_pos = plane_frac_ * L;
		}
		// Same speed scale as the arrows/legend so a streamline's colour reads on the |u| ramp; falls
		// back to a fixed 1.8U when auto-range is off.
		view.speed_scale = auto_range_ ? std::max(auto_speed_max_, 1e-3f)
			: 1.8f * (float)info_.U;
		view.density = tracer_grid_density_;
		view.step_ds = (float)info_.h; // one cell per integration step (each step adds exactly h of arc)
		// Integration-step cap. The boring filter (below) thresholds a streamline's ARC LENGTH against a
		// multiple of Lx (slider top = 1.1·Lx), but a streamline integrates at most max_points·h of arc.
		// On a FINE grid (small h) a fixed step count spans < Lx, so NO streamline can reach the threshold
		// and the filter collapses to all-or-nothing after a resize. Floor the cap at a couple of domain
		// lengths of arc (steps to cross the domain × kTracerSpan) so the filter keeps a working gradient
		// at any resolution; a larger user "Tracer length" still wins. At the nominal grid (nx≈200) the
		// floor (500) is below the default trail (600), so nominal behaviour is unchanged.
		const float kTracerSpan = 2.5f; // guaranteed arc reach ≥ 2.5·Lx (well above the 1.1·Lx filter top)
		const int span_floor = (int)std::ceil(kTracerSpan * (float)info_.Lx / std::max(1e-6f, (float)info_.h));
		view.max_points = std::max(tracer_trail_, span_floor);
		// "Boring" filter: hide streamlines whose arc is below Lx·boring. The slider minimum (0.99) disables
		// the filter entirely (show all) — a fixed multiple can't reliably clear a straight crosser's arc
		// (≈ Lx − h) across domain sizes, so treat the bottom notch as "off". The span_floor above keeps
		// the threshold reachable regardless of resolution (else refining the grid makes it binary).
		view.min_length = (tracer_boring_ <= 0.995f) ? 0.0f : tracer_boring_ * (float)info_.Lx;
		view.dt = dt;
		view.hold_seconds = 1.0f; // keep a tracer for 1 s after it was last interesting (anti-flicker)
		view.instant = tracer_boring_instant_; // dragging the slider ⇒ bypass the hold (live filter)

		ensureFabricBvh();
		auto consume = [&](const FlowField& f) { FlowField ff = f; ff.fabric = fabric_bvh_.get(); tracers_.advance(view, ff); };
		const bool got = paraglider_worker_ ? paraglider_worker_->withFlowField(consume)
			: worker_->withFlowField(consume);
		if (!got) { tracer_strips_ = 0; return; }

		tracer_strips_ = tracers_.strips();
		tracer_firsts_ = tracers_.firsts();
		tracer_counts_ = tracers_.counts();
		if (tracer_strips_ == 0) return;

		if (!tracers_logged_)
		{
			tracers_logged_ = true;
			std::fprintf(stderr, "[G1] flow tracers live: %d inlet streamlines, %s mode\n",
				tracer_strips_, tracer_3d_ ? "3D" : "2D");
		}

		const std::vector<float>& v = tracers_.vertices();
		glBindVertexArray(tracer_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, tracer_vbo_);
		if ((int)v.size() > tracer_vbo_capacity_)
		{
			glBufferData(GL_ARRAY_BUFFER, v.size() * sizeof(float), v.data(), GL_DYNAMIC_DRAW);
			tracer_vbo_capacity_ = (int)v.size();
		}
		else glBufferSubData(GL_ARRAY_BUFFER, 0, v.size() * sizeof(float), v.data());
		glBindVertexArray(0);
	}

	void SliceViewer::drawTracers(const QMatrix4x4& mvp)
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!show_tracers_ || tracer_strips_ <= 0) return;
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDepthMask(GL_FALSE); // blend over the scene; still depth-TESTED against solids
		// Pull the ribbons slightly toward the camera in depth so a 2D streamline lying ON the slice
		// plane (or a 3D one grazing the bed) wins the z-fight instead of flickering; a genuine solid
		// in front still occludes (the bias is tiny).
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(-1.0f, -1.0f);
		tracer_prog_.bind();
		tracer_prog_.setUniformValue("uMVP", mvp);
		// Screen-space ribbon: constant PIXEL thickness. Viewport in framebuffer pixels (logical size
		// × device-pixel ratio, matching the GL viewport Qt set).
		const float dpr = (float)devicePixelRatioF();
		tracer_prog_.setUniformValue("uViewport", QVector2D((float)width() * dpr, (float)height() * dpr));
		tracer_prog_.setUniformValue("uHalfPx", 0.5f * tracer_width_);
		glBindVertexArray(tracer_vao_);
		glMultiDrawArrays(GL_TRIANGLE_STRIP, tracer_firsts_.data(), tracer_counts_.data(), tracer_strips_);
		glBindVertexArray(0);
		tracer_prog_.release();
		glPolygonOffset(0.0f, 0.0f);
		glDisable(GL_POLYGON_OFFSET_FILL);
		glDepthMask(GL_TRUE);
		glDisable(GL_BLEND);
	}

	bool SliceViewer::projectPoint(const QMatrix4x4& mvp, const QVector3D& w, QPointF& px) const
	{
		QVector4D clip = mvp * QVector4D(w, 1.0f);
		if (!(clip.w() > 1e-3f)) return false; // behind / on the camera plane (also rejects NaN)
		float nx = clip.x() / clip.w(), ny = clip.y() / clip.w();
		// Reject non-finite or comfortably off-screen points. Feeding an EXTREME pixel coordinate to
		// QPainter::drawText corrupts the GL paint engine's batched-glyph flush, which silently drops
		// the WHOLE text batch — including the legend numbers drawn earlier this frame (the "labels
		// flash for one frame" bug). Clamp the domain of what we ever hand to drawText.
		if (!std::isfinite(nx) || !std::isfinite(ny)) return false;
		if (nx < -1.3f || nx > 1.3f || ny < -1.3f || ny > 1.3f) return false;
		px.setX((nx * 0.5f + 0.5f) * width());
		px.setY((1.0f - (ny * 0.5f + 0.5f)) * height());
		return true;
	}

	void SliceViewer::drawAxes(const QMatrix4x4& mvp)
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!show_axes_) return;
		if (axes_dirty_) buildAxesGeometry();
		const QVector4D col[3] = {
			QVector4D(0.90f, 0.28f, 0.24f, 1.0f), // X red
			QVector4D(0.30f, 0.78f, 0.32f, 1.0f), // Y green
			QVector4D(0.32f, 0.55f, 0.95f, 1.0f)  // Z blue
		};
		glLineWidth(2.0f);
		prog_.bind();
		prog_.setUniformValue("uFlat", 1);

		// World-origin triad: depth-tested (part of the scene) so it sits correctly among the geometry.
		prog_.setUniformValue("uMVP", mvp);
		glBindVertexArray(axis_vao_);
		for (int a = 0; a < 3; ++a)
		{
			prog_.setUniformValue("uColor", col[a]);
			glDrawArrays(GL_LINES, axis_vert_off_[a], axis_vert_cnt_[a]);
		}

		// Camera-aligned orientation gizmo (bottom-left), always on top (depth test off). A
		// rotation-only view + orthographic anchor keeps it fixed on screen while it spins with the view.
		QMatrix4x4 rot = camera_.view();
		rot.setColumn(3, QVector4D(0, 0, 0, 1)); // strip translation → rotation only
		float aspect = (float)width() / (float)std::max(1, height());
		QMatrix4x4 cproj;
		cproj.translate(-0.84f, -0.80f, 0.0f);      // bottom-left NDC anchor
		cproj.scale(0.15f / aspect, 0.15f, 0.001f); // square (aspect-corrected) + flatten z so it never clips
		corner_mvp_ = cproj * rot;
		glDisable(GL_DEPTH_TEST);
		prog_.setUniformValue("uMVP", corner_mvp_);
		glBindVertexArray(corner_vao_);
		for (int a = 0; a < 3; ++a)
		{
			prog_.setUniformValue("uColor", col[a]);
			glDrawArrays(GL_LINES, a * 2, 2);
		}
		glBindVertexArray(0);
		prog_.release();
		glEnable(GL_DEPTH_TEST);
		glLineWidth(1.0f);
	}

	void SliceViewer::drawAxesLabels(QPainter& p)
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!show_axes_ || !have_info_) return;
		p.setRenderHint(QPainter::TextAntialiasing, true);
		const QColor cx(232, 96, 84), cy(110, 205, 120), cz(120, 160, 245), cw(235, 237, 240);
		const float Lx = (float)info_.Lx, Ly = (float)info_.Ly, Lz = (float)info_.Lz;
		const QMatrix4x4 mvp = camera_.viewProjection();

		auto text = [&](const QMatrix4x4& m, const QVector3D& w, const QString& s, const QColor& c, QPointF off)
		{
			QPointF px;
			if (projectPoint(m, w, px)) { p.setPen(c); p.drawText(px + off, s); }
		};

		// Axis-end letters + the origin marker (bold).
		QFont f = p.font(); f.setPointSizeF(10.0); f.setBold(true); p.setFont(f);
		text(mvp, QVector3D(Lx, 0, 0), "X", cx, QPointF(4, -2));
		text(mvp, QVector3D(0, Ly, 0), "Y", cy, QPointF(4, -2));
		text(mvp, QVector3D(0, 0, Lz), "Z", cz, QPointF(4, -2));
		text(mvp, QVector3D(0, 0, 0), "0", cw, QPointF(-12, 13));

		// Metre tick labels (smaller, non-bold) so the user can eyeball world coordinates.
		f.setBold(false); f.setPointSizeF(8.0); p.setFont(f);
		auto ticks = [&](char ax, float L)
		{
			for (float t = axis_tick_step_; t <= L + 1e-4f; t += axis_tick_step_)
			{
				QVector3D w = (ax == 'x') ? QVector3D(t, 0, 0) : (ax == 'y') ? QVector3D(0, t, 0) : QVector3D(0, 0, t);
				text(mvp, w, QString::number(t, 'g', 3), (ax == 'x') ? cx : (ax == 'y') ? cy : cz, QPointF(3, 11));
			}
		};
		ticks('x', Lx); ticks('y', Ly); ticks('z', Lz);

		// Corner-gizmo letters (uses the transform drawAxes() set this same frame).
		f.setBold(true); f.setPointSizeF(9.0); p.setFont(f);
		text(corner_mvp_, QVector3D(1.18f, 0, 0), "X", cx, QPointF(-4, 4));
		text(corner_mvp_, QVector3D(0, 1.18f, 0), "Y", cy, QPointF(-4, 4));
		text(corner_mvp_, QVector3D(0, 0, 1.18f), "Z", cz, QPointF(-4, 4));
	}

	void SliceViewer::drawLegendWith(QPainter& p)
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!have_info_) return;
		p.setRenderHint(QPainter::Antialiasing, true);
		p.setRenderHint(QPainter::TextAntialiasing, true);

		// One colorbar: value→colour ramp with a title (field + units) and min/mid/max tick labels.
		// `sample(t)` returns the fill colour for the normalised bar position t in [0,1] (bottom→top).
		auto draw_bar = [&](int x, int y, int w, int h, const QString& title,
			double vmin, double vmax, const QString& units, auto sample)
		{
			const int steps = 64;
			for (int s = 0; s < steps; ++s)
			{
				float t0 = (float)s / steps;
				int yy = y + h - (int)((s + 1) * (float)h / steps);
				int hh = (int)std::ceil((float)h / steps) + 1;
				p.fillRect(QRect(x, yy, w, hh), sample(t0));
			}
			p.setPen(QColor(255, 255, 255)); // white bar outline
			p.drawRect(QRect(x, y, w, h));
			QFont f = p.font(); f.setPointSizeF(9.0); p.setFont(f);
			// Title + tick numbers in WHITE on a subtle dark backing so they read over the coloured bar,
			// the scene and the arrows (and never spill past the widget's right edge).
			QRect titleRect(x - 52, y - 20, w + 104, 16);
			p.fillRect(titleRect, QColor(18, 20, 24, 160));
			p.setPen(QColor(255, 255, 255));
			p.drawText(titleRect, Qt::AlignCenter, title + "  [" + units + "]");
			auto label = [&](double val, int yy)
			{
				QRect r(x + w + 5, yy - 8, 46, 16);
				p.fillRect(r, QColor(18, 20, 24, 160));
				p.setPen(QColor(255, 255, 255));
				p.drawText(r, Qt::AlignLeft | Qt::AlignVCenter, QString::number(val, 'g', 3));
			};
			label(vmax, y);
			label(0.5 * (vmin + vmax), y + h / 2);
			label(vmin, y + h);
		};

		const int W = width();
		int barW = 16, barH = std::min(220, height() - 90);
		int barX = W - barW - 66;
		int barY = 40;

		// --- Active slice field ramp (the shared 5-stop colormap over [vmin_,vmax_]) --------------
		QString title, units;
		switch (field_)
		{
		case Field::SpeedMag: title = "|u|"; units = "m/s"; break;
		case Field::VelU: title = "u"; units = "m/s"; break;
		case Field::VelV: title = "v"; units = "m/s"; break;
		case Field::VelW: title = "w"; units = "m/s"; break;
		case Field::Pressure: title = "p"; units = "Pa"; break;
		}
		draw_bar(barX, barY, barW, barH, title, vmin_, vmax_, units,
			[](float t) { float r, g, b; scour_colormap(t, r, g, b); return QColor::fromRgbF(r, g, b); });

		// --- Cp wind-load ramp (shown when "Colour building by Cp" is ON and the building is Cp-coloured) --
		// Same 5-stop ramp the overlay uses, over the current SYMMETRIC Cp range [vox_cp_lo_, vox_cp_hi_], so
		// it's clear what surface pressures the building colours mean (blue = suction, red = pressure). Ticks
		// read min / 0 / max (0 falls on the mid, since the range is symmetric about Cp = 0). Placed to the
		// LEFT of the slice bar so the two colorbars + their labels don't collide.
		if (colour_by_cp_ && has_vox_cp_)
		{
			const int cpX = barX - 130;
			draw_bar(cpX, barY, barW, barH, "Cp", vox_cp_lo_, vox_cp_hi_, "pressure coeff.",
				[](float t) { float r, g, b; scour_colormap(t, r, g, b); return QColor::fromRgbF(r, g, b); });
		}

		// External-aero convention is deliberately fixed: prescribed freestream enters at
		// X-min and travels along +X. Keep this visible even when no arrows cross the camera.
		if (paraglider_worker_)
		{
			QFont f = p.font(); f.setPointSizeF(10.0); f.setBold(true); p.setFont(f);
			const QRect cue(18, 16, 210, 28);
			p.fillRect(cue, QColor(18, 20, 24, 190));
			p.setPen(QColor(255, 120, 90));
			p.drawText(cue.adjusted(9, 0, -5, 0), Qt::AlignLeft | Qt::AlignVCenter,
				QString::fromUtf8("FREESTREAM  +X  →"));
		}
	}

	void SliceViewer::buildSliceGeometry()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_ || !have_info_) return;
		SliceParams sp = currentParams();
		std::vector<float> pos((size_t)NRES * NRES * 3);
		for (int b = 0; b < NRES; ++b)
			for (int a = 0; a < NRES; ++a)
			{
				float x, y, z;
				slice_vertex_world(sp, a, b, x, y, z);
				size_t o = (size_t)(b * NRES + a) * 3;
				pos[o + 0] = x; pos[o + 1] = y; pos[o + 2] = z;
			}
		glBindBuffer(GL_ARRAY_BUFFER, pos_vbo_);
		glBufferSubData(GL_ARRAY_BUFFER, 0, pos.size() * sizeof(float), pos.data());
		geometry_dirty_ = false;
	}

	void SliceViewer::buildBoxGeometry()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_) return;
		float Lx = (float)info_.Lx, Ly = (float)info_.Ly, Lz = (float)info_.Lz;
		if (Lx <= 0) { Lx = 10; Ly = 10; Lz = 5; }
		float c[8][3] = {
			{ 0, 0, 0 }, { Lx, 0, 0 }, { Lx, Ly, 0 }, { 0, Ly, 0 },
			{ 0, 0, Lz }, { Lx, 0, Lz }, { Lx, Ly, Lz }, { 0, Ly, Lz }
		};
		int edges[12][2] = {
			{ 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
			{ 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 },
			{ 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }
		};
		std::vector<float> v;
		v.reserve(12 * 2 * 3);
		for (auto& e : edges)
			for (int p = 0; p < 2; ++p)
			{
				v.push_back(c[e[p]][0]); v.push_back(c[e[p]][1]); v.push_back(c[e[p]][2]);
			}
		box_vertex_count_ = (int)(v.size() / 3);
		glBindVertexArray(box_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, box_vbo_);
		glBufferData(GL_ARRAY_BUFFER, v.size() * sizeof(float), v.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);
		glBindVertexArray(0);
	}

	// Cell-boundary lines where the grid meets the CURRENT slice plane. The lines are drawn at the per-axis
	// cell-face coordinates (the graded metric arrays from setGridLines, or uniform i·h if none) so a graded
	// mesh reads directly: lines bunch up in the fine h_fine core and spread out in the coarse far field.
	void SliceViewer::buildGridGeometry()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_) return;
		grid_dirty_ = false;
		float Lx = (float)info_.Lx, Ly = (float)info_.Ly, Lz = (float)info_.Lz;
		if (Lx <= 0) { Lx = 10; Ly = 10; Lz = 5; }
		// Per-axis face coordinates: the provided metric arrays, or a uniform i·h fallback (nx cells over Lx).
		auto faces = [this](const std::vector<float>& arr, int n, float L)
		{
			if (!arr.empty()) return arr;
			std::vector<float> f;
			f.reserve((std::size_t)n + 1);
			const float h = n > 0 ? L / (float)n : L;
			for (int i = 0; i <= n; ++i) f.push_back((float)i * h);
			return f;
		};
		const std::vector<float> xf = faces(grid_xf_, info_.nx, Lx);
		const std::vector<float> yf = faces(grid_yf_, info_.ny, Ly);
		const std::vector<float> zf = faces(grid_zf_, info_.nz, Lz);

		std::vector<float> v;
		auto seg = [&](float ax, float ay, float az, float bx, float by, float bz)
		{ v.push_back(ax); v.push_back(ay); v.push_back(az); v.push_back(bx); v.push_back(by); v.push_back(bz); };

		if (axis_ == Axis::Z) // z-normal plane: draw the x-y cell grid at z = z0
		{
			const float z0 = plane_frac_ * Lz;
			for (float x : xf) seg(x, 0.0f, z0, x, Ly, z0);
			for (float y : yf) seg(0.0f, y, z0, Lx, y, z0);
		}
		else if (axis_ == Axis::X) // x-normal plane: y-z grid at x = x0
		{
			const float x0 = plane_frac_ * Lx;
			for (float y : yf) seg(x0, y, 0.0f, x0, y, Lz);
			for (float z : zf) seg(x0, 0.0f, z, x0, Ly, z);
		}
		else // y-normal plane: x-z grid at y = y0
		{
			const float y0 = plane_frac_ * Ly;
			for (float x : xf) seg(x, y0, 0.0f, x, y0, Lz);
			for (float z : zf) seg(0.0f, y0, z, Lx, y0, z);
		}

		grid_vertex_count_ = (int)(v.size() / 3);
		glBindVertexArray(grid_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, grid_vbo_);
		glBufferData(GL_ARRAY_BUFFER, v.size() * sizeof(float), v.data(), GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);
		glBindVertexArray(0);
	}

	void SliceViewer::drawGrid(const QMatrix4x4& mvp)
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!show_grid_ || !gl_ready_) return;
		if (grid_dirty_) buildGridGeometry();
		if (grid_vertex_count_ <= 0) return;
		prog_.bind();
		prog_.setUniformValue("uMVP", mvp);
		prog_.setUniformValue("uClipPlane", computeClipPlane());
		glDisable(GL_CLIP_DISTANCE0); // grid lines are a reference overlay — never clipped
		prog_.setUniformValue("uFlat", 1);
		prog_.setUniformValue("uColor", QVector4D(0.35f, 0.85f, 0.55f, 0.7f)); // translucent green
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDepthFunc(GL_LEQUAL); // the lines are coplanar with the slice mesh — let equal-depth fragments win
		glBindVertexArray(grid_vao_);
		glDrawArrays(GL_LINES, 0, grid_vertex_count_);
		glBindVertexArray(0);
		glDepthFunc(GL_LESS);
		prog_.release();
	}

	namespace
	{
		// "Nice" tick spacing so each axis carries ~5 marks (10 m ⇒ 2 m, 5 m ⇒ 2 m, etc.).
		inline float axis_tick_step(float maxL)
		{
			if (maxL <= 6.0f) return 1.0f;
			if (maxL <= 15.0f) return 2.0f;
			if (maxL <= 40.0f) return 5.0f;
			return 10.0f;
		}
	}

	void SliceViewer::buildAxesGeometry()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_) return;
		float Lx = (float)info_.Lx, Ly = (float)info_.Ly, Lz = (float)info_.Lz;
		if (Lx <= 0) { Lx = 10; Ly = 10; Lz = 5; }
		const float maxL = std::max({ Lx, Ly, Lz });
		const float step = axis_tick_step(maxL);
		axis_tick_step_ = step;
		const float tk = std::clamp(0.03f * maxL, 0.05f, 0.5f); // tick-mark length [m]

		std::vector<float> v;
		auto seg = [&](float ax, float ay, float az, float bx, float by, float bz)
		{ v.push_back(ax); v.push_back(ay); v.push_back(az); v.push_back(bx); v.push_back(by); v.push_back(bz); };

		// X axis (+ ticks marked in +y and +z so they read from any view).
		axis_vert_off_[0] = 0;
		seg(0, 0, 0, Lx, 0, 0);
		for (float t = step; t <= Lx + 1e-4f; t += step) { seg(t, 0, 0, t, tk, 0); seg(t, 0, 0, t, 0, tk); }
		axis_vert_cnt_[0] = (int)v.size() / 3 - axis_vert_off_[0];
		// Y axis (ticks in +x and +z).
		axis_vert_off_[1] = (int)v.size() / 3;
		seg(0, 0, 0, 0, Ly, 0);
		for (float t = step; t <= Ly + 1e-4f; t += step) { seg(0, t, 0, tk, t, 0); seg(0, t, 0, 0, t, tk); }
		axis_vert_cnt_[1] = (int)v.size() / 3 - axis_vert_off_[1];
		// Z axis (ticks in +x and +y).
		axis_vert_off_[2] = (int)v.size() / 3;
		seg(0, 0, 0, 0, 0, Lz);
		for (float t = step; t <= Lz + 1e-4f; t += step) { seg(0, 0, t, tk, 0, t); seg(0, 0, t, 0, tk, t); }
		axis_vert_cnt_[2] = (int)v.size() / 3 - axis_vert_off_[2];

		glBindVertexArray(axis_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, axis_vbo_);
		glBufferData(GL_ARRAY_BUFFER, v.size() * sizeof(float), v.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);
		glBindVertexArray(0);
		axes_dirty_ = false;
	}

	void SliceViewer::buildCornerGizmo()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_) return;
		// Three unit axes from the origin; drawn with a camera-rotation-only transform so it always
		// shows the world X/Y/Z directions as seen from the current view.
		const float g[] = {
			0, 0, 0, 1, 0, 0, // X
			0, 0, 0, 0, 1, 0, // Y
			0, 0, 0, 0, 0, 1  // Z
		};
		glBindVertexArray(corner_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, corner_vbo_);
		glBufferData(GL_ARRAY_BUFFER, sizeof(g), g, GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);
		glBindVertexArray(0);
	}

	QVector4D SliceViewer::computeClipPlane() const
	{
		// Disabled / no domain yet ⇒ an all-keep plane: dot(pos, 0) + 1 = 1 ≥ 0 everywhere.
		if (!clip_enabled_ || !have_info_) return QVector4D(0.0f, 0.0f, 0.0f, 1.0f);

		const float Lx = (float)info_.Lx, Ly = (float)info_.Ly, Lz = (float)info_.Lz;
		QVector3D P, n; // a point on the plane + the KEPT-side normal (kept where dot(x-P, n) ≥ 0)

		if (clip_mode_ == 3) // face-camera: normal = view direction, swept about the target along it
		{
			const QVector3D fwd = camera_.forward(); // eye → scene (unit)
			const QVector3D ctr = camera_.targetPoint();
			const float span = std::sqrt(Lx * Lx + Ly * Ly + Lz * Lz); // full domain depth range
			P = ctr + fwd * ((clip_frac_ - 0.5f) * span);
			n = fwd; // keep the far side ⇒ hide everything between the camera and the plane
		}
		else // axis-aligned X / Y / Z
		{
			const int a = clip_mode_;
			const float L = (a == 0) ? Lx : (a == 1) ? Ly : Lz;
			const float p = clip_frac_ * L;
			const QVector3D e(a == 0 ? 1.0f : 0.0f, a == 1 ? 1.0f : 0.0f, a == 2 ? 1.0f : 0.0f);
			P = e * p;
			const float eye_a = camera_.eye()[a];
			n = e * ((eye_a > p) ? -1.0f : 1.0f); // keep the side away from the camera ⇒ hide the near side
		}

		if (clip_flip_) n = -n;
		const float d = -QVector3D::dotProduct(P, n);
		return QVector4D(n.x(), n.y(), n.z(), d);
	}

	void SliceViewer::drawClipPlaneViz(const QMatrix4x4& mvp, const QVector4D& plane)
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_ || !clip_enabled_ || !have_info_) return;
		QVector3D n(plane.x(), plane.y(), plane.z());
		if (n.lengthSquared() < 1e-8f) return;
		n.normalize();

		// Centre the visualisation quad on the domain centre projected onto the plane (so it covers the
		// scene wherever the cut crosses it), sized to the domain diagonal. Two in-plane axes u,v span it.
		const float Lx = (float)info_.Lx, Ly = (float)info_.Ly, Lz = (float)info_.Lz;
		const QVector3D dctr(0.5f * Lx, 0.5f * Ly, 0.5f * Lz);
		const float sd = QVector3D::dotProduct(dctr, n) + plane.w(); // signed distance centre → plane
		const QVector3D P = dctr - n * sd;                           // projection of the centre onto the plane
		const QVector3D ref = (std::fabs(n.z()) < 0.9f) ? QVector3D(0, 0, 1) : QVector3D(1, 0, 0);
		const QVector3D u = QVector3D::crossProduct(n, ref).normalized();
		const QVector3D v = QVector3D::crossProduct(n, u).normalized();
		const float R = 0.5f * std::sqrt(Lx * Lx + Ly * Ly + Lz * Lz);
		const QVector3D c[4] = { P - u * R - v * R, P + u * R - v * R, P + u * R + v * R, P - u * R + v * R };
		float verts[12];
		for (int i = 0; i < 4; ++i) { verts[i * 3] = c[i].x(); verts[i * 3 + 1] = c[i].y(); verts[i * 3 + 2] = c[i].z(); }

		glBindVertexArray(clip_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, clip_vbo_);
		glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);

		prog_.bind();
		prog_.setUniformValue("uMVP", mvp);
		prog_.setUniformValue("uFlat", 1);
		// Faint translucent fill (no depth write, so it never occludes the exposed interior) + bright outline.
		// Drawn with GL_CLIP_DISTANCE0 OFF (the caller disables it) so the plane itself is not clipped.
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDepthMask(GL_FALSE);
		prog_.setUniformValue("uColor", QVector4D(0.20f, 0.85f, 0.95f, 0.12f));
		glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
		glLineWidth(2.0f);
		prog_.setUniformValue("uColor", QVector4D(0.30f, 0.90f, 1.0f, 0.9f));
		glDrawArrays(GL_LINE_LOOP, 0, 4);
		glLineWidth(1.0f);
		glDepthMask(GL_TRUE);
		glDisable(GL_BLEND);
		glBindVertexArray(0);
		prog_.release();
	}

	void SliceViewer::uploadMesh()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_ || !mesh_upload_pending_) return;
		mesh_upload_pending_ = false;

		const paracfd::core::TriMesh& m = pending_mesh_;
		if (m.empty())
		{
			has_mesh_ = false;
			mesh_index_count_ = 0;
			pending_mesh_ = paracfd::core::TriMesh{};
			return;
		}

		glBindVertexArray(mesh_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, mesh_pos_vbo_);
		glBufferData(GL_ARRAY_BUFFER, m.positions.size() * sizeof(float), m.positions.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, mesh_norm_vbo_);
		glBufferData(GL_ARRAY_BUFFER, m.normals.size() * sizeof(float), m.normals.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_idx_ebo_);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, m.indices.size() * sizeof(unsigned int), m.indices.data(), GL_STATIC_DRAW);
		glBindVertexArray(0);
		mesh_index_count_ = (int)m.indices.size();

		// The model transform is the gizmo TRS (seeded centre-on-bed in setMesh) or an explicit override
		// matrix (setMeshTranslate/Placement) — both are already set, so uploadMesh only pushes the
		// geometry. modelMatrix() supplies the placement at draw time.

		has_mesh_ = true;
		mesh_triangle_colour_upload_pending_ = true; // revalidate colours against the new triangle count
		pending_mesh_ = paracfd::core::TriMesh{}; // free CPU copy; it lives in GL now
	}

	void SliceViewer::uploadTriangleSurfaceColours()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_ || !mesh_triangle_colour_upload_pending_) return;
		mesh_triangle_colour_upload_pending_ = false;
		has_mesh_triangle_colours_ = has_mesh_ && mesh_index_count_ > 0 &&
			mesh_triangle_colours_.size() == static_cast<std::size_t>(mesh_index_count_ / 3) * 8;
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, mesh_triangle_colour_ssbo_);
		if (has_mesh_triangle_colours_)
			glBufferData(GL_SHADER_STORAGE_BUFFER, mesh_triangle_colours_.size() * sizeof(float), mesh_triangle_colours_.data(), GL_DYNAMIC_DRAW);
		else glBufferData(GL_SHADER_STORAGE_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	}

	void SliceViewer::uploadParagliderDebugBoxes()
	{
		PARACFD_ASSERT_GL_THREAD();if(!gl_ready_||!paraglider_debug_upload_pending_)return;paraglider_debug_upload_pending_=false;
		auto upload=[&](unsigned vao,unsigned vbo,const std::vector<float>& lines,int& count){glBindVertexArray(vao);glBindBuffer(GL_ARRAY_BUFFER,vbo);glBufferData(GL_ARRAY_BUFFER,lines.size()*sizeof(float),lines.empty()?nullptr:lines.data(),GL_DYNAMIC_DRAW);glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(float),(void*)0);glEnableVertexAttribArray(0);glBindVertexArray(0);count=static_cast<int>(lines.size()/3);};
		upload(amr_debug_vao_,amr_debug_vbo_,pending_amr_debug_lines_,amr_debug_vertex_count_);upload(eb_debug_vao_,eb_debug_vbo_,pending_eb_debug_lines_,eb_debug_vertex_count_);
	}

	void SliceViewer::uploadVoxelOverlay()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_ || !vox_upload_pending_) return;
		vox_upload_pending_ = false;

		const paracfd::core::MacGrid g = pending_vox_grid_;
		const std::vector<unsigned char>& s = pending_vox_solid_;
		if (s.empty() || (int)s.size() != g.p_count() || g.h <= 0.0)
		{
			has_vox_ = false;
			vox_vertex_count_ = 0;
			pending_vox_solid_.clear();
			return;
		}

		// The mask is the flow solid mask (0 = fluid, nonzero = solid). Emit the exposed voxel-surface
		// faces — a face wherever a solid cell's across-neighbour is fluid or outside the domain — i.e.
		// the solid staircase skin, not its interior. Each face = 2 triangles, each vertex = 3 pos + 3 nrm.
		auto solid_at = [&](int i, int j, int k) -> bool
		{
			if (i < 0 || i >= g.nx || j < 0 || j >= g.ny || k < 0 || k >= g.nz) return false; // out of domain = fluid
			return s[(std::size_t)g.pidx(i, j, k)] != 0;
		};

		// Colour by Cp when enabled AND the worker has published a Cp field matching this grid. Each solid
		// cell's exposed faces take the cell's mean Cp, mapped through the shared colormap ramp over the
		// symmetric range [vox_cp_lo_,vox_cp_hi_] (blue = suction … red = pressure).
		const bool want_cp = colour_by_cp_ && (int)vox_cell_cp_.size() == g.p_count();
		const float cp_span = (vox_cp_hi_ > vox_cp_lo_) ? (vox_cp_hi_ - vox_cp_lo_) : 1.0f;

		std::vector<float> pos, nrm, col;
		float cr = 0.5f, cg = 0.5f, cb = 0.5f; // current cell's Cp colour (set per solid cell below)
		auto quad = [&](float ox, float oy, float oz, float ux, float uy, float uz, float vx, float vy, float vz, float nx, float ny, float nz)
		{
			float p0[3] = { ox, oy, oz };
			float p1[3] = { ox + ux, oy + uy, oz + uz };
			float p2[3] = { ox + ux + vx, oy + uy + vy, oz + uz + vz };
			float p3[3] = { ox + vx, oy + vy, oz + vz };
			const float* tris[6] = { p0, p1, p2, p0, p2, p3 };
			for (const float* p : tris)
			{
				pos.push_back(p[0]); pos.push_back(p[1]); pos.push_back(p[2]);
				nrm.push_back(nx); nrm.push_back(ny); nrm.push_back(nz);
				if (want_cp) { col.push_back(cr); col.push_back(cg); col.push_back(cb); }
			}
		};
		for (int k = 0; k < g.nz; ++k)
			for (int j = 0; j < g.ny; ++j)
				for (int i = 0; i < g.nx; ++i)
				{
					if (!solid_at(i, j, k)) continue;
					if (want_cp)
					{
						const float cp = vox_cell_cp_[(std::size_t)g.pidx(i, j, k)];
						if (cp == cp) // finite (NaN off the surface → neutral grey)
							scour_colormap((cp - vox_cp_lo_) / cp_span, cr, cg, cb);
						else { cr = cg = cb = 0.6f; }
					}
					// Cell's minimum-corner face coords + per-axis widths from the grid metrics — on a GRADED
					// grid these are the true (non-uniform) cell bounds xf[i]/dx[i]; the accessors fall back to
					// exact i·h / h on a uniform grid, so the staircase lands on the same voxels the flow uses.
					const float x0 = (float)g.xf(i), y0 = (float)g.yf(j), z0 = (float)g.zf(k);
					const float hx = (float)g.dx(i), hy = (float)g.dy(j), hz = (float)g.dz(k);
					if (!solid_at(i - 1, j, k)) quad(x0, y0, z0, 0, hy, 0, 0, 0, hz, -1, 0, 0); // -x
					if (!solid_at(i + 1, j, k)) quad(x0 + hx, y0, z0, 0, 0, hz, 0, hy, 0, 1, 0, 0); // +x
					if (!solid_at(i, j - 1, k)) quad(x0, y0, z0, 0, 0, hz, hx, 0, 0, 0, -1, 0); // -y
					if (!solid_at(i, j + 1, k)) quad(x0, y0 + hy, z0, hx, 0, 0, 0, 0, hz, 0, 1, 0); // +y
					if (!solid_at(i, j, k - 1)) quad(x0, y0, z0, hx, 0, 0, 0, hy, 0, 0, 0, -1); // -z
					if (!solid_at(i, j, k + 1)) quad(x0, y0, z0 + hz, 0, hy, 0, hx, 0, 0, 0, 0, 1); // +z
				}

		vox_vertex_count_ = (int)(pos.size() / 3);
		if (vox_vertex_count_ == 0) { has_vox_ = false; has_vox_cp_ = false; return; }

		glBindVertexArray(vox_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, vox_pos_vbo_);
		glBufferData(GL_ARRAY_BUFFER, pos.size() * sizeof(float), pos.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vox_norm_vbo_);
		glBufferData(GL_ARRAY_BUFFER, nrm.size() * sizeof(float), nrm.data(), GL_STATIC_DRAW);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(1);
		has_vox_cp_ = want_cp && col.size() == pos.size();
		if (has_vox_cp_)
		{
			glBindBuffer(GL_ARRAY_BUFFER, vox_color_vbo_);
			glBufferData(GL_ARRAY_BUFFER, col.size() * sizeof(float), col.data(), GL_STATIC_DRAW);
			glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
			glEnableVertexAttribArray(2);
		}
		else
			glDisableVertexAttribArray(2); // fall back to the flat uBaseColor path
		glBindVertexArray(0);

		has_vox_ = true;
		// NB: pending_vox_solid_ is RETAINED (not cleared) so a fresh Cp field — or a colour_by_cp_ toggle —
		// can re-extract + re-colour the overlay without a mask republish from the worker.
		std::fprintf(stderr, "[G1] voxel overlay: %d solid exposed-face triangles%s\n",
			vox_vertex_count_ / 3, has_vox_cp_ ? " (Cp-coloured)" : "");
	}

	void SliceViewer::resizeGL(int w, int h)
	{
		PARACFD_ASSERT_GL_THREAD();
		glViewport(0, 0, w, h);
		camera_.setViewport(w, h);
	}

	void SliceViewer::paintGL()
	{
		PARACFD_ASSERT_GL_THREAD();

		// QPainter drives the 2D legend overlay (feature 3) AFTER the raw GL scene. All raw GL is
		// fenced inside begin/endNativePainting so Qt's paint engine restores its own GL state.
		QPainter painter(this);
		painter.beginNativePainting();

		if (geometry_dirty_) { buildSliceGeometry(); buildBoxGeometry(); } // both track info_.Lx/Ly/Lz — rebuild together on a domain resize

		// Sync the voxel overlay to the live flow solid mask (obstacle ∪ any runtime-marked solid) so
		// EVERY cell the fluid treats as solid is drawn — a mis-placed solid
		// voxel becomes visible. The worker bumps its mask generation only on a mask change (load /
		// obstacle rebuild / bed re-mask), so this copies the host mask + rebuilds the exposed-face
		// geometry ONLY then (no per-frame mask D2H). uploadVoxelOverlay() (below) does the extraction.
		if (worker_ && worker_->maskGeneration() != vox_mask_gen_)
		{
			std::vector<unsigned char> mask;
			paracfd::core::MacGrid mg;
			if (worker_->copyMask(mask, mg))
			{
				vox_mask_gen_ = worker_->maskGeneration();
				setVoxelOverlay(mask, mg); // queues uploadVoxelOverlay() for this same paint
			}
		}

		// Sync the per-cell Cp field for the building overlay. The worker republishes it (bumping its load
		// generation) only when it re-integrates the loads (~every 30 steps), so this pulls the fresh Cp +
		// its symmetric colour range and queues a recolour ONLY then — never per-frame. Runs after the mask
		// sync so both feed the SAME uploadVoxelOverlay() below.
		if (worker_ && worker_->loadGeneration() != vox_load_gen_)
		{
			std::vector<float> cp;
			float lo = 0.0f, hi = 0.0f;
			paracfd::core::MacGrid lg;
			// The worker hands back the TIME-AVERAGED Cp (+ its symmetric range) once an averaging window is
			// active, else the live instantaneous Cp — so the building recolours to the converged field the
			// moment the user starts averaging, with no extra plumbing here.
			if (worker_->copyDisplayCp(cp, lo, hi, lg))
			{
				vox_load_gen_ = worker_->loadGeneration();
				vox_cell_cp_.swap(cp);
				vox_cp_lo_ = lo;
				vox_cp_hi_ = hi;
				if (colour_by_cp_ && !pending_vox_solid_.empty()) vox_upload_pending_ = true; // recolour
			}
		}

		if (mesh_upload_pending_) uploadMesh();
		if (mesh_triangle_colour_upload_pending_) uploadTriangleSurfaceColours();
		if (paraglider_debug_upload_pending_) uploadParagliderDebugBoxes();
		if (vox_upload_pending_) uploadVoxelOverlay();

		// CUDA writes the slice colours directly into the registered VBO (zero copy), sampling the
		// worker's post-step device snapshots under display_mutex(). The interop runs on its own
		// stream, so this never blocks on the worker's step(). The same locked scope runs the
		// throttled auto-range reduction (fluid cells only, 3 scalars back — never the whole field)
		// FIRST, so the fill + legend + arrows all use the freshly-tracked [vmin_,vmax_].
		if (worker_ && worker_->display_ready() && cuda_res_ && have_info_)
		{
			std::lock_guard<std::mutex> lk(worker_->display_mutex());
			// Sampler grid = the core's DEVICE-view metrics (grid_fx maps each vertex to the correct GRADED
			// cell on the device). currentParams() alone carries only nx/ny/nz/h + the true extent (sp.Lx…);
			// displayGrid() adds the device metric ARRAYS with identical dims. Uniform ⇒ null metrics (x/h).
			SliceParams sp = currentParams();
			sp.grid = worker_->displayGrid();
			if (auto_range_ && (range_ctr_++ % kRangeEvery == 0))
			{
				FieldRange fr;
				if (slice_gl_reduce(cuda_res_, worker_->disp_u(), worker_->disp_v(), worker_->disp_w(),
					worker_->disp_p(), worker_->disp_solid(), sp.grid, field_, &fr))
					applyAutoRange(fr);
			}
			if (show_slice_) // slice fill skipped when hidden (no fill, no draw)
				slice_gl_fill(cuda_res_, worker_->disp_u(), worker_->disp_v(), worker_->disp_w(), worker_->disp_p(), sp);
		}
		else if (paraglider_worker_ && have_info_)
		{
			// The production field is block-AMR/FP32. The worker publishes a throttled,
			// coarse uniform host resampling exclusively for visualization, which lets the
			// mature slice renderer remain useful while AMR-native GL interop is developed.
			paraglider_worker_->withFlowField([&](const FlowField& field)
			{
				SliceParams sp = currentParams();
				sp.grid = field.grid;
				if (auto_range_ && (range_ctr_++ % kRangeEvery == 0))
				{
					float reduced[3]{};
					slice_reduce_cpu(field.u, field.v, field.w, field.p, nullptr, field.grid, field_, reduced);
					FieldRange range;
					range.field_min = reduced[0]; range.field_max = reduced[1]; range.speed_max = reduced[2];
					range.valid = std::isfinite(reduced[0]) && std::isfinite(reduced[1]) && reduced[0] <= reduced[1];
					applyAutoRange(range);
					sp.vmin = vmin_; sp.vmax = vmax_;
				}
				if (show_slice_)
				{
					host_slice_colours_.resize(static_cast<std::size_t>(sp.nu) * sp.nv);
					slice_fill_cpu(field.u, field.v, field.w, field.p, sp, host_slice_colours_.data());
					glBindBuffer(GL_ARRAY_BUFFER, color_vbo_);
					glBufferSubData(GL_ARRAY_BUFFER, 0,
						host_slice_colours_.size() * sizeof(float4), host_slice_colours_.data());
					glBindBuffer(GL_ARRAY_BUFFER, 0);
				}
			});
		}

		updateArrows(); // advect the arrow tracers + upload instance data (reads the worker's host flow)
		updateTracers(); // advect the streaklines + upload line geometry (same host-flow snapshot)

		// QPainter::beginNativePainting() resets GL to DEFAULT state (depth test OFF, depth mask
		// undefined, clear colour black, blending on) — so the scene's depth buffer must be
		// re-established every frame here, else geometry draws in submission order (voxels on top,
		// polygons unsorted). Restore the exact state initializeGL() used to set once.
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);
		glDepthMask(GL_TRUE);
		glDisable(GL_BLEND);
		glDisable(GL_CULL_FACE);
		glClearColor(0.09f, 0.10f, 0.12f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		prog_.bind();
		QMatrix4x4 mvp = camera_.viewProjection();
		prog_.setUniformValue("uMVP", mvp);

		// Clip plane (see inside hollow structures): the KEPT half-space is dot(pos,n)+d ≥ 0. Every scene
		// vertex shader writes gl_ClipDistance[0] = dot(pos, uClipPlane), but it only cuts geometry while
		// GL_CLIP_DISTANCE0 is enabled — enabled per SOLID draw below (mesh, voxels) and left OFF for
		// the flow slice + arrows + reference geometry, so the flow field inside the cavity stays visible.
		const QVector4D clipPlane = computeClipPlane();
		prog_.setUniformValue("uClipPlane", clipPlane);
		glDisable(GL_CLIP_DISTANCE0);

		// Slice mesh (per-vertex CUDA colours). Hidden by the "Show slice" toggle. Also suppressed while the
		// grid overlay is on, so "Show grid on slice" shows the bare cell wireframe (not the filled colours).
		// Never clipped.
		if (show_slice_ && !show_grid_)
		{
			prog_.setUniformValue("uFlat", 0);
			glBindVertexArray(slice_vao_);
			glDrawElements(GL_TRIANGLES, index_count_, GL_UNSIGNED_INT, (void*)0);
		}

		// Domain wireframe (flat grey).
		prog_.setUniformValue("uFlat", 1);
		prog_.setUniformValue("uColor", QVector4D(0.55f, 0.58f, 0.62f, 1.0f));
		glBindVertexArray(box_vao_);
		glDrawArrays(GL_LINES, 0, box_vertex_count_);
		if(amr_debug_vertex_count_>0){prog_.setUniformValue("uColor",QVector4D(0.10f,0.85f,1.0f,1.0f));glBindVertexArray(amr_debug_vao_);glDrawArrays(GL_LINES,0,amr_debug_vertex_count_);}
		if(eb_debug_vertex_count_>0){glLineWidth(1.5f);prog_.setUniformValue("uColor",QVector4D(1.0f,0.25f,0.08f,1.0f));glBindVertexArray(eb_debug_vao_);glDrawArrays(GL_LINES,0,eb_debug_vertex_count_);glLineWidth(1.0f);}

		glBindVertexArray(0);
		prog_.release();

		// Loaded STEP model (lit, its own shader + model transform). Hidden by the "Show model" toggle.
		if (show_model_ && has_mesh_ && mesh_index_count_ > 0)
		{
			mesh_prog_.bind();
			mesh_prog_.setUniformValue("uEye", camera_.eye());
			mesh_prog_.setUniformValue("uClipPlane", clipPlane);
			mesh_prog_.setUniformValue("uBaseColor", QVector4D(0.74f, 0.71f, 0.66f, 1.0f));
			mesh_prog_.setUniformValue("uUseVertexColor", 0); // STEP model: flat base colour
			mesh_prog_.setUniformValue("uUseTriangleColor", has_mesh_triangle_colours_ ? 1 : 0);
			if (has_mesh_triangle_colours_) glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, mesh_triangle_colour_ssbo_);
			if (clip_enabled_) glEnable(GL_CLIP_DISTANCE0); // the STEP model is a SOLID ⇒ clipped
			glBindVertexArray(mesh_vao_);
			const QMatrix4x4 model = modelMatrix();
			mesh_prog_.setUniformValue("uMVP", mvp * model);
			mesh_prog_.setUniformValue("uModel", model);
			glDrawElements(GL_TRIANGLES, mesh_index_count_, GL_UNSIGNED_INT, (void*)0);
			if (has_mesh_triangle_colours_) glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, 0);
			glBindVertexArray(0);
			glDisable(GL_CLIP_DISTANCE0);
			mesh_prog_.release();
		}

		// Voxelization overlay (exposed staircase faces) of the rigid solid mask (obstacle, dark-blue).
		// Reuses the lit mesh shader (identity model).
		if (has_vox_ && show_solid_vox_ && vox_vertex_count_ > 0)
		{
			QMatrix4x4 ident;
			mesh_prog_.bind();
			mesh_prog_.setUniformValue("uMVP", mvp);
			mesh_prog_.setUniformValue("uModel", ident);
			mesh_prog_.setUniformValue("uEye", camera_.eye());
			mesh_prog_.setUniformValue("uClipPlane", clipPlane);
			if (clip_enabled_) glEnable(GL_CLIP_DISTANCE0); // voxel solids ⇒ clipped
			glBindVertexArray(vox_vao_);
			// Cp colouring: per-vertex ramp colours when built (checkbox ON + a Cp field published); else the
			// uniform dark-blue solid. The overlay is re-extracted on toggle so has_vox_cp_ tracks the choice.
			const bool cp_shade = colour_by_cp_ && has_vox_cp_;
			mesh_prog_.setUniformValue("uUseVertexColor", cp_shade ? 1 : 0);
			mesh_prog_.setUniformValue("uUseTriangleColor", 0);
			mesh_prog_.setUniformValue("uBaseColor", QVector4D(0.13f, 0.28f, 0.68f, 1.0f)); // dark-blue solid (uUseVertexColor==0)
			glDrawArrays(GL_TRIANGLES, 0, vox_vertex_count_);
			glBindVertexArray(0);
			glDisable(GL_CLIP_DISTANCE0);
			mesh_prog_.release();
		}

		// World-origin XYZ triad + camera-aligned corner gizmo (depth-tested world axes; gizmo on top).
		drawAxes(mvp);

		// Grid overlay: cell-boundary lines where the grid meets the slice plane (shows the graded mesh).
		drawGrid(mvp);

		// Clip-plane visualisation (translucent quad + outline showing where the cut is). Drawn with the
		// clip test OFF so the plane itself is not clipped; only shown while the feature is enabled.
		if (clip_enabled_) drawClipPlaneViz(mvp, clipPlane);

		// Animated flow arrows (drawn last: blended, no depth write). Never clipped (flow stays visible).
		drawArrows(mvp);

		// Grid-seeded streaklines (blended line strips, no depth write; depth-tested against solids so
		// they read as 3D). Never clipped, like the arrows. Drawn after the arrows so both flow layers
		// blend over the scene.
		drawTracers(mvp);

		// Model-placement gizmo (fluid-viewer model): the manipulator glyph, always on top (depth test
		// off) so it stays grabbable. No-op unless a mode is selected and a gizmo-editable model exists.
		drawGizmo(mvp);

		// Hand a CLEAN GL state back to QPainter. Its text draws via a glyph-cache FBO that renders
		// nothing if depth test / stray bindings are left enabled — so the legend numbers + axis labels
		// silently vanish (solid fillRects still show). begin/endNativePainting does NOT reliably restore
		// this here (the depth-state gotcha), so reset it explicitly before the 2D overlay.
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_BLEND);
		glDisable(GL_CLIP_DISTANCE0); // never leave the clip test on for Qt's paint engine
		glDepthMask(GL_TRUE);
		glActiveTexture(GL_TEXTURE0);
		glBindVertexArray(0);

		painter.endNativePainting();
		drawLegendWith(painter);
		drawAxesLabels(painter); // X/Y/Z + metre tick labels (same QPainter-after-GL path as the legend)

		// fps bookkeeping.
		++fps_frames_;
		++frames_total_;
		double el = fps_timer_.elapsed() / 1000.0;
		if (el >= 0.5)
		{
			avg_fps_ = fps_frames_ / el;
			fps_frames_ = 0;
			fps_timer_.restart();
			emit fpsUpdated(avg_fps_);
		}
	}

	void SliceViewer::mousePressEvent(QMouseEvent* e)
	{
		drag_btn_ = e->button();
		last_mouse_ = e->position().toPoint();
		// Left-click on any gizmo handle grabs it (and suppresses the camera orbit for this drag). The
		// grabbed handle picks the operation dynamically — translate / rotate / scale / uniform-scale.
		if (e->button() == Qt::LeftButton && gizmoEditable() && gizmo_on_)
		{
			int op = -1, axis = -1;
			if (pickGizmo(e->position(), op, axis))
			{
				beginGizmoDrag(op, axis, e->position());
				drag_btn_ = Qt::NoButton; // don't orbit while dragging the manipulator
				update();
			}
		}
	}

	void SliceViewer::mouseMoveEvent(QMouseEvent* e)
	{
		const QPointF posf = e->position();
		const QPoint p = posf.toPoint();

		if (gz_drag_op_ >= 0) // dragging a gizmo handle
		{
			dragGizmo(posf, QPointF(last_mouse_));
			last_mouse_ = p;
			update();
			return;
		}
		if (drag_btn_ == Qt::NoButton) // hover: highlight the handle under the cursor (mouse tracking)
		{
			if (gizmoEditable() && gizmo_on_)
			{
				int op = -1, axis = -1;
				pickGizmo(posf, op, axis);
				if (op != gz_hover_op_ || axis != gz_hover_axis_) { gz_hover_op_ = op; gz_hover_axis_ = axis; update(); }
			}
			last_mouse_ = p;
			return;
		}

		const QPoint d = p - last_mouse_;
		last_mouse_ = p;
		if (drag_btn_ == Qt::LeftButton)
			camera_.orbit(-d.x() * 0.4f, d.y() * 0.4f);
		else if (drag_btn_ == Qt::RightButton || drag_btn_ == Qt::MiddleButton)
			camera_.pan((float)d.x(), (float)d.y());
		update();
	}

	void SliceViewer::mouseReleaseEvent(QMouseEvent* e)
	{
		if (gz_drag_op_ >= 0)
		{
			gz_drag_op_ = -1;
			gz_drag_axis_ = -1;
			fabric_bvh_dirty_ = true;
			emit modelPlacementChanged(); // re-voxelization happens on the next Apply, not per-drag
			update();
		}
		drag_btn_ = Qt::NoButton;
		Q_UNUSED(e);
	}

	void SliceViewer::wheelEvent(QWheelEvent* e)
	{
		float steps = e->angleDelta().y() / 120.0f;
		camera_.zoom(std::pow(1.15f, steps));
		update();
	}
}
