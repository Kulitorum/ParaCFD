// slice_viewer.cpp — see slice_viewer.h. GL 4.3 core rendering of a CUDA-filled slice.
#include "gui/slice_viewer.h"

#include "core/geometry/model_placement.h"
#include "gui/colormap.h"
#include "gui/gl_thread_check.h"
#include "gui/paraglider_sim_worker.h"

#include <QMatrix3x3>
#include <QFontMetrics>
#include <QMouseEvent>
#include <QPainter>
#include <QVector2D>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <set>
#include <vector>

namespace paracfd::gui
{
	namespace
	{
		constexpr std::size_t issue_kind_index(GeometryIssueKind kind)
		{
			return static_cast<std::size_t>(kind);
		}

		bool valid_issue_kind(GeometryIssueKind kind)
		{
			return issue_kind_index(kind) < static_cast<std::size_t>(GeometryIssueKind::Count);
		}

		QVector4D issue_colour(GeometryIssueKind kind)
		{
			switch (kind)
			{
			case GeometryIssueKind::UnexpectedGap:       return {1.00f, 0.15f, 0.08f, 1.00f};
			case GeometryIssueKind::RepairedSeam:        return {1.00f, 0.48f, 0.05f, 0.95f};
			case GeometryIssueKind::UnsharedSeam:        return {1.00f, 0.68f, 0.08f, 0.95f};
			case GeometryIssueKind::IsolatedFabricComponent:return {1.00f, 0.08f, 0.32f, 1.00f};
			case GeometryIssueKind::IntentionalOpening:  return {0.05f, 0.88f, 1.00f, 0.95f};
			case GeometryIssueKind::UnclassifiedOpening:return {0.18f, 0.66f, 0.92f, 0.95f};
			default:                                      return {1.00f, 0.15f, 0.08f, 1.00f};
			}
		}

		float issue_width_pixels(GeometryIssueKind kind)
		{
			switch (kind)
			{
			case GeometryIssueKind::UnexpectedGap:
			case GeometryIssueKind::IsolatedFabricComponent: return 5.0f;
			case GeometryIssueKind::RepairedSeam:
			case GeometryIssueKind::UnsharedSeam: return 3.5f;
			default: return 3.0f;
			}
		}

		QString issue_label(GeometryIssueKind kind)
		{
			switch (kind)
			{
			case GeometryIssueKind::UnexpectedGap:        return QStringLiteral("Unexpected gap");
			case GeometryIssueKind::RepairedSeam:         return QStringLiteral("Repaired seam");
			case GeometryIssueKind::UnsharedSeam:         return QStringLiteral("Unshared seam");
			case GeometryIssueKind::IsolatedFabricComponent:return QStringLiteral("Detached boundaries");
			case GeometryIssueKind::IntentionalOpening:   return QStringLiteral("Intentional opening");
			case GeometryIssueKind::UnclassifiedOpening: return QStringLiteral("Unclassified opening");
			default:                                      return QStringLiteral("Geometry diagnostic");
			}
		}

		std::vector<float> debug_box_lines(const std::vector<std::array<float, 6>>& boxes)
		{
			std::vector<float> lines;lines.reserve(boxes.size()*12*2*3);
			std::set<std::array<float,6>> unique_lines;
			static constexpr int edges[12][2]={{0,1},{0,2},{0,4},{1,3},{1,5},{2,3},{2,6},{3,7},{4,5},{4,6},{5,7},{6,7}};
			for(const auto& b:boxes)
			{
				const float p[8][3]={{b[0],b[1],b[2]},{b[3],b[1],b[2]},{b[0],b[4],b[2]},{b[3],b[4],b[2]},{b[0],b[1],b[5]},{b[3],b[1],b[5]},{b[0],b[4],b[5]},{b[3],b[4],b[5]}};
				for(const auto& edge:edges)
				{
					std::array<float,3> a{p[edge[0]][0],p[edge[0]][1],p[edge[0]][2]};
					std::array<float,3> c{p[edge[1]][0],p[edge[1]][1],p[edge[1]][2]};
					for(float& coordinate:a)if(coordinate==0.0f)coordinate=0.0f;
					for(float& coordinate:c)if(coordinate==0.0f)coordinate=0.0f;
					if(c<a)std::swap(a,c);
					unique_lines.insert({a[0],a[1],a[2],c[0],c[1],c[2]});
				}
			}
			for(const auto& line:unique_lines)lines.insert(lines.end(),line.begin(),line.end());
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
layout(location=2) in vec3 aColor; // optional per-vertex debug colour
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
uniform int  uUseVertexColor; // 0 = flat uBaseColor; nonzero = per-vertex colour
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
uniform float uSize;                 // uniform glyph scale (both length and lateral width)
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
	float L = uLen * (0.35 + 0.9 * aSpeed) * uSize; // length grows with speed and the uniform size control
	vec3 world = aPos + T * (aGlyph.x * L) + N * (aGlyph.y * L);
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

		// Static CAD-diagnostic ribbons. Input points/tangents are model-local so the exact same
		// placement matrix as the STEP surface keeps them registered while the wing is manipulated.
		const char* kGeometryIssueVert = R"(#version 430 core
layout(location=0) in vec4 aPos;   // xyz model-local, w = ribbon side (-1 / +1)
layout(location=1) in vec3 aTan;   // model-local polyline tangent
layout(location=2) in vec4 aColor;
uniform mat4 uVP;
uniform mat4 uModel;
uniform vec2 uViewport;
uniform float uHalfPx;
uniform vec4 uClipPlane;
out vec4 vColor;
void main()
{
	vColor = aColor;
	vec4 world = uModel * vec4(aPos.xyz, 1.0);
	vec4 worldT = uModel * vec4(aPos.xyz + aTan, 1.0);
	vec4 clip = uVP * world;
	vec4 clipT = uVP * worldT;
	vec2 s0 = clip.xy / clip.w;
	vec2 s1 = clipT.xy / clipT.w;
	vec2 dir = (s1 - s0) * uViewport;
	dir = (dot(dir, dir) > 1e-12) ? normalize(dir) : vec2(1.0, 0.0);
	vec2 nrm = vec2(-dir.y, dir.x);
	vec2 offNdc = nrm * (aPos.w * uHalfPx) * 2.0 / uViewport;
	clip.xy += offNdc * clip.w;
	gl_ClipDistance[0] = dot(world, uClipPlane);
	gl_Position = clip;
}
)";
		const char* kGeometryIssueFrag = R"(#version 430 core
in vec4 vColor;
out vec4 fragColor;
void main() { fragColor = vColor; }
)";

		const char* kVolumeVert = R"(#version 430 core
out vec2 vNdc;
void main()
{
	vec2 p = gl_VertexID == 0 ? vec2(-1.0,-1.0) : gl_VertexID == 1 ? vec2(3.0,-1.0) : vec2(-1.0,3.0);
	vNdc = p;
	gl_Position = vec4(p,0.0,1.0);
}
)";
		const char* kVolumeFrag = R"(#version 430 core
in vec2 vNdc;
out vec4 fragColor;
uniform sampler3D uVolume;
uniform mat4 uInvVP;
uniform vec3 uBoxMax;
uniform vec4 uClipPlane;
uniform float uThreshold;
uniform float uFocus;
uniform float uOpacity;
uniform float uCellSize;
uniform int uSteps;

bool rayBox(vec3 ro,vec3 rd,out float t0,out float t1)
{
	vec3 safeDir=vec3(rd.x<0.0?-max(abs(rd.x),1e-7):max(abs(rd.x),1e-7),
		rd.y<0.0?-max(abs(rd.y),1e-7):max(abs(rd.y),1e-7),
		rd.z<0.0?-max(abs(rd.z),1e-7):max(abs(rd.z),1e-7));
	vec3 a=(vec3(0.0)-ro)/safeDir,b=(uBoxMax-ro)/safeDir;
	vec3 lo=min(a,b),hi=max(a,b);
	t0=max(max(lo.x,lo.y),lo.z);t1=min(min(hi.x,hi.y),hi.z);
	t0=max(t0,0.0);return t1>t0;
}
void main()
{
	vec4 na=uInvVP*vec4(vNdc,-1.0,1.0),fa=uInvVP*vec4(vNdc,1.0,1.0);
	vec3 ro=na.xyz/na.w,farPoint=fa.xyz/fa.w,rd=normalize(farPoint-ro);
	float t0,t1;if(!rayBox(ro,rd,t0,t1))discard;
	float dt=(t1-t0)/float(max(uSteps,1));vec4 accum=vec4(0.0);
	for(int step=0;step<512;++step)
	{
		if(step>=uSteps||accum.a>0.985)break;
		vec3 world=ro+rd*(t0+(float(step)+0.5)*dt);
		if(dot(vec4(world,1.0),uClipPlane)<0.0)continue;
		float cp=texture(uVolume,clamp(world/uBoxMax,vec3(0.0),vec3(1.0))).r;
		float magnitude=abs(cp);if(magnitude<uThreshold)continue;
		float strength=smoothstep(uThreshold,max(uThreshold+1e-6,uFocus),magnitude);
		vec3 cold=mix(vec3(0.08,0.22,0.72),vec3(0.12,0.78,1.0),0.35+0.65*strength);
		vec3 hot=mix(vec3(0.88,0.18,0.06),vec3(1.0,0.78,0.12),0.25*(1.0-strength));
		vec3 color=cp<0.0?cold:hot;
		float alpha=1.0-exp(-0.045*uOpacity*strength*dt/max(uCellSize,1e-6));
		accum.rgb+=(1.0-accum.a)*alpha*color;accum.a+=(1.0-accum.a)*alpha;
	}
	if(accum.a<0.005)discard;fragColor=accum;
}
)";
	}

	SliceViewer::SliceViewer(QWidget* parent) : QOpenGLWidget(parent)
	{
		setFocusPolicy(Qt::StrongFocus);
		setMouseTracking(true); // hover-highlight the gizmo handle under the cursor (no button held)
		arrows_.set_count(arrow_density_);
	}

	void SliceViewer::setParagliderWorker(ParagliderSimWorker* w)
	{
		paraglider_worker_ = w;
		if(!w){simulation_running_=simulation_auto_paused_=simulation_settling_ready_=false;simulation_pause_label_.clear();simulation_error_label_.clear();}
		range_valid_ = false;
		resetFlowAnimation();
		update();
	}

	void SliceViewer::resetFlowAnimation()
	{
		arrows_.reset();tracers_.reset();arrow_draw_count_=0;tracer_strips_=0;
		tracer_firsts_.clear();tracer_counts_.clear();
		arrow_clock_started_=false;tracer_clock_started_=false;
		update();
	}

	void SliceViewer::setSimulationState(bool running,bool auto_paused,bool auto_pause_enabled,
		bool settling_ready,double settling_score,double flow_throughs,
		const QString& automatic_pause_label,const QString& error_label,bool mean_force_ready,
		double mean_force_drift,double mean_force_tolerance,double maximum_flow_throughs)
	{
		simulation_running_=running;simulation_auto_paused_=auto_paused;simulation_auto_pause_enabled_=auto_pause_enabled;simulation_settling_ready_=settling_ready;simulation_settling_score_=settling_score;simulation_flow_throughs_=flow_throughs;simulation_pause_label_=automatic_pause_label;simulation_error_label_=error_label;simulation_mean_force_ready_=mean_force_ready;simulation_mean_force_drift_=mean_force_drift;simulation_mean_force_tolerance_=mean_force_tolerance;simulation_maximum_flow_throughs_=maximum_flow_throughs;update();
	}

	void SliceViewer::setPlaybackState(bool active,std::size_t frame,std::size_t frame_count,
		double physical_time,bool discontinuity)
	{
		const bool changed=playback_active_!=active;playback_active_=active;
		playback_frame_=frame;playback_frame_count_=frame_count;
		playback_physical_time_=std::max(0.0,physical_time);
		slice_sample_dirty_=true;iso_surface_dirty_=true;range_valid_=false;
		if(changed||discontinuity){resetFlowAnimation();probe_valid_=false;}
		update();
	}

	void SliceViewer::setSimulationProgress(long long step,double physical_time,double wall_time)
	{
		simulation_step_=step;simulation_physical_time_=std::max(0.0,physical_time);simulation_wall_time_=std::max(0.0,wall_time);update();
	}

	void SliceViewer::setShowArrows(bool on)
	{
		show_arrows_ = on;
		if (on) arrows_.reset(); // re-seed so they don't reappear frozen
		update();
	}

	void SliceViewer::setShowTracers(bool on)
	{
		show_tracers_ = on;
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
		setArrowMode(three_d?0:1);
	}

	void SliceViewer::setArrowMode(int mode)
	{
		mode=std::clamp(mode,0,2);if(arrow_mode_==mode)return;arrow_mode_=mode;arrow_3d_=mode==0;
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
		if (geometry_issue_vbo_) glDeleteBuffers(1, &geometry_issue_vbo_);
		if (geometry_error_vbo_) glDeleteBuffers(1, &geometry_error_vbo_);
		if (arrow_glyph_vbo_) glDeleteBuffers(1, &arrow_glyph_vbo_);
		if (arrow_inst_vbo_) glDeleteBuffers(1, &arrow_inst_vbo_);
		if (tracer_vbo_) glDeleteBuffers(1, &tracer_vbo_);
		if (axis_vbo_) glDeleteBuffers(1, &axis_vbo_);
		if (corner_vbo_) glDeleteBuffers(1, &corner_vbo_);
		if (clip_vbo_) glDeleteBuffers(1, &clip_vbo_);
		if (gizmo_vbo_) glDeleteBuffers(1, &gizmo_vbo_);
		if (iso_vbo_) glDeleteBuffers(1,&iso_vbo_);
		if (pressure_force_inst_vbo_) glDeleteBuffers(1,&pressure_force_inst_vbo_);
		if (volume_texture_) glDeleteTextures(1,&volume_texture_);
		if (slice_vao_) glDeleteVertexArrays(1, &slice_vao_);
		if (box_vao_) glDeleteVertexArrays(1, &box_vao_);
		if (grid_vao_) glDeleteVertexArrays(1, &grid_vao_);
		if (mesh_vao_) glDeleteVertexArrays(1, &mesh_vao_);
		if (amr_debug_vao_) glDeleteVertexArrays(1, &amr_debug_vao_);
		if (eb_debug_vao_) glDeleteVertexArrays(1, &eb_debug_vao_);
		if (geometry_issue_vao_) glDeleteVertexArrays(1, &geometry_issue_vao_);
		if (geometry_error_vao_) glDeleteVertexArrays(1, &geometry_error_vao_);
		if (arrow_vao_) glDeleteVertexArrays(1, &arrow_vao_);
		if (tracer_vao_) glDeleteVertexArrays(1, &tracer_vao_);
		if (axis_vao_) glDeleteVertexArrays(1, &axis_vao_);
		if (corner_vao_) glDeleteVertexArrays(1, &corner_vao_);
		if (clip_vao_) glDeleteVertexArrays(1, &clip_vao_);
		if (gizmo_vao_) glDeleteVertexArrays(1, &gizmo_vao_);
		if (iso_vao_) glDeleteVertexArrays(1,&iso_vao_);
		if (volume_vao_) glDeleteVertexArrays(1,&volume_vao_);
		if (pressure_force_vao_) glDeleteVertexArrays(1,&pressure_force_vao_);
		doneCurrent();
	}

	void SliceViewer::setInfo(const SimInfo& info)
	{
		info_ = info;
		if (!(info_.finest_h > 0.0)) info_.finest_h = info_.h;
		have_info_ = true;
		updateSliceResolution();
		camera_.frameDomain((float)info.Lx, (float)info.Ly, (float)info.Lz);
		setDefaultPlane();
		updateRange();
		range_valid_ = false; // re-snap the auto-range on the new domain
		slice_sample_dirty_=true;slice_sample_generation_=0;
		auto_speed_max_ = 1.8f * (float)info.U;
		geometry_dirty_ = true;
		axes_dirty_ = true; // rebuild the world triad + ticks for the new domain extents
		grid_dirty_ = true;
		// Base arrow length ~2% of the domain diagonal (a few cells) — visible over both the wake
		// slice and the model without cluttering. Re-seed the tracers for the new domain.
		float diag = std::sqrt((float)(info.Lx * info.Lx + info.Ly * info.Ly + info.Lz * info.Lz));
		arrow_len_ = std::max(4.0f * (float)info.h, 0.02f * diag);
		resetFlowAnimation(); // re-seed all animated flow state for the new domain
		iso_volume_.reset();cp_volume_.reset();volume_texture_ready_=false;volume_texture_generation_=0;iso_surface_dirty_=true;
	}

	void SliceViewer::frameDomainView()
	{
		if(!have_info_)return;
		camera_.frameDomain((float)info_.Lx,(float)info_.Ly,(float)info_.Lz);
		update();
	}

	void SliceViewer::frameThinYDebugView()
	{
		if(!have_info_)return;camera_.frameDomain((float)info_.Lx,(float)info_.Ly,(float)info_.Lz);camera_.setOrientation(90.0f,0.0f);update();
	}

	void SliceViewer::frameSliceView()
	{
		if(!have_info_)return;
		QVector3D lo(0,0,0),hi((float)info_.Lx,(float)info_.Ly,(float)info_.Lz);
		const int a=(int)axis_;const float L=a==0?(float)info_.Lx:a==1?(float)info_.Ly:(float)info_.Lz;
		lo[a]=hi[a]=plane_frac_*L;camera_.frameBounds(lo,hi,1.15f);
		if(axis_==Axis::X)camera_.setOrientation(0,0);
		else if(axis_==Axis::Y)camera_.setOrientation(90,0);
		else camera_.setOrientation(0,90);
		update();
	}

	bool SliceViewer::frameWingView()
	{
		if(!gz_valid_)return false;
		const QMatrix4x4 transform=modelMatrix();
		QVector3D lo,hi;
		bool first=true;
		for(int ix=0;ix<2;++ix)for(int iy=0;iy<2;++iy)for(int iz=0;iz<2;++iz)
		{
			const QVector3D source(ix?gz_bbox_max_.x():gz_bbox_min_.x(),
				iy?gz_bbox_max_.y():gz_bbox_min_.y(),iz?gz_bbox_max_.z():gz_bbox_min_.z());
			const QVector3D world=transform.map(source);
			if(first){lo=hi=world;first=false;}
			else
			{
				lo.setX(std::min(lo.x(),world.x()));lo.setY(std::min(lo.y(),world.y()));lo.setZ(std::min(lo.z(),world.z()));
				hi.setX(std::max(hi.x(),world.x()));hi.setY(std::max(hi.y(),world.y()));hi.setZ(std::max(hi.z(),world.z()));
			}
		}
		if(first)return false;
		camera_.frameBounds(lo,hi);
		update();
		return true;
	}

	void SliceViewer::setReferenceU(double U)
	{
		if (!have_info_) return;
		info_.U = U;
		if (!auto_range_) updateRange(); // refresh the fixed per-field range to the new current
		slice_sample_dirty_=true;
		update();
	}

	void SliceViewer::updateSliceResolution()
	{
		if (!have_info_) return;
		const double h = info_.finest_h > 0.0 ? info_.finest_h : info_.h;
		const double extent_u = axis_ == Axis::X ? info_.Ly : info_.Lx;
		const double extent_v = axis_ == Axis::Z ? info_.Ly : info_.Lz;
		auto samples = [&](double extent)
		{
			if (!(extent > 0.0) || !(h > 0.0)) return 2;
			return std::clamp(static_cast<int>(std::ceil(extent / h)) + 1, 2, kMaxSliceAxis);
		};
		int nu = samples(extent_u), nv = samples(extent_v);
		const double vertex_count = static_cast<double>(nu) * nv;
		if (vertex_count > kMaxSliceVertices)
		{
			const double scale = std::sqrt(static_cast<double>(kMaxSliceVertices) / vertex_count);
			nu = std::max(2, static_cast<int>(std::floor((nu - 1) * scale)) + 1);
			nv = std::max(2, static_cast<int>(std::floor((nv - 1) * scale)) + 1);
		}
		if (nu == slice_nu_ && nv == slice_nv_) return;
		slice_nu_ = nu;
		slice_nv_ = nv;
		geometry_dirty_ = true;
		std::fprintf(stderr, "[viewer] AMR slice %dx%d at finest h=%.6g m\n", slice_nu_, slice_nv_, h);
	}

	void SliceViewer::setDefaultPlane()
	{
		plane_frac_ = 0.5f;
		if(clip_follows_slice_){clip_mode_=(int)axis_;clip_frac_=plane_frac_;}
	}

	void SliceViewer::setPlaneFraction(float frac)
	{
		const float next=std::min(1.0f,std::max(0.0f,frac));if(std::abs(next-plane_frac_)<1e-7f)return;plane_frac_=next;
		if(clip_follows_slice_)clip_frac_=plane_frac_;
		geometry_dirty_ = true;
		grid_dirty_ = true; // the grid overlay lives on the plane → moves with it
		slice_sample_dirty_=true;
		emit planeFractionChanged(plane_frac_);
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
		surface_cp_plus_.clear();surface_cp_minus_.clear();surface_delta_cp_.clear();clearTrianglePressureForces();probe_valid_=false;
		// Seed host-side placement before the deferred GL upload.
		gz_bbox_min_ = QVector3D(mesh.bbox_min[0], mesh.bbox_min[1], mesh.bbox_min[2]);
		gz_bbox_max_ = QVector3D(mesh.bbox_max[0], mesh.bbox_max[1], mesh.bbox_max[2]);
		computeDefaultXform();
		gz_drag_op_ = -1; gz_drag_axis_ = -1;
		fabric_bvh_dirty_ = true;
		pressure_force_dirty_=true;
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
		surface_cp_plus_.clear();surface_cp_minus_.clear();surface_delta_cp_.clear();clearTrianglePressureForces();probe_valid_=false;
		fabric_bvh_.reset();
		fabric_bvh_dirty_ = false;
		mesh_index_count_ = 0;
		mesh_override_ = false; // a fresh model (fluid viewer) is gizmo-editable again
		gz_valid_ = false;
		gz_drag_op_ = -1; gz_drag_axis_ = -1;
		gz_hover_op_ = -1; gz_hover_axis_ = -1;
		update();
	}


	void SliceViewer::setMeshTranslate(double tx, double ty, double tz)
	{
		mesh_override_mat_ = QMatrix4x4();
		mesh_override_mat_.translate((float)tx, (float)ty, (float)tz);
		mesh_override_ = true;
		fabric_bvh_dirty_ = true;
		pressure_force_dirty_ = true;
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
		pressure_force_dirty_=true;
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
		// Generic viewer fallback placement. The paraglider window normally supplies an
		// explicit aerodynamic placement immediately after loading the mesh.
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
		pressure_force_dirty_ = true;
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
		pressure_force_dirty_ = true;
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
		pressure_force_dirty_=true;
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

	void SliceViewer::setGeometryIssuePolylines(std::vector<GeometryIssuePolyline> polylines)
	{
		constexpr int stride = 11; // pos.xyz+side, tangent.xyz, rgba
		geometry_issue_vertices_.clear();
		geometry_issue_polyline_counts_.fill(0);
		geometry_issue_max_gap_m_.fill(std::nullopt);
		for (auto& firsts : geometry_issue_firsts_) firsts.clear();
		for (auto& counts : geometry_issue_counts_) counts.clear();

		// Group vertices by kind so each warning class can use a distinct screen-space width without
		// issuing a draw call for every CAD edge. Invalid coordinates split a producer polyline rather
		// than accidentally bridging across missing diagnostic data.
		for (std::size_t kind_index = 0; kind_index < kGeometryIssueKindCount; ++kind_index)
		{
			const auto kind = static_cast<GeometryIssueKind>(kind_index);
			const QVector4D colour = issue_colour(kind);
			for (const GeometryIssuePolyline& polyline : polylines)
			{
				if (!valid_issue_kind(polyline.kind) || issue_kind_index(polyline.kind) != kind_index) continue;
				bool emitted = false;
				std::vector<QVector3D> run;
				run.reserve(polyline.points.size());
				auto emit_run = [&]
				{
					if (run.size() < 2) { run.clear(); return; }
					const int first = static_cast<int>(geometry_issue_vertices_.size() / stride);
					for (std::size_t point = 0; point < run.size(); ++point)
					{
						const QVector3D& p = run[point];
						QVector3D tangent = run[std::min(point + 1, run.size() - 1)]
							- run[point == 0 ? 0 : point - 1];
						if (tangent.lengthSquared() <= 1e-20f)
							tangent = point + 1 < run.size() ? run[point + 1] - p : p - run[point - 1];
						if (tangent.lengthSquared() <= 1e-20f) tangent = QVector3D(1, 0, 0);
						else tangent.normalize();
						for (float side : {-1.0f, 1.0f})
							geometry_issue_vertices_.insert(geometry_issue_vertices_.end(),
								{p.x(), p.y(), p.z(), side, tangent.x(), tangent.y(), tangent.z(),
								 colour.x(), colour.y(), colour.z(), colour.w()});
					}
					geometry_issue_firsts_[kind_index].push_back(first);
					geometry_issue_counts_[kind_index].push_back(static_cast<int>(2 * run.size()));
					emitted = true;
					run.clear();
				};

				for (const auto& point : polyline.points)
				{
					if (!std::isfinite(point[0]) || !std::isfinite(point[1]) || !std::isfinite(point[2]))
					{
						emit_run();
						continue;
					}
					const QVector3D p(point[0], point[1], point[2]);
					if (!run.empty() && (p - run.back()).lengthSquared() <= 1e-20f) continue;
					run.push_back(p);
				}
				emit_run();
				if (!emitted) continue;
				++geometry_issue_polyline_counts_[kind_index];
				if (polyline.measured_gap_m && std::isfinite(*polyline.measured_gap_m)
					&& *polyline.measured_gap_m >= 0.0)
				{
					auto& maximum = geometry_issue_max_gap_m_[kind_index];
					maximum = maximum ? std::max(*maximum, *polyline.measured_gap_m)
						: *polyline.measured_gap_m;
				}
			}
		}
		geometry_issue_upload_pending_ = true;
		update();
	}

	void SliceViewer::clearGeometryIssuePolylines()
	{
		setGeometryIssuePolylines({});
	}

	void SliceViewer::setGeometryErrorDiagnosticSegments(
		std::vector<std::array<float, 6>> world_segments,const QString& overlay_text,
		bool frame_segments)
	{
		constexpr int stride=11;geometry_error_vertices_.clear();geometry_error_firsts_.clear();
		geometry_error_counts_.clear();geometry_error_bounds_valid_=false;
		const QVector4D colour(1.0f,0.04f,0.02f,1.0f);
		for(const auto& segment:world_segments)
		{
			bool finite=true;for(float coordinate:segment)finite=finite&&std::isfinite(coordinate);if(!finite)continue;
			const QVector3D a(segment[0],segment[1],segment[2]),b(segment[3],segment[4],segment[5]);
			QVector3D tangent=b-a;if(tangent.lengthSquared()<=1e-20f)continue;tangent.normalize();
			const int first=static_cast<int>(geometry_error_vertices_.size()/stride);
			for(const QVector3D& point:{a,b})for(float side:{-1.0f,1.0f})
				geometry_error_vertices_.insert(geometry_error_vertices_.end(),
					{point.x(),point.y(),point.z(),side,tangent.x(),tangent.y(),tangent.z(),
					 colour.x(),colour.y(),colour.z(),colour.w()});
			geometry_error_firsts_.push_back(first);geometry_error_counts_.push_back(4);
			if(!geometry_error_bounds_valid_){geometry_error_bounds_min_=geometry_error_bounds_max_=a;geometry_error_bounds_valid_=true;}
			for(const QVector3D& point:{a,b})
			{
				geometry_error_bounds_min_.setX(std::min(geometry_error_bounds_min_.x(),point.x()));
				geometry_error_bounds_min_.setY(std::min(geometry_error_bounds_min_.y(),point.y()));
				geometry_error_bounds_min_.setZ(std::min(geometry_error_bounds_min_.z(),point.z()));
				geometry_error_bounds_max_.setX(std::max(geometry_error_bounds_max_.x(),point.x()));
				geometry_error_bounds_max_.setY(std::max(geometry_error_bounds_max_.y(),point.y()));
				geometry_error_bounds_max_.setZ(std::max(geometry_error_bounds_max_.z(),point.z()));
			}
		}
		geometry_error_overlay_text_=overlay_text.trimmed();
		// A failed-cell diagnostic remains active even when preprocessing could only
		// identify broad local context.  In that case the UI deliberately draws no
		// fabric ribbons instead of claiming the whole candidate region is defective.
		geometry_error_active_=!geometry_error_firsts_.empty()||!geometry_error_overlay_text_.isEmpty();
		if(geometry_error_active_&&geometry_error_overlay_text_.isEmpty())
			geometry_error_overlay_text_=QString("GEOMETRY ERROR DIAGNOSTIC — %1 SEGMENT%2")
				.arg(geometry_error_firsts_.size()).arg(geometry_error_firsts_.size()==1?"":"S");
		geometry_error_upload_pending_=true;
		if(frame_segments&&geometry_error_active_)frameGeometryErrorDiagnosticView();else update();
	}

	void SliceViewer::clearGeometryErrorDiagnostic()
	{
		geometry_error_vertices_.clear();geometry_error_firsts_.clear();geometry_error_counts_.clear();
		geometry_error_overlay_text_.clear();geometry_error_bounds_valid_=false;
		geometry_error_active_=false;geometry_error_upload_pending_=true;update();
	}

	bool SliceViewer::frameGeometryErrorDiagnosticView()
	{
		if(!geometry_error_active_||!geometry_error_bounds_valid_)return false;
		camera_.frameBounds(geometry_error_bounds_min_,geometry_error_bounds_max_,1.35f);update();return true;
	}

	std::size_t SliceViewer::geometryIssueCount(GeometryIssueKind kind) const
	{
		return valid_issue_kind(kind) ? geometry_issue_polyline_counts_[issue_kind_index(kind)] : 0;
	}

	std::optional<double> SliceViewer::geometryIssueMaximumGap(GeometryIssueKind kind) const
	{
		return valid_issue_kind(kind) ? geometry_issue_max_gap_m_[issue_kind_index(kind)] : std::nullopt;
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
		float h=have_info_?std::max(1e-6f,(float)info_.finest_h):1.0f;
		switch (field_)
		{
		case Field::SpeedMag: vmin_ = 0.0f; vmax_ = 1.8f * U; break;
		case Field::VelU: vmin_ = -0.5f * U; vmax_ = 1.7f * U; break;
		case Field::VelV: vmin_ = -0.9f * U; vmax_ = 0.9f * U; break;
		case Field::VelW: vmin_ = -0.9f * U; vmax_ = 0.9f * U; break;
		case Field::PressureDelta: vmin_ = -0.7f * rhoU2; vmax_ = 0.7f * rhoU2; break;
		case Field::PressureCoefficient:vmin_=-2.0f;vmax_=1.0f;break;
		case Field::PressureGradient:vmin_=0;vmax_=rhoU2/h;break;
		case Field::VorticityMagnitude:vmin_=0;vmax_=2.0f*U/h;break;
		case Field::QCriterion:{const float scale=U*U/(h*h);vmin_=-scale;vmax_=scale;break;}
		}
		if (vmax_ <= vmin_) vmax_ = vmin_ + 1.0f;
	}

	void SliceViewer::setAutoRange(bool on)
	{
		auto_range_ = on;
		range_valid_ = false;     // re-snap when toggled back on
		if (!on) updateRange();   // restore the fixed per-field colour range
		slice_sample_dirty_=true;
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

		float tmin = field_is_magnitude(field_) ? 0.0f : fr.field_min;
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
			const char* fn=field_short_name(field_);
			std::fprintf(stderr, "[viewer] auto-range %s: vmin=%.4g vmax=%.4g (raw[%.4g,%.4g]) speed_max=%.4g\n",
				fn, vmin_, vmax_, fr.field_min, fr.field_max, auto_speed_max_);
		}
	}

	SliceParams SliceViewer::currentParams() const
	{
		SliceParams sp;
		sp.grid.nx = info_.nx; sp.grid.ny = info_.ny; sp.grid.nz = info_.nz; sp.grid.h = info_.h;
		// True AMR domain extent supplied by the worker's host visualization snapshot.
		sp.Lx = (float)info_.Lx; sp.Ly = (float)info_.Ly; sp.Lz = (float)info_.Lz;
		sp.axis = axis_;
		float L = (axis_ == Axis::X) ? (float)info_.Lx : (axis_ == Axis::Y) ? (float)info_.Ly : (float)info_.Lz;
		sp.plane_pos = plane_frac_ * L;
		sp.nu = slice_nu_; sp.nv = slice_nv_;
		sp.field = field_;
		sp.vmin = vmin_; sp.vmax = vmax_;
		sp.rho=(float)info_.rho;sp.reference_speed=(float)info_.U;
		return sp;
	}

	void SliceViewer::initializeGL()
	{
		PARACFD_ASSERT_GL_THREAD();
		initializeOpenGLFunctions();

		const char* ver = reinterpret_cast<const char*>(glGetString(GL_VERSION));
		std::fprintf(stderr, "[viewer] GL context: %s\n", ver ? ver : "(null)");

		glClearColor(0.09f, 0.10f, 0.12f, 1.0f);
		glEnable(GL_DEPTH_TEST);

		if (!prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kVert)
			|| !prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kFrag)
			|| !prog_.link())
		{
			std::fprintf(stderr, "[viewer] shader error: %s\n", prog_.log().toUtf8().constData());
		}

		if (!mesh_prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kMeshVert)
			|| !mesh_prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kMeshFrag)
			|| !mesh_prog_.link())
		{
			std::fprintf(stderr, "[viewer] mesh shader error: %s\n", mesh_prog_.log().toUtf8().constData());
		}

		if (!arrow_prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kArrowVert)
			|| !arrow_prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kArrowFrag)
			|| !arrow_prog_.link())
		{
			std::fprintf(stderr, "[viewer] arrow shader error: %s\n", arrow_prog_.log().toUtf8().constData());
		}

		if (!tracer_prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kTracerVert)
			|| !tracer_prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kTracerFrag)
			|| !tracer_prog_.link())
		{
			std::fprintf(stderr, "[viewer] tracer shader error: %s\n", tracer_prog_.log().toUtf8().constData());
		}
		if (!geometry_issue_prog_.addShaderFromSourceCode(QOpenGLShader::Vertex, kGeometryIssueVert)
			|| !geometry_issue_prog_.addShaderFromSourceCode(QOpenGLShader::Fragment, kGeometryIssueFrag)
			|| !geometry_issue_prog_.link())
		{
			std::fprintf(stderr, "[viewer] geometry-issue shader error: %s\n",
				geometry_issue_prog_.log().toUtf8().constData());
		}
		if(!volume_prog_.addShaderFromSourceCode(QOpenGLShader::Vertex,kVolumeVert)
			||!volume_prog_.addShaderFromSourceCode(QOpenGLShader::Fragment,kVolumeFrag)
			||!volume_prog_.link())
			std::fprintf(stderr,"[viewer] volume shader error: %s\n",volume_prog_.log().toUtf8().constData());

		glGenVertexArrays(1, &slice_vao_);
		glGenVertexArrays(1, &box_vao_);
		glGenVertexArrays(1, &grid_vao_);
		glGenVertexArrays(1, &mesh_vao_);
		glGenVertexArrays(1, &amr_debug_vao_);
		glGenVertexArrays(1, &eb_debug_vao_);
		glGenVertexArrays(1, &geometry_issue_vao_);
		glGenVertexArrays(1, &geometry_error_vao_);
		glGenVertexArrays(1, &arrow_vao_);
		glGenVertexArrays(1, &tracer_vao_);
		glGenVertexArrays(1, &axis_vao_);
		glGenVertexArrays(1, &corner_vao_);
		glGenVertexArrays(1, &clip_vao_);
		glGenVertexArrays(1,&iso_vao_);
		glGenVertexArrays(1,&volume_vao_);
		glGenVertexArrays(1,&pressure_force_vao_);
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
		glGenBuffers(1, &geometry_issue_vbo_);
		glGenBuffers(1, &geometry_error_vbo_);
		glGenBuffers(1, &arrow_glyph_vbo_);
		glGenBuffers(1, &arrow_inst_vbo_);
		glGenBuffers(1, &tracer_vbo_);
		glGenBuffers(1, &axis_vbo_);
		glGenBuffers(1, &corner_vbo_);
		glGenBuffers(1, &clip_vbo_);
		glGenBuffers(1,&iso_vbo_);
		glGenBuffers(1,&pressure_force_inst_vbo_);
		glGenTextures(1,&volume_texture_);

		glBindVertexArray(slice_vao_);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, idx_ebo_);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);

		// The axis-dependent slice buffers are sized in buildSliceGeometry() from the
		// domain extent and finest active AMR spacing.
		glBindBuffer(GL_ARRAY_BUFFER, color_vbo_);
		glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ARRAY_BUFFER, pos_vbo_);
		glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);
		glBindVertexArray(0);

		glBindVertexArray(iso_vao_);
		glBindBuffer(GL_ARRAY_BUFFER,iso_vbo_);
		glBufferData(GL_ARRAY_BUFFER,0,nullptr,GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,6*sizeof(float),(void*)0);
		glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,6*sizeof(float),(void*)(3*sizeof(float)));
		glEnableVertexAttribArray(0);glEnableVertexAttribArray(1);glBindVertexArray(0);

		gl_ready_ = true;
		buildBoxGeometry();
		buildArrowGlyph();
		buildTracerBuffers();
		buildGeometryIssueBuffers();
		buildGeometryErrorDiagnosticBuffers();
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

		// Pressure-force glyphs share the immutable arrow mesh but have their own per-instance data.
		glBindVertexArray(pressure_force_vao_);
		glBindBuffer(GL_ARRAY_BUFFER,arrow_glyph_vbo_);
		glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,2*sizeof(float),(void*)0);glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER,pressure_force_inst_vbo_);
		for(int attribute=1;attribute<=4;++attribute)glEnableVertexAttribArray(attribute);
		glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,stride,(void*)0);
		glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,stride,(void*)(3*sizeof(float)));
		glVertexAttribPointer(3,1,GL_FLOAT,GL_FALSE,stride,(void*)(6*sizeof(float)));
		glVertexAttribPointer(4,1,GL_FLOAT,GL_FALSE,stride,(void*)(7*sizeof(float)));
		for(int attribute=1;attribute<=4;++attribute)glVertexAttribDivisor(attribute,1);
		glBindVertexArray(0);
	}

	void SliceViewer::setTriangleSurfaceProbeData(const std::vector<float>& plus_cp,
		const std::vector<float>& minus_cp,const std::vector<float>& delta_cp)
	{
		if(plus_cp.size()!=minus_cp.size()||plus_cp.size()!=delta_cp.size())
		{
			surface_cp_plus_.clear();surface_cp_minus_.clear();surface_delta_cp_.clear();return;
		}
		surface_cp_plus_=plus_cp;surface_cp_minus_=minus_cp;surface_delta_cp_=delta_cp;
	}

	void SliceViewer::setTrianglePressureForces(const std::vector<float>& xyz)
	{
		triangle_pressure_forces_=xyz;pressure_force_dirty_=true;update();
	}

	void SliceViewer::clearTrianglePressureForces()
	{
		triangle_pressure_forces_.clear();pressure_force_instances_.clear();pressure_force_draw_count_=0;
		pressure_force_dirty_=true;update();
	}

	void SliceViewer::buildPressureForceInstances()
	{
		pressure_force_dirty_=false;pressure_force_instances_.clear();pressure_force_draw_count_=0;
		const std::size_t triangles=std::min(fabric_mesh_.triangle_count(),triangle_pressure_forces_.size()/3);
		if(triangles==0)return;
		std::vector<float> magnitude(triangles,0.0f),finite_magnitude;
		for(std::size_t triangle=0;triangle<triangles;++triangle)
		{
			const float x=triangle_pressure_forces_[3*triangle],y=triangle_pressure_forces_[3*triangle+1],z=triangle_pressure_forces_[3*triangle+2];
			const float m=std::sqrt(x*x+y*y+z*z);magnitude[triangle]=std::isfinite(m)?m:0.0f;if(magnitude[triangle]>0)finite_magnitude.push_back(magnitude[triangle]);
		}
		if(finite_magnitude.empty())return;
		const std::size_t percentile=static_cast<std::size_t>(0.95*(finite_magnitude.size()-1));
		std::nth_element(finite_magnitude.begin(),finite_magnitude.begin()+percentile,finite_magnitude.end());
		const float focus=std::max(finite_magnitude[percentile],1e-12f);
		const std::size_t stride=std::max<std::size_t>(1,(triangles+static_cast<std::size_t>(pressure_force_density_)-1)/static_cast<std::size_t>(pressure_force_density_));
		const QMatrix4x4 model=modelMatrix();pressure_force_instances_.reserve(std::min<std::size_t>(triangles,pressure_force_density_)*8);
		for(std::size_t triangle=0;triangle<triangles;triangle+=stride)
		{
			const float m=magnitude[triangle];if(!(m>0))continue;
			const auto i0=fabric_mesh_.indices[3*triangle],i1=fabric_mesh_.indices[3*triangle+1],i2=fabric_mesh_.indices[3*triangle+2];
			auto vertex=[&](std::uint32_t i){return QVector3D(fabric_mesh_.positions[3*i],fabric_mesh_.positions[3*i+1],fabric_mesh_.positions[3*i+2]);};
			QVector3D position=model.map((vertex(i0)+vertex(i1)+vertex(i2))/3.0f);
			QVector3D direction(triangle_pressure_forces_[3*triangle],triangle_pressure_forces_[3*triangle+1],triangle_pressure_forces_[3*triangle+2]);direction/=m;
			position+=direction*std::max(1e-5f,0.001f*arrow_len_);
			const float normalized=std::clamp(std::sqrt(m/focus),0.12f,1.0f);
			pressure_force_instances_.insert(pressure_force_instances_.end(),{position.x(),position.y(),position.z(),direction.x(),direction.y(),direction.z(),normalized,0.84f});
		}
		pressure_force_draw_count_=static_cast<int>(pressure_force_instances_.size()/8);
		if(gl_ready_&&pressure_force_draw_count_>0)
		{
			glBindBuffer(GL_ARRAY_BUFFER,pressure_force_inst_vbo_);
			glBufferData(GL_ARRAY_BUFFER,pressure_force_instances_.size()*sizeof(float),pressure_force_instances_.data(),GL_DYNAMIC_DRAW);
			glBindBuffer(GL_ARRAY_BUFFER,0);
		}
	}

	void SliceViewer::drawPressureForces(const QMatrix4x4& mvp)
	{
		if(!show_pressure_forces_||pressure_force_draw_count_<=0)return;
		arrow_prog_.bind();arrow_prog_.setUniformValue("uMVP",mvp);arrow_prog_.setUniformValue("uEye",camera_.eye());
		arrow_prog_.setUniformValue("uLen",0.45f*arrow_len_);arrow_prog_.setUniformValue("uSize",pressure_force_size_);
		arrow_prog_.setUniformValue("uColorMode",1);arrow_prog_.setUniformValue("uSolidColor",QVector3D(1.0f,0.70f,0.16f));
		glEnable(GL_BLEND);glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);glDepthMask(GL_FALSE);glDepthFunc(GL_LEQUAL);
		glBindVertexArray(pressure_force_vao_);glDrawArraysInstanced(GL_TRIANGLES,0,arrow_glyph_verts_,pressure_force_draw_count_);glBindVertexArray(0);
		glDepthFunc(GL_LESS);glDepthMask(GL_TRUE);glDisable(GL_BLEND);arrow_prog_.release();
	}

	void SliceViewer::updateArrows()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_ || !show_arrows_ || !paraglider_worker_ || !have_info_)
		{
			arrow_draw_count_ = 0;
			return;
		}
		float dt = 0.0f;
		if (!arrow_clock_started_) { arrow_clock_.start(); arrow_clock_started_ = true; }
		else dt = (float)(arrow_clock_.restart() / 1000.0);

		ArrowView view;
		view.three_d = arrow_mode_==0;
		view.static_grid = arrow_mode_==2;
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
		view.age_rate = arrow_speed_mult_;     // preserve travel distance before a particle ages out

		ensureFabricBvh();
		auto consume = [&](const FlowField& f) { FlowField ff = f; ff.fabric = fabric_bvh_.get(); arrows_.advance(dt, view, ff); };
		const bool got = paraglider_worker_->withFlowField(consume);
		if (!got) { arrow_draw_count_ = 0; return; }

		const std::vector<float>& inst = arrows_.instance_data();
		int n = arrows_.count();
		arrow_draw_count_ = n;
		if (n == 0) return;

		if (!arrows_logged_)
		{
			arrows_logged_ = true;
			const char* mode=arrow_mode_==0?"3D":arrow_mode_==1?"2D":"static-grid";
			std::fprintf(stderr, "[viewer] flow arrows live: %d particles, %s mode\n", n, mode);
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
		// The scalar slice never writes depth, so arrows remain visible over it while ordinary depth
		// testing still lets the zero-thickness STEP canopy occlude glyphs on its far side.
		glDepthMask(GL_FALSE);
		arrow_prog_.bind();
		arrow_prog_.setUniformValue("uMVP", mvp);
		arrow_prog_.setUniformValue("uEye", camera_.eye());
		arrow_prog_.setUniformValue("uLen", arrow_len_);
		arrow_prog_.setUniformValue("uSize", arrow_size_);
		// Over a visible colour slice the speed-coloured arrows blend into the field, so draw them
		// solid WHITE for contrast (reads over any colormap stop); when the slice is hidden, colour
		// them by speed (matches the legend) so they stay informative against the dark background.
		arrow_prog_.setUniformValue("uColorMode", show_slice_ ? 1 : 0);
		arrow_prog_.setUniformValue("uSolidColor", arrow_mode_==2?QVector3D(0.03f,0.04f,0.05f):QVector3D(1.0f,1.0f,1.0f));
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
		if (!gl_ready_ || !show_tracers_ || !paraglider_worker_ || !have_info_)
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
		// Integration-step cap. Ensure a fine display grid can still follow recirculating paths for a
		// couple of domain lengths; a larger user "Tracer length" remains authoritative.
		const float kTracerSpan = 2.5f; // enough arc reach for recirculating paths before the cap
		const int span_floor = (int)std::ceil(kTracerSpan * (float)info_.Lx / std::max(1e-6f, (float)info_.h));
		view.max_points = std::max(tracer_trail_, span_floor);
		// Rank the current field by path tortuosity and hide this fraction from the boring end.
		// Quantiles use the full slider range even when most paths are exactly straight.
		view.boring_hide_fraction = tracer_boring_;
		view.dt = dt;
		view.hold_seconds = 1.0f; // keep a tracer for 1 s after it was last interesting (anti-flicker)
		view.instant = tracer_boring_instant_; // dragging the slider ⇒ bypass the hold (live filter)

		ensureFabricBvh();
		auto consume = [&](const FlowField& f) { FlowField ff = f; ff.fabric = fabric_bvh_.get(); tracers_.advance(view, ff); };
		const bool got = paraglider_worker_->withFlowField(consume);
		if (!got) { tracer_strips_ = 0; return; }

		tracer_strips_ = tracers_.strips();
		tracer_firsts_ = tracers_.firsts();
		tracer_counts_ = tracers_.counts();
		if (tracer_strips_ == 0) return;

		if (!tracers_logged_)
		{
			tracers_logged_ = true;
			std::fprintf(stderr, "[viewer] flow tracers live: %d inlet streamlines, %s mode\n",
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

	void SliceViewer::drawGeometryIssuePolylines(const QMatrix4x4& mvp,
		const QVector4D& clip_plane)
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!show_geometry_issues_ || geometry_issue_vertices_.empty()) return;
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LEQUAL);
		glDepthMask(GL_FALSE);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		// A small toward-camera polygon offset makes a diagnostic curve on the source face readable.
		// Depth testing stays enabled, so a genuinely nearer skin still occludes a back-side finding.
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(-2.0f, -2.0f);
		if (clip_enabled_) glEnable(GL_CLIP_DISTANCE0);
		geometry_issue_prog_.bind();
		geometry_issue_prog_.setUniformValue("uVP", mvp);
		geometry_issue_prog_.setUniformValue("uModel", modelMatrix());
		geometry_issue_prog_.setUniformValue("uClipPlane", clip_plane);
		const float dpr = static_cast<float>(devicePixelRatioF());
		geometry_issue_prog_.setUniformValue("uViewport",
			QVector2D(static_cast<float>(width()) * dpr, static_cast<float>(height()) * dpr));
		glBindVertexArray(geometry_issue_vao_);
		for (std::size_t kind_index = 0; kind_index < kGeometryIssueKindCount; ++kind_index)
		{
			const auto& firsts = geometry_issue_firsts_[kind_index];
			if (firsts.empty()) continue;
			geometry_issue_prog_.setUniformValue("uHalfPx",
				0.5f * issue_width_pixels(static_cast<GeometryIssueKind>(kind_index)) * dpr);
			glMultiDrawArrays(GL_TRIANGLE_STRIP, firsts.data(), geometry_issue_counts_[kind_index].data(),
				static_cast<GLsizei>(firsts.size()));
		}
		glBindVertexArray(0);
		geometry_issue_prog_.release();
		glDisable(GL_CLIP_DISTANCE0);
		glPolygonOffset(0.0f, 0.0f);
		glDisable(GL_POLYGON_OFFSET_FILL);
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);
		glDepthFunc(GL_LESS);
	}

	void SliceViewer::drawGeometryErrorDiagnosticSegments(const QMatrix4x4& mvp)
	{
		PARACFD_ASSERT_GL_THREAD();
		if(!geometry_error_active_||geometry_error_firsts_.empty())return;
		// Focused errors are intentionally an always-readable overlay: the diagnostic canopy is
		// transparent and these opaque ribbons do not disappear behind either canopy skin.
		glDisable(GL_DEPTH_TEST);glDepthMask(GL_FALSE);glDisable(GL_BLEND);glDisable(GL_CLIP_DISTANCE0);
		geometry_issue_prog_.bind();QMatrix4x4 identity;identity.setToIdentity();
		geometry_issue_prog_.setUniformValue("uVP",mvp);
		geometry_issue_prog_.setUniformValue("uModel",identity);
		geometry_issue_prog_.setUniformValue("uClipPlane",QVector4D(0,0,0,1));
		const float dpr=static_cast<float>(devicePixelRatioF());
		geometry_issue_prog_.setUniformValue("uViewport",QVector2D(static_cast<float>(width())*dpr,static_cast<float>(height())*dpr));
		geometry_issue_prog_.setUniformValue("uHalfPx",3.5f*dpr);
		glBindVertexArray(geometry_error_vao_);
		glMultiDrawArrays(GL_TRIANGLE_STRIP,geometry_error_firsts_.data(),geometry_error_counts_.data(),
			static_cast<GLsizei>(geometry_error_firsts_.size()));
		glBindVertexArray(0);geometry_issue_prog_.release();
		glDepthMask(GL_TRUE);glEnable(GL_DEPTH_TEST);glDepthFunc(GL_LESS);
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
			p.drawText(titleRect,Qt::AlignCenter,units.isEmpty()?title:title+"  ["+units+"]");
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
		const QString title=QString::fromLatin1(field_short_name(field_));
		const QString units=QString::fromLatin1(field_units(field_));
		draw_bar(barX, barY, barW, barH, title, vmin_, vmax_, units,
			[](float t) { float r, g, b; scour_colormap(t, r, g, b); return QColor::fromRgbF(r, g, b); });

		// External-aero convention is deliberately fixed: prescribed freestream enters at
		// X-min and travels along +X. Keep this visible even when no arrows cross the camera.
		if (paraglider_worker_)
		{
			QFont f = p.font(); f.setPointSizeF(10.0); f.setBold(true); p.setFont(f);
			const bool stopped_with_error=!playback_active_&&!simulation_error_label_.isEmpty();
			const QColor state_colour=playback_active_?QColor(61,205,224):
				(stopped_with_error?QColor(255,105,90):(simulation_running_?QColor(75,220,120):QColor(255,185,65)));
			const QRect state_cue(std::max(238,(width()-410)/2),16,410,70);p.fillRect(state_cue,QColor(18,20,24,220));p.setPen(state_colour);p.drawRect(state_cue.adjusted(0,0,-1,-1));p.setBrush(state_colour);p.setPen(Qt::NoPen);p.drawEllipse(QRect(state_cue.x()+10,state_cue.y()+10,12,12));p.setBrush(Qt::NoBrush);p.setPen(state_colour);
			QString state_text;
			if(playback_active_)state_text="PLAYBACK";
			else if(stopped_with_error)state_text="STOPPED — ERROR";
			else if(simulation_running_)state_text="RUNNING";
			else if(simulation_auto_paused_)state_text=QString("PAUSED — %1").arg(simulation_pause_label_.isEmpty()?QString("AUTOMATIC"):simulation_pause_label_);
			else state_text="PAUSED";
			if(!simulation_case_label_.isEmpty())state_text+=QString(" — %1").arg(simulation_case_label_);p.drawText(state_cue.adjusted(30,2,-6,-44),Qt::AlignLeft|Qt::AlignVCenter,state_text);
			QFont state_detail=f;state_detail.setBold(false);state_detail.setPointSizeF(8.5);p.setFont(state_detail);p.setPen(QColor(225,228,234));QString detail;
			if(playback_active_)detail=QString("recorded frame %1 / %2 · simulation timeline")
				.arg(playback_frame_+1).arg(playback_frame_count_);
			else if(stopped_with_error)detail=QFontMetrics(p.font()).elidedText(simulation_error_label_,Qt::ElideRight,state_cue.width()-20);
			else if(!simulation_auto_pause_enabled_)detail="auto-pause disabled";
			else if(simulation_settling_ready_&&simulation_flow_throughs_<kAutoPauseMinimumFlowThroughs)
				detail=QString("score %1 · warm-up %2/%3 · max %4 flow").arg(simulation_settling_score_,0,'f',2).arg(simulation_flow_throughs_,0,'f',2).arg(kAutoPauseMinimumFlowThroughs,0,'f',2).arg(simulation_maximum_flow_throughs_,0,'g',3);
			else if(simulation_settling_ready_&&simulation_mean_force_ready_)
				detail=QString("settle %1 · mean %2/%3% · flow %4/%5").arg(simulation_settling_score_,0,'f',2).arg(100.0*simulation_mean_force_drift_,0,'g',3).arg(100.0*simulation_mean_force_tolerance_,0,'g',3).arg(simulation_flow_throughs_,0,'f',2).arg(simulation_maximum_flow_throughs_,0,'g',3);
			else if(simulation_settling_ready_)
				detail=QString("settle %1 · mean warming · flow %2/%3").arg(simulation_settling_score_,0,'f',2).arg(simulation_flow_throughs_,0,'f',2).arg(simulation_maximum_flow_throughs_,0,'g',3);
			else detail=QString("observing %1/%2 flow-throughs").arg(simulation_flow_throughs_,0,'f',2).arg(simulation_maximum_flow_throughs_,0,'g',3);
			p.drawText(state_cue.adjusted(10,24,-6,-24),Qt::AlignLeft|Qt::AlignVCenter,detail);
			const QString time_detail=playback_active_?
				QString("simulation t=%1 s  ·  source step %2").arg(playback_physical_time_,0,'f',3).arg(simulation_step_):
				QString("step %1  ·  simulation t=%2 s  ·  wall=%3 s").arg(simulation_step_).arg(simulation_physical_time_,0,'f',3).arg(simulation_wall_time_,0,'f',2);
			p.drawText(state_cue.adjusted(10,46,-6,-3),Qt::AlignLeft|Qt::AlignVCenter,time_detail);p.setFont(f);
			const QRect cue(18, 16, 210, 28);
			p.fillRect(cue, QColor(18, 20, 24, 190));
			p.setPen(QColor(255, 120, 90));
			p.drawText(cue.adjusted(9, 0, -5, 0), Qt::AlignLeft | Qt::AlignVCenter,
				QString::fromUtf8("FREESTREAM  +X  →"));
			if(show_slice_)
			{
				QFont detail=f;detail.setBold(false);detail.setPointSizeF(9.0);p.setFont(detail);
				const QRect resolutionCue(18,46,270,22);p.fillRect(resolutionCue,QColor(18,20,24,170));p.setPen(QColor(220,225,232));
				p.drawText(resolutionCue.adjusted(9,0,-5,0),Qt::AlignLeft|Qt::AlignVCenter,
					QString("SLICE  %1 x %2   finest h=%3 m").arg(slice_nu_).arg(slice_nv_).arg(info_.finest_h,0,'g',4));
			}
			if(thin_debug_)
			{
				QFont debugFont=f;debugFont.setPointSizeF(11.0);debugFont.setBold(true);p.setFont(debugFont);const QRect debugCue(18,74,340,30);p.fillRect(debugCue,QColor(115,25,20,215));p.setPen(QColor(255,235,220));p.drawText(debugCue.adjusted(9,0,-5,0),Qt::AlignLeft|Qt::AlignVCenter,QString("CROPPED-Y  %1 CELLS   STEP %2").arg(thin_debug_layers_).arg(simulation_step_));
			}
		}
	}

	void SliceViewer::drawGeometryIssueLegend(QPainter& painter)
	{
		std::vector<std::size_t> visible_kinds;
		std::size_t total = 0;
		if (show_geometry_issues_)
		for (std::size_t kind = 0; kind < kGeometryIssueKindCount; ++kind)
		{
			if (geometry_issue_polyline_counts_[kind] == 0) continue;
			visible_kinds.push_back(kind);
			total += geometry_issue_polyline_counts_[kind];
		}
		if (visible_kinds.empty() && geometry_quality_status_.isEmpty()) return;

		const bool has_gap = geometry_issue_polyline_counts_[issue_kind_index(GeometryIssueKind::UnexpectedGap)] > 0;
		const bool has_warning = has_gap
			|| geometry_issue_polyline_counts_[issue_kind_index(GeometryIssueKind::RepairedSeam)] > 0
			|| geometry_issue_polyline_counts_[issue_kind_index(GeometryIssueKind::UnsharedSeam)] > 0
			|| geometry_issue_polyline_counts_[issue_kind_index(GeometryIssueKind::IsolatedFabricComponent)] > 0
			|| geometry_issue_polyline_counts_[issue_kind_index(GeometryIssueKind::UnclassifiedOpening)] > 0;
		const bool informational_status = !geometry_quality_status_.isEmpty()
			&& !geometry_quality_warning_;
		const QColor border = (has_gap || geometry_quality_warning_) ? QColor(255, 65, 45)
			: (has_warning || informational_status) ? QColor(255, 174, 42) : QColor(30, 210, 235);
		QFont font = painter.font();
		font.setFamily(QStringLiteral("Consolas"));
		font.setPointSizeF(8.5);
		QFont title_font = font; title_font.setBold(true);
		const int row_height = 19;
		const int box_width = std::min(520,std::max(180,width()-32));
		const int content_width=box_width-22;
		const QFontMetrics metrics(font);
		const int status_height=geometry_quality_status_.isEmpty()?0:std::max(row_height,
			metrics.boundingRect(QRect(0,0,content_width,1000),Qt::TextWordWrap,
				geometry_quality_status_).height()+4);
		const int box_height = 30 + row_height * static_cast<int>(visible_kinds.size())
			+ status_height + 8;
		const QRect box(std::max(8, width() - box_width - 16),
			std::max(8, height() - box_height - 16), box_width, box_height);

		painter.save();
		painter.setRenderHint(QPainter::Antialiasing, true);
		painter.fillRect(box, QColor(13, 18, 25, 225));
		painter.setPen(QPen(border, 1.5));
		painter.drawRect(box.adjusted(0, 0, -1, -1));
		painter.setFont(title_font);
		painter.setPen(border);
		const QString title = total ? QString("MODEL PROBLEMS  %1").arg(total)
			: QStringLiteral("SOLID CHECK");
		painter.drawText(box.adjusted(10, 4, -8, -box.height() + 27), Qt::AlignVCenter,title);
		painter.setFont(font);

		int y = box.y() + 31;
		if(!geometry_quality_status_.isEmpty())
		{
			painter.setPen(border);painter.drawText(QRect(box.x()+11,y,content_width,status_height),
				Qt::AlignLeft|Qt::AlignVCenter|Qt::TextWordWrap,geometry_quality_status_);
			y+=status_height;
		}
		for (const std::size_t kind_index : visible_kinds)
		{
			const auto kind = static_cast<GeometryIssueKind>(kind_index);
			const QVector4D c = issue_colour(kind);
			const QColor colour = QColor::fromRgbF(c.x(), c.y(), c.z(), c.w());
			painter.setPen(QPen(colour, std::min(5.0f, issue_width_pixels(kind)), Qt::SolidLine,
				Qt::RoundCap, Qt::RoundJoin));
			painter.drawLine(box.x() + 11, y + row_height / 2, box.x() + 35, y + row_height / 2);
			QString text = QString("%1  ×%2").arg(issue_label(kind)).arg(geometry_issue_polyline_counts_[kind_index]);
			if (geometry_issue_max_gap_m_[kind_index])
			{
				const double millimetres = *geometry_issue_max_gap_m_[kind_index] * 1000.0;
				text += millimetres < 0.01
					? QString("  max %1 µm").arg(millimetres * 1000.0, 0, 'g', 3)
					: QString("  max %1 mm").arg(millimetres, 0, 'g', 3);
			}
			painter.setPen(QColor(235, 240, 246));
			painter.drawText(QRect(box.x() + 43, y, box.width() - 51, row_height),
				Qt::AlignLeft | Qt::AlignVCenter, text);
			y += row_height;
		}
		painter.restore();
	}

	void SliceViewer::drawGeometryErrorDiagnosticOverlay(QPainter& painter)
	{
		if(!geometry_error_active_)return;
		const int box_width=std::max(220,std::min(560,width()-36)),box_height=64;
		const QRect box(18,std::max(18,height()-box_height-18),box_width,box_height);
		painter.save();painter.setRenderHint(QPainter::Antialiasing,true);
		painter.fillRect(box,QColor(24,8,8,232));painter.setPen(QPen(QColor(255,42,28),2.0));
		painter.drawRect(box.adjusted(0,0,-1,-1));QFont font=painter.font();font.setFamily(QStringLiteral("Consolas"));
		font.setPointSizeF(10.0);font.setBold(true);painter.setFont(font);painter.setPen(QColor(255,72,54));
		painter.drawText(box.adjusted(12,5,-10,-box.height()+27),Qt::AlignLeft|Qt::AlignVCenter,
			QStringLiteral("CFD PREPROCESSING NEEDS ATTENTION"));
		font.setPointSizeF(8.5);font.setBold(false);painter.setFont(font);painter.setPen(QColor(255,232,228));
		painter.drawText(box.adjusted(12,29,-10,-6),Qt::AlignLeft|Qt::AlignVCenter|Qt::TextWordWrap,
			geometry_error_overlay_text_);painter.restore();
	}

	void SliceViewer::buildSliceGeometry()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_ || !have_info_) return;
		SliceParams sp = currentParams();
		const std::size_t vertices = static_cast<std::size_t>(sp.nu) * sp.nv;
		if (slice_buffer_nu_ != sp.nu || slice_buffer_nv_ != sp.nv)
		{
			std::vector<unsigned int> idx;
			idx.reserve(static_cast<std::size_t>(sp.nu - 1) * (sp.nv - 1) * 6);
			for (int b = 0; b < sp.nv - 1; ++b)
				for (int a = 0; a < sp.nu - 1; ++a)
				{
					const unsigned int v00 = static_cast<unsigned int>(b * sp.nu + a);
					const unsigned int v10 = v00 + 1, v01 = v00 + sp.nu, v11 = v01 + 1;
					idx.push_back(v00); idx.push_back(v10); idx.push_back(v11);
					idx.push_back(v00); idx.push_back(v11); idx.push_back(v01);
				}
			glBindVertexArray(slice_vao_);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, idx_ebo_);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx.size() * sizeof(unsigned int), idx.data(), GL_STATIC_DRAW);
			glBindBuffer(GL_ARRAY_BUFFER, color_vbo_);
			glBufferData(GL_ARRAY_BUFFER, vertices * 4 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
			glBindBuffer(GL_ARRAY_BUFFER, pos_vbo_);
			glBufferData(GL_ARRAY_BUFFER, vertices * 3 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
			glBindVertexArray(0);
			index_count_ = static_cast<int>(idx.size());
			slice_buffer_nu_ = sp.nu;
			slice_buffer_nv_ = sp.nv;
		}
		std::vector<float> pos(vertices * 3);
		for (int b = 0; b < sp.nv; ++b)
			for (int a = 0; a < sp.nu; ++a)
			{
				float x, y, z;
				slice_vertex_world(sp, a, b, x, y, z);
				size_t o = static_cast<std::size_t>(b * sp.nu + a) * 3;
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

	void SliceViewer::buildGeometryIssueBuffers()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_) return;
		constexpr int stride = 11 * sizeof(float);
		glBindVertexArray(geometry_issue_vao_);
		glBindBuffer(GL_ARRAY_BUFFER, geometry_issue_vbo_);
		glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, stride, (void*)0);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)(4 * sizeof(float)));
		glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, stride, (void*)(7 * sizeof(float)));
		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);
		glEnableVertexAttribArray(2);
		glBindVertexArray(0);
	}

	void SliceViewer::uploadGeometryIssuePolylines()
	{
		PARACFD_ASSERT_GL_THREAD();
		if (!gl_ready_ || !geometry_issue_upload_pending_) return;
		geometry_issue_upload_pending_ = false;
		glBindBuffer(GL_ARRAY_BUFFER, geometry_issue_vbo_);
		glBufferData(GL_ARRAY_BUFFER, geometry_issue_vertices_.size() * sizeof(float),
			geometry_issue_vertices_.empty() ? nullptr : geometry_issue_vertices_.data(), GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	void SliceViewer::buildGeometryErrorDiagnosticBuffers()
	{
		PARACFD_ASSERT_GL_THREAD();if(!gl_ready_)return;constexpr int stride=11*sizeof(float);
		glBindVertexArray(geometry_error_vao_);glBindBuffer(GL_ARRAY_BUFFER,geometry_error_vbo_);
		glBufferData(GL_ARRAY_BUFFER,0,nullptr,GL_DYNAMIC_DRAW);
		glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,stride,(void*)0);
		glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,stride,(void*)(4*sizeof(float)));
		glVertexAttribPointer(2,4,GL_FLOAT,GL_FALSE,stride,(void*)(7*sizeof(float)));
		glEnableVertexAttribArray(0);glEnableVertexAttribArray(1);glEnableVertexAttribArray(2);
		glBindVertexArray(0);
	}

	void SliceViewer::uploadGeometryErrorDiagnosticSegments()
	{
		PARACFD_ASSERT_GL_THREAD();if(!gl_ready_||!geometry_error_upload_pending_)return;
		geometry_error_upload_pending_=false;glBindBuffer(GL_ARRAY_BUFFER,geometry_error_vbo_);
		glBufferData(GL_ARRAY_BUFFER,geometry_error_vertices_.size()*sizeof(float),
			geometry_error_vertices_.empty()?nullptr:geometry_error_vertices_.data(),GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER,0);
	}

	void SliceViewer::buildIsoSurface(const ScalarVolume& volume)
	{
		iso_vertices_.clear();iso_negative_vertices_=iso_positive_vertices_=0;
		if(!volume.valid())return;
		const std::size_t cubes=static_cast<std::size_t>(volume.nx-1)*(volume.ny-1)*(volume.nz-1);
		// Marching tetrahedra expands every sampled cube into six tetrahedra. Keeping this
		// near 200k cubes gives an interactive preview while retaining the large-scale wake.
		const int stride=std::max(1,static_cast<int>(std::ceil(std::cbrt(static_cast<double>(cubes)/200000.0))));
		static constexpr int tetrahedra[6][4]={{0,5,1,6},{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6}};
		static constexpr int edges[6][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
		auto value=[&](int i,int j,int k){return volume.values[(static_cast<std::size_t>(k)*volume.ny+j)*volume.nx+i];};
		auto append_surface=[&](float level)
		{
			for(int k=0;k<volume.nz-1;k+=stride)for(int j=0;j<volume.ny-1;j+=stride)for(int i=0;i<volume.nx-1;i+=stride)
			{
				const int i1=std::min(i+stride,volume.nx-1),j1=std::min(j+stride,volume.ny-1),k1=std::min(k+stride,volume.nz-1);
				const int ci[8]={i,i1,i1,i,i,i1,i1,i},cj[8]={j,j,j1,j1,j,j,j1,j1},ck[8]={k,k,k,k,k1,k1,k1,k1};
				QVector3D position[8];float scalar[8];
				for(int corner=0;corner<8;++corner){position[corner]=QVector3D((ci[corner]+0.5f)*volume.h,(cj[corner]+0.5f)*volume.h,(ck[corner]+0.5f)*volume.h);scalar[corner]=value(ci[corner],cj[corner],ck[corner]);}
				for(const auto& tetra:tetrahedra)
				{
					QVector3D crossing[4];int count=0;
					for(const auto& edge:edges)
					{
						const int a=tetra[edge[0]],b=tetra[edge[1]];const bool sa=scalar[a]>=level,sb=scalar[b]>=level;
						if(sa==sb)continue;const float denominator=scalar[b]-scalar[a];
						const float t=std::abs(denominator)>1e-20f?std::clamp((level-scalar[a])/denominator,0.0f,1.0f):0.5f;
						if(count<4)crossing[count++]=position[a]+t*(position[b]-position[a]);
					}
					if(count<3)continue;
					QVector3D centre;for(int q=0;q<count;++q)centre+=crossing[q];centre/=static_cast<float>(count);
					QVector3D normal=QVector3D::crossProduct(crossing[1]-crossing[0],crossing[2]-crossing[0]);
					if(normal.lengthSquared()<1e-16f)continue;normal.normalize();
					QVector3D basis=(crossing[0]-centre).normalized(),side=QVector3D::crossProduct(normal,basis).normalized();
					std::array<int,4> order{0,1,2,3};
					std::sort(order.begin(),order.begin()+count,[&](int a,int b){const QVector3D pa=crossing[a]-centre,pb=crossing[b]-centre;return std::atan2(QVector3D::dotProduct(pa,side),QVector3D::dotProduct(pa,basis))<std::atan2(QVector3D::dotProduct(pb,side),QVector3D::dotProduct(pb,basis));});
					for(int triangle=1;triangle<count-1;++triangle)
					{
						const QVector3D p0=crossing[order[0]],p1=crossing[order[triangle]],p2=crossing[order[triangle+1]];
						QVector3D n=QVector3D::crossProduct(p1-p0,p2-p0);if(n.lengthSquared()<1e-16f)continue;n.normalize();
						for(const QVector3D& p:{p0,p1,p2})iso_vertices_.insert(iso_vertices_.end(),{p.x(),p.y(),p.z(),n.x(),n.y(),n.z()});
					}
				}
			}
		};
		if(field_uses_signed_iso_pair(volume.field))
		{
			iso_threshold_value_=(0.04f+0.96f*iso_level_fraction_)*volume.focus_abs;
			append_surface(-iso_threshold_value_);iso_negative_vertices_=static_cast<int>(iso_vertices_.size()/6);
			append_surface(iso_threshold_value_);iso_positive_vertices_=static_cast<int>(iso_vertices_.size()/6)-iso_negative_vertices_;
		}
		else
		{
			iso_threshold_value_=(0.04f+0.92f*iso_level_fraction_)*volume.focus_positive;
			append_surface(iso_threshold_value_);iso_positive_vertices_=static_cast<int>(iso_vertices_.size()/6);
		}
		iso_upload_pending_=true;iso_surface_dirty_=false;
	}

	void SliceViewer::uploadIsoSurface()
	{
		if(!iso_upload_pending_)return;iso_upload_pending_=false;
		glBindBuffer(GL_ARRAY_BUFFER,iso_vbo_);glBufferData(GL_ARRAY_BUFFER,iso_vertices_.size()*sizeof(float),iso_vertices_.empty()?nullptr:iso_vertices_.data(),GL_DYNAMIC_DRAW);glBindBuffer(GL_ARRAY_BUFFER,0);
	}

	void SliceViewer::drawIsoSurface(const QMatrix4x4& mvp,const QVector4D& clip_plane)
	{
		if(!show_iso_surface_||iso_positive_vertices_+iso_negative_vertices_<=0)return;
		mesh_prog_.bind();mesh_prog_.setUniformValue("uEye",camera_.eye());mesh_prog_.setUniformValue("uClipPlane",clip_plane);
		mesh_prog_.setUniformValue("uUseVertexColor",0);mesh_prog_.setUniformValue("uUseTriangleColor",0);mesh_prog_.setUniformValue("uMVP",mvp);mesh_prog_.setUniformValue("uModel",QMatrix4x4{});
		glEnable(GL_BLEND);glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);glDepthMask(GL_FALSE);if(clip_enabled_)glEnable(GL_CLIP_DISTANCE0);glBindVertexArray(iso_vao_);
		if(iso_negative_vertices_>0){mesh_prog_.setUniformValue("uBaseColor",QVector4D(0.10f,0.48f,1.0f,0.48f));glDrawArrays(GL_TRIANGLES,0,iso_negative_vertices_);}
		if(iso_positive_vertices_>0){mesh_prog_.setUniformValue("uBaseColor",QVector4D(field_==Field::QCriterion?0.22f:1.0f,field_==Field::QCriterion?0.88f:0.34f,field_==Field::QCriterion?0.82f:0.10f,0.50f));glDrawArrays(GL_TRIANGLES,iso_negative_vertices_,iso_positive_vertices_);}
		glBindVertexArray(0);glDisable(GL_CLIP_DISTANCE0);glDepthMask(GL_TRUE);glDisable(GL_BLEND);mesh_prog_.release();
	}

	void SliceViewer::uploadVolumeTexture(const ScalarVolume& volume)
	{
		if(!volume.valid())return;glActiveTexture(GL_TEXTURE0);glBindTexture(GL_TEXTURE_3D,volume_texture_);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);
		glPixelStorei(GL_UNPACK_ALIGNMENT,1);glTexImage3D(GL_TEXTURE_3D,0,GL_R32F,volume.nx,volume.ny,volume.nz,0,GL_RED,GL_FLOAT,volume.values.data());
		// QPainter's glyph atlas assumes the OpenGL default row alignment. Leaving this at one
		// shears newly uploaded font rows into the diagonal fragments seen in the probe/legend text.
		glPixelStorei(GL_UNPACK_ALIGNMENT,4);glBindTexture(GL_TEXTURE_3D,0);
		volume_texture_generation_=volume.generation;volume_texture_ready_=true;
	}

	void SliceViewer::updateScalarRepresentations()
	{
		if(!paraglider_worker_||!have_info_)return;
		if(show_iso_surface_)
		{
			// Rebuilding a marching-tetrahedra mesh every solver publication can stall the UI and
			// hold the snapshot mutex. Derived wake structures remain readable at this cadence.
			const bool update_due=iso_surface_dirty_||!iso_update_clock_.isValid()||iso_update_clock_.elapsed()>=750;
			if(update_due)
			{
				if(!iso_volume_)iso_volume_=std::make_unique<ScalarVolume>();const auto old_generation=iso_volume_->generation;const Field old_field=iso_volume_->field;
				if(paraglider_worker_->scalarVolume(field_,*iso_volume_)&&(iso_surface_dirty_||old_generation!=iso_volume_->generation||old_field!=iso_volume_->field))buildIsoSurface(*iso_volume_);
				iso_update_clock_.restart();
			}
		}
		if(show_volume_)
		{
			if(!cp_volume_)cp_volume_=std::make_unique<ScalarVolume>();
			if(paraglider_worker_->scalarVolume(Field::PressureCoefficient,*cp_volume_)&&(!volume_texture_ready_||volume_texture_generation_!=cp_volume_->generation))uploadVolumeTexture(*cp_volume_);
		}
		if(iso_upload_pending_)uploadIsoSurface();
	}

	void SliceViewer::drawVolume(const QMatrix4x4& mvp,const QVector4D& clip_plane)
	{
		if(!show_volume_||!volume_texture_ready_||!cp_volume_||!cp_volume_->valid())return;
		bool invertible=false;const QMatrix4x4 inverse=mvp.inverted(&invertible);if(!invertible)return;
		volume_prog_.bind();volume_prog_.setUniformValue("uInvVP",inverse);volume_prog_.setUniformValue("uBoxMax",QVector3D((float)info_.Lx,(float)info_.Ly,(float)info_.Lz));
		volume_prog_.setUniformValue("uClipPlane",clip_plane);volume_prog_.setUniformValue("uThreshold",cp_volume_->focus_abs*(0.01f+0.99f*volume_threshold_));
		volume_prog_.setUniformValue("uFocus",cp_volume_->focus_abs);volume_prog_.setUniformValue("uOpacity",volume_opacity_);volume_prog_.setUniformValue("uCellSize",cp_volume_->h);
		volume_prog_.setUniformValue("uSteps",96+static_cast<int>(volume_detail_*320.0f));volume_prog_.setUniformValue("uVolume",0);
		glActiveTexture(GL_TEXTURE0);glBindTexture(GL_TEXTURE_3D,volume_texture_);glEnable(GL_BLEND);glBlendFunc(GL_ONE,GL_ONE_MINUS_SRC_ALPHA);glDepthMask(GL_FALSE);glDisable(GL_DEPTH_TEST);
		glBindVertexArray(volume_vao_);glDrawArrays(GL_TRIANGLES,0,3);glBindVertexArray(0);glEnable(GL_DEPTH_TEST);glDepthMask(GL_TRUE);glDisable(GL_BLEND);glBindTexture(GL_TEXTURE_3D,0);volume_prog_.release();
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

		if (mesh_upload_pending_) uploadMesh();
		if (mesh_triangle_colour_upload_pending_) uploadTriangleSurfaceColours();
		if (paraglider_debug_upload_pending_) uploadParagliderDebugBoxes();
		if (geometry_issue_upload_pending_) uploadGeometryIssuePolylines();
		if (geometry_error_upload_pending_) uploadGeometryErrorDiagnosticSegments();
		updateScalarRepresentations();
		if(show_pressure_forces_&&pressure_force_dirty_)buildPressureForceInstances();

		if (paraglider_worker_ && have_info_ && show_slice_)
		{
			const std::uint64_t generation=paraglider_worker_->visualizationGeneration();
			if(slice_sample_dirty_||generation!=slice_sample_generation_)
			{
			// Sample this 2-D plane directly from the finest active AMR brick at each display
			// vertex. Arrows and tracers retain a cheaper coarse 3-D snapshot, but the scalar
			// plane now exposes the actual refinement around the wing.
			SliceParams sp = currentParams();
			FieldRange range;
			if (paraglider_worker_->sampleAmrSlice(sp, host_slice_values_, range))
			{
				if (auto_range_ && (range_ctr_++ % kRangeEvery == 0))
				{
					applyAutoRange(range);
				}
				if (show_slice_)
				{
					host_slice_colours_.resize(host_slice_values_.size());
					const float span = vmax_ > vmin_ ? vmax_ - vmin_ : 1.0f;
					for (std::size_t q = 0; q < host_slice_values_.size(); ++q)
					{
						float r, g, b;
						scour_colormap((host_slice_values_[q] - vmin_) / span, r, g, b);
						host_slice_colours_[q] = {r, g, b, 0.72f};
					}
					glBindBuffer(GL_ARRAY_BUFFER, color_vbo_);
					glBufferSubData(GL_ARRAY_BUFFER, 0,
						host_slice_colours_.size() * sizeof(float4), host_slice_colours_.data());
					glBindBuffer(GL_ARRAY_BUFFER, 0);
				}
				slice_sample_dirty_=false;slice_sample_generation_=generation;
			}
			}
		}

		updateArrows(); // advect the arrow tracers + upload instance data (reads the worker's host flow)
		updateTracers(); // advect the streaklines + upload line geometry (same host-flow snapshot)

		// QPainter::beginNativePainting() resets GL to DEFAULT state (depth test OFF, depth mask
		// undefined, clear colour black, blending on) — so the scene's depth buffer must be
		// re-established every frame here, otherwise geometry draws in submission order.
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);
		glDepthMask(GL_TRUE);
		glDisable(GL_BLEND);
		glDisable(GL_CULL_FACE);
		glClearColor(0.09f, 0.10f, 0.12f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		QMatrix4x4 mvp = camera_.viewProjection();

		// Clip plane (see inside hollow structures): the KEPT half-space is dot(pos,n)+d ≥ 0. Every scene
		// vertex shader writes gl_ClipDistance[0] = dot(pos, uClipPlane), but it only cuts geometry while
		// GL_CLIP_DISTANCE0 is enabled per STEP-surface draw and left off for
		// the flow slice + arrows + reference geometry, so the flow field inside the cavity stays visible.
		const QVector4D clipPlane = computeClipPlane();
		drawVolume(mvp,clipPlane);
		prog_.bind();
		prog_.setUniformValue("uMVP", mvp);
		prog_.setUniformValue("uClipPlane", clipPlane);
		glDisable(GL_CLIP_DISTANCE0);

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
			const bool error_diagnostic=geometry_error_active_;
			mesh_prog_.setUniformValue("uBaseColor", QVector4D(0.74f, 0.71f, 0.66f,error_diagnostic?0.22f:1.0f));
			mesh_prog_.setUniformValue("uUseVertexColor", 0); // STEP model: flat base colour
			mesh_prog_.setUniformValue("uUseTriangleColor", has_mesh_triangle_colours_ ? 1 : 0);
			if (has_mesh_triangle_colours_) glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, mesh_triangle_colour_ssbo_);
			if (clip_enabled_) glEnable(GL_CLIP_DISTANCE0); // visual section only; CFD remains two-sided fabric
			if(error_diagnostic){glEnable(GL_BLEND);glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);glDepthMask(GL_FALSE);}
			glBindVertexArray(mesh_vao_);
			const QMatrix4x4 model = modelMatrix();
			mesh_prog_.setUniformValue("uMVP", mvp * model);
			mesh_prog_.setUniformValue("uModel", model);
			glDrawElements(GL_TRIANGLES, mesh_index_count_, GL_UNSIGNED_INT, (void*)0);
			if (has_mesh_triangle_colours_) glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, 0);
			glBindVertexArray(0);
			glDisable(GL_CLIP_DISTANCE0);
			if(error_diagnostic){glDepthMask(GL_TRUE);glDisable(GL_BLEND);}
			mesh_prog_.release();
		}
		drawIsoSurface(mvp,clipPlane);

		// Draw the field slice after the opaque canopy so the translucent plane does not
		// write depth before the STEP surface. This keeps both the field and wing readable.
		if (show_slice_ && !show_grid_)
		{
			prog_.bind();
			prog_.setUniformValue("uMVP",mvp);
			prog_.setUniformValue("uClipPlane",clipPlane);
			prog_.setUniformValue("uFlat",0);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			glDepthMask(GL_FALSE);
			glBindVertexArray(slice_vao_);
			glDrawElements(GL_TRIANGLES,index_count_,GL_UNSIGNED_INT,(void*)0);
			glBindVertexArray(0);
			glDepthMask(GL_TRUE);
			glDisable(GL_BLEND);
			prog_.release();
		}

		// World-origin XYZ triad + camera-aligned corner gizmo (depth-tested world axes; gizmo on top).
		drawAxes(mvp);

		// Grid overlay: cell-boundary lines where the grid meets the slice plane (shows the graded mesh).
		drawGrid(mvp);

		// Clip-plane visualisation (translucent quad + outline showing where the cut is). Drawn with the
		// clip test OFF so the plane itself is not clipped; only shown while the feature is enabled.
		if (clip_enabled_&&!clip_follows_slice_) drawClipPlaneViz(mvp, clipPlane);
		// CAD diagnostics are drawn after the surface and slice. Their slight depth bias keeps a curve
		// on the fabric legible, while normal depth testing still hides back-side findings.
		drawGeometryIssuePolylines(mvp, clipPlane);
		drawGeometryErrorDiagnosticSegments(mvp);

		// Flow arrows are blended without writing depth. The scalar plane cannot occlude them, but the
		// opaque STEP surface already in the depth buffer can.
		drawArrows(mvp);
		drawPressureForces(mvp);

		// Grid-seeded streaklines use the same rule: the scalar plane does not occlude them, while the
		// STEP surface does. Drawn after arrows so both flow layers blend over the scene.
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
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
		glDisable(GL_CLIP_DISTANCE0); // never leave the clip test on for Qt's paint engine
		glDepthMask(GL_TRUE);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D,0);
		glBindTexture(GL_TEXTURE_3D,0);
		glUseProgram(0);
		glBindBuffer(GL_ARRAY_BUFFER,0);
		glBindVertexArray(0);
		glPixelStorei(GL_UNPACK_ALIGNMENT,4);

		painter.endNativePainting();
		drawLegendWith(painter);
		drawAxesLabels(painter); // X/Y/Z + metre tick labels (same QPainter-after-GL path as the legend)
		drawGeometryIssueLegend(painter);
		drawGeometryErrorDiagnosticOverlay(painter);
		drawProbe(painter);

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

	void SliceViewer::updateProbeAt(const QPointF& screen)
	{
		probe_valid_=false;probe_triangle_=-1;if(!have_info_)return;
		QVector3D origin,direction;mouseRay(screen,origin,direction);
		const float diagonal=std::sqrt((float)(info_.Lx*info_.Lx+info_.Ly*info_.Ly+info_.Lz*info_.Lz));
		const float ray_length=std::max(4.0f*diagonal,4.0f*camera_.distance());
		float best_distance=std::numeric_limits<float>::infinity();QVector3D hit_world;
		if(show_model_&&!fabric_mesh_.empty())
		{
			ensureFabricBvh();if(fabric_bvh_)
			{
				const paracfd::core::Vec3d ray_a{origin.x(),origin.y(),origin.z()};
				const paracfd::core::Vec3d ray_b{origin.x()+direction.x()*ray_length,origin.y()+direction.y()*ray_length,origin.z()+direction.z()*ray_length};
				double t_min=1e-8;
				for(int crossing=0;crossing<32;++crossing)
				{
					const auto hit=fabric_bvh_->intersect_segment(ray_a,ray_b,t_min,1.0);if(!hit.hit)break;
					const QVector3D candidate((float)hit.position.x,(float)hit.position.y,(float)hit.position.z);const QVector4D clip=computeClipPlane();
					if(QVector3D::dotProduct(candidate,QVector3D(clip.x(),clip.y(),clip.z()))+clip.w()>=0)
					{best_distance=static_cast<float>(hit.t)*ray_length;hit_world=candidate;probe_triangle_=static_cast<int>(hit.triangle_id);break;}
					t_min=std::min(1.0,hit.t+1e-7);
				}
			}
		}
		if(show_slice_)
		{
			const int axis=static_cast<int>(axis_);const float extent=axis==0?(float)info_.Lx:axis==1?(float)info_.Ly:(float)info_.Lz;
			if(std::abs(direction[axis])>1e-7f)
			{
				const float distance=(plane_frac_*extent-origin[axis])/direction[axis];const QVector3D candidate=origin+distance*direction;
				if(distance>=0&&distance<best_distance&&candidate.x()>=0&&candidate.x()<=(float)info_.Lx&&candidate.y()>=0&&candidate.y()<=(float)info_.Ly&&candidate.z()>=0&&candidate.z()<=(float)info_.Lz)
				{best_distance=distance;hit_world=candidate;probe_triangle_=-1;}
			}
		}
		if(!std::isfinite(best_distance))return;
		probe_world_=hit_world;QStringList lines;
		lines<<QString("x %1   y %2   z %3 m").arg(hit_world.x(),0,'f',4).arg(hit_world.y(),0,'f',4).arg(hit_world.z(),0,'f',4);
		FieldProbeSample sample;
		if(paraglider_worker_&&paraglider_worker_->samplePoint(hit_world.x(),hit_world.y(),hit_world.z(),sample))
		{
			const double speed=std::sqrt(sample.u*sample.u+sample.v*sample.v+sample.w*sample.w);
			lines<<QString("u [%1, %2, %3]  |u| %4 m/s").arg(sample.u,0,'g',5).arg(sample.v,0,'g',5).arg(sample.w,0,'g',5).arg(speed,0,'g',5);
			lines<<QString::fromUtf8("Δp %1 Pa   Cp %2   |∇p| %3 Pa/m").arg(sample.pressure_delta,0,'g',6).arg(sample.pressure_coefficient,0,'g',5).arg(sample.pressure_gradient,0,'g',5);
			lines<<QString::fromUtf8("|ω| %1 1/s   Q %2 1/s²").arg(sample.vorticity_magnitude,0,'g',5).arg(sample.q_criterion,0,'g',5);
		}
		if(probe_triangle_>=0&&probe_triangle_<static_cast<int>(surface_cp_plus_.size()))
		{
			const float plus=surface_cp_plus_[probe_triangle_],minus=surface_cp_minus_[probe_triangle_],delta=surface_delta_cp_[probe_triangle_];
			lines<<QString::fromUtf8("surface #%1   Cp+ %2   Cp− %3   ΔCp %4").arg(probe_triangle_).arg(plus,0,'g',5).arg(minus,0,'g',5).arg(delta,0,'g',5);
		}
		probe_text_=lines.join('\n');probe_valid_=true;update();
	}

	void SliceViewer::drawProbe(QPainter& painter)
	{
		if(!probe_valid_)return;QPointF screen;if(!projectPoint(camera_.viewProjection(),probe_world_,screen))return;
		painter.save();painter.setRenderHint(QPainter::Antialiasing,true);QPen marker(QColor(255,224,118));marker.setWidthF(2.0);painter.setPen(marker);painter.drawEllipse(screen,5,5);painter.drawLine(screen+QPointF(-10,0),screen+QPointF(10,0));painter.drawLine(screen+QPointF(0,-10),screen+QPointF(0,10));
		QFont font=painter.font();font.setFamily("Consolas");font.setPointSizeF(9.0);painter.setFont(font);QFontMetrics metrics(font);
		QRect text_rect=metrics.boundingRect(QRect(0,0,430,200),Qt::TextWordWrap,probe_text_).adjusted(-10,-8,10,8);
		int x=static_cast<int>(screen.x()+16),y=static_cast<int>(screen.y()+16);if(x+text_rect.width()>width()-8)x=static_cast<int>(screen.x()-text_rect.width()-16);if(y+text_rect.height()>height()-8)y=height()-text_rect.height()-8;x=std::max(8,x);y=std::max(8,y);text_rect.moveTopLeft(QPoint(x,y));
		painter.fillRect(text_rect,QColor(13,19,27,232));painter.setPen(QColor(255,224,118));painter.drawRoundedRect(text_rect,3,3);painter.setPen(QColor(235,242,248));painter.drawText(text_rect.adjusted(10,8,-10,-8),Qt::TextWordWrap,probe_text_);painter.restore();
	}

	void SliceViewer::beginSliceDrag(const QPointF& screen)
	{
		slice_dragging_=true;slice_drag_start_=screen;slice_drag_start_fraction_=plane_frac_;drag_btn_=Qt::NoButton;
	}

	void SliceViewer::dragSlice(const QPointF& screen)
	{
		const int axis=static_cast<int>(axis_);const float extent=axis==0?(float)info_.Lx:axis==1?(float)info_.Ly:(float)info_.Lz;
		QVector3D centre(0.5f*(float)info_.Lx,0.5f*(float)info_.Ly,0.5f*(float)info_.Lz);centre[axis]=slice_drag_start_fraction_*extent;
		QVector3D shifted=centre;shifted[axis]+=0.1f*extent;QPointF a,b;float fraction_delta=0.0f;
		if(projectPoint(camera_.viewProjection(),centre,a)&&projectPoint(camera_.viewProjection(),shifted,b))
		{
			const QPointF projected=b-a,mouse=screen-slice_drag_start_;const double length2=projected.x()*projected.x()+projected.y()*projected.y();
			if(length2>9.0)fraction_delta=static_cast<float>((mouse.x()*projected.x()+mouse.y()*projected.y())/length2*0.1);
			else fraction_delta=static_cast<float>(-(screen.y()-slice_drag_start_.y())/std::max(120,height()));
		}
		setPlaneFraction(slice_drag_start_fraction_+fraction_delta);
	}

	void SliceViewer::mousePressEvent(QMouseEvent* e)
	{
		drag_btn_ = e->button();
		last_mouse_ = e->position().toPoint();
		press_mouse_=last_mouse_;
		if(e->button()==Qt::LeftButton&&(e->modifiers()&Qt::AltModifier)&&show_slice_){beginSliceDrag(e->position());return;}
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
		if(slice_dragging_){dragSlice(posf);last_mouse_=p;return;}

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
		if(slice_dragging_){slice_dragging_=false;drag_btn_=Qt::NoButton;update();return;}
		if (gz_drag_op_ >= 0)
		{
			gz_drag_op_ = -1;
			gz_drag_axis_ = -1;
			fabric_bvh_dirty_ = true;
			pressure_force_dirty_ = true;
			emit modelPlacementChanged();
			update();
		}
		const bool click=e->button()==Qt::LeftButton&&drag_btn_==Qt::LeftButton&&(e->position().toPoint()-press_mouse_).manhattanLength()<5;
		drag_btn_ = Qt::NoButton;if(click){updateProbeAt(e->position());update();}
		Q_UNUSED(e);
	}

	void SliceViewer::wheelEvent(QWheelEvent* e)
	{
		float steps = e->angleDelta().y() / 120.0f;
		if((e->modifiers()&Qt::ShiftModifier)&&show_slice_)setPlaneFraction(plane_frac_+0.02f*steps);
		else camera_.zoom(std::pow(1.15f, steps));
		update();
	}
}
