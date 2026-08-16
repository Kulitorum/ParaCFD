// Inlet-seeded streamlines for the paraglider viewer.
//
// Companion to the arrow field (flow_particles.h): where the arrows are short instanced glyphs, the
// tracers are LONG lines. A grid of seeds on the INLET face each integrate a streamline through the
// current velocity field — stepping along the flow all the way across the domain until the line
// leaves an edge, meets a solid, stalls, or curls up in an eddy (a length cap bounds the spiral).
// The lines are re-integrated every frame, so they warp as the wake sheds, and each vertex is
// coloured by the LOCAL flow speed through the shared colormap (no animation, no highlight). Pure
// CPU state (Qt-free, no GL); the DRAW is camera-facing RIBBONS — a
// glMultiDrawArrays(GL_TRIANGLE_STRIP) over the packed vertex buffer (slice_viewer.cpp) — so the
// thickness is real geometry (glLineWidth is capped by most core-profile drivers).
//
// Two modes (RESEARCH-agnostic — pure visualisation):
//   3D : seed the inlet plane (y,z grid), integrate with the full {u,v,w}.
//   2D : seed the inlet edge of the current slice plane, integrate the in-plane velocity (the
//        streamlines live on the slice — the classic top-down wake view).
//
// Units: SI (m, m/s). Positions are world metres in [0,Lx]x[0,Ly]x[0,Lz]. The velocity sampling +
// solid test mirror FlowParticles (cell-centred trilinear reconstruction from the MAC faces).
#pragma once

#include "gui/flow_particles.h" // FlowField + MacGrid

#include <array>
#include <vector>

namespace paracfd::gui
{
	// Seed + integration settings. The inlet seed lattice is sized from `density` (count along the
	// larger inlet dimension) so the spacing is isotropic. speed_scale maps |vel| to the hot end of
	// the shared colormap so tracer colour reads on the SAME scale as the |u| slice/legend.
	struct TracerView
	{
		bool three_d = true;      // false ⇒ seed + integrate on the current slice plane
		int axis = 2;             // slice normal (0=X,1=Y,2=Z) — used only in 2D
		float plane_pos = 0.0f;   // world coord [m] of the plane along `axis` (2D)
		float speed_scale = 1.0f; // |vel| that maps to the hot end of the ramp [m/s]
		int density = 12;         // inlet seed count along the larger inlet dimension
		int max_points = 600;     // max streamline length in integration steps (bounds eddies)
		float step_ds = 0.05f;    // arc-length step [m] (caller sets to ~h)
		// Fraction [0,1] of the least-interesting current paths to hide. Interest is deterministic
		// arc/chord tortuosity; the empirical rank adapts the whole control range to each flow field.
		float boring_hide_fraction = 0.0f;
		// Temporal hold: a tracer keeps being drawn for this long [s] after it was last "interesting"
		// (outside the hidden quantile), so paths hovering at the threshold don't flicker.
		// `dt` is the
		// seconds elapsed since the previous advance() (drives the hold countdown).
		float dt = 0.0f;
		float hold_seconds = 1.0f;
		// Instant mode (set while the boring slider is being dragged): bypass the retention — a tracer
		// that stops qualifying drops immediately (hold forced to 0) so the filter tracks the drag
		// frame-by-frame. Interesting tracers still refresh their hold, so retention resumes cleanly on
		// release with no pop.
		bool instant = false;
	};

	class FlowTracers
	{
	public:
		// Integrate once per published CFD field and cache the paths/scores. Slider changes only rank
		// and repack the cached geometry; a null/empty field clears the buffers.
		void advance(const TracerView& view, const FlowField& f);
		// Force the seed lattice to rebuild + the temporal hold to reset (call on a domain/axis change).
		void reset() { cfg_valid_ = false; cache_valid_ = false; }

		// Packed ribbon geometry for glMultiDrawArrays(GL_TRIANGLE_STRIP, firsts, counts, strips).
		// Interleaved stride 11 floats per vertex: [pos.x pos.y pos.z side  tan.x tan.y tan.z  r g b a]
		// (`side` = ±1 selects the ribbon edge, expanded in the shader; `tan` is the unit line
		// tangent). Two vertices per streamline point ⇒ counts[i] = 2*npoints_i.
		const std::vector<float>& vertices() const { return vbo_; }
		const std::vector<int>& firsts() const { return firsts_; }
		const std::vector<int>& counts() const { return counts_; }
		int strips() const { return (int)counts_.size(); }
		int vertexCount() const { return (int)(vbo_.size() / kStride); }

		static constexpr int kStride = 11; // floats per vertex

	private:
		struct Pt { float x, y, z, spd; }; // world position [m] + speed [m/s]
		void sample(const FlowField& f, float x, float y, float z, double& uu, double& vv, double& ww) const;
		bool crosses_fabric(const FlowField& f, float ax, float ay, float az, float bx, float by, float bz) const;
		// (Re)build the persistent inlet seed lattice when the config changes; resets the holds.
		void build_seeds(const TracerView& view, const FlowField& f);
		// Integrate one streamline from the seed into line_; returns the point count.
		int integrate_line(const TracerView& view, const FlowField& f, float sx, float sy, float sz);
		// Build a ribbon strip from a cached line into the packed buffers.
		void emit_ribbon(const TracerView& view, const Pt* line, int n);

		std::vector<std::array<float, 3>> seeds_; // persistent inlet seed points [m]
		std::vector<float> hold_;                 // per-seed visibility remaining [s] (temporal hold)
		std::vector<Pt> line_;                    // scratch: the current streamline's points
		std::vector<Pt> cached_points_;           // flattened paths for the current flow generation
		std::vector<std::size_t> cached_offsets_; // one offset per seed plus a terminal offset
		std::vector<float> cached_score_;         // arc/chord tortuosity for each seed
		std::vector<int> cached_rank_;            // ascending-interest empirical rank; -1 if too short
		std::vector<float> vbo_;                  // interleaved ribbon vertices (see vertices())
		std::vector<int> firsts_, counts_;

		// Seed-config signature: rebuild seeds + reset the holds when any of these change.
		bool cfg_valid_ = false,cache_valid_ = false;
		std::uint64_t cached_generation_ = 0;
		bool s_three_d_ = true;
		int s_axis_ = 2, s_density_ = 0, s_nx_ = 0, s_ny_ = 0, s_nz_ = 0,s_max_points_=0;
		float s_plane_ = 0.0f,s_step_ds_=0.0f;
	};
}
