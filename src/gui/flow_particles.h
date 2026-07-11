// flow_particles.h — variomigo-style animated particle-arrow field for the G1 viewer.
//
// A CPU particle system (Qt-free, no GL): N massless tracers seeded at random fluid points
// and advected each frame along the bilinearly/trilinearly sampled local MAC velocity. A
// particle respawns when it ages out, leaves the domain, or enters a solid cell — so arrows
// stream along the flow and die into obstacles / the bed. Ported from variomigo's Dart
// `FlowParticles` (app/lib/features/ridge_airflow/ridge_flow_painters.dart); here the DRAW is
// GL 4.3 instancing (slice_viewer.cpp), and this class only owns the CPU state + per-frame
// integration, emitting an interleaved instance array for the arrow shader.
//
// Two modes (RESEARCH-agnostic — pure visualisation):
//   3D : seed anywhere in the volume, advect with the full {u,v,w}.
//   2D : seed + advect ON the current slice plane using only the in-plane velocity (arrows
//        live on the visualisation plane; the out-of-plane component is projected out).
//
// Units: SI (m, m/s). Positions are world metres in [0,Lx]x[0,Ly]x[0,Lz].
#pragma once

#include "core/fluid/mac_grid.h"

#include <cstdint>
#include <random>
#include <vector>

namespace scour::gui
{
	// A borrowed, read-only view of the live host velocity + solid snapshot (owned by the
	// SimWorker, valid only under its flow lock). Cell-centred sampling reconstructs velocity
	// from the MAC faces; `solid` (nx*ny*nz) marks obstacle/bed cells (treated as zero flow).
	struct FlowField
	{
		const double* u = nullptr;
		const double* v = nullptr;
		const double* w = nullptr;
		const unsigned char* solid = nullptr;
		scour::core::MacGrid grid;
		// Current flow direction (ChannelBC::flow_sign): +1 ⇒ inlet at x-min (default), -1 ⇒ inlet at
		// x-max (reversed tide). The tracers seed from whichever face is the inlet.
		int flow_sign = 1;
	};

	// What plane / mode the arrows are advected in, plus the colour scale (so arrow speed maps
	// through the shared colormap over the SAME range as the slice's |u| legend).
	struct ArrowView
	{
		bool three_d = true;     // false ⇒ constrain to the slice plane (2D)
		int axis = 2;            // slice normal (Axis: 0=X,1=Y,2=Z) — used only in 2D
		float plane_pos = 0.0f;  // world coord [m] of the plane along `axis` (2D)
		float speed_scale = 1.0f; // |vel| that maps to the hot end of the ramp [m/s]
		float gain = 1.0f;       // advection speed multiplier (1 = real time)
	};

	class FlowParticles
	{
	public:
		// (Re)allocate to `n` particles and respawn them all. Applied live by the density slider.
		void set_count(int n);
		int count() const { return count_; }

		// Respawn every particle (call on a plane/axis/mode change so none linger off-plane).
		void reset() { needs_reset_ = true; }

		// Advance all particles by `dt` seconds through `f` under `view`. Fills instance_data().
		// Seeding needs the field (to reject solid spawn points); a null/empty field is a no-op.
		void advance(float dt, const ArrowView& view, const FlowField& f);

		// Interleaved per-particle instance data for the arrow shader, stride 8 floats:
		//   [pos.x pos.y pos.z  dir.x dir.y dir.z  speedNorm  alpha]
		// Only the first count()*8 entries are valid. dir is a unit vector (0 if stalled).
		const std::vector<float>& instance_data() const { return inst_; }

	private:
		void spawn(int k, const ArrowView& view, const FlowField& f, bool initial);
		// Cell-centred trilinear velocity at world (x,y,z); solids contribute zero. Returns the
		// three components; caller projects for 2D.
		void sample(const FlowField& f, float x, float y, float z, double& uu, double& vv, double& ww) const;
		bool is_solid(const FlowField& f, float x, float y, float z) const;

		int count_ = 0;
		bool needs_reset_ = false;

		std::vector<float> px_, py_, pz_;   // world position [m]
		std::vector<float> age_, life_;     // age + lifetime [s]
		std::vector<float> dirx_, diry_, dirz_; // unit velocity direction (world)
		std::vector<float> spd_;            // normalised speed [0,1] (for colour + length)
		std::vector<float> alpha_;          // fade-in/out
		std::vector<float> inst_;           // interleaved output (stride 8)

		std::mt19937 rng_{ 0xC0B0Du };
	};
}
