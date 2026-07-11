// slice_field.h — GUI G1 slice colour-fill kernel (pure CUDA, no Qt/GL).
//
// Samples a scalar field of the MAC state on an axis-aligned slice plane and maps it
// through a colour map into a dense float4 (RGBA) vertex-colour array. This is the
// numerical half of the CUDA-GL interop path: the SAME device kernel writes either a
// plain device buffer (unit test / headless) or a buffer mapped from a registered GL
// VBO (slice_gl.cu) — zero copy in the GL case.
//
// Every device path here has an identical __host__ __device__ implementation, so the
// GPU kernel and the CPU reference are bit-for-bit the same arithmetic; the GPU-vs-CPU
// gtest (tests/test_slice_field.cu) only has to catch launch/index/memory mistakes.
// Units: SI (m, s, m/s, Pa). Fields are double (matching the M1/M2 core).
#pragma once

#include "core/fluid/mac_grid.h"

#include <cuda_runtime.h> // float4, cudaStream_t

namespace windcfd::gui
{
	// Which scalar of the {u,v,w,p} state to colour.
	enum class Field : int
	{
		SpeedMag = 0, // |(u,v,w)| at the cell centre
		VelU = 1,     // x-velocity (cell centre)
		VelV = 2,     // y-velocity
		VelW = 3,     // z-velocity
		Pressure = 4  // cell pressure
	};

	// Plane normal axis.
	enum class Axis : int
	{
		X = 0, // plane x = plane_pos; in-plane axes (y,z)
		Y = 1, // plane y = plane_pos; in-plane axes (x,z)
		Z = 2  // plane z = plane_pos; in-plane axes (x,y)
	};

	// Description of the slice mesh + colour mapping. POD so it is passed by value into
	// the device kernel. The mesh is a regular (nu x nv) grid of vertices; vertex (a,b)
	// (a in [0,nu), b in [0,nv)) maps to a world point on the plane, is looked up in the
	// nearest MAC cell, reduced to the chosen scalar, and colour-mapped over [vmin,vmax].
	struct SliceParams
	{
		windcfd::core::MacGrid grid; // MAC geometry (nx,ny,nz,h)
		Axis axis = Axis::Y;       // plane normal
		float plane_pos = 0.0f;    // world coordinate [m] of the plane along `axis`
		int nu = 128;              // mesh resolution along in-plane axis 1
		int nv = 128;              // mesh resolution along in-plane axis 2
		Field field = Field::SpeedMag;
		float vmin = 0.0f; // colour-map lower bound (field units)
		float vmax = 1.0f; // colour-map upper bound
	};

	// World position of mesh vertex (a,b) on the slice plane. Shared by the colour path
	// and the host that builds the static GL position buffer, so geometry and sampling
	// agree exactly. Returns SI metres.
	WINDCFD_HD inline void slice_vertex_world(const SliceParams& sp, int a, int b, float& x, float& y, float& z)
	{
		const windcfd::core::MacGrid& g = sp.grid;
		float Lx = (float)(g.nx * g.h), Ly = (float)(g.ny * g.h), Lz = (float)(g.nz * g.h);
		float s = (sp.nu > 1) ? (float)a / (float)(sp.nu - 1) : 0.0f;
		float t = (sp.nv > 1) ? (float)b / (float)(sp.nv - 1) : 0.0f;
		if (sp.axis == Axis::X) { x = sp.plane_pos; y = s * Ly; z = t * Lz; }       // in-plane (y,z)
		else if (sp.axis == Axis::Y) { x = s * Lx; y = sp.plane_pos; z = t * Lz; } // in-plane (x,z)
		else { x = s * Lx; y = t * Ly; z = sp.plane_pos; }                         // Axis::Z: (x,y)
	}

	// Fill `out` (device, nu*nv float4 RGBA in [0,1]) from the device MAC fields. `p` may
	// be null (then Pressure renders as 0). Launches on `stream` (0 = default).
	void slice_fill_gpu(const double* u, const double* v, const double* w, const double* p,
		const SliceParams& sp, float4* out, cudaStream_t stream);

	// CPU reference: identical arithmetic on HOST copies of the fields. `out` is host,
	// nu*nv float4. Used by the GPU-vs-CPU parity test.
	void slice_fill_cpu(const double* u, const double* v, const double* w, const double* p,
		const SliceParams& sp, float4* out);

	// --- Auto-range reduction (GUI colour-scale follows the live data) -----------------------
	// Reduced live range of the displayed field, computed by a GPU reduction over the FULL 3-D
	// domain (fluid cells only). Drives the slice colormap, the legend and the arrow speed scale
	// so all three agree. Only 3 scalars ever cross PCIe (never the whole field).
	struct FieldRange
	{
		float field_min = 0.0f; // min of the displayed field's cell-centred scalar over fluid cells
		float field_max = 0.0f; // max of the displayed field
		float speed_max = 0.0f; // max |(u,v,w)| over fluid cells (arrow colour scale)
		bool valid = false;     // false if no fluid cell contributed (field_min > field_max)
	};

	// GPU block-reduction: fills `out3_dev` (device, 3 floats) = [field_min, field_max, speed_max]
	// over cells where `solid` is 0 (null ⇒ all fluid); non-finite cells are skipped so a NaN blow-up
	// never poisons the scale. Initialises out3_dev itself; launches on `stream`. `field` selects the
	// scalar (same cell-centred sampling as the slice fill), so the reduced range matches the display.
	// `p` may be null (⇒ Pressure samples as 0), matching the slice fill.
	void slice_reduce_gpu(const double* u, const double* v, const double* w, const double* p,
		const unsigned char* solid, windcfd::core::MacGrid g, Field field, float* out3_dev, cudaStream_t stream);

	// CPU reference (host fields → out3[3]) with identical arithmetic. Used by the GPU-vs-CPU test.
	void slice_reduce_cpu(const double* u, const double* v, const double* w, const double* p,
		const unsigned char* solid, windcfd::core::MacGrid g, Field field, float out3[3]);
}
