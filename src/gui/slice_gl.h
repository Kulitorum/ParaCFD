// slice_gl.h — CUDA-GL interop wrapper for the G1 slice viewer. Registers a GL vertex
// buffer with CUDA once (cudaGraphicsGLRegisterBuffer), then each frame maps it, lets
// the slice kernel write RGBA vertex colours DIRECTLY into GL memory (zero copy), and
// unmaps it. The opaque handle keeps GL headers out of the Qt widget translation unit.
#pragma once

#include "gui/slice_field.h"

namespace windcfd::gui
{
	// Register GL buffer `vbo` (raw name) for CUDA write-discard access. Must be called
	// with the GL context current (main thread). Returns an opaque handle, or nullptr on
	// failure (see stderr). `vbo` must already be sized to nu*nv*sizeof(float4).
	void* slice_gl_register(unsigned int vbo);

	// Unregister a handle from slice_gl_register (GL context current, main thread).
	void slice_gl_unregister(void* handle);

	// Map the registered VBO, run the slice colour kernel reading the device MAC fields,
	// then unmap. `p` (pressure) may be null (⇒ Pressure renders as 0). Must run on the main
	// (GL) thread; the caller is responsible for serialising against the simulation worker's
	// step(). Returns false on a CUDA/interop error.
	bool slice_gl_fill(void* handle, const double* u, const double* v, const double* w, const double* p,
		const SliceParams& sp);

	// GPU-reduce the live [min,max] of `field` (+ max speed) over the device snapshot, fluid cells
	// only (`solid` may be null). Runs on the resource's own stream + a persistent 3-float scratch;
	// copies back ONLY the 3 scalars (never the whole field), so the auto colour-range stays off the
	// PCIe/critical path. `p` may be null (⇒ Pressure samples as 0). Main (GL) thread; serialise vs
	// the worker like slice_gl_fill. Returns false on error or when no fluid cell contributed
	// (out->valid is set accordingly).
	bool slice_gl_reduce(void* handle, const double* u, const double* v, const double* w, const double* p,
		const unsigned char* solid, windcfd::core::MacGrid g, Field field, FieldRange* out);
}
