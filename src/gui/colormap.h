// colormap.h — the ONE source of truth for the G1 viewer's perceptual speed→colour ramp.
//
// A 5-stop blue→cyan→green→yellow→red gradient over t in [0,1]. Shared verbatim by:
//   - the slice colour kernel (gui/slice_field.cu, __host__ __device__),
//   - the flow-arrow fragment shader (a GLSL MIRROR of these exact stops), and
//   - the on-screen legend (gui/slice_viewer.cpp, QPainter swatches),
// so the slice, the arrows and the colorbar all agree. Plain float lerps ⇒ the GPU kernel
// and its CPU reference stay bit-for-bit identical (the slice_field parity test).
//
// If you change a stop here, change the GLSL mirror `cmap()` in slice_viewer.cpp's arrow
// shader too (there is a comment marking it).
#pragma once

#ifdef __CUDACC__
#	define SCOUR_CM_HD __host__ __device__
#else
#	define SCOUR_CM_HD
#endif

namespace scour::gui
{
	// Perceptual 5-stop gradient blue->cyan->green->yellow->red. `t` is clamped to [0,1].
	SCOUR_CM_HD inline void scour_colormap(float t, float& r, float& g, float& b)
	{
		t = t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t);
		const float rs[5] = { 0.23f, 0.10f, 0.20f, 0.98f, 0.85f };
		const float gs[5] = { 0.30f, 0.65f, 0.78f, 0.90f, 0.18f };
		const float bs[5] = { 0.75f, 0.90f, 0.30f, 0.20f, 0.15f };
		float x = t * 4.0f; // [0,4]
		int i = (int)x;
		if (i > 3) i = 3; // segment [i,i+1]
		float f = x - (float)i;
		r = rs[i] + (rs[i + 1] - rs[i]) * f;
		g = gs[i] + (gs[i + 1] - gs[i]) * f;
		b = bs[i] + (bs[i + 1] - bs[i]) * f;
	}
}
