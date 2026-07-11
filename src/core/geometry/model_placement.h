// model_placement.h — the ONE definition of where a loaded STEP model sits in the domain.
//
// The display (slice_viewer) and the voxelizer MUST agree on the model's placement or the
// solid mask lands somewhere other than where the mesh is drawn. This header is that single
// shared source of truth: an AFFINE placement world(v) = M·v + t, where M is a 3×3 linear
// part (rotation·scale) and t a translation (metres). The default M = identity reduces it to
// the increment-1 translation-only "centre on x/y, drop onto the bed" transform, byte-for-byte.
// The interactive placement gizmo (scour-gui) fills M with a rotation·scale so a loaded model
// can be moved/rotated/scaled and then voxelized exactly where the user placed it.
//
// OCC-free + Qt-free (pure inline maths on a TriMesh bbox), so it lives in libscour's
// include path and is included by both the core voxelizer and the GUI viewer.
// Units: SI METRES.
#pragma once

#include "core/geometry/tri_mesh.h"

namespace scour::core
{
	// Affine placement applied to every mesh coordinate: world(v) = M·v + t. M is row-major
	// (m[0..2] = row 0, m[3..5] = row 1, m[6..8] = row 2); the default is the identity so a
	// placement with only tx/ty/tz set behaves exactly like the old translation-only struct.
	struct ModelPlacement
	{
		double tx = 0.0, ty = 0.0, tz = 0.0;
		double m[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 }; // linear part (rotation·scale), row-major

		// world = M·v + t.
		void apply(double x, double y, double z, double& ox, double& oy, double& oz) const
		{
			ox = m[0] * x + m[1] * y + m[2] * z + tx;
			oy = m[3] * x + m[4] * y + m[5] * z + ty;
			oz = m[6] * x + m[7] * y + m[8] * z + tz;
		}

		// Determinant of the linear part = the factor by which M scales volume (1 for a pure
		// rotation/translation). Used to place-correct the mesh volume so the voxelizer's
		// mesh-vs-voxel volume gate stays meaningful when the model is scaled.
		double linear_det() const
		{
			return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
		}
	};

	// Centre the model in x/y over a domain of extent (Lx,Ly) and drop it onto the bed so its
	// minimum-z vertex rests on z = 0. Identity linear part ⇒ identical to the slice viewer's
	// default display transform (translation only).
	inline ModelPlacement place_model_on_bed(const TriMesh& mesh, double Lx, double Ly)
	{
		ModelPlacement p;
		if (mesh.empty()) return p;
		const double cx = 0.5 * ((double)mesh.bbox_min[0] + (double)mesh.bbox_max[0]);
		const double cy = 0.5 * ((double)mesh.bbox_min[1] + (double)mesh.bbox_max[1]);
		p.tx = 0.5 * Lx - cx;
		p.ty = 0.5 * Ly - cy;
		p.tz = -(double)mesh.bbox_min[2];
		return p;
	}
}
