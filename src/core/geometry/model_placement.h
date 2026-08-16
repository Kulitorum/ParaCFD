// model_placement.h — the ONE definition of where a loaded STEP model sits in the domain.
//
// The display, BVH, AMR, and embedded-boundary preprocessing must agree on the STEP placement.
// This header is their single shared source of truth: world(v) = M·v + t, where M is a 3×3
// linear part and t a translation in metres. The interactive viewer fills this transform so
// the fabric is preprocessed exactly where it is drawn.
//
// OCC-free + Qt-free (pure inline maths on a TriMesh bbox).
// Units: SI METRES.
#pragma once

#include "core/geometry/tri_mesh.h"

#include <cmath>

namespace paracfd::core
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

		// Determinant of the linear part (1 for a proper rotation). Its sign is needed
		// when transforming face normals through a reflected placement.
		double linear_det() const
		{
			return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
		}
	};

	// Face-only horizontal bbox inference deliberately determines axes, not polarity.
	// ParaCFD keeps Z as up, treats the longer of X/Y as span, and maps the assumed
	// negative chord direction into its fixed +X freestream. `yaw_degrees` is therefore
	// either +90 degrees for span-X/chord-Y data or 180 degrees for span-Y/chord-X data.
	// The caller must still expose a leading/trailing-edge flip.
	struct HorizontalWingAxes
	{
		bool valid = false;
		int span_axis = -1;
		int chord_axis = -1;
		double aspect_ratio = 1.0;
		double yaw_degrees = 0.0;
	};

	inline HorizontalWingAxes infer_horizontal_wing_axes(const TriMesh& mesh, double minimum_ratio = 1.25)
	{
		HorizontalWingAxes result;
		const double dx = static_cast<double>(mesh.bbox_max[0]) - mesh.bbox_min[0];
		const double dy = static_cast<double>(mesh.bbox_max[1]) - mesh.bbox_min[1];
		const double shorter = dx < dy ? dx : dy, longer = dx > dy ? dx : dy;
		if (!(shorter > 0.0) || !(longer > 0.0) || !(minimum_ratio > 1.0)) return result;
		result.aspect_ratio = longer / shorter;
		if (result.aspect_ratio < minimum_ratio) return result;
		result.valid = true;
		result.span_axis = dx > dy ? 0 : 1;
		result.chord_axis = 1 - result.span_axis;
		result.yaw_degrees = result.span_axis == 0 ? 90.0 : 180.0;
		return result;
	}

	// Return a copy of `mesh` with the placement applied to every vertex position and the bbox
	// recomputed. Indices and CAD provenance are preserved.
	inline TriMesh placed_mesh(const TriMesh& mesh, const ModelPlacement& p)
	{
		TriMesh out = mesh;
		const std::size_t nv = mesh.vertex_count();
		double bmin[3] = { 1e300, 1e300, 1e300 }, bmax[3] = { -1e300, -1e300, -1e300 };
		for (std::size_t i = 0; i < nv; ++i)
		{
			double ox, oy, oz;
			p.apply((double)mesh.positions[3 * i], (double)mesh.positions[3 * i + 1], (double)mesh.positions[3 * i + 2], ox, oy, oz);
			out.positions[3 * i] = (float)ox;
			out.positions[3 * i + 1] = (float)oy;
			out.positions[3 * i + 2] = (float)oz;
			const double v[3] = { ox, oy, oz };
			for (int c = 0; c < 3; ++c) { if (v[c] < bmin[c]) bmin[c] = v[c]; if (v[c] > bmax[c]) bmax[c] = v[c]; }
		}
		// Normals transform by inverse-transpose(M). The old rotation/uniform-scale path happened
		// to remain visually plausible when normals were copied unchanged, but EB plus/minus sides
		// are geometric and must agree with the displayed surface after any affine placement.
		if (mesh.normals.size() == mesh.positions.size())
		{
			const double det = p.linear_det();
			if (std::abs(det) > 1e-30)
			{
				// Triangle indices are not reversed under a reflecting placement, so their
				// geometric cross product transforms as det(M) M^-T n. Preserve that local
				// winding convention rather than silently flipping the displayed side.
				const double orientation = det < 0.0 ? -1.0 : 1.0;
				const double invt[9] = {
					(p.m[4] * p.m[8] - p.m[5] * p.m[7]) / det,
					(p.m[5] * p.m[6] - p.m[3] * p.m[8]) / det,
					(p.m[3] * p.m[7] - p.m[4] * p.m[6]) / det,
					(p.m[2] * p.m[7] - p.m[1] * p.m[8]) / det,
					(p.m[0] * p.m[8] - p.m[2] * p.m[6]) / det,
					(p.m[1] * p.m[6] - p.m[0] * p.m[7]) / det,
					(p.m[1] * p.m[5] - p.m[2] * p.m[4]) / det,
					(p.m[2] * p.m[3] - p.m[0] * p.m[5]) / det,
					(p.m[0] * p.m[4] - p.m[1] * p.m[3]) / det };
				for (std::size_t i = 0; i < nv; ++i)
				{
					const double x = mesh.normals[3 * i], y = mesh.normals[3 * i + 1], z = mesh.normals[3 * i + 2];
					double nx = orientation * (invt[0] * x + invt[1] * y + invt[2] * z);
					double ny = orientation * (invt[3] * x + invt[4] * y + invt[5] * z);
					double nz = orientation * (invt[6] * x + invt[7] * y + invt[8] * z);
					const double len = std::sqrt(nx * nx + ny * ny + nz * nz);
					if (len > 1e-30) { nx /= len; ny /= len; nz /= len; }
					out.normals[3 * i] = (float)nx; out.normals[3 * i + 1] = (float)ny; out.normals[3 * i + 2] = (float)nz;
				}
			}
		}
		if (nv)
			for (int c = 0; c < 3; ++c) { out.bbox_min[c] = (float)bmin[c]; out.bbox_max[c] = (float)bmax[c]; }
		return out;
	}
}
