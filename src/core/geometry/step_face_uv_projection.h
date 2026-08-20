// Internal OCCT-facing helper for constructing one deterministic, trim-valid
// parameter-space chart for a STEP contact sample chain.  Unlike step_import.h,
// this header intentionally exposes OpenCascade types and is only consumed by
// the isolated geometry library and its focused OCCT regression probe.
#pragma once

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>

#include <array>
#include <string>
#include <vector>

namespace paracfd::core::detail
{
	// Boundary uses prefer the exact oriented face-edge pcurve(s). For this path,
	// Each entry in `sample_geometric_tolerances_mm` is the native source-chain
	// side of that sample's exact CAD contact certificate. The implementation
	// combines only that entry with the actual target edge/face native tolerances
	// and independently verifies the target 3D curve and pcurve-on-surface round
	// trips. If no supplied edge has a pcurve on `face`, the same trim-valid
	// surface-candidate path used for an interior use is applied. The operation
	// fails closed when a periodic chart has no continuous branch or has more than
	// one equally valid branch.
	bool project_trim_valid_face_uv(
		const std::vector<gp_Pnt>& world_points,
		const TopoDS_Face& face,
		const std::vector<TopoDS_Edge>& exact_boundary_edges,
		bool boundary_use,
		const std::vector<double>& sample_geometric_tolerances_mm,
		std::vector<std::array<double, 2>>& sample_uv,
		std::string& error,
		bool* used_exact_boundary_pcurve = nullptr);

	// Backwards-compatible convenience wrapper for chains whose samples all have
	// the same source-side geometric certificate.
	bool project_trim_valid_face_uv(
		const std::vector<gp_Pnt>& world_points,
		const TopoDS_Face& face,
		const std::vector<TopoDS_Edge>& exact_boundary_edges,
		bool boundary_use,
		double geometric_tolerance_mm,
		std::vector<std::array<double, 2>>& sample_uv,
		std::string& error,
		bool* used_exact_boundary_pcurve = nullptr);
}
