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
	// Boundary uses prefer the exact edge pcurve(s).  If no supplied edge has a
	// pcurve on `face`, the same trim-valid surface-candidate path used for an
	// interior use is applied.  The operation fails closed when a periodic chart
	// has no continuous branch or has more than one equally valid branch.
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
