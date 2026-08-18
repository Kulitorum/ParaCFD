// geometry_quality.h — OCCT-free, model-local diagnostics produced by CAD preprocessing.
//
// Geometry inspection is static preprocessing.  These records deliberately retain enough
// source provenance and line geometry for a GUI or command-line client to explain degraded
// input without exposing OpenCascade types outside paracfd_geometry.  Coordinates are SI
// metres in the loaded STEP model's local frame; the renderer applies the same placement as
// the corresponding TriMesh.
#pragma once

#include <array>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace paracfd::core
{
	enum class GeometryIssueStage : std::uint8_t
	{
		cad_input = 0,
		normalization = 1,
		tessellation = 2,
		embedded_boundary = 3
	};

	// Keep the vocabulary broader than the current producer. STEP import currently emits exact
	// disconnected components and explicitly unsupported named surfaces; later seam/gap/EB
	// audits append to the same report without changing client APIs.
	enum class GeometryIssueKind : std::uint8_t
	{
		intentional_opening_boundary = 0,
		unclassified_free_boundary = 1,
		matched_unshared_seam = 2,
		near_miss_gap = 3,
		disconnected_fabric_component = 4,
		proper_surface_intersection = 5,
		coincident_surface_overlap = 6,
		invalid_self_intersection = 7,
		nonconforming_tessellation = 8,
		missing_boundary_provenance = 9,
		unresolved_embedded_boundary = 10,
		// A source representation explicitly names a surface role that ParaCFD knows is
		// unsupported for this simulation.  Unlike disconnected-component detection this
		// remains valid if an exporter stretches the surface until it touches another face.
		unsupported_named_surface = 11
	};

	enum class GeometryIssueSeverity : std::uint8_t
	{
		information = 0,
		warning_run = 1,
		blocker = 2
	};

	enum class GeometryIssueState : std::uint8_t
	{
		detected = 0,
		accepted_intentional = 1,
		repaired = 2,
		represented_as_opening = 3,
		unresolved = 4
	};

	enum class GeometryFlowEligibility : std::uint8_t
	{
		clean = 0,
		ready_with_warnings = 1,
		blocked = 2
	};

	struct GeometryIssuePolyline
	{
		std::uint32_t source_edge_id = ~std::uint32_t{0};
		std::vector<std::array<float, 3>> points;
		// Absent for a purely topological issue such as a disconnected component.
		std::optional<double> measured_gap_m;
	};

	struct GeometryIssue
	{
		std::uint64_t id = 0; // deterministic within one imported STEP topology
		GeometryIssueStage stage = GeometryIssueStage::cad_input;
		GeometryIssueKind kind = GeometryIssueKind::disconnected_fabric_component;
		GeometryIssueSeverity severity = GeometryIssueSeverity::warning_run;
		GeometryIssueState state = GeometryIssueState::detected;
		std::string summary;

		// Surface area represented by the issue's triangles, not a guessed leak area.
		double surface_area_m2 = 0.0;
		double boundary_length_m = 0.0;
		std::vector<std::uint32_t> source_face_ids;
		std::vector<std::uint32_t> source_triangle_ids;
		// Exact (trimmed but otherwise unmodified) STEP RepresentationItem names that caused
		// a semantic classification. Empty for purely topological findings.
		std::vector<std::string> source_representation_names;
		std::vector<GeometryIssuePolyline> boundary_polylines;
	};

	struct GeometryQualityReport
	{
		// Number of face components found using only exact shared TopoDS edges and exact
		// full edge-on-face contact certificates. Free boundaries/openings never split a face
		// component and proximity alone never joins two components.
		std::size_t exact_connected_component_count = 0;
		double largest_component_area_m2 = 0.0;
		std::vector<GeometryIssue> issues;

		GeometryFlowEligibility flow_eligibility() const
		{
			bool warning = false;
			for (const GeometryIssue& issue : issues)
			{
				if (issue.severity == GeometryIssueSeverity::blocker)
					return GeometryFlowEligibility::blocked;
				warning = warning || issue.severity == GeometryIssueSeverity::warning_run;
			}
			return warning ? GeometryFlowEligibility::ready_with_warnings
				: GeometryFlowEligibility::clean;
		}

		std::size_t warning_count() const
		{
			std::size_t count = 0;
			for (const GeometryIssue& issue : issues)
				count += issue.severity == GeometryIssueSeverity::warning_run;
			return count;
		}

		std::size_t blocker_count() const
		{
			std::size_t count = 0;
			for (const GeometryIssue& issue : issues)
				count += issue.severity == GeometryIssueSeverity::blocker;
			return count;
		}
	};

	// These findings are reversible source-selection suggestions, not geometry repairs.  Keep
	// the policy in the OCC-free core so GUI and command-line solvers cannot silently disagree.
	inline bool geometry_issue_is_default_exclusion(const GeometryIssue& issue)
	{
		return issue.kind == GeometryIssueKind::disconnected_fabric_component
			|| issue.kind == GeometryIssueKind::unsupported_named_surface;
	}

	inline std::vector<std::uint32_t> default_excluded_triangle_ids(
		const GeometryQualityReport& report)
	{
		std::vector<std::uint32_t> result;
		for (const GeometryIssue& issue : report.issues)
			if (geometry_issue_is_default_exclusion(issue))
				result.insert(result.end(), issue.source_triangle_ids.begin(),
					issue.source_triangle_ids.end());
		std::sort(result.begin(), result.end());
		result.erase(std::unique(result.begin(), result.end()), result.end());
		return result;
	}
}
