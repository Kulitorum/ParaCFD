#include "core/geometry/step_import.h"
#include "core/geometry/triangle_bvh.h"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
#include <string>

using namespace paracfd::core;

namespace
{
	const char* reason_name(FabricEdgeUnknownReason reason)
	{
		switch (reason)
		{
		case FabricEdgeUnknownReason::none: return "none";
		case FabricEdgeUnknownReason::invalid_triangle: return "invalid triangle";
		case FabricEdgeUnknownReason::source_provenance: return "source provenance";
		case FabricEdgeUnknownReason::indexed_fan_mismatch: return "indexed fan mismatch";
		case FabricEdgeUnknownReason::proximity_ambiguity: return "proximity ambiguity";
		case FabricEdgeUnknownReason::mixed_open_attachment: return "mixed open/attachment";
		case FabricEdgeUnknownReason::cad_fan_mismatch: return "CAD fan mismatch";
		case FabricEdgeUnknownReason::no_geometric_evidence: return "no geometric evidence";
		case FabricEdgeUnknownReason::nonreciprocal_attachment: return "nonreciprocal attachment";
		}
		return "invalid reason";
	}
}

int main(int argc, char** argv)
{
	const std::string path = argc > 1 ? argv[1] : "Test-Data/PlanBParakite.step";
	const double deflection_mm = argc > 2 ? std::strtod(argv[2], nullptr) : 2.0;
	std::string error;
	const StepGeometry imported = load_step_geometry(path, deflection_mm, &error);
	const TriMesh& mesh = imported.mesh;
	if (mesh.empty())
	{
		std::fprintf(stderr, "STEP import failed: %s\n", error.c_str());
		return 1;
	}

	const TriangleBvh bvh(mesh);
	std::array<std::uint64_t, 3> kinds{};
	std::array<std::uint64_t, 5> roles{};
	std::array<std::uint64_t, 9> reasons{};
	std::map<std::uint32_t, std::uint64_t> unknown_by_cad_edge;
	std::map<std::uint32_t, std::uint64_t> unknown_by_source_face;
	std::map<std::uint32_t, std::uint64_t> unknown_by_declared_incidence;
	std::uint64_t unknown_without_cad_edge = 0;
	std::uint64_t unknown_periodic = 0;
	double minimum_unknown_cad_tolerance = std::numeric_limits<double>::infinity();
	double maximum_unknown_cad_tolerance = 0.0;
	std::uint64_t unknown_cad_tolerance_covers_clearance = 0;
	unsigned printed_unknown_examples = 0;

	for (std::size_t triangle = 0; triangle < mesh.triangle_count(); ++triangle)
		for (unsigned edge = 0; edge < 3; ++edge)
		{
			const FabricEdgeCertificate certificate =
				bvh.edge_certificate(static_cast<std::uint32_t>(triangle), static_cast<int>(edge));
			++kinds[static_cast<std::size_t>(certificate.kind)];
			++roles[static_cast<std::size_t>(certificate.role)];
			if (certificate.kind != FabricEdgeKind::unknown) continue;
			++reasons[static_cast<std::size_t>(certificate.unknown_reason)];
			const std::uint32_t cad_edge = mesh.cad_edge_id(triangle, edge);
			if (cad_edge == TriMesh::kNoCadEdgeId) ++unknown_without_cad_edge;
			else
			{
				++unknown_by_cad_edge[cad_edge];
				const double tolerance = mesh.cad_edge_tolerance(triangle, edge);
				minimum_unknown_cad_tolerance = std::min(minimum_unknown_cad_tolerance, tolerance);
				maximum_unknown_cad_tolerance = std::max(maximum_unknown_cad_tolerance, tolerance);
				if (tolerance >= bvh.edge_clearance_tolerance())
					++unknown_cad_tolerance_covers_clearance;
			}
			if (mesh.has_face_provenance()) ++unknown_by_source_face[mesh.source_face_ids[triangle]];
			++unknown_by_declared_incidence[mesh.cad_edge_incident_face_count(triangle, edge)];
			if (mesh.cad_edge_is_periodic_seam(triangle, edge)) ++unknown_periodic;
			if (printed_unknown_examples < 24 && certificate.unknown_reason !=
				FabricEdgeUnknownReason::invalid_triangle)
			{
				const auto a = mesh.vertex_position_double(mesh.indices[3 * triangle + edge]);
				const auto b = mesh.vertex_position_double(mesh.indices[3 * triangle + (edge + 1) % 3]);
				std::printf("example t=%zu e=%u face=%u CAD=%u incidence=%u reason=%s "
					"a=[%.9g %.9g %.9g] b=[%.9g %.9g %.9g]\n", triangle, edge,
					mesh.has_face_provenance() ? mesh.source_face_ids[triangle] : 0u, cad_edge,
					mesh.cad_edge_incident_face_count(triangle, edge),
					reason_name(certificate.unknown_reason), a[0], a[1], a[2], b[0], b[1], b[2]);
				++printed_unknown_examples;
			}
		}

	std::printf("STEP %s\n", path.c_str());
	std::printf("triangles %zu; BVH triangles %zu; half-edges %zu\n", mesh.triangle_count(),
		bvh.triangle_count(), 3 * mesh.triangle_count());
	std::printf("kind attached/free/unknown: %llu / %llu / %llu\n",
		static_cast<unsigned long long>(kinds[0]), static_cast<unsigned long long>(kinds[1]),
		static_cast<unsigned long long>(kinds[2]));
	std::printf("role tessellation/seam/junction/unknown: %llu / %llu / %llu / %llu\n",
		static_cast<unsigned long long>(roles[1]), static_cast<unsigned long long>(roles[2]),
		static_cast<unsigned long long>(roles[3]), static_cast<unsigned long long>(roles[4]));
	for (std::size_t reason = 1; reason < reasons.size(); ++reason)
		if (reasons[reason] != 0)
			std::printf("unknown reason %-26s %llu\n",
				reason_name(static_cast<FabricEdgeUnknownReason>(reason)),
				static_cast<unsigned long long>(reasons[reason]));
	std::printf("unknown with/without CAD edge: %llu / %llu; periodic %llu\n",
		static_cast<unsigned long long>(kinds[2] - unknown_without_cad_edge),
		static_cast<unsigned long long>(unknown_without_cad_edge),
		static_cast<unsigned long long>(unknown_periodic));
	std::printf("unknown CAD tolerance min/max %.9g / %.9g m; covers BVH clearance %llu/%llu\n",
		minimum_unknown_cad_tolerance, maximum_unknown_cad_tolerance,
		static_cast<unsigned long long>(unknown_cad_tolerance_covers_clearance),
		static_cast<unsigned long long>(kinds[2] - unknown_without_cad_edge));
	std::printf("unknown by declared CAD incidence:");
	for (const auto& [incidence, count] : unknown_by_declared_incidence)
		std::printf(" %u=%llu", incidence, static_cast<unsigned long long>(count));
	std::printf("\n");

	auto print_top = [](const char* label, const auto& counts)
	{
		std::multimap<std::uint64_t, std::uint32_t, std::greater<>> ordered;
		for (const auto& [id, count] : counts) ordered.emplace(count, id);
		std::printf("top %s:", label);
		unsigned emitted = 0;
		for (const auto& [count, id] : ordered)
		{
			if (emitted++ == 12) break;
			std::printf(" %u=%llu", id, static_cast<unsigned long long>(count));
		}
		std::printf("\n");
	};
	print_top("unknown CAD edges", unknown_by_cad_edge);
	print_top("unknown source faces", unknown_by_source_face);
	std::printf("exact fabric components: %zu; largest area %.9g m^2; quality warnings/blockers: %zu / %zu\n",
		imported.quality.exact_connected_component_count,
		imported.quality.largest_component_area_m2, imported.quality.warning_count(),
		imported.quality.blocker_count());
	for (const GeometryIssue& issue : imported.quality.issues)
	{
		std::printf("  disconnected component faces/triangles/area/boundary-lines/length: "
			"%zu / %zu / %.9g m^2 / %zu / %.9g m\n", issue.source_face_ids.size(),
			issue.source_triangle_ids.size(), issue.surface_area_m2,
			issue.boundary_polylines.size(), issue.boundary_length_m);
		std::printf("    source faces:");
		for (std::uint32_t face : issue.source_face_ids) std::printf(" %u", face);
		std::printf("\n");
	}
	return 0;
}
