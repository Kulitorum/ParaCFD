// tri_mesh.h — a plain, dependency-free triangle mesh in SI METRES.
//
// This struct is the hand-off contract between the OpenCascade-backed STEP importer and
// the OCC-free paraglider geometry/CFD core. It intentionally contains no OpenCascade or Qt
// types, so BVH/AMR/EB preprocessing can stay isolated from the CAD dependency.
//
// ⚠ Units: METRES. The STEP importer scales OCC's native millimetres by 0.001 before
// filling this; every consumer (display, BVH, AMR, and EB) works in metres.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace paracfd::core
{
	// Per-half-edge provenance state.  `none` is affirmative: this half-edge is an
	// ordinary tessellation/synthetic edge rather than a mapped CAD boundary.
	// `unknown_boundary` is deliberately distinct; it means the source face boundary
	// could not be mapped unambiguously and must never be inferred to be a fluid opening.
	enum class CadEdgeProvenanceState : std::uint8_t
	{
		none = 0,
		known = 1,
		unknown_boundary = 2
	};

	// A lit triangle mesh in SI METRES with face-oriented per-vertex normals. No global
	// watertight/outward shell orientation is required. Flat SoA
	// arrays: positions/normals are 3 floats per vertex (x,y,z); indices are triangle vertex
	// indices, 3 per triangle (0-based). bbox_min/max are the axis-aligned bounds in metres.
	struct TriMesh
	{
		// OpenCascade import retains its transformed SI coordinates here before the
		// display/GPU copy is rounded to float.  Programmatic meshes may omit this
		// optional CPU geometry sidecar; vertex_position_double() then promotes the
		// ordinary float coordinate instead.
		std::vector<double> positions_fp64;   // optional 3 * vertex_count, metres
		std::vector<float> positions;       // 3 * vertex_count, metres
		std::vector<float> normals;         // 3 * vertex_count, unit, area-weighted per vertex
		// Optional face-local CAD parameters. STEP import duplicates vertices per source face,
		// therefore one (u,v) pair per mesh vertex retains the BRep tessellation provenance
		// without exposing OpenCascade types to the CFD core. Missing UV data is stored as NaN.
		std::vector<float> vertex_uv;       // 2 * vertex_count, source-face parameter coordinates
		// Optional discrete-topology identity, parallel to the face-local render vertices.
		// Vertices at a conforming CAD seam/contact deliberately remain duplicated because
		// their UV coordinates and shading normals are face-local, but every copy receives
		// the same topology ID. Distinct IDs may occupy the same position (for example two
		// intentional fabric boundaries); consumers must never weld them by proximity.
		// Programmatic/legacy meshes may omit this sidecar, in which case the ordinary
		// vertex index is the topology identity.
		std::vector<std::uint32_t> topology_vertex_ids; // optional vertex_count
		std::vector<std::uint32_t> indices; // 3 * triangle_count
		// Stable within one imported STEP shape: the zero-based TopoDS face traversal index that
		// produced each triangle. This is the key used to accumulate CFD patches back to CAD faces.
		std::vector<std::uint32_t> source_face_ids; // triangle_count

		// OCC-free CAD-boundary provenance, parallel to `indices`. Half-edge h of triangle t is
		// the directed final-mesh edge indices[3*t+h] -> indices[3*t+(h+1)%3]. CAD edge IDs are
		// zero-based TopoDS edge traversal IDs and therefore stay stable when the same imported
		// shape is re-tessellated. Internal tessellation and newly introduced clipping edges use
		// state `none`; unavailable/ambiguous source-face boundary mappings use
		// `unknown_boundary`. Both carry kNoCadEdgeId with incident-face count and tolerance zero,
		// but only the former may be inferred to be an opening. Valid tolerances are the source
		// BRep edge tolerances in metres.
		// The incident count is the number of incident TopoDS_Face uses (half-sheet sectors),
		// not the number of unique face objects: a periodic seam used twice by one face has
		// count two. It is not limited to the manifold value two.
		static constexpr std::uint32_t kNoCadEdgeId = std::numeric_limits<std::uint32_t>::max();
		std::vector<std::uint8_t> triangle_cad_edge_provenance_states; // 3 * triangle_count
		std::vector<std::uint32_t> triangle_cad_edge_ids; // 3 * triangle_count
		std::vector<std::uint32_t> triangle_cad_edge_incident_face_counts; // 3 * triangle_count
		std::vector<double> triangle_cad_edge_tolerances; // 3 * triangle_count, metres; 0 if unavailable
		std::vector<std::uint8_t> triangle_cad_edge_is_periodic_seam; // 3 * triangle_count, 0/1
		// Optional exact BRep edge/face contact certificate.  OpenCascade determines this
		// before tessellation, so an independently tessellated rib/skin T-junction does not
		// depend on polygonal coincidence.  `contact_ids` identify the source CAD contact
		// curve (not a tessellation-edge equality class); `certified_fan_degrees` count the
		// physical half-sheet sectors along it.  Zero/no_contact means no extra certificate,
		// never an opening.  Ordinary programmatic meshes may omit both arrays.
		static constexpr std::uint64_t kNoCadContactId = std::numeric_limits<std::uint64_t>::max();
		std::vector<std::uint64_t> triangle_cad_edge_contact_ids; // 3 * triangle_count
		std::vector<std::uint32_t> triangle_cad_edge_certified_fan_degrees; // 3 * triangle_count
		// Exact conformer-atom identity, including fan-one/open-boundary atoms.  This
		// sidecar is deliberately separate from the public contact certificate above:
		// fan-one atoms are not contacts, but two different exact CAD atoms must not be
		// collapsed merely because they have the same endpoint topology IDs.
		static constexpr std::uint64_t kNoCadEdgeAtomId = std::numeric_limits<std::uint64_t>::max();
		std::vector<std::uint64_t> triangle_cad_edge_atom_ids; // optional 3 * triangle_count

		std::array<float, 3> bbox_min{ { 0.0f, 0.0f, 0.0f } };
		std::array<float, 3> bbox_max{ { 0.0f, 0.0f, 0.0f } };

		std::size_t vertex_count() const { return positions.size() / 3; }
		std::size_t triangle_count() const { return indices.size() / 3; }
		bool empty() const { return indices.empty(); }
		bool has_fp64_positions() const { return positions_fp64.size() == 3 * vertex_count(); }
		bool has_uv() const { return vertex_uv.size() == 2 * vertex_count(); }
		bool has_topology_vertex_ids() const
		{
			return !positions.empty() && topology_vertex_ids.size() == vertex_count();
		}
		bool has_malformed_topology_vertex_ids() const
		{
			return !topology_vertex_ids.empty()
				&& topology_vertex_ids.size() != vertex_count();
		}
		bool has_face_provenance() const { return source_face_ids.size() == triangle_count(); }
		bool has_cad_edge_provenance() const
		{
			return !indices.empty() && triangle_cad_edge_provenance_states.size() == indices.size()
				&& triangle_cad_edge_ids.size() == indices.size()
				&& triangle_cad_edge_incident_face_counts.size() == indices.size()
				&& triangle_cad_edge_tolerances.size() == indices.size()
				&& triangle_cad_edge_is_periodic_seam.size() == indices.size();
		}
		bool has_cad_edge_contact_provenance() const
		{
			return !indices.empty() && triangle_cad_edge_contact_ids.size() == indices.size()
				&& triangle_cad_edge_certified_fan_degrees.size() == indices.size();
		}
		bool has_cad_edge_atom_provenance() const
		{
			return !indices.empty() && triangle_cad_edge_atom_ids.size() == indices.size();
		}
		bool has_malformed_cad_edge_atom_provenance() const
		{
			return !triangle_cad_edge_atom_ids.empty()
				&& triangle_cad_edge_atom_ids.size() != indices.size();
		}

		CadEdgeProvenanceState cad_edge_provenance_state(std::size_t triangle,
			unsigned half_edge) const
		{
			const std::size_t offset = 3 * triangle + half_edge;
			return has_cad_edge_provenance() && half_edge < 3
				&& offset < triangle_cad_edge_provenance_states.size()
				? static_cast<CadEdgeProvenanceState>(triangle_cad_edge_provenance_states[offset])
				: CadEdgeProvenanceState::none;
		}

		bool cad_edge_provenance_unknown(std::size_t triangle, unsigned half_edge) const
		{
			return cad_edge_provenance_state(triangle, half_edge)
				== CadEdgeProvenanceState::unknown_boundary;
		}

		std::array<double, 3> vertex_position_double(std::size_t vertex) const
		{
			const std::size_t offset = 3 * vertex;
			if (has_fp64_positions())
				return { { positions_fp64[offset], positions_fp64[offset + 1], positions_fp64[offset + 2] } };
			return { { static_cast<double>(positions[offset]), static_cast<double>(positions[offset + 1]),
				static_cast<double>(positions[offset + 2]) } };
		}

		std::uint32_t topology_vertex_id(std::size_t vertex) const
		{
			return has_topology_vertex_ids() && vertex < topology_vertex_ids.size()
				? topology_vertex_ids[vertex] : static_cast<std::uint32_t>(vertex);
		}

		std::uint32_t cad_edge_id(std::size_t triangle, unsigned half_edge) const
		{
			const std::size_t offset = 3 * triangle + half_edge;
			return has_cad_edge_provenance() && half_edge < 3 && offset < triangle_cad_edge_ids.size()
				? triangle_cad_edge_ids[offset] : kNoCadEdgeId;
		}

		std::uint32_t cad_edge_incident_face_count(std::size_t triangle, unsigned half_edge) const
		{
			const std::size_t offset = 3 * triangle + half_edge;
			return has_cad_edge_provenance() && half_edge < 3
				&& offset < triangle_cad_edge_incident_face_counts.size()
				? triangle_cad_edge_incident_face_counts[offset] : 0u;
		}

		double cad_edge_tolerance(std::size_t triangle, unsigned half_edge) const
		{
			const std::size_t offset = 3 * triangle + half_edge;
			return has_cad_edge_provenance() && half_edge < 3
				&& offset < triangle_cad_edge_tolerances.size()
				? triangle_cad_edge_tolerances[offset] : 0.0;
		}

		bool cad_edge_is_periodic_seam(std::size_t triangle, unsigned half_edge) const
		{
			const std::size_t offset = 3 * triangle + half_edge;
			return has_cad_edge_provenance() && half_edge < 3
				&& offset < triangle_cad_edge_is_periodic_seam.size()
				&& triangle_cad_edge_is_periodic_seam[offset] != 0;
		}

		std::uint64_t cad_edge_contact_id(std::size_t triangle, unsigned half_edge) const
		{
			const std::size_t offset = 3 * triangle + half_edge;
			return has_cad_edge_contact_provenance() && half_edge < 3
				&& offset < triangle_cad_edge_contact_ids.size()
				? triangle_cad_edge_contact_ids[offset] : kNoCadContactId;
		}

		std::uint32_t cad_edge_certified_fan_degree(std::size_t triangle,
			unsigned half_edge) const
		{
			const std::size_t offset = 3 * triangle + half_edge;
			return has_cad_edge_contact_provenance() && half_edge < 3
				&& offset < triangle_cad_edge_certified_fan_degrees.size()
				? triangle_cad_edge_certified_fan_degrees[offset] : 0u;
		}

		std::uint64_t cad_edge_atom_id(std::size_t triangle, unsigned half_edge) const
		{
			const std::size_t offset = 3 * triangle + half_edge;
			return has_cad_edge_atom_provenance() && half_edge < 3
				&& offset < triangle_cad_edge_atom_ids.size()
				? triangle_cad_edge_atom_ids[offset] : kNoCadEdgeAtomId;
		}

		std::array<float, 3> bbox_size() const
		{
			return { { bbox_max[0] - bbox_min[0], bbox_max[1] - bbox_min[1], bbox_max[2] - bbox_min[2] } };
		}
	};
}
