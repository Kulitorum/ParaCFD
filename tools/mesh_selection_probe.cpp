#include "core/geometry/mesh_selection.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

using namespace paracfd::core;

namespace
{
	void require(bool condition, const char* message)
	{
		if (!condition) throw std::runtime_error(message);
	}
}

int main()
{
	try
	{
		TriMesh source;
		source.positions = {0,0,0, 1,0,0, 0,1,0, 100,0,0, 101,0,0, 100,1,0};
		source.positions_fp64.assign(source.positions.begin(), source.positions.end());
		source.normals = {0,0,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1};
		source.vertex_uv = {0,0, 1,0, 0,1, 0,0, 1,0, 0,1};
		source.topology_vertex_ids = {101,102,103,201,202,203};
		source.indices = {0,1,2, 3,4,5};
		source.source_face_ids = {7,11};
		source.triangle_cad_edge_provenance_states.assign(6,
			static_cast<std::uint8_t>(CadEdgeProvenanceState::known));
		source.triangle_cad_edge_ids = {1,2,3,4,5,6};
		source.triangle_cad_edge_incident_face_counts.assign(6,1);
		source.triangle_cad_edge_tolerances.assign(6,1e-7);
		source.triangle_cad_edge_is_periodic_seam.assign(6,0);
		source.triangle_cad_edge_contact_ids = {10,20,30,40,50,60};
		source.triangle_cad_edge_certified_fan_degrees.assign(6,2);
		source.bbox_min = {0,0,0}; source.bbox_max = {101,1,0};

		const std::vector<std::uint32_t> excluded{1};
		const TriMeshSelection selected = exclude_triangles(source, excluded);
		require(selected.mesh.triangle_count() == 1, "wrong selected triangle count");
		require(selected.mesh.vertex_count() == 3, "unused vertices were retained");
		require(selected.source_triangle_ids == std::vector<std::uint32_t>{0},
			"source triangle map was not retained");
		require(selected.mesh.source_face_ids == std::vector<std::uint32_t>{7},
			"CAD face provenance was not retained");
		require(selected.mesh.has_topology_vertex_ids()
			&& selected.mesh.topology_vertex_ids == std::vector<std::uint32_t>({101,102,103}),
			"discrete topology vertex IDs were not retained");
		require(selected.mesh.has_cad_edge_provenance(), "CAD edge provenance was not retained");
		require(selected.mesh.has_cad_edge_contact_provenance(),
			"CAD contact provenance was not retained");
		require(std::abs(selected.mesh.bbox_max[0] - 1.0f) < 1e-7f,
			"excluded geometry still influenced bbox");
		TriMesh legacy;
		legacy.positions = {0,0,0, 1,0,0};
		require(!legacy.has_topology_vertex_ids() && legacy.topology_vertex_id(1) == 1,
			"legacy meshes do not fall back to their render vertex indices");
		std::puts("mesh selection: PASS");
		return 0;
	}
	catch (const std::exception& error)
	{
		std::fprintf(stderr, "mesh selection: FAIL: %s\n", error.what());
		return 1;
	}
}
