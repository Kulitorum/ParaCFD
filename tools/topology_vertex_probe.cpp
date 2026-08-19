#include "core/geometry/embedded_boundary.h"
#include "core/geometry/mesh_clip.h"
#include "core/geometry/triangle_bvh.h"

#include <array>
#include <cstdint>
#include <cstdio>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace paracfd::core;

namespace
{
	void require(bool condition, const char* message)
	{
		if (!condition) throw std::runtime_error(message);
	}

	void append_triangle(TriMesh& mesh, const std::array<Vec3d, 3>& points,
		const std::array<std::uint32_t, 3>& topology, std::uint32_t source_face)
	{
		const std::uint32_t base = static_cast<std::uint32_t>(mesh.vertex_count());
		for (int corner = 0; corner < 3; ++corner)
		{
			const Vec3d point = points[corner];
			mesh.positions.insert(mesh.positions.end(), {static_cast<float>(point.x),
				static_cast<float>(point.y), static_cast<float>(point.z)});
			mesh.positions_fp64.insert(mesh.positions_fp64.end(), {point.x, point.y, point.z});
			mesh.topology_vertex_ids.push_back(topology[corner]);
		}
		mesh.indices.insert(mesh.indices.end(), {base, base + 1, base + 2});
		mesh.source_face_ids.push_back(source_face);
	}

	void initialize_cad_sidecars(TriMesh& mesh)
	{
		const std::size_t half_edges=mesh.indices.size();
		mesh.triangle_cad_edge_provenance_states.assign(half_edges,
			static_cast<std::uint8_t>(CadEdgeProvenanceState::none));
		mesh.triangle_cad_edge_ids.assign(half_edges,TriMesh::kNoCadEdgeId);
		mesh.triangle_cad_edge_incident_face_counts.assign(half_edges,0u);
		mesh.triangle_cad_edge_tolerances.assign(half_edges,0.0);
		mesh.triangle_cad_edge_is_periodic_seam.assign(half_edges,0u);
		mesh.triangle_cad_edge_contact_ids.assign(half_edges,TriMesh::kNoCadContactId);
		mesh.triangle_cad_edge_certified_fan_degrees.assign(half_edges,0u);
	}

	TriMesh three_sector_junction(double edge_length)
	{
		TriMesh mesh;
		const Vec3d a{0.2, 0.2, 0.5};
		const Vec3d b{0.2 + edge_length, 0.2, 0.5};
		append_triangle(mesh, {a, b, {0.5, 0.8, 0.5}}, {10, 11, 20}, 0);
		append_triangle(mesh, {a, b, {0.5, 0.2, 0.9}}, {10, 11, 21}, 1);
		append_triangle(mesh, {a, b, {0.5, 0.2, 0.1}}, {10, 11, 22}, 2);
		mesh.bbox_min = {0.2f, 0.2f, 0.1f};
		mesh.bbox_max = {0.8f, 0.8f, 0.9f};
		return mesh;
	}

	TriMesh face_local_square(bool reverse=false)
	{
		TriMesh mesh;
		auto first=[&]{append_triangle(mesh, {{{0,0,0},{1,0,0},{1,1,0}}}, {10,11,12}, 0);};
		auto second=[&]{append_triangle(mesh, {{{0,0,0},{1,1,0},{0,1,0}}}, {10,12,13}, 1);};
		if(reverse){second();first();}else{first();second();}
		mesh.bbox_min={0,0,0};
		mesh.bbox_max={1,1,0};
		return mesh;
	}

	std::vector<std::uint32_t> topology_ids_at(const TriMesh& mesh,Vec3d point)
	{
		std::vector<std::uint32_t> result;
		for(std::size_t vertex=0;vertex<mesh.vertex_count();++vertex)
		{
			const auto p=mesh.vertex_position_double(vertex);
			if(p[0]==point.x&&p[1]==point.y&&p[2]==point.z)
				result.push_back(mesh.topology_vertex_id(vertex));
		}
		return result;
	}

	ExactCellInput captured_input;
	int capture_calls = 0;

	ExactCellDecomposition capture_exact_input(const ExactCellInput& input)
	{
		captured_input = input;
		++capture_calls;
		ExactCellDecomposition result;
		result.errors.push_back("intentional topology-input capture");
		return result;
	}

	std::size_t capture_exact_vertex_count(const TriMesh& mesh)
	{
		captured_input = {};
		capture_calls = 0;
		const TriangleBvh bvh(mesh);
		EmbeddedBoundaryBuildOptions options;
		options.complex_subdivisions = 2;
		options.exact_cell_decomposer = capture_exact_input;
		const EmbeddedBoundary boundary = build_embedded_boundary(mesh, bvh,
			{{0, 0, 0}, 1, 1, 1, 1.0}, options);
		require(!boundary.ready_for_flow(),
			"capture decomposer did not leave the intentionally rejected cell unresolved");
		require(capture_calls == 1 && captured_input.triangles.size() == 3,
			"complex junction did not reach the exact-cell callback exactly once");
		return captured_input.vertices.size();
	}
}

int main()
{
	try
	{
		// The edge is intentionally shorter than the FP64 clearance band. Geometric
		// recovery cannot establish it, so this specifically exercises the indexed
		// topology sidecar rather than the legacy proximity path.
		const TriMesh indexed = three_sector_junction(1.0e-13);
		TriangleBvh indexed_bvh(indexed);
		std::uint64_t attachment = FabricEdgeCertificate::no_attachment;
		for (std::uint32_t triangle = 0; triangle < 3; ++triangle)
		{
			const FabricEdgeCertificate certificate = indexed_bvh.edge_certificate(triangle, 0);
			require(certificate.kind == FabricEdgeKind::attached
				&& certificate.role == FabricEdgeRole::junction
				&& certificate.incident_fan_degree == 3,
				"canonical topology IDs did not certify the three-sector junction");
			if (triangle == 0) attachment = certificate.attachment_id;
			else require(certificate.attachment_id == attachment,
				"junction half-edges did not receive one reciprocal attachment token");
		}
		require(attachment != FabricEdgeCertificate::no_attachment,
			"indexed junction did not receive an attachment token");

		TriMesh legacy = indexed;
		legacy.topology_vertex_ids.clear();
		TriangleBvh legacy_bvh(legacy);
		require(legacy_bvh.edge_certificate(0, 0).kind == FabricEdgeKind::unknown,
			"short face-local edge unexpectedly bypassed the legacy ambiguity path");

		const TriMesh conforming = three_sector_junction(0.6);
		const std::size_t conforming_vertices = capture_exact_vertex_count(conforming);
		TriMesh face_local = conforming;
		face_local.topology_vertex_ids.clear();
		const std::size_t face_local_vertices = capture_exact_vertex_count(face_local);
		require(conforming_vertices == 5 && face_local_vertices == 9,
			"exact-cell input did not reuse canonical topology vertices with legacy fallback");
		for (const ExactCellTriangle& triangle : captured_input.triangles)
			for (std::uint32_t vertex : triangle.vertices)
				require(vertex < captured_input.vertices.size(),
					"exact-cell fallback emitted an out-of-range vertex index");

		TriMesh coincident_but_distinct;
		const Vec3d edge_a{0.1,0.1,0.5},edge_b{0.9,0.1,0.5};
		append_triangle(coincident_but_distinct,{edge_a,edge_b,{0.5,0.8,0.5}},
			{100,101,102},0);
		append_triangle(coincident_but_distinct,{edge_a,edge_b,{0.5,0.1,0.9}},
			{200,201,202},1);
		TriangleBvh distinct_bvh(coincident_but_distinct);
		require(distinct_bvh.edge_certificate(0,0).kind==FabricEdgeKind::confirmed_free
			&&distinct_bvh.edge_certificate(1,0).kind==FabricEdgeKind::confirmed_free,
			"distinct explicit topology IDs were proximity-welded at identical coordinates");

		TriMesh incomplete_fan;
		append_triangle(incomplete_fan,{edge_a,edge_b,{0.5,0.8,0.5}},
			{300,301,302},0);
		append_triangle(incomplete_fan,{edge_a,edge_b,{0.5,0.1,0.9}},
			{300,301,303},1);
		initialize_cad_sidecars(incomplete_fan);
		for(const std::size_t half_edge:{std::size_t{0},std::size_t{3}})
		{
			incomplete_fan.triangle_cad_edge_contact_ids[half_edge]=77;
			incomplete_fan.triangle_cad_edge_certified_fan_degrees[half_edge]=3;
		}
		TriangleBvh incomplete_bvh(incomplete_fan);
		for(std::uint32_t triangle=0;triangle<2;++triangle)
		{
			const FabricEdgeCertificate certificate=incomplete_bvh.edge_certificate(triangle,0);
			require(certificate.kind==FabricEdgeKind::unknown
				&&certificate.unknown_reason==FabricEdgeUnknownReason::indexed_fan_mismatch,
				"indexed edge group ignored a larger exact CAD contact fan");
		}

		const TriMesh clipped=clip_mesh_to_axis_slab(face_local_square(),0,0.25,1.0);
		require(clipped.has_topology_vertex_ids(),
			"axis-slab clipping discarded the topology vertex sidecar");
		const auto retained_bottom=topology_ids_at(clipped,{1,0,0});
		const auto retained_top=topology_ids_at(clipped,{1,1,0});
		require(!retained_bottom.empty()&&std::set<std::uint32_t>(retained_bottom.begin(),retained_bottom.end())==std::set<std::uint32_t>{11}
			&&!retained_top.empty()&&std::set<std::uint32_t>(retained_top.begin(),retained_top.end())==std::set<std::uint32_t>{12},
			"axis-slab clipping did not preserve retained source topology IDs");
		const auto bottom_cut=topology_ids_at(clipped,{0.25,0,0});
		const auto shared_cut=topology_ids_at(clipped,{0.25,0.25,0});
		const auto top_cut=topology_ids_at(clipped,{0.25,1,0});
		require(bottom_cut.size()==1&&shared_cut.size()>=2&&top_cut.size()==1,
			"axis-slab clipping did not emit the expected cut-edge vertices");
		require(std::set<std::uint32_t>(shared_cut.begin(),shared_cut.end()).size()==1,
			"adjacent triangles did not reuse one clipped topological edge point");
		require(bottom_cut.front()!=shared_cut.front()&&bottom_cut.front()!=top_cut.front()
			&&shared_cut.front()!=top_cut.front(),
			"distinct cut-edge points received the same generated topology ID");
		const TriMesh clipped_replay=clip_mesh_to_axis_slab(face_local_square(),0,0.25,1.0);
		require(clipped.topology_vertex_ids==clipped_replay.topology_vertex_ids,
			"axis-slab generated topology IDs are not deterministic");
		const TriMesh clipped_reordered=clip_mesh_to_axis_slab(face_local_square(true),0,0.25,1.0);
		const auto reordered_bottom=topology_ids_at(clipped_reordered,{0.25,0,0});
		const auto reordered_shared=topology_ids_at(clipped_reordered,{0.25,0.25,0});
		const auto reordered_top=topology_ids_at(clipped_reordered,{0.25,1,0});
		require(!reordered_bottom.empty()&&!reordered_shared.empty()&&!reordered_top.empty()
			&&reordered_bottom.front()==bottom_cut.front()
			&&reordered_shared.front()==shared_cut.front()
			&&reordered_top.front()==top_cut.front(),
			"generated topology IDs depend on source triangle order");
		TriangleBvh clipped_bvh(clipped);
		require(!clipped_bvh.empty(),"valid clipped topology was rejected by the BVH boundary");

		TriMesh contact_clip_source;
		append_triangle(contact_clip_source,{{{0,0,0},{1,0,0},{0,1,0}}},
			{400,401,402},0);
		initialize_cad_sidecars(contact_clip_source);
		contact_clip_source.triangle_cad_edge_tolerances[0]=2.0e-6;
		contact_clip_source.triangle_cad_edge_contact_ids[0]=91;
		contact_clip_source.triangle_cad_edge_certified_fan_degrees[0]=3;
		const TriMesh contact_clipped=clip_mesh_to_axis_slab(contact_clip_source,0,0.25,1.0);
		bool retained_interior_contact=false;
		for(std::size_t triangle=0;triangle<contact_clipped.triangle_count();++triangle)
			for(unsigned half_edge=0;half_edge<3;++half_edge)
				retained_interior_contact=retained_interior_contact||(
					contact_clipped.cad_edge_id(triangle,half_edge)==TriMesh::kNoCadEdgeId
					&&contact_clipped.cad_edge_contact_id(triangle,half_edge)==91
					&&contact_clipped.cad_edge_certified_fan_degree(triangle,half_edge)==3
					&&contact_clipped.cad_edge_tolerance(triangle,half_edge)==2.0e-6);
		require(retained_interior_contact,
			"axis-slab clipping dropped a face-interior CAD contact certificate");

		TriMesh inconsistent=three_sector_junction(0.6);
		inconsistent.positions_fp64[3*3]+=1e-6;
		bool rejected=false;
		try{TriangleBvh invalid_bvh(inconsistent);(void)invalid_bvh;}
		catch(const std::invalid_argument& error)
		{
			rejected=std::string(error.what())==
				"TriangleBvh topology vertex ID 10 has inconsistent coordinates at vertices 0 and 3";
		}
		require(rejected,
			"BVH construction did not deterministically reject one topology ID at two coordinates");

		TriMesh malformed=three_sector_junction(0.6);
		malformed.topology_vertex_ids.pop_back();
		TriangleBvh reusable(conforming);
		bool malformed_rejected=false;
		try{reusable.build(malformed);}
		catch(const std::invalid_argument& error)
		{
			malformed_rejected=std::string(error.what()).find(
				"topology vertex sidecar")!=std::string::npos;
		}
		require(malformed_rejected&&reusable.empty(),
			"malformed topology sidecar failed open or left a stale reusable BVH");
		bool clip_rejected=false;
		try{(void)clip_mesh_to_axis_slab(malformed,0,0.25,1.0);}
		catch(const std::invalid_argument&){clip_rejected=true;}
		require(clip_rejected,
			"axis-slab clipping silently discarded a malformed topology sidecar");

		std::puts("topology vertex IDs: PASS");
		return 0;
	}
	catch (const std::exception& error)
	{
		std::fprintf(stderr, "topology vertex IDs: FAIL: %s\n", error.what());
		return 1;
	}
}
