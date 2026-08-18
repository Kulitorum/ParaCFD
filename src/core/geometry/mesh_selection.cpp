#include "core/geometry/mesh_selection.h"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace paracfd::core
{
	namespace
	{
		void append_vertex(TriMesh& destination, const TriMesh& source, std::size_t vertex)
		{
			if (vertex >= source.vertex_count())
				throw std::invalid_argument("TriMesh selection contains an out-of-range vertex index");
			const std::size_t xyz = 3 * vertex;
			destination.positions.insert(destination.positions.end(),
				source.positions.begin() + xyz, source.positions.begin() + xyz + 3);
			if (source.has_fp64_positions())
				destination.positions_fp64.insert(destination.positions_fp64.end(),
					source.positions_fp64.begin() + xyz, source.positions_fp64.begin() + xyz + 3);
			if (source.normals.size() == source.positions.size())
				destination.normals.insert(destination.normals.end(),
					source.normals.begin() + xyz, source.normals.begin() + xyz + 3);
			if (source.has_uv())
			{
				const std::size_t uv = 2 * vertex;
				destination.vertex_uv.insert(destination.vertex_uv.end(),
					source.vertex_uv.begin() + uv, source.vertex_uv.begin() + uv + 2);
			}
		}

		void append_half_edge_metadata(TriMesh& destination, const TriMesh& source,
			std::size_t triangle)
		{
			const std::size_t first = 3 * triangle;
			if (source.has_cad_edge_provenance())
			{
				destination.triangle_cad_edge_provenance_states.insert(
					destination.triangle_cad_edge_provenance_states.end(),
					source.triangle_cad_edge_provenance_states.begin() + first,
					source.triangle_cad_edge_provenance_states.begin() + first + 3);
				destination.triangle_cad_edge_ids.insert(destination.triangle_cad_edge_ids.end(),
					source.triangle_cad_edge_ids.begin() + first,
					source.triangle_cad_edge_ids.begin() + first + 3);
				destination.triangle_cad_edge_incident_face_counts.insert(
					destination.triangle_cad_edge_incident_face_counts.end(),
					source.triangle_cad_edge_incident_face_counts.begin() + first,
					source.triangle_cad_edge_incident_face_counts.begin() + first + 3);
				destination.triangle_cad_edge_tolerances.insert(
					destination.triangle_cad_edge_tolerances.end(),
					source.triangle_cad_edge_tolerances.begin() + first,
					source.triangle_cad_edge_tolerances.begin() + first + 3);
				destination.triangle_cad_edge_is_periodic_seam.insert(
					destination.triangle_cad_edge_is_periodic_seam.end(),
					source.triangle_cad_edge_is_periodic_seam.begin() + first,
					source.triangle_cad_edge_is_periodic_seam.begin() + first + 3);
			}
			if (source.has_cad_edge_contact_provenance())
			{
				destination.triangle_cad_edge_contact_ids.insert(
					destination.triangle_cad_edge_contact_ids.end(),
					source.triangle_cad_edge_contact_ids.begin() + first,
					source.triangle_cad_edge_contact_ids.begin() + first + 3);
				destination.triangle_cad_edge_certified_fan_degrees.insert(
					destination.triangle_cad_edge_certified_fan_degrees.end(),
					source.triangle_cad_edge_certified_fan_degrees.begin() + first,
					source.triangle_cad_edge_certified_fan_degrees.begin() + first + 3);
			}
		}

		void compute_referenced_bbox(TriMesh& mesh)
		{
			if (mesh.empty()) return;
			double lo[3] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
				std::numeric_limits<double>::max()};
			double hi[3] = {-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
				-std::numeric_limits<double>::max()};
			for (std::size_t vertex = 0; vertex < mesh.vertex_count(); ++vertex)
			{
				const auto point = mesh.vertex_position_double(vertex);
				for (int axis = 0; axis < 3; ++axis)
				{
					lo[axis] = std::min(lo[axis], point[axis]);
					hi[axis] = std::max(hi[axis], point[axis]);
				}
			}
			for (int axis = 0; axis < 3; ++axis)
			{
				mesh.bbox_min[axis] = static_cast<float>(lo[axis]);
				mesh.bbox_max[axis] = static_cast<float>(hi[axis]);
			}
		}
	}

	TriMeshSelection select_triangles(const TriMesh& source,
		std::span<const std::uint8_t> keep_triangle)
	{
		if (keep_triangle.size() != source.triangle_count())
			throw std::invalid_argument("TriMesh selection mask size does not match triangle count");

		TriMeshSelection result;
		std::vector<std::uint32_t> vertex_map(source.vertex_count(), TriMesh::kNoCadEdgeId);
		for (std::size_t triangle = 0; triangle < source.triangle_count(); ++triangle)
		{
			if (!keep_triangle[triangle]) continue;
			for (unsigned corner = 0; corner < 3; ++corner)
			{
				const std::uint32_t old_vertex = source.indices[3 * triangle + corner];
				if (old_vertex >= vertex_map.size())
					throw std::invalid_argument("TriMesh selection contains an out-of-range vertex index");
				std::uint32_t& new_vertex = vertex_map[old_vertex];
				if (new_vertex == TriMesh::kNoCadEdgeId)
				{
					new_vertex = static_cast<std::uint32_t>(result.mesh.vertex_count());
					append_vertex(result.mesh, source, old_vertex);
				}
				result.mesh.indices.push_back(new_vertex);
			}
			if (source.has_face_provenance())
				result.mesh.source_face_ids.push_back(source.source_face_ids[triangle]);
			append_half_edge_metadata(result.mesh, source, triangle);
			result.source_triangle_ids.push_back(static_cast<std::uint32_t>(triangle));
		}
		compute_referenced_bbox(result.mesh);
		return result;
	}

	TriMeshSelection exclude_triangles(const TriMesh& source,
		std::span<const std::uint32_t> excluded_triangle_ids)
	{
		std::vector<std::uint8_t> keep(source.triangle_count(), 1);
		for (std::uint32_t triangle : excluded_triangle_ids)
		{
			if (triangle >= keep.size())
				throw std::invalid_argument("TriMesh exclusion contains an out-of-range triangle ID");
			keep[triangle] = 0;
		}
		return select_triangles(source, keep);
	}
}
