#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace paracfd::core
{
	struct UvPoint
	{
		double u = 0.0;
		double v = 0.0;
	};

	inline bool operator==(const UvPoint& a, const UvPoint& b)
	{
		return a.u == b.u && a.v == b.v;
	}

	struct UvVertex
	{
		static constexpr std::uint32_t no_topology_id =
			std::numeric_limits<std::uint32_t>::max();

		UvPoint uv;
		std::uint32_t topology_id = no_topology_id;
	};

	struct UvTriangle
	{
		std::array<std::uint32_t, 3> vertices{};
	};

	struct UvConstraintSample
	{
		UvPoint uv;
		std::uint32_t topology_id = UvVertex::no_topology_id;
	};

	struct UvConstraintOptions
	{
		// Absolute tolerance in the face's native UV units. It is used only to classify
		// a supplied sample on an existing vertex/edge; it is never used to weld mesh
		// vertices to one another. Callers should derive this from the face/pcurve audit.
		double coordinate_tolerance = 1.0e-12;
	};

	struct UvEdgeView
	{
		std::array<std::uint32_t, 2> vertices{};
		bool boundary = false;
		std::vector<std::uint64_t> constraint_ids;
	};

	// OCC-free deterministic planar straight-line-graph insertion layer.
	// Input triangles must be a consistently CCW, two-manifold planar triangulation.
	// Validation rejects nonconforming geometric edge intersections and positive-area
	// triangle overlap in addition to indexed incidence errors. Domain boundaries are
	// inferred from one-sided edges, so inner boundary loops (holes) are retained without
	// an inside/outside fill convention.
	class UvConstraintMesh
	{
	public:
		UvConstraintMesh(std::vector<UvVertex> vertices, std::vector<UvTriangle> triangles,
			UvConstraintOptions options = {});

		// Transactional: on failure neither triangles, vertices, topology IDs nor existing
		// constraint tags are changed. Samples are inserted first, then every prescribed
		// segment cuts a crossed triangle strip and retriangulates its two cavities. No
		// triangulation-dependent Steiner point is added on the prescribed chain.
		bool insert_polyline(std::span<const UvConstraintSample> samples,
			std::uint64_t constraint_id, std::string* error = nullptr);

		// Assigns deterministic face-local IDs to vertices not claimed by a shared CAD
		// sample. `next_id` is a caller-owned global allocator initialized above the
		// reserved shared-ID range; it is advanced on success and reused across faces.
		// UINT32_MAX is the unassigned sentinel and is never allocated.
		bool assign_unclaimed_topology_ids(std::uint32_t& next_id,
			std::string* error = nullptr);

		bool validate(std::string* error = nullptr) const;
		bool contains(const UvPoint& point, bool include_boundary = true) const;
		bool has_edge(std::uint32_t a, std::uint32_t b) const;
		bool edge_has_constraint(std::uint32_t a, std::uint32_t b,
			std::uint64_t constraint_id) const;
		std::optional<std::uint32_t> vertex_for_topology_id(std::uint32_t topology_id) const;
		std::vector<UvEdgeView> edges() const;
		std::size_t boundary_loop_count() const;
		double area() const;

		const std::vector<UvVertex>& vertices() const { return vertices_; }
		const std::vector<UvTriangle>& triangles() const { return triangles_; }

	private:
		using Edge = std::array<std::uint32_t, 2>;

		struct EdgeLess
		{
			bool operator()(const Edge& a, const Edge& b) const
			{
				return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
			}
		};

		static Edge edge_key(std::uint32_t a, std::uint32_t b);
		bool insert_sample(const UvConstraintSample& sample, std::uint32_t& vertex,
			std::string& error);
		bool split_edge(const Edge& edge, const UvConstraintSample& sample,
			std::uint32_t& vertex, std::string& error);
		bool split_triangle(std::size_t triangle, const UvConstraintSample& sample,
			std::uint32_t& vertex, std::string& error);
		bool insert_atomic_segment(std::uint32_t a, std::uint32_t b,
			std::uint64_t constraint_id, std::string& error);
		bool triangulate_polygon(const std::vector<std::uint32_t>& polygon,
			std::vector<UvTriangle>& output, std::string& error) const;
		bool claim_topology_id(std::uint32_t vertex, std::uint32_t topology_id,
			std::string& error);
		void add_constraint_tag(const Edge& edge, std::uint64_t constraint_id);

		std::vector<UvVertex> vertices_;
		std::vector<UvTriangle> triangles_;
		UvConstraintOptions options_;
		std::map<Edge, std::vector<std::uint64_t>, EdgeLess> constraint_tags_;
	};
}
