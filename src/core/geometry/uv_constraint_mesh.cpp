#include "core/geometry/uv_constraint_mesh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace paracfd::core
{
	namespace
	{
		using Edge = std::array<std::uint32_t, 2>;
		// UINT32_MAX is reserved as the missing topology-ID sentinel. Keeping the
		// indexable vertex count below it also makes every size_t -> uint32_t cast
		// explicit and prevents uint32_t iteration from wrapping.
		constexpr std::size_t maximum_uv_vertex_count =
			static_cast<std::size_t>(UvVertex::no_topology_id);

		long double orient(const UvPoint& a, const UvPoint& b, const UvPoint& c)
		{
			return (static_cast<long double>(b.u) - a.u)
				* (static_cast<long double>(c.v) - a.v)
				- (static_cast<long double>(b.v) - a.v)
				* (static_cast<long double>(c.u) - a.u);
		}

		double squared_distance(const UvPoint& a, const UvPoint& b)
		{
			const double du = a.u - b.u;
			const double dv = a.v - b.v;
			return du * du + dv * dv;
		}

		bool finite(const UvPoint& point)
		{
			return std::isfinite(point.u) && std::isfinite(point.v);
		}

		double segment_parameter_tolerance(const UvPoint& a, const UvPoint& b,
			double coordinate_tolerance)
		{
			const double length = std::hypot(b.u - a.u, b.v - a.v);
			if (!(length > 0.0) || !std::isfinite(length)) return 1.0;
			return std::min(1.0, coordinate_tolerance / length);
		}

		bool point_on_closed_segment_exact(const UvPoint& point, const UvPoint& a,
			const UvPoint& b)
		{
			return orient(a, b, point) == 0.0L
				&& point.u >= std::min(a.u, b.u) && point.u <= std::max(a.u, b.u)
				&& point.v >= std::min(a.v, b.v) && point.v <= std::max(a.v, b.v);
		}

		bool nonconforming_edge_intersection(const Edge& first, const Edge& second,
			const std::vector<UvVertex>& vertices)
		{
			const UvPoint& a = vertices[first[0]].uv;
			const UvPoint& b = vertices[first[1]].uv;
			const UvPoint& c = vertices[second[0]].uv;
			const UvPoint& d = vertices[second[1]].uv;
			const long double ab_c = orient(a, b, c), ab_d = orient(a, b, d);
			const long double cd_a = orient(c, d, a), cd_b = orient(c, d, b);
			if (((ab_c > 0.0L && ab_d < 0.0L) || (ab_c < 0.0L && ab_d > 0.0L))
				&& ((cd_a > 0.0L && cd_b < 0.0L) || (cd_a < 0.0L && cd_b > 0.0L)))
				return true;

			const auto endpoint_is_shared = [](std::uint32_t endpoint, const Edge& edge)
			{
				return endpoint == edge[0] || endpoint == edge[1];
			};
			if (point_on_closed_segment_exact(a, c, d)
				&& !endpoint_is_shared(first[0], second)) return true;
			if (point_on_closed_segment_exact(b, c, d)
				&& !endpoint_is_shared(first[1], second)) return true;
			if (point_on_closed_segment_exact(c, a, b)
				&& !endpoint_is_shared(second[0], first)) return true;
			if (point_on_closed_segment_exact(d, a, b)
				&& !endpoint_is_shared(second[1], first)) return true;
			return false;
		}

		bool separated_by_supporting_edge(const UvTriangle& axes,
			const UvTriangle& candidate, const std::vector<UvVertex>& vertices)
		{
			for (unsigned edge = 0; edge < 3; ++edge)
			{
				const UvPoint& a = vertices[axes.vertices[edge]].uv;
				const UvPoint& b = vertices[axes.vertices[(edge + 1) % 3]].uv;
				bool candidate_reaches_left_half_plane = false;
				for (std::uint32_t vertex : candidate.vertices)
					if (orient(a, b, vertices[vertex].uv) > 0.0L)
					{
						candidate_reaches_left_half_plane = true;
						break;
					}
				if (!candidate_reaches_left_half_plane) return true;
			}
			return false;
		}

		bool triangles_have_positive_area_overlap(const UvTriangle& a,
			const UvTriangle& b, const std::vector<UvVertex>& vertices)
		{
			// Both inputs are strictly CCW here. The separating-axis theorem for two
			// convex triangles says their interiors overlap iff neither triangle has a
			// supporting edge with the other triangle wholly on or outside its right side.
			return !separated_by_supporting_edge(a, b, vertices)
				&& !separated_by_supporting_edge(b, a, vertices);
		}

		bool point_on_segment(const UvPoint& point, const UvPoint& a, const UvPoint& b,
			double tolerance, bool include_ends)
		{
			const double du = b.u - a.u;
			const double dv = b.v - a.v;
			const double length = std::hypot(du, dv);
			if (!(length > tolerance)) return squared_distance(point, a) <= tolerance * tolerance;

			const long double distance_numerator = std::abs(orient(a, b, point));
			if (distance_numerator > static_cast<long double>(tolerance) * length) return false;
			const double projection = ((point.u - a.u) * du + (point.v - a.v) * dv)
				/ (length * length);
			const double parameter_tolerance = tolerance / length;
			return include_ends
				? projection >= -parameter_tolerance && projection <= 1.0 + parameter_tolerance
				: projection > parameter_tolerance && projection < 1.0 - parameter_tolerance;
		}

		bool point_in_triangle(const UvPoint& point, const UvPoint& a, const UvPoint& b,
			const UvPoint& c, double tolerance, bool strict)
		{
			const long double ab = orient(a, b, point);
			const long double bc = orient(b, c, point);
			const long double ca = orient(c, a, point);
			const long double scale = std::max({std::hypot(b.u - a.u, b.v - a.v),
				std::hypot(c.u - b.u, c.v - b.v), std::hypot(a.u - c.u, a.v - c.v)});
			if (!(scale > 0.0L)) return false;
			const long double threshold = static_cast<long double>(tolerance) * scale;
			return strict ? ab > threshold && bc > threshold && ca > threshold
				: ab >= -threshold && bc >= -threshold && ca >= -threshold;
		}

		bool proper_intersection(const UvPoint& a, const UvPoint& b, const UvPoint& c,
			const UvPoint& d, double tolerance, double& segment_parameter)
		{
			const long double ab_c = orient(a, b, c);
			const long double ab_d = orient(a, b, d);
			const long double cd_a = orient(c, d, a);
			const long double cd_b = orient(c, d, b);
			const long double ab_length = std::hypot(b.u - a.u, b.v - a.v);
			const long double cd_length = std::hypot(d.u - c.u, d.v - c.v);
			if (!(ab_length > tolerance) || !(cd_length > tolerance)) return false;
			const long double ab_threshold = static_cast<long double>(tolerance) * ab_length;
			const long double cd_threshold = static_cast<long double>(tolerance) * cd_length;
			if (!((ab_c > ab_threshold && ab_d < -ab_threshold)
				|| (ab_c < -ab_threshold && ab_d > ab_threshold))) return false;
			if (!((cd_a > cd_threshold && cd_b < -cd_threshold)
				|| (cd_a < -cd_threshold && cd_b > cd_threshold))) return false;

			const long double r_u = static_cast<long double>(b.u) - a.u;
			const long double r_v = static_cast<long double>(b.v) - a.v;
			const long double s_u = static_cast<long double>(d.u) - c.u;
			const long double s_v = static_cast<long double>(d.v) - c.v;
			const long double denominator = r_u * s_v - r_v * s_u;
			const long double q_u = static_cast<long double>(c.u) - a.u;
			const long double q_v = static_cast<long double>(c.v) - a.v;
			segment_parameter = static_cast<double>((q_u * s_v - q_v * s_u) / denominator);
			const double other_parameter = static_cast<double>(
				(q_u * r_v - q_v * r_u) / denominator);
			const double parameter_tolerance = segment_parameter_tolerance(a, b, tolerance);
			const double other_parameter_tolerance = segment_parameter_tolerance(c, d, tolerance);
			return segment_parameter > parameter_tolerance
				&& segment_parameter < 1.0 - parameter_tolerance
				&& other_parameter > other_parameter_tolerance
				&& other_parameter < 1.0 - other_parameter_tolerance;
		}

		std::string edge_text(const Edge& edge)
		{
			return "(" + std::to_string(edge[0]) + "," + std::to_string(edge[1]) + ")";
		}
	}

	UvConstraintMesh::UvConstraintMesh(std::vector<UvVertex> vertices,
		std::vector<UvTriangle> triangles, UvConstraintOptions options)
		: vertices_(std::move(vertices)), triangles_(std::move(triangles)), options_(options)
	{
	}

	UvConstraintMesh::Edge UvConstraintMesh::edge_key(std::uint32_t a, std::uint32_t b)
	{
		return a < b ? Edge{a, b} : Edge{b, a};
	}

	void UvConstraintMesh::add_constraint_tag(const Edge& edge, std::uint64_t constraint_id)
	{
		auto& tags = constraint_tags_[edge];
		const auto position = std::lower_bound(tags.begin(), tags.end(), constraint_id);
		if (position == tags.end() || *position != constraint_id) tags.insert(position, constraint_id);
	}

	bool UvConstraintMesh::claim_topology_id(std::uint32_t vertex,
		std::uint32_t topology_id, std::string& error)
	{
		if (vertex >= vertices_.size())
		{
			error = "canonical topology ID claim has an out-of-range UV vertex";
			return false;
		}
		if (topology_id == UvVertex::no_topology_id)
		{
			error = "a constraint sample has no canonical topology ID";
			return false;
		}
		for (std::size_t other = 0; other < vertices_.size(); ++other)
			if (other != static_cast<std::size_t>(vertex)
				&& vertices_[other].topology_id == topology_id)
			{
				error = "canonical topology ID " + std::to_string(topology_id)
					+ " is already assigned to a different UV vertex";
				return false;
			}
		std::uint32_t& current = vertices_[vertex].topology_id;
		if (current != UvVertex::no_topology_id && current != topology_id)
		{
			error = "UV vertex " + std::to_string(vertex) + " already has canonical topology ID "
				+ std::to_string(current) + ", not " + std::to_string(topology_id);
			return false;
		}
		current = topology_id;
		return true;
	}

	bool UvConstraintMesh::split_edge(const Edge& edge, const UvConstraintSample& sample,
		std::uint32_t& vertex, std::string& error)
	{
		std::vector<std::size_t> incident;
		for (std::size_t triangle = 0; triangle < triangles_.size(); ++triangle)
		{
			const auto& indices = triangles_[triangle].vertices;
			for (unsigned i = 0; i < 3; ++i)
				if (edge_key(indices[i], indices[(i + 1) % 3]) == edge)
				{
					incident.push_back(triangle);
					break;
				}
		}
		if (incident.empty() || incident.size() > 2)
		{
			error = "cannot split non-manifold or missing edge " + edge_text(edge);
			return false;
		}
		if (vertices_.size() >= maximum_uv_vertex_count)
		{
			error = "cannot split UV edge: uint32 vertex index capacity exhausted";
			return false;
		}

		vertex = static_cast<std::uint32_t>(vertices_.size());
		vertices_.push_back({sample.uv, sample.topology_id});
		std::vector<UvTriangle> appended;
		appended.reserve(incident.size());
		for (std::size_t triangle : incident)
		{
			const auto old = triangles_[triangle].vertices;
			bool found = false;
			for (unsigned i = 0; i < 3; ++i)
			{
				const std::uint32_t from = old[i];
				const std::uint32_t to = old[(i + 1) % 3];
				if (edge_key(from, to) != edge) continue;
				const std::uint32_t opposite = old[(i + 2) % 3];
				triangles_[triangle].vertices = {from, vertex, opposite};
				appended.push_back({{vertex, to, opposite}});
				found = true;
				break;
			}
			if (!found)
			{
				error = "edge incidence changed during deterministic split";
				return false;
			}
		}
		triangles_.insert(triangles_.end(), appended.begin(), appended.end());

		const auto tag = constraint_tags_.find(edge);
		if (tag != constraint_tags_.end())
		{
			const std::vector<std::uint64_t> inherited = tag->second;
			constraint_tags_.erase(tag);
			constraint_tags_[edge_key(edge[0], vertex)] = inherited;
			constraint_tags_[edge_key(vertex, edge[1])] = inherited;
		}
		return true;
	}

	bool UvConstraintMesh::split_triangle(std::size_t triangle,
		const UvConstraintSample& sample, std::uint32_t& vertex, std::string& error)
	{
		if (triangle >= triangles_.size())
		{
			error = "triangle location became invalid during point insertion";
			return false;
		}
		if (vertices_.size() >= maximum_uv_vertex_count)
		{
			error = "cannot split UV triangle: uint32 vertex index capacity exhausted";
			return false;
		}
		const auto old = triangles_[triangle].vertices;
		vertex = static_cast<std::uint32_t>(vertices_.size());
		vertices_.push_back({sample.uv, sample.topology_id});
		triangles_[triangle].vertices = {old[0], old[1], vertex};
		triangles_.push_back({{old[1], old[2], vertex}});
		triangles_.push_back({{old[2], old[0], vertex}});
		return true;
	}

	bool UvConstraintMesh::insert_sample(const UvConstraintSample& sample,
		std::uint32_t& vertex, std::string& error)
	{
		if (!finite(sample.uv))
		{
			error = "constraint sample has a non-finite UV coordinate";
			return false;
		}
		if (sample.topology_id == UvVertex::no_topology_id)
		{
			error = "constraint sample has no canonical topology ID";
			return false;
		}

		const double tolerance2 = options_.coordinate_tolerance * options_.coordinate_tolerance;
		std::optional<std::size_t> matching_vertex;
		for (std::size_t candidate = 0; candidate < vertices_.size(); ++candidate)
		{
			const double distance = squared_distance(sample.uv, vertices_[candidate].uv);
			if (distance <= tolerance2)
			{
				if (matching_vertex)
				{
					error = "sample lies within the UV tolerance of multiple mesh vertices";
					return false;
				}
				matching_vertex = candidate;
			}
		}
		if (matching_vertex)
		{
			vertex = static_cast<std::uint32_t>(*matching_vertex);
			return claim_topology_id(vertex, sample.topology_id, error);
		}

		std::map<Edge, unsigned, EdgeLess> incidence;
		for (const UvTriangle& triangle : triangles_)
			for (unsigned i = 0; i < 3; ++i)
				++incidence[edge_key(triangle.vertices[i], triangle.vertices[(i + 1) % 3])];
		std::optional<Edge> containing_edge;
		for (const auto& [edge, count] : incidence)
		{
			(void)count;
			if (!point_on_segment(sample.uv, vertices_[edge[0]].uv, vertices_[edge[1]].uv,
				options_.coordinate_tolerance, false)) continue;
			if (containing_edge)
			{
				error = "sample lies on multiple distinct mesh edges; input is not a planar triangulation";
				return false;
			}
			containing_edge = edge;
		}
		if (containing_edge) return split_edge(*containing_edge, sample, vertex, error);

		for (std::size_t triangle = 0; triangle < triangles_.size(); ++triangle)
		{
			const auto& indices = triangles_[triangle].vertices;
			if (point_in_triangle(sample.uv, vertices_[indices[0]].uv, vertices_[indices[1]].uv,
				vertices_[indices[2]].uv, options_.coordinate_tolerance, true))
				return split_triangle(triangle, sample, vertex, error);
		}
		error = "constraint sample is outside the triangulated face or inside one of its holes";
		return false;
	}

	bool UvConstraintMesh::triangulate_polygon(const std::vector<std::uint32_t>& polygon,
		std::vector<UvTriangle>& output, std::string& error) const
	{
		if (polygon.size() < 3)
		{
			error = "constraint cavity side has fewer than three vertices";
			return false;
		}
		std::vector<std::uint32_t> remaining = polygon;
		long double polygon_area2 = 0.0L;
		for (std::size_t i = 0; i < remaining.size(); ++i)
		{
			const UvPoint& a = vertices_[remaining[i]].uv;
			const UvPoint& b = vertices_[remaining[(i + 1) % remaining.size()]].uv;
			polygon_area2 += static_cast<long double>(a.u) * b.v
				- static_cast<long double>(a.v) * b.u;
		}
		if (!(polygon_area2 > 0.0L))
		{
			error = "constraint cavity side is not a positive-area CCW polygon";
			return false;
		}

		while (remaining.size() > 3)
		{
			std::optional<std::size_t> selected;
			std::uint32_t selected_vertex = std::numeric_limits<std::uint32_t>::max();
			for (std::size_t i = 0; i < remaining.size(); ++i)
			{
				const std::uint32_t previous = remaining[(i + remaining.size() - 1) % remaining.size()];
				const std::uint32_t current = remaining[i];
				const std::uint32_t next = remaining[(i + 1) % remaining.size()];
				const UvPoint& a = vertices_[previous].uv;
				const UvPoint& b = vertices_[current].uv;
				const UvPoint& c = vertices_[next].uv;
				if (!(orient(a, b, c) > 0.0L)) continue;

				bool blocked = false;
				for (std::uint32_t candidate : remaining)
				{
					if (candidate == previous || candidate == current || candidate == next) continue;
					if (point_in_triangle(vertices_[candidate].uv, a, b, c,
						options_.coordinate_tolerance, false))
					{
						blocked = true;
						break;
					}
				}
				if (!blocked && current < selected_vertex)
				{
					selected = i;
					selected_vertex = current;
				}
			}
			if (!selected)
			{
				error = "deterministic ear clipping found no valid constraint-cavity ear";
				return false;
			}
			const std::size_t i = *selected;
			output.push_back({{remaining[(i + remaining.size() - 1) % remaining.size()],
				remaining[i], remaining[(i + 1) % remaining.size()]}});
			remaining.erase(remaining.begin() + static_cast<std::ptrdiff_t>(i));
		}
		if (!(orient(vertices_[remaining[0]].uv, vertices_[remaining[1]].uv,
			vertices_[remaining[2]].uv) > 0.0L))
		{
			error = "constraint cavity ended with a degenerate triangle";
			return false;
		}
		output.push_back({{remaining[0], remaining[1], remaining[2]}});
		return true;
	}

	bool UvConstraintMesh::insert_atomic_segment(std::uint32_t a, std::uint32_t b,
		std::uint64_t constraint_id, std::string& error)
	{
		if (a == b)
		{
			error = "constraint contains a zero-length segment";
			return false;
		}
		const Edge requested = edge_key(a, b);
		if (has_edge(a, b))
		{
			add_constraint_tag(requested, constraint_id);
			return true;
		}

		const UvPoint& start = vertices_[a].uv;
		const UvPoint& finish = vertices_[b].uv;
		for (std::size_t vertex = 0; vertex < vertices_.size(); ++vertex)
			if (vertex != static_cast<std::size_t>(a)
				&& vertex != static_cast<std::size_t>(b) && point_on_segment(vertices_[vertex].uv,
				start, finish, options_.coordinate_tolerance, false))
			{
				error = "constraint segment passes through unsampled UV vertex "
					+ std::to_string(vertex) + "; include it in the canonical sample chain";
				return false;
			}

		struct EdgeInfo
		{
			unsigned incidence = 0;
		};
		std::map<Edge, EdgeInfo, EdgeLess> edge_info;
		for (const UvTriangle& triangle : triangles_)
			for (unsigned i = 0; i < 3; ++i)
				++edge_info[edge_key(triangle.vertices[i], triangle.vertices[(i + 1) % 3])].incidence;

		struct Crossing
		{
			double parameter = 0.0;
			Edge edge{};
		};
		std::vector<Crossing> crossings;
		for (const auto& [edge, info] : edge_info)
		{
			double parameter = 0.0;
			if (!proper_intersection(start, finish, vertices_[edge[0]].uv,
				vertices_[edge[1]].uv, options_.coordinate_tolerance, parameter)) continue;
			if (info.incidence != 2)
			{
				error = "constraint segment crosses face boundary edge " + edge_text(edge)
					+ " and would leave the trimmed domain";
				return false;
			}
			if (constraint_tags_.contains(edge))
			{
				error = "constraint segment crosses an existing constraint at edge " + edge_text(edge)
					+ "; provide their intersection as a canonical sample first";
				return false;
			}
			crossings.push_back({parameter, edge});
		}
		std::sort(crossings.begin(), crossings.end(), [](const Crossing& lhs, const Crossing& rhs)
		{
			return lhs.parameter < rhs.parameter
				|| (lhs.parameter == rhs.parameter && lhs.edge < rhs.edge);
		});
		const double crossing_parameter_tolerance = segment_parameter_tolerance(start, finish,
			options_.coordinate_tolerance);
		for (std::size_t i = 1; i < crossings.size(); ++i)
			if (std::abs(crossings[i].parameter - crossings[i - 1].parameter)
				<= crossing_parameter_tolerance)
			{
				error = "constraint trace meets more than one mesh edge at the same parameter";
				return false;
			}

		std::vector<double> parameters{0.0};
		for (const Crossing& crossing : crossings) parameters.push_back(crossing.parameter);
		parameters.push_back(1.0);
		std::vector<std::size_t> cavity;
		for (std::size_t interval = 0; interval + 1 < parameters.size(); ++interval)
		{
			const double parameter = 0.5 * (parameters[interval] + parameters[interval + 1]);
			const UvPoint point{start.u + parameter * (finish.u - start.u),
				start.v + parameter * (finish.v - start.v)};
			std::optional<std::size_t> containing;
			for (std::size_t triangle = 0; triangle < triangles_.size(); ++triangle)
			{
				const auto& indices = triangles_[triangle].vertices;
				if (!point_in_triangle(point, vertices_[indices[0]].uv, vertices_[indices[1]].uv,
					vertices_[indices[2]].uv, options_.coordinate_tolerance, true)) continue;
				if (containing)
				{
					error = "constraint trace lies in overlapping input triangles";
					return false;
				}
				containing = triangle;
			}
			if (!containing)
			{
				error = "constraint trace leaves the trimmed face or enters a hole";
				return false;
			}
			if (cavity.empty() || cavity.back() != *containing) cavity.push_back(*containing);
		}
		if (cavity.empty())
		{
			error = "constraint trace did not cross any triangle interior";
			return false;
		}
		std::set<std::size_t> cavity_set(cavity.begin(), cavity.end());
		if (cavity_set.size() != cavity.size())
		{
			error = "constraint triangle strip is self-intersecting";
			return false;
		}

		struct DirectedBoundary
		{
			unsigned count = 0;
			std::uint32_t from = 0;
			std::uint32_t to = 0;
		};
		std::map<Edge, DirectedBoundary, EdgeLess> strip_edges;
		for (std::size_t triangle : cavity)
		{
			const auto& indices = triangles_[triangle].vertices;
			for (unsigned i = 0; i < 3; ++i)
			{
				const std::uint32_t from = indices[i];
				const std::uint32_t to = indices[(i + 1) % 3];
				auto& record = strip_edges[edge_key(from, to)];
				++record.count;
				if (record.count == 1)
				{
					record.from = from;
					record.to = to;
				}
			}
		}
		std::map<std::uint32_t, std::uint32_t> successor;
		std::size_t boundary_edge_count = 0;
		for (const auto& [edge, record] : strip_edges)
		{
			(void)edge;
			if (record.count == 2) continue;
			if (record.count != 1 || successor.contains(record.from))
			{
				error = "constraint triangle strip is not a two-manifold disk";
				return false;
			}
			successor[record.from] = record.to;
			++boundary_edge_count;
		}
		std::vector<std::uint32_t> cycle;
		std::uint32_t current = a;
		for (std::size_t step = 0; step < boundary_edge_count; ++step)
		{
			cycle.push_back(current);
			const auto next = successor.find(current);
			if (next == successor.end())
			{
				error = "constraint cavity boundary does not contain the segment start";
				return false;
			}
			current = next->second;
			if (current == a && step + 1 != boundary_edge_count)
			{
				error = "constraint cavity boundary contains multiple loops";
				return false;
			}
		}
		if (current != a || cycle.size() != boundary_edge_count)
		{
			error = "constraint cavity boundary is not one closed loop";
			return false;
		}
		const auto b_position = std::find(cycle.begin(), cycle.end(), b);
		if (b_position == cycle.end())
		{
			error = "constraint cavity boundary does not contain the segment finish";
			return false;
		}
		const std::size_t split = static_cast<std::size_t>(b_position - cycle.begin());
		std::vector<std::uint32_t> side_a(cycle.begin(), cycle.begin() + split + 1);
		std::vector<std::uint32_t> side_b(cycle.begin() + split, cycle.end());
		side_b.push_back(a);

		std::vector<UvTriangle> replacement;
		if (!triangulate_polygon(side_a, replacement, error)
			|| !triangulate_polygon(side_b, replacement, error)) return false;

		std::vector<UvTriangle> retained;
		retained.reserve(triangles_.size() - cavity.size() + replacement.size());
		for (std::size_t triangle = 0; triangle < triangles_.size(); ++triangle)
			if (!cavity_set.contains(triangle)) retained.push_back(triangles_[triangle]);
		retained.insert(retained.end(), replacement.begin(), replacement.end());
		triangles_ = std::move(retained);
		add_constraint_tag(requested, constraint_id);
		return true;
	}

	bool UvConstraintMesh::insert_polyline(std::span<const UvConstraintSample> samples,
		std::uint64_t constraint_id, std::string* error)
	{
		auto fail = [&](const std::string& message)
		{
			if (error) *error = message;
			return false;
		};
		if (samples.size() < 2) return fail("a constraint polyline needs at least two samples");
		if (samples.size() > maximum_uv_vertex_count)
			return fail("constraint sample chain exceeds uint32 vertex index capacity");
		std::string validation_error;
		if (!validate(&validation_error)) return fail("invalid input UV mesh: " + validation_error);

		UvConstraintMesh working = *this;
		std::vector<std::uint32_t> sample_vertices;
		sample_vertices.reserve(samples.size());
		for (const UvConstraintSample& sample : samples)
		{
			std::uint32_t vertex = 0;
			if (!working.insert_sample(sample, vertex, validation_error)) return fail(validation_error);
			sample_vertices.push_back(vertex);
		}
		for (std::size_t i = 1; i < sample_vertices.size(); ++i)
			if (!working.insert_atomic_segment(sample_vertices[i - 1], sample_vertices[i],
				constraint_id, validation_error)) return fail(validation_error);
		if (!working.validate(&validation_error))
			return fail("constraint insertion produced an invalid UV mesh: " + validation_error);
		*this = std::move(working);
		if (error) error->clear();
		return true;
	}

	bool UvConstraintMesh::assign_unclaimed_topology_ids(std::uint32_t& next_id,
		std::string* error)
	{
		if (vertices_.size() > maximum_uv_vertex_count)
		{
			if (error) *error = "UV mesh exceeds uint32 vertex index capacity";
			return false;
		}
		// Keep allocator exhaustion transactional just like constraint insertion: neither
		// the caller's cursor nor a prefix of this face may change on failure.
		UvConstraintMesh working = *this;
		const std::uint64_t sentinel = UvVertex::no_topology_id;
		std::uint64_t candidate_id = next_id;
		std::unordered_set<std::uint32_t> occupied;
		for (const UvVertex& vertex : working.vertices_)
			if (vertex.topology_id != UvVertex::no_topology_id)
				occupied.insert(vertex.topology_id);
		for (UvVertex& vertex : working.vertices_)
		{
			if (vertex.topology_id != UvVertex::no_topology_id) continue;
			while (candidate_id < sentinel
				&& occupied.contains(static_cast<std::uint32_t>(candidate_id))) ++candidate_id;
			if (candidate_id >= sentinel)
			{
				if (error) *error = "canonical topology ID allocator exhausted";
				return false;
			}
			vertex.topology_id = static_cast<std::uint32_t>(candidate_id);
			occupied.insert(vertex.topology_id);
			++candidate_id;
		}
		*this = std::move(working);
		next_id = static_cast<std::uint32_t>(candidate_id);
		if (error) error->clear();
		return true;
	}

	std::vector<UvEdgeView> UvConstraintMesh::edges() const
	{
		std::map<Edge, unsigned, EdgeLess> incidence;
		for (const UvTriangle& triangle : triangles_)
			for (unsigned i = 0; i < 3; ++i)
				++incidence[edge_key(triangle.vertices[i], triangle.vertices[(i + 1) % 3])];
		std::vector<UvEdgeView> result;
		result.reserve(incidence.size());
		for (const auto& [edge, count] : incidence)
		{
			UvEdgeView view;
			view.vertices = edge;
			view.boundary = count == 1;
			const auto tag = constraint_tags_.find(edge);
			if (tag != constraint_tags_.end()) view.constraint_ids = tag->second;
			result.push_back(std::move(view));
		}
		return result;
	}

	bool UvConstraintMesh::has_edge(std::uint32_t a, std::uint32_t b) const
	{
		const Edge target = edge_key(a, b);
		for (const UvTriangle& triangle : triangles_)
			for (unsigned i = 0; i < 3; ++i)
				if (edge_key(triangle.vertices[i], triangle.vertices[(i + 1) % 3]) == target)
					return true;
		return false;
	}

	bool UvConstraintMesh::edge_has_constraint(std::uint32_t a, std::uint32_t b,
		std::uint64_t constraint_id) const
	{
		const auto found = constraint_tags_.find(edge_key(a, b));
		return found != constraint_tags_.end()
			&& std::binary_search(found->second.begin(), found->second.end(), constraint_id);
	}

	std::optional<std::uint32_t> UvConstraintMesh::vertex_for_topology_id(
		std::uint32_t topology_id) const
	{
		for (std::size_t vertex = 0; vertex < vertices_.size(); ++vertex)
			if (vertices_[vertex].topology_id == topology_id)
			{
				if (vertex >= maximum_uv_vertex_count) return std::nullopt;
				return static_cast<std::uint32_t>(vertex);
			}
		return std::nullopt;
	}

	double UvConstraintMesh::area() const
	{
		long double sum = 0.0L;
		for (const UvTriangle& triangle : triangles_)
			sum += orient(vertices_[triangle.vertices[0]].uv, vertices_[triangle.vertices[1]].uv,
				vertices_[triangle.vertices[2]].uv);
		return static_cast<double>(0.5L * sum);
	}

	bool UvConstraintMesh::contains(const UvPoint& point, bool include_boundary) const
	{
		for (const UvTriangle& triangle : triangles_)
			if (point_in_triangle(point, vertices_[triangle.vertices[0]].uv,
				vertices_[triangle.vertices[1]].uv, vertices_[triangle.vertices[2]].uv,
				options_.coordinate_tolerance, !include_boundary)) return true;
		return false;
	}

	std::size_t UvConstraintMesh::boundary_loop_count() const
	{
		std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> adjacency;
		for (const UvEdgeView& edge : edges())
			if (edge.boundary)
			{
				adjacency[edge.vertices[0]].push_back(edge.vertices[1]);
				adjacency[edge.vertices[1]].push_back(edge.vertices[0]);
			}
		std::unordered_set<std::uint32_t> visited;
		std::size_t components = 0;
		for (const auto& [start, neighbours] : adjacency)
		{
			(void)neighbours;
			if (visited.contains(start)) continue;
			++components;
			std::vector<std::uint32_t> stack{start};
			visited.insert(start);
			while (!stack.empty())
			{
				const std::uint32_t vertex = stack.back();
				stack.pop_back();
				for (std::uint32_t neighbour : adjacency[vertex])
					if (visited.insert(neighbour).second) stack.push_back(neighbour);
			}
		}
		return components;
	}

	bool UvConstraintMesh::validate(std::string* error) const
	{
		auto fail = [&](const std::string& message)
		{
			if (error) *error = message;
			return false;
		};
		if (!(options_.coordinate_tolerance >= 0.0)
			|| !std::isfinite(options_.coordinate_tolerance))
			return fail("coordinate tolerance is not finite and non-negative");
		if (vertices_.size() > maximum_uv_vertex_count)
			return fail("UV mesh exceeds uint32 vertex index capacity");
		for (std::size_t vertex = 0; vertex < vertices_.size(); ++vertex)
			if (!finite(vertices_[vertex].uv))
				return fail("vertex " + std::to_string(vertex) + " has non-finite UV coordinates");
		std::unordered_map<std::uint32_t, std::size_t> topology_owner;
		for (std::size_t vertex = 0; vertex < vertices_.size(); ++vertex)
		{
			const std::uint32_t topology = vertices_[vertex].topology_id;
			if (topology == UvVertex::no_topology_id) continue;
			if (!topology_owner.emplace(topology, vertex).second)
				return fail("canonical topology ID " + std::to_string(topology)
					+ " occurs on more than one UV vertex");
		}

		struct EdgeIncidence
		{
			unsigned count = 0;
			int directed_sum = 0;
		};
		std::map<Edge, EdgeIncidence, EdgeLess> incidence;
		for (std::size_t triangle = 0; triangle < triangles_.size(); ++triangle)
		{
			const auto& indices = triangles_[triangle].vertices;
			if (indices[0] >= vertices_.size() || indices[1] >= vertices_.size()
				|| indices[2] >= vertices_.size())
				return fail("triangle " + std::to_string(triangle) + " has an invalid vertex index");
			if (indices[0] == indices[1] || indices[1] == indices[2] || indices[2] == indices[0])
				return fail("triangle " + std::to_string(triangle) + " repeats a vertex");
			if (!(orient(vertices_[indices[0]].uv, vertices_[indices[1]].uv,
				vertices_[indices[2]].uv) > 0.0L))
				return fail("triangle " + std::to_string(triangle) + " is not strictly CCW");
			for (unsigned i = 0; i < 3; ++i)
			{
				const std::uint32_t from = indices[i];
				const std::uint32_t to = indices[(i + 1) % 3];
				const Edge edge = edge_key(from, to);
				EdgeIncidence& record = incidence[edge];
				++record.count;
				record.directed_sum += from < to ? 1 : -1;
				if (record.count > 2)
					return fail("edge " + edge_text(edge) + " is non-manifold");
			}
		}
		for (const auto& [edge, record] : incidence)
			if (record.count == 2 && record.directed_sum != 0)
				return fail("edge " + edge_text(edge)
					+ " has inconsistent adjacent triangle winding");

		std::vector<Edge> geometric_edges;
		geometric_edges.reserve(incidence.size());
		for (const auto& [edge, record] : incidence)
		{
			(void)record;
			geometric_edges.push_back(edge);
		}
		for (std::size_t first = 0; first < geometric_edges.size(); ++first)
			for (std::size_t second = first + 1; second < geometric_edges.size(); ++second)
				if (nonconforming_edge_intersection(geometric_edges[first],
					geometric_edges[second], vertices_))
					return fail("geometric mesh edges " + edge_text(geometric_edges[first])
						+ " and " + edge_text(geometric_edges[second])
						+ " cross, overlap, or meet without conforming vertex incidence");

		for (std::size_t first = 0; first < triangles_.size(); ++first)
			for (std::size_t second = first + 1; second < triangles_.size(); ++second)
				if (triangles_have_positive_area_overlap(triangles_[first], triangles_[second],
					vertices_))
					return fail("triangles " + std::to_string(first) + " and "
						+ std::to_string(second) + " have positive-area overlap");
		std::unordered_map<std::uint32_t, unsigned> boundary_degree;
		for (const auto& [edge, record] : incidence)
			if (record.count == 1)
			{
				++boundary_degree[edge[0]];
				++boundary_degree[edge[1]];
			}
		for (const auto& [vertex, degree] : boundary_degree)
			if (degree != 2)
				return fail("boundary vertex " + std::to_string(vertex)
					+ " does not belong to one closed boundary loop");
		for (const auto& [edge, tags] : constraint_tags_)
		{
			if (!incidence.contains(edge)) return fail("constraint tag refers to missing edge " + edge_text(edge));
			if (tags.empty() || !std::is_sorted(tags.begin(), tags.end())
				|| std::adjacent_find(tags.begin(), tags.end()) != tags.end())
				return fail("constraint edge tags are empty, unsorted, or duplicated");
		}
		if (error) error->clear();
		return true;
	}
}
