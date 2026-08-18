// Diagnostic only: ask OCCT to imprint the imported BRep faces before tessellation.
//
// This tool deliberately does not call ParaCFD's triangle/embedded-boundary pipeline.
// It compares native BRep topology before and after General Fuse, and reports sewing
// and ShapeFix as separate (non-production) gap-repair experiments.

#include <BOPAlgo_Splitter.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAbs_CurveType.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Geom_Surface.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IntTools_CommonPrt.hxx>
#include <IntTools_EdgeFace.hxx>
#include <IntTools_Range.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <ShapeFix_Shape.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
	using Clock = std::chrono::steady_clock;
	using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

	struct Options
	{
		std::string step = "Test-Data/PlanBParakite.step";
		double fuzzy_mm = 0.01;
		double sew_mm = 0.01;
		double deflection_mm = 2.0;
		double gap_search_mm = 5.0;
		unsigned gap_samples = 9;
		std::size_t gap_max_edges = 32;
		std::size_t subset_faces = 0;
		bool run_split = false;
		bool run_sewing = false;
		bool run_shape_fix = false;
		bool split_after_sew = false;
		bool verbose_gaps = false;
	};

	void usage(const char* exe)
	{
		std::printf(
			"Usage: %s [--step FILE] [--fuzzy-mm T] [--sew-mm T] "
			"[--deflection-mm D] [--gap-search-mm D] [--gap-samples N] [--gap-max-edges N] "
			"[--subset N] [--run-split] [--run-sew] [--run-shape-fix] "
			"[--split-after-sew] [--verbose-gaps]\n", exe);
	}

	bool parse_options(int argc, char** argv, Options& options)
	{
		for (int i = 1; i < argc; ++i)
		{
			const std::string arg = argv[i];
			auto value = [&](const char* name) -> const char*
			{
				if (i + 1 >= argc)
				{
					std::fprintf(stderr, "%s requires a value\n", name);
					return nullptr;
				}
				return argv[++i];
			};
			if (arg == "--step")
			{
				const char* v = value("--step"); if (!v) return false; options.step = v;
			}
			else if (arg == "--fuzzy-mm")
			{
				const char* v = value("--fuzzy-mm"); if (!v) return false;
				options.fuzzy_mm = std::strtod(v, nullptr);
			}
			else if (arg == "--sew-mm")
			{
				const char* v = value("--sew-mm"); if (!v) return false;
				options.sew_mm = std::strtod(v, nullptr);
			}
			else if (arg == "--deflection-mm")
			{
				const char* v = value("--deflection-mm"); if (!v) return false;
				options.deflection_mm = std::strtod(v, nullptr);
			}
			else if (arg == "--gap-search-mm")
			{
				const char* v = value("--gap-search-mm"); if (!v) return false;
				options.gap_search_mm = std::strtod(v, nullptr);
			}
			else if (arg == "--gap-samples")
			{
				const char* v = value("--gap-samples"); if (!v) return false;
				options.gap_samples = static_cast<unsigned>(std::strtoul(v, nullptr, 10));
			}
			else if (arg == "--gap-max-edges")
			{
				const char* v = value("--gap-max-edges"); if (!v) return false;
				options.gap_max_edges = static_cast<std::size_t>(std::strtoull(v, nullptr, 10));
			}
			else if (arg == "--subset")
			{
				const char* v = value("--subset"); if (!v) return false;
				options.subset_faces = static_cast<std::size_t>(std::strtoull(v, nullptr, 10));
			}
			else if (arg == "--run-split") options.run_split = true;
			else if (arg == "--run-sew") options.run_sewing = true;
			else if (arg == "--run-shape-fix") options.run_shape_fix = true;
			else if (arg == "--verbose-gaps") options.verbose_gaps = true;
			else if (arg == "--split-after-sew")
			{
				options.split_after_sew = true;
				options.run_split = true;
				options.run_sewing = true;
			}
			else if (arg == "--help" || arg == "-h") { usage(argv[0]); std::exit(0); }
			else
			{
				std::fprintf(stderr, "unknown option: %s\n", arg.c_str());
				return false;
			}
		}
		return options.fuzzy_mm >= 0.0 && options.sew_mm >= 0.0
			&& options.deflection_mm > 0.0 && options.gap_search_mm > 0.0
			&& options.gap_samples >= 3;
	}

	TopoDS_Shape read_step(const std::string& path)
	{
		STEPControl_Reader reader;
		if (reader.ReadFile(path.c_str()) != IFSelect_RetDone)
			throw std::runtime_error("STEPControl_Reader::ReadFile failed");
		if (reader.TransferRoots() <= 0)
			throw std::runtime_error("STEPControl_Reader::TransferRoots transferred no roots");
		TopoDS_Shape shape = reader.OneShape();
		if (shape.IsNull()) throw std::runtime_error("STEP reader returned a null shape");
		return shape;
	}

	std::vector<TopoDS_Face> unique_faces(const TopoDS_Shape& shape)
	{
		ShapeMap map;
		TopExp::MapShapes(shape, TopAbs_FACE, map);
		std::vector<TopoDS_Face> result;
		result.reserve(static_cast<std::size_t>(map.Extent()));
		for (Standard_Integer i = 1; i <= map.Extent(); ++i)
			result.push_back(TopoDS::Face(map(i)));
		return result;
	}

	TopoDS_Shape face_compound(const std::vector<TopoDS_Face>& faces)
	{
		TopoDS_Compound compound;
		BRep_Builder builder;
		builder.MakeCompound(compound);
		for (const TopoDS_Face& face : faces) builder.Add(compound, face);
		return compound;
	}

	std::array<double, 6> bounds(const TopoDS_Shape& shape)
	{
		Bnd_Box box;
		BRepBndLib::AddOptimal(shape, box, false, true);
		std::array<double, 6> result{};
		if (box.IsVoid())
		{
			result.fill(std::numeric_limits<double>::quiet_NaN());
			return result;
		}
		box.Get(result[0], result[1], result[2], result[3], result[4], result[5]);
		return result;
	}

	std::vector<TopoDS_Face> central_subset(const std::vector<TopoDS_Face>& all,
		std::size_t requested)
	{
		if (requested == 0 || requested >= all.size()) return all;
		const TopoDS_Shape all_shape = face_compound(all);
		const auto all_box = bounds(all_shape);
		const std::array<double, 3> centre{
			0.5 * (all_box[0] + all_box[3]), 0.5 * (all_box[1] + all_box[4]),
			0.5 * (all_box[2] + all_box[5]) };
		std::vector<std::pair<double, std::size_t>> order;
		order.reserve(all.size());
		for (std::size_t i = 0; i < all.size(); ++i)
		{
			const auto box = bounds(all[i]);
			const double dx = 0.5 * (box[0] + box[3]) - centre[0];
			const double dy = 0.5 * (box[1] + box[4]) - centre[1];
			const double dz = 0.5 * (box[2] + box[5]) - centre[2];
			order.emplace_back(dx * dx + dy * dy + dz * dz, i);
		}
		std::nth_element(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(requested),
			order.end());
		order.resize(requested);
		std::sort(order.begin(), order.end(), [](const auto& a, const auto& b)
			{ return a.second < b.second; });
		std::vector<TopoDS_Face> subset;
		subset.reserve(requested);
		for (const auto& [distance, index] : order)
		{
			(void)distance;
			subset.push_back(all[index]);
		}
		return subset;
	}

	struct DisjointSet
	{
		explicit DisjointSet(std::size_t size) : parent(size), rank(size, 0)
		{
			std::iota(parent.begin(), parent.end(), 0);
		}
		std::size_t find(std::size_t x)
		{
			while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
			return x;
		}
		void unite(std::size_t a, std::size_t b)
		{
			a = find(a); b = find(b); if (a == b) return;
			if (rank[a] < rank[b]) std::swap(a, b);
			parent[b] = a;
			if (rank[a] == rank[b]) ++rank[a];
		}
		std::vector<std::size_t> parent;
		std::vector<unsigned char> rank;
	};

	struct BoundaryComponent
	{
		double length_mm = 0.0;
		std::size_t edges = 0;
		Bnd_Box box;
	};

	struct TopologyStats
	{
		std::size_t faces = 0;
		std::size_t edges = 0;
		std::size_t vertices = 0;
		std::size_t free_edges = 0;
		std::size_t manifold_edges = 0;
		std::size_t junction_edges = 0;
		std::size_t unused_edges = 0;
		std::size_t degenerate_edges = 0;
		double min_edge_tolerance_mm = std::numeric_limits<double>::infinity();
		double max_edge_tolerance_mm = 0.0;
		double free_length_mm = 0.0;
		std::vector<BoundaryComponent> free_components;
	};

	TopologyStats topology_stats(const TopoDS_Shape& shape)
	{
		TopologyStats stats;
		ShapeMap faces, edges, vertices;
		TopExp::MapShapes(shape, TopAbs_FACE, faces);
		TopExp::MapShapes(shape, TopAbs_EDGE, edges);
		TopExp::MapShapes(shape, TopAbs_VERTEX, vertices);
		stats.faces = static_cast<std::size_t>(faces.Extent());
		stats.edges = static_cast<std::size_t>(edges.Extent());
		stats.vertices = static_cast<std::size_t>(vertices.Extent());
		std::vector<std::size_t> uses(stats.edges, 0);
		for (Standard_Integer face_id = 1; face_id <= faces.Extent(); ++face_id)
			for (TopExp_Explorer edge_exp(faces(face_id), TopAbs_EDGE); edge_exp.More(); edge_exp.Next())
			{
				const Standard_Integer edge_id = edges.FindIndex(edge_exp.Current());
				if (edge_id > 0) ++uses[static_cast<std::size_t>(edge_id - 1)];
			}

		std::vector<std::size_t> free_ids;
		for (Standard_Integer edge_id = 1; edge_id <= edges.Extent(); ++edge_id)
		{
			const TopoDS_Edge edge = TopoDS::Edge(edges(edge_id));
			const double tolerance = BRep_Tool::Tolerance(edge);
			stats.min_edge_tolerance_mm = std::min(stats.min_edge_tolerance_mm, tolerance);
			stats.max_edge_tolerance_mm = std::max(stats.max_edge_tolerance_mm, tolerance);
			if (BRep_Tool::Degenerated(edge)) ++stats.degenerate_edges;
			const std::size_t use_count = uses[static_cast<std::size_t>(edge_id - 1)];
			if (use_count == 0) ++stats.unused_edges;
			else if (use_count == 1)
			{
				++stats.free_edges;
				free_ids.push_back(static_cast<std::size_t>(edge_id - 1));
			}
			else if (use_count == 2) ++stats.manifold_edges;
			else ++stats.junction_edges;
		}
		if (!std::isfinite(stats.min_edge_tolerance_mm)) stats.min_edge_tolerance_mm = 0.0;

		DisjointSet dsu(static_cast<std::size_t>(std::max(1, vertices.Extent())));
		for (std::size_t edge_index : free_ids)
		{
			TopoDS_Vertex first, last;
			TopExp::Vertices(TopoDS::Edge(edges(static_cast<Standard_Integer>(edge_index + 1))),
				first, last, false);
			const Standard_Integer a = first.IsNull() ? 0 : vertices.FindIndex(first);
			const Standard_Integer b = last.IsNull() ? 0 : vertices.FindIndex(last);
			if (a > 0 && b > 0) dsu.unite(static_cast<std::size_t>(a - 1),
				static_cast<std::size_t>(b - 1));
		}
		std::map<std::size_t, BoundaryComponent> components;
		std::size_t orphan = static_cast<std::size_t>(vertices.Extent());
		for (std::size_t edge_index : free_ids)
		{
			const TopoDS_Edge edge = TopoDS::Edge(edges(static_cast<Standard_Integer>(edge_index + 1)));
			TopoDS_Vertex first, last;
			TopExp::Vertices(edge, first, last, false);
			const Standard_Integer vertex_id = !first.IsNull() ? vertices.FindIndex(first)
				: (!last.IsNull() ? vertices.FindIndex(last) : 0);
			const std::size_t root = vertex_id > 0
				? dsu.find(static_cast<std::size_t>(vertex_id - 1)) : orphan++;
			BoundaryComponent& component = components[root];
			GProp_GProps properties;
			BRepGProp::LinearProperties(edge, properties);
			component.length_mm += properties.Mass();
			++component.edges;
			BRepBndLib::Add(edge, component.box);
		}
		for (auto& [root, component] : components)
		{
			(void)root;
			stats.free_length_mm += component.length_mm;
			stats.free_components.push_back(std::move(component));
		}
		std::sort(stats.free_components.begin(), stats.free_components.end(),
			[](const BoundaryComponent& a, const BoundaryComponent& b)
			{ return a.length_mm > b.length_mm; });
		return stats;
	}

	void print_topology(const char* label, const TopoDS_Shape& shape, unsigned component_limit = 8)
	{
		const TopologyStats stats = topology_stats(shape);
		const auto box = bounds(shape);
		std::printf("\n[%s]\n", label);
		std::printf("  valid BRep: %s\n", BRepCheck_Analyzer(shape).IsValid() ? "yes" : "NO");
		std::printf("  faces/edges/vertices: %zu / %zu / %zu\n",
			stats.faces, stats.edges, stats.vertices);
		std::printf("  edge face-use fan free/manifold/junction/unused: %zu / %zu / %zu / %zu\n",
			stats.free_edges, stats.manifold_edges, stats.junction_edges, stats.unused_edges);
		std::printf("  degenerate edges: %zu; edge tolerance min/max: %.9g / %.9g mm\n",
			stats.degenerate_edges, stats.min_edge_tolerance_mm, stats.max_edge_tolerance_mm);
		std::printf("  bbox mm: [%.9g %.9g %.9g] -> [%.9g %.9g %.9g]\n",
			box[0], box[1], box[2], box[3], box[4], box[5]);
		std::printf("  free boundary: %zu components, %zu edges, %.9g mm total\n",
			stats.free_components.size(), stats.free_edges, stats.free_length_mm);
		for (std::size_t i = 0; i < stats.free_components.size() && i < component_limit; ++i)
		{
			const BoundaryComponent& component = stats.free_components[i];
			std::array<double, 6> b{};
			if (!component.box.IsVoid()) component.box.Get(b[0], b[1], b[2], b[3], b[4], b[5]);
			std::printf("    #%zu edges=%zu length=%.9g bbox=[%.6g %.6g %.6g]-[%.6g %.6g %.6g]\n",
				i + 1, component.edges, component.length_mm,
				b[0], b[1], b[2], b[3], b[4], b[5]);
		}
	}

	bool exact_full_edge_on_face(const TopoDS_Edge& edge, const TopoDS_Face& face)
	{
		BRepAdaptor_Curve curve(edge);
		const double first = curve.FirstParameter();
		const double last = curve.LastParameter();
		if (!std::isfinite(first) || !std::isfinite(last) || !(last > first)) return false;
		IntTools_EdgeFace intersection;
		intersection.SetEdge(edge);
		intersection.SetFace(face);
		intersection.SetRange(first, last);
		intersection.SetFuzzyValue(0.0);
		intersection.Perform();
		if (!intersection.IsDone()) return false;
		std::vector<std::pair<double, double>> common_ranges;
		for (NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
			intersection.CommonParts()); common.More(); common.Next())
		{
			if (common.Value().Type() != TopAbs_EDGE) continue;
			double lo = 0.0, hi = 0.0;
			common.Value().Range1(lo, hi);
			if (hi < lo) std::swap(lo, hi);
			lo = std::max(lo, first);
			hi = std::min(hi, last);
			if (hi > lo) common_ranges.emplace_back(lo, hi);
		}
		if (common_ranges.empty()) return false;
		std::sort(common_ranges.begin(), common_ranges.end());
		const double parameter_tolerance = std::max(Precision::PConfusion(),
			64.0 * std::numeric_limits<double>::epsilon()
			* std::max({ 1.0, std::abs(first), std::abs(last) }));
		double covered = first;
		for (const auto& range : common_ranges)
		{
			if (range.first > covered + parameter_tolerance) return false;
			covered = std::max(covered, range.second);
		}
		return covered >= last - parameter_tolerance;
	}

	std::string curve_description(const BRepAdaptor_Curve& curve)
	{
		const char* name = "Other";
		switch (curve.GetType())
		{
		case GeomAbs_Line: name = "Line"; break;
		case GeomAbs_Circle: name = "Circle"; break;
		case GeomAbs_Ellipse: name = "Ellipse"; break;
		case GeomAbs_Hyperbola: name = "Hyperbola"; break;
		case GeomAbs_Parabola: name = "Parabola"; break;
		case GeomAbs_BezierCurve: name = "Bezier"; break;
		case GeomAbs_BSplineCurve: name = "BSpline"; break;
		case GeomAbs_OffsetCurve: name = "Offset"; break;
		default: break;
		}
		std::ostringstream result;
		result << name;
		if (curve.GetType() == GeomAbs_BezierCurve || curve.GetType() == GeomAbs_BSplineCurve)
			result << " degree=" << curve.Degree();
		return result.str();
	}

	std::string surface_description(const TopoDS_Face& face)
	{
		BRepAdaptor_Surface surface(face, true);
		const char* name = "Other";
		switch (surface.GetType())
		{
		case GeomAbs_Plane: name = "Plane"; break;
		case GeomAbs_Cylinder: name = "Cylinder"; break;
		case GeomAbs_Cone: name = "Cone"; break;
		case GeomAbs_Sphere: name = "Sphere"; break;
		case GeomAbs_Torus: name = "Torus"; break;
		case GeomAbs_BezierSurface: name = "Bezier"; break;
		case GeomAbs_BSplineSurface: name = "BSpline"; break;
		case GeomAbs_SurfaceOfRevolution: name = "Revolution"; break;
		case GeomAbs_SurfaceOfExtrusion: name = "Extrusion"; break;
		case GeomAbs_OffsetSurface: name = "Offset"; break;
		default: break;
		}
		std::ostringstream result;
		result << name;
		if (surface.GetType() == GeomAbs_BezierSurface
			|| surface.GetType() == GeomAbs_BSplineSurface)
			result << " degree=" << surface.UDegree() << 'x' << surface.VDegree();
		return result.str();
	}

	struct GapExample
	{
		std::size_t edge = 0;
		double minimum_mm = 0.0;
		double union_max_mm = 0.0;
		double single_face_max_mm = 0.0;
		std::array<double, 6> box{};
		std::array<double, 3> first{};
		std::array<double, 3> last{};
		double length_mm = 0.0;
		std::string curve;
		std::vector<int> owner_faces;
		std::vector<int> best_faces;
		std::vector<int> union_faces;
		std::vector<double> union_profile;
		std::vector<double> best_profile;
		bool partial = false;
		bool near_noncoincident = false;
	};

	void print_edge_face_gap_report(const TopoDS_Shape& shape, double search_mm,
		unsigned sample_count, std::size_t max_edges, bool verbose)
	{
		const auto start = Clock::now();
		ShapeMap faces, edges;
		TopExp::MapShapes(shape, TopAbs_FACE, faces);
		TopExp::MapShapes(shape, TopAbs_EDGE, edges);
		std::vector<std::set<Standard_Integer>> owners(static_cast<std::size_t>(edges.Extent()));
		std::vector<std::size_t> uses(static_cast<std::size_t>(edges.Extent()), 0);
		struct FaceProjection
		{
			TopoDS_Face face;
			Bnd_Box box;
			Handle(Geom_Surface) surface;
			gp_Trsf inverse_location;
		};
		std::vector<FaceProjection> projections(static_cast<std::size_t>(faces.Extent()));
		for (Standard_Integer face_id = 1; face_id <= faces.Extent(); ++face_id)
		{
			FaceProjection& projection = projections[static_cast<std::size_t>(face_id - 1)];
			projection.face = TopoDS::Face(faces(face_id));
			BRepBndLib::AddOptimal(faces(face_id), projection.box,
				false, true);
			TopLoc_Location location;
			projection.surface = BRep_Tool::Surface(projection.face, location);
			projection.inverse_location = location.Transformation().Inverted();
			for (TopExp_Explorer edge_exp(faces(face_id), TopAbs_EDGE); edge_exp.More(); edge_exp.Next())
			{
				const Standard_Integer edge_id = edges.FindIndex(edge_exp.Current());
				if (edge_id <= 0) continue;
				++uses[static_cast<std::size_t>(edge_id - 1)];
				owners[static_cast<std::size_t>(edge_id - 1)].insert(face_id);
			}
		}
		if (verbose)
		{
			std::printf("\n[Source face inventory -- zero-based TopExp IDs]\n");
			for (Standard_Integer face_id = 1; face_id <= faces.Extent(); ++face_id)
			{
				const TopoDS_Face face = TopoDS::Face(faces(face_id));
				const auto face_box = bounds(face);
				std::size_t edge_count = 0;
				for (TopExp_Explorer edge_exp(face, TopAbs_EDGE); edge_exp.More(); edge_exp.Next())
					++edge_count;
				GProp_GProps properties;
				BRepGProp::SurfaceProperties(face, properties);
				const char* orientation = face.Orientation() == TopAbs_REVERSED ? "REVERSED"
					: (face.Orientation() == TopAbs_FORWARD ? "FORWARD" : "OTHER");
				std::printf("  face%d %-24s orient=%-8s edges=%zu area=%.9g mm^2 "
					"bbox=[%.5g %.5g %.5g]-[%.5g %.5g %.5g]\n", face_id - 1,
					surface_description(face).c_str(), orientation, edge_count, properties.Mass(),
					face_box[0], face_box[1], face_box[2], face_box[3], face_box[4], face_box[5]);
			}
			std::printf("  STEPControl_Reader exposes no persistent labels here; IDs above are "
				"the provenance mapping used by this diagnostic.\n");
		}

		std::vector<Standard_Integer> free_edge_ids;
		for (Standard_Integer edge_id = 1; edge_id <= edges.Extent(); ++edge_id)
			if (uses[static_cast<std::size_t>(edge_id - 1)] == 1
				&& !BRep_Tool::Degenerated(TopoDS::Edge(edges(edge_id))))
				free_edge_ids.push_back(edge_id);
		const std::size_t total_free_edges = free_edge_ids.size();
		if (max_edges != 0 && free_edge_ids.size() > max_edges)
		{
			std::vector<Standard_Integer> sampled;
			sampled.reserve(max_edges);
			for (std::size_t sample = 0; sample < max_edges; ++sample)
				sampled.push_back(free_edge_ids[sample * free_edge_ids.size() / max_edges]);
			free_edge_ids.swap(sampled);
		}
		const std::size_t free_edges = free_edge_ids.size();
		std::size_t exact_single_face_support = 0;
		std::size_t exact_union_support = 0;
		std::size_t near_but_noncoincident = 0;
		std::size_t partial_contact_only = 0;
		std::size_t no_nearby_face = 0;
		std::size_t failed_distance = 0;
		std::vector<double> finite_union_gaps;
		std::vector<GapExample> examples;
		const std::array<double, 9> thresholds{
			1.0e-6, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 0.5, 1.0, 2.0, 5.0 };
		std::array<std::size_t, thresholds.size() + 1> histogram{};

		for (Standard_Integer edge_id : free_edge_ids)
		{
			const std::size_t index = static_cast<std::size_t>(edge_id - 1);
			const TopoDS_Edge edge = TopoDS::Edge(edges(edge_id));
			BRepAdaptor_Curve curve(edge);
			const double first = curve.FirstParameter(), last = curve.LastParameter();
			if (!std::isfinite(first) || !std::isfinite(last) || !(last > first))
			{
				++failed_distance;
				continue;
			}
			Bnd_Box edge_box;
			BRepBndLib::AddOptimal(edge, edge_box, false, true);
			edge_box.Enlarge(search_mm);
			std::vector<Standard_Integer> candidates;
			for (Standard_Integer face_id = 1; face_id <= faces.Extent(); ++face_id)
			{
				if (owners[index].contains(face_id)
					|| edge_box.IsOut(projections[static_cast<std::size_t>(face_id - 1)].box)) continue;
				candidates.push_back(face_id);
			}
			if (candidates.empty())
			{
				++no_nearby_face;
				++histogram.back();
				continue;
			}

			std::vector<gp_Pnt> sample_points;
			sample_points.reserve(sample_count);
			for (unsigned sample = 0; sample < sample_count; ++sample)
			{
				const double fraction = static_cast<double>(sample) / (sample_count - 1);
				sample_points.push_back(curve.Value(first + fraction * (last - first)));
			}
			std::vector<double> union_distance(sample_count,
				std::numeric_limits<double>::infinity());
			std::vector<int> union_face(sample_count, -1);
			double best_single_face_max = std::numeric_limits<double>::infinity();
			std::vector<Standard_Integer> best_single_faces;
			std::vector<double> best_single_profile;
			double minimum_edge_face = std::numeric_limits<double>::infinity();
			double intrinsic_tolerance = std::max(Precision::Confusion(), BRep_Tool::Tolerance(edge));
			bool exact_on_one_face = false;
			bool any_distance = false;
			for (Standard_Integer face_id : candidates)
			{
				const FaceProjection& projection = projections[static_cast<std::size_t>(face_id - 1)];
				if (projection.surface.IsNull()) continue;
				double candidate_max = 0.0;
				double candidate_min = std::numeric_limits<double>::infinity();
				bool candidate_valid = true;
				bool candidate_has_sample = false;
				std::vector<double> candidate_profile(sample_count,
					std::numeric_limits<double>::infinity());
				for (unsigned sample = 0; sample < sample_count; ++sample)
				{
					gp_Pnt local_point = sample_points[sample];
					local_point.Transform(projection.inverse_location);
					GeomAPI_ProjectPointOnSurf point_projection(local_point, projection.surface);
					if (point_projection.NbPoints() <= 0)
					{
						candidate_valid = false;
						continue;
					}
					double u = 0.0, v = 0.0;
					point_projection.LowerDistanceParameters(u, v);
					const double face_tolerance = std::max(Precision::PConfusion(),
						BRep_Tool::Tolerance(projection.face));
					const BRepClass_FaceClassifier classifier(projection.face, gp_Pnt2d(u, v),
						face_tolerance, true);
					if (classifier.State() != TopAbs_IN && classifier.State() != TopAbs_ON)
					{
						candidate_valid = false;
						continue;
					}
					const double distance = point_projection.LowerDistance();
					candidate_has_sample = true;
					candidate_profile[sample] = distance;
					candidate_max = std::max(candidate_max, distance);
					candidate_min = std::min(candidate_min, distance);
					if (distance < union_distance[sample])
					{
						union_distance[sample] = distance;
						union_face[sample] = face_id - 1;
					}
				}
				if (!candidate_has_sample || !std::isfinite(candidate_min)) continue;
				minimum_edge_face = std::min(minimum_edge_face, candidate_min);
				if (candidate_min > search_mm) continue;
				any_distance = true;
				intrinsic_tolerance = std::max(intrinsic_tolerance,
					BRep_Tool::Tolerance(projection.face));
				if (candidate_valid
					&& candidate_max <= std::max(1.0e-5, 16.0 * intrinsic_tolerance))
					exact_on_one_face = exact_on_one_face
						|| exact_full_edge_on_face(edge, projection.face);
				if (candidate_valid)
				{
					const double comparison_tolerance = std::max(1.0e-12,
						1.0e-9 * std::max(1.0, std::abs(best_single_face_max)));
					if (!std::isfinite(best_single_face_max)
						|| candidate_max < best_single_face_max - comparison_tolerance)
					{
						best_single_face_max = candidate_max;
						best_single_faces.assign(1, face_id - 1);
						best_single_profile = candidate_profile;
					}
					else if (std::abs(candidate_max - best_single_face_max)
						<= comparison_tolerance)
						best_single_faces.push_back(face_id - 1);
				}
			}
			if (!any_distance)
			{
				++no_nearby_face;
				++histogram.back();
				continue;
			}
			double union_max = 0.0;
			for (double distance : union_distance) union_max = std::max(union_max, distance);
			if (!std::isfinite(union_max))
			{
				++failed_distance;
				++histogram.back();
				continue;
			}
			if (exact_on_one_face) ++exact_single_face_support;
			if (union_max <= intrinsic_tolerance) ++exact_union_support;
			else if (union_max <= search_mm) ++near_but_noncoincident;
			const bool partial = minimum_edge_face <= intrinsic_tolerance
				&& union_max > intrinsic_tolerance;
			if (partial)
				++partial_contact_only;
			finite_union_gaps.push_back(union_max);
			std::size_t bin = 0;
			while (bin < thresholds.size() && union_max > thresholds[bin]) ++bin;
			++histogram[bin];
			const bool near_noncoincident = union_max > intrinsic_tolerance
				&& union_max <= search_mm;
			if (near_noncoincident || (verbose && partial))
			{
				GapExample example;
				example.edge = index;
				example.minimum_mm = minimum_edge_face;
				example.union_max_mm = union_max;
				example.single_face_max_mm = best_single_face_max;
				example.box = bounds(edge);
				example.first = { sample_points.front().X(), sample_points.front().Y(),
					sample_points.front().Z() };
				example.last = { sample_points.back().X(), sample_points.back().Y(),
					sample_points.back().Z() };
				GProp_GProps edge_properties;
				BRepGProp::LinearProperties(edge, edge_properties);
				example.length_mm = edge_properties.Mass();
				example.curve = curve_description(curve);
				for (Standard_Integer owner : owners[index])
					example.owner_faces.push_back(owner - 1);
				for (Standard_Integer candidate : best_single_faces)
					example.best_faces.push_back(candidate);
				example.union_faces = union_face;
				example.union_profile = union_distance;
				example.best_profile = best_single_profile;
				example.partial = partial;
				example.near_noncoincident = near_noncoincident;
				examples.push_back(std::move(example));
			}
		}

		std::sort(finite_union_gaps.begin(), finite_union_gaps.end());
		auto percentile = [&](double p)
		{
			if (finite_union_gaps.empty()) return std::numeric_limits<double>::quiet_NaN();
			const std::size_t index = static_cast<std::size_t>(std::llround(
				p * static_cast<double>(finite_union_gaps.size() - 1)));
			return finite_union_gaps[index];
		};
		std::printf("\n[Native BRep free-edge support gap audit]\n");
		std::printf("  samples/edge=%u; nearby-face search=%.6g mm; runtime %.3f s\n",
			sample_count, search_mm, std::chrono::duration<double>(Clock::now() - start).count());
		std::printf("  non-degenerate free edges analyzed/total: %zu / %zu%s\n", free_edges,
			total_free_edges, free_edges == total_free_edges ? "" : " (uniform ID sample)");
		std::printf("  exact full edge-on-one-face contacts (OCCT curve/surface): %zu\n",
			exact_single_face_support);
		std::printf("  fully supported by union of faces within native tolerance: %zu\n",
			exact_union_support);
		std::printf("  near but genuinely noncoincident (sample max <= %.6g mm): %zu\n",
			search_mm, near_but_noncoincident);
		std::printf("  partial/intersection contact only (sampled min native, rest separated): %zu\n",
			partial_contact_only);
		std::printf("  no nearby supporting face / distance failures: %zu / %zu\n",
			no_nearby_face, failed_distance);
		std::printf("  union-support max-gap p50/p90/p99: %.9g / %.9g / %.9g mm\n",
			percentile(0.50), percentile(0.90), percentile(0.99));
		std::printf("  max-gap histogram:");
		for (std::size_t bin = 0; bin < thresholds.size(); ++bin)
			std::printf(" <=%.4g:%zu", thresholds[bin], histogram[bin]);
		std::printf(" >%.4g/no-support:%zu\n", thresholds.back(), histogram.back());

		if (verbose)
			std::sort(examples.begin(), examples.end(), [](const GapExample& a, const GapExample& b)
				{ return a.edge < b.edge; });
		else
			std::sort(examples.begin(), examples.end(), [](const GapExample& a, const GapExample& b)
				{ return a.union_max_mm > b.union_max_mm; });
		const std::size_t example_limit = verbose ? examples.size()
			: std::min<std::size_t>(examples.size(), 12);
		for (std::size_t i = 0; i < example_limit; ++i)
		{
			const GapExample& example = examples[i];
			std::printf("    edge %zu min/max-union/best-single %.6g / %.6g / %.6g mm "
				"bbox=[%.4g %.4g %.4g]-[%.4g %.4g %.4g]\n", example.edge,
				example.minimum_mm, example.union_max_mm, example.single_face_max_mm,
				example.box[0], example.box[1], example.box[2], example.box[3],
				example.box[4], example.box[5]);
			if (!verbose) continue;
			const char* role = example.partial && example.near_noncoincident
				? "near noncoincident + partial/contact"
				: (example.partial ? "partial/intersection contact" : "near noncoincident seam");
			std::printf("      role=%s; curve=%s; length=%.9g mm\n", role,
				example.curve.c_str(), example.length_mm);
			std::printf("      endpoints=[%.9g %.9g %.9g] -> [%.9g %.9g %.9g] mm\n",
				example.first[0], example.first[1], example.first[2],
				example.last[0], example.last[1], example.last[2]);
			auto print_faces = [&](const char* label, std::vector<int> ids)
			{
				ids.erase(std::remove(ids.begin(), ids.end(), -1), ids.end());
				std::sort(ids.begin(), ids.end());
				ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
				std::printf("      %s:", label);
				if (ids.empty()) std::printf(" none");
				for (int id : ids)
				{
					if (id < 0 || id >= faces.Extent()) continue;
					std::printf(" face%d(%s)", id,
						surface_description(TopoDS::Face(faces(id + 1))).c_str());
				}
				std::printf("\n");
			};
			print_faces("owner source face(s)", example.owner_faces);
			print_faces("best single supporting face(s)", example.best_faces);
			print_faces("union supporting face(s)", example.union_faces);
			std::printf("      union gap profile [fraction:gap@face]:");
			for (std::size_t sample = 0; sample < example.union_profile.size(); ++sample)
			{
				const double fraction = example.union_profile.size() > 1
					? static_cast<double>(sample) / (example.union_profile.size() - 1) : 0.0;
				if (std::isfinite(example.union_profile[sample]))
					std::printf(" %.3g:%.6g@%d", fraction, example.union_profile[sample],
						example.union_faces[sample]);
				else std::printf(" %.3g:out", fraction);
			}
			std::printf(" mm\n");
			if (!example.best_profile.empty())
			{
				std::printf("      best-single gap profile:");
				for (double gap : example.best_profile)
					if (std::isfinite(gap)) std::printf(" %.6g", gap);
					else std::printf(" out");
				std::printf(" mm\n");
			}
		}
		if (verbose)
			std::printf("  face IDs are zero-based TopExp traversal IDs; no STEP/XCAF labels "
				"were available through STEPControl_Reader.\n");
		std::printf("  interpretation: exact/union-supported edges lack topology only; "
			"nonzero max gaps are geometric mismatch, not a tessellation crack.\n");
	}

	struct Chain
	{
		std::vector<gp_Pnt> points;
	};

	bool chains_match(const Chain& a, const Chain& b, double tolerance_mm)
	{
		if (a.points.size() != b.points.size()) return false;
		auto match = [&](bool reverse)
		{
			for (std::size_t i = 0; i < a.points.size(); ++i)
			{
				const std::size_t j = reverse ? b.points.size() - 1 - i : i;
				if (a.points[i].Distance(b.points[j]) > tolerance_mm) return false;
			}
			return true;
		};
		return match(false) || match(true);
	}

	struct QuantizedPoint
	{
		std::int64_t x = 0, y = 0, z = 0;
		friend bool operator<(const QuantizedPoint& a, const QuantizedPoint& b)
		{
			return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z);
		}
	};

	struct SegmentKey
	{
		QuantizedPoint a, b;
		friend bool operator<(const SegmentKey& x, const SegmentKey& y)
		{
			return std::tie(x.a, x.b) < std::tie(y.a, y.b);
		}
	};

	QuantizedPoint quantize(const gp_Pnt& p, double spacing)
	{
		return { static_cast<std::int64_t>(std::llround(p.X() / spacing)),
			static_cast<std::int64_t>(std::llround(p.Y() / spacing)),
			static_cast<std::int64_t>(std::llround(p.Z() / spacing)) };
	}

	struct TessellationStats
	{
		std::size_t missing_faces = 0;
		std::size_t missing_edge_chains = 0;
		std::size_t shared_edges = 0;
		std::size_t conforming_shared_edges = 0;
		std::size_t nonconforming_shared_edges = 0;
		std::size_t free_boundary_segments = 0;
		std::size_t coincident_free_segment_uses = 0;
		std::size_t unique_clear_free_segments = 0;
		std::size_t triangles = 0;
	};

	TessellationStats tessellation_stats(TopoDS_Shape shape, double deflection_mm,
		double geometric_match_mm)
	{
		BRepTools::Clean(shape);
		BRepMesh_IncrementalMesh mesher(shape, deflection_mm);
		mesher.Perform();
		TessellationStats stats;
		ShapeMap faces, edges;
		TopExp::MapShapes(shape, TopAbs_FACE, faces);
		TopExp::MapShapes(shape, TopAbs_EDGE, edges);
		std::vector<std::vector<Chain>> chains(static_cast<std::size_t>(edges.Extent()));
		std::vector<std::size_t> edge_uses(static_cast<std::size_t>(edges.Extent()), 0);
		for (Standard_Integer face_id = 1; face_id <= faces.Extent(); ++face_id)
		{
			const TopoDS_Face face = TopoDS::Face(faces(face_id));
			TopLoc_Location location;
			const Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, location);
			if (triangulation.IsNull()) { ++stats.missing_faces; continue; }
			stats.triangles += static_cast<std::size_t>(triangulation->NbTriangles());
			for (TopExp_Explorer edge_exp(face, TopAbs_EDGE); edge_exp.More(); edge_exp.Next())
			{
				const TopoDS_Edge edge = TopoDS::Edge(edge_exp.Current());
				const Standard_Integer edge_id = edges.FindIndex(edge);
				if (edge_id <= 0) continue;
				++edge_uses[static_cast<std::size_t>(edge_id - 1)];
				const Handle(Poly_PolygonOnTriangulation) polygon =
					BRep_Tool::PolygonOnTriangulation(edge, triangulation, location);
				if (polygon.IsNull() || polygon->NbNodes() < 2)
				{
					++stats.missing_edge_chains;
					continue;
				}
				Chain chain;
				chain.points.reserve(static_cast<std::size_t>(polygon->NbNodes()));
				const gp_Trsf transform = location.Transformation();
				for (Standard_Integer node = 1; node <= polygon->NbNodes(); ++node)
				{
					gp_Pnt p = triangulation->Node(polygon->Node(node));
					p.Transform(transform);
					chain.points.push_back(p);
				}
				chains[static_cast<std::size_t>(edge_id - 1)].push_back(std::move(chain));
			}
		}

		std::map<SegmentKey, std::size_t> free_segments;
		const double quantization = std::max(1.0e-9, geometric_match_mm);
		for (std::size_t edge = 0; edge < chains.size(); ++edge)
		{
			if (edge_uses[edge] > 1)
			{
				++stats.shared_edges;
				bool conforming = chains[edge].size() == edge_uses[edge] && !chains[edge].empty();
				const double tolerance = std::max(geometric_match_mm,
					BRep_Tool::Tolerance(TopoDS::Edge(edges(static_cast<Standard_Integer>(edge + 1)))));
				for (std::size_t i = 1; conforming && i < chains[edge].size(); ++i)
					conforming = chains_match(chains[edge][0], chains[edge][i], tolerance);
				if (conforming) ++stats.conforming_shared_edges;
				else ++stats.nonconforming_shared_edges;
			}
			else if (edge_uses[edge] == 1 && chains[edge].size() == 1)
			{
				const Chain& chain = chains[edge][0];
				for (std::size_t i = 1; i < chain.points.size(); ++i)
				{
					SegmentKey key{ quantize(chain.points[i - 1], quantization),
						quantize(chain.points[i], quantization) };
					if (key.b < key.a) std::swap(key.a, key.b);
					++free_segments[key];
					++stats.free_boundary_segments;
				}
			}
		}
		for (const auto& [segment, uses] : free_segments)
		{
			(void)segment;
			if (uses > 1) stats.coincident_free_segment_uses += uses;
			else ++stats.unique_clear_free_segments;
		}
		return stats;
	}

	void print_tessellation(const char* label, TopoDS_Shape shape, double deflection_mm,
		double match_mm)
	{
		const auto start = Clock::now();
		const TessellationStats stats = tessellation_stats(shape, deflection_mm, match_mm);
		const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
		std::printf("  %s tessellation at %.6g mm: %zu triangles in %.3f s\n", label,
			deflection_mm, stats.triangles, seconds);
		std::printf("    missing faces/edge chains: %zu / %zu\n",
			stats.missing_faces, stats.missing_edge_chains);
		std::printf("    shared CAD edges conforming/nonconforming: %zu / %zu (total %zu)\n",
			stats.conforming_shared_edges, stats.nonconforming_shared_edges, stats.shared_edges);
		std::printf("    topologically-free tess segments: %zu; geometrically paired uses at %.6g mm: "
			"%zu; unique/open: %zu\n", stats.free_boundary_segments, match_mm,
			stats.coincident_free_segment_uses, stats.unique_clear_free_segments);
	}

	struct ProvenanceStats
	{
		std::size_t source_faces = 0;
		std::size_t unchanged_sources = 0;
		std::size_t modified_sources = 0;
		std::size_t deleted_sources = 0;
		std::size_t sources_without_image = 0;
		std::size_t result_faces = 0;
		std::size_t result_with_one_parent = 0;
		std::size_t result_with_multiple_parents = 0;
		std::size_t orphan_result_faces = 0;
	};

	ProvenanceStats provenance_stats(BOPAlgo_Splitter& splitter,
		const std::vector<TopoDS_Face>& sources, const TopoDS_Shape& result)
	{
		ShapeMap result_faces;
		TopExp::MapShapes(result, TopAbs_FACE, result_faces);
		std::vector<std::set<std::size_t>> parents(static_cast<std::size_t>(result_faces.Extent()));
		ProvenanceStats stats;
		stats.source_faces = sources.size();
		stats.result_faces = static_cast<std::size_t>(result_faces.Extent());
		for (std::size_t source_id = 0; source_id < sources.size(); ++source_id)
		{
			const TopoDS_Face& source = sources[source_id];
			const NCollection_List<TopoDS_Shape>& modified = splitter.Modified(source);
			bool has_image = false;
			for (NCollection_List<TopoDS_Shape>::Iterator it(modified); it.More(); it.Next())
				for (TopExp_Explorer face_exp(it.Value(), TopAbs_FACE); face_exp.More(); face_exp.Next())
				{
					const Standard_Integer result_id = result_faces.FindIndex(face_exp.Current());
					if (result_id > 0)
					{
						parents[static_cast<std::size_t>(result_id - 1)].insert(source_id);
						has_image = true;
					}
				}
			if (!has_image)
			{
				const Standard_Integer result_id = result_faces.FindIndex(source);
				if (result_id > 0)
				{
					parents[static_cast<std::size_t>(result_id - 1)].insert(source_id);
					has_image = true;
					++stats.unchanged_sources;
				}
			}
			else ++stats.modified_sources;
			if (splitter.IsDeleted(source)) ++stats.deleted_sources;
			if (!has_image) ++stats.sources_without_image;
		}
		for (const auto& source_set : parents)
		{
			if (source_set.empty()) ++stats.orphan_result_faces;
			else if (source_set.size() == 1) ++stats.result_with_one_parent;
			else ++stats.result_with_multiple_parents;
		}
		return stats;
	}

	void run_split(const std::vector<TopoDS_Face>& faces, const Options& options)
	{
		NCollection_List<TopoDS_Shape> arguments;
		for (const TopoDS_Face& face : faces) arguments.Append(face);
		BOPAlgo_Splitter splitter;
		splitter.SetArguments(arguments);
		splitter.SetFuzzyValue(options.fuzzy_mm);
		splitter.SetRunParallel(true);
		splitter.SetUseOBB(true);
		splitter.SetNonDestructive(true);
		splitter.SetToFillHistory(true);
		std::printf("\n[General Fuse / Splitter]\n");
		std::printf("  arguments=%zu, fuzzy=%.9g mm, parallel=yes, non-destructive=yes\n",
			faces.size(), options.fuzzy_mm);
		const auto start = Clock::now();
		splitter.Perform();
		const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
		std::printf("  runtime: %.3f s; errors=%s; warnings=%s\n", seconds,
			splitter.HasErrors() ? "YES" : "no", splitter.HasWarnings() ? "yes" : "no");
		if (splitter.HasErrors()) splitter.DumpErrors(std::cerr);
		if (splitter.HasWarnings()) splitter.DumpWarnings(std::cerr);
		if (splitter.HasErrors() || splitter.Shape().IsNull()) return;
		print_topology("General Fuse result", splitter.Shape());
		print_tessellation("General Fuse result", splitter.Shape(), options.deflection_mm,
			std::max(options.fuzzy_mm, 1.0e-7));
		const ProvenanceStats history = provenance_stats(splitter, faces, splitter.Shape());
		std::printf("  history source unchanged/modified/deleted/no-image: %zu / %zu / %zu / %zu\n",
			history.unchanged_sources, history.modified_sources, history.deleted_sources,
			history.sources_without_image);
		std::printf("  history result one-parent/multi-parent/orphan: %zu / %zu / %zu\n",
			history.result_with_one_parent, history.result_with_multiple_parents,
			history.orphan_result_faces);
		std::printf("  source-face provenance feasibility: %s\n",
			history.sources_without_image == 0 && history.result_with_multiple_parents == 0
				&& history.orphan_result_faces == 0 ? "YES" : "NO / needs geometric fallback");
	}

	TopoDS_Shape run_sewing(const std::vector<TopoDS_Face>& faces, const Options& options)
	{
		std::printf("\n[Optional gap repair: BRepBuilderAPI_Sewing]\n");
		std::printf("  tolerance=%.9g mm; cutting=yes; non-manifold=yes\n", options.sew_mm);
		BRepBuilderAPI_Sewing sewing(options.sew_mm, true, true, true, true);
		sewing.SetLocalTolerancesMode(false);
		for (const TopoDS_Face& face : faces) sewing.Add(face);
		const auto start = Clock::now();
		sewing.Perform();
		const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
		std::printf("  runtime %.3f s; tool free/contiguous/multiple/deleted: %d / %d / %d / %d\n",
			seconds, sewing.NbFreeEdges(), sewing.NbContigousEdges(), sewing.NbMultipleEdges(),
			sewing.NbDeletedFaces());
		if (sewing.SewedShape().IsNull())
		{
			std::printf("  sewing returned a null shape\n");
			return {};
		}
		print_topology("Sewing result (diagnostic only)", sewing.SewedShape());
		print_tessellation("Sewing result", sewing.SewedShape(), options.deflection_mm,
			std::max(options.sew_mm, 1.0e-7));
		return sewing.SewedShape();
	}

	void run_shape_fix(const TopoDS_Shape& selected, const Options& options)
	{
		std::printf("\n[Optional local repair: ShapeFix_Shape]\n");
		std::printf("  precision/max tolerance %.9g mm; no sewing or surface filling\n",
			options.sew_mm);
		ShapeFix_Shape fix(selected);
		fix.SetPrecision(options.sew_mm);
		fix.SetMinTolerance(1.0e-7);
		fix.SetMaxTolerance(std::max(options.sew_mm, 1.0e-7));
		const auto start = Clock::now();
		const bool changed = fix.Perform();
		const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
		const TopoDS_Shape result = fix.Shape();
		std::printf("  runtime %.3f s; reported changed=%s\n", seconds, changed ? "yes" : "no");
		if (!result.IsNull()) print_topology("ShapeFix result (diagnostic only)", result);
	}
}

int main(int argc, char** argv)
{
	Options options;
	if (!parse_options(argc, argv, options)) { usage(argv[0]); return 2; }
	try
	{
		std::printf("OCCT BRep normalization diagnostic -- production EB is not used or modified\n");
		std::printf("STEP: %s\n", options.step.c_str());
		const auto read_start = Clock::now();
		const TopoDS_Shape input = read_step(options.step);
		const double read_seconds = std::chrono::duration<double>(Clock::now() - read_start).count();
		const std::vector<TopoDS_Face> all_faces = unique_faces(input);
		const std::vector<TopoDS_Face> selected = central_subset(all_faces, options.subset_faces);
		const TopoDS_Shape selected_shape = face_compound(selected);
		std::printf("STEP read/transfer %.3f s; source faces %zu; selected %zu%s\n", read_seconds,
			all_faces.size(), selected.size(), selected.size() == all_faces.size() ? " (all)" : " (central subset)");
		print_topology("Original selected BRep faces", selected_shape);
		print_edge_face_gap_report(selected_shape, options.gap_search_mm, options.gap_samples,
			options.gap_max_edges, options.verbose_gaps);
		print_tessellation("Original", selected_shape, options.deflection_mm,
			std::max(options.fuzzy_mm, 1.0e-7));
		TopoDS_Shape sewed;
		if (options.run_sewing || options.split_after_sew) sewed = run_sewing(selected, options);
		if (options.run_split)
		{
			if (options.split_after_sew && !sewed.IsNull())
			{
				const std::vector<TopoDS_Face> sewed_faces = unique_faces(sewed);
				std::printf("\nGeneral Fuse will use %zu faces from the sewing result; "
					"history is relative to those faces.\n", sewed_faces.size());
				run_split(sewed_faces, options);
			}
			else run_split(selected, options);
		}
		if (options.run_shape_fix) run_shape_fix(selected_shape, options);
		std::printf("\nNOTE: General Fuse only imprints/splits existing faces; it does not fill openings.\n");
		std::printf("      Sewing and ShapeFix results above are separate experiments, not inputs to Fuse.\n");
		return 0;
	}
	catch (const Standard_Failure& failure)
	{
		std::fprintf(stderr, "OCCT exception: %s\n", failure.GetMessageString());
	}
	catch (const std::exception& exception)
	{
		std::fprintf(stderr, "error: %s\n", exception.what());
	}
	return 1;
}
