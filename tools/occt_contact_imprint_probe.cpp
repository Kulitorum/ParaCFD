// Diagnostic only: evaluate exact, local CAD-contact imprinting before tessellation.
//
// This deliberately stays outside the production STEP/EB path.  It compares local
// OCCT transactions for a certified free-edge-on-face contact:
//   1. a two-face General Fuse;
//   2. a target-face split by the exact owner edge;
//   3. the target split followed by native-tolerance topological sewing;
//   4. a Cartesian-cell General Fuse using the original BRep faces (not facets).
// There is no fuzzy/gap repair, tolerance inflation, DRAWEXE, or GUI code here.

#include "core/geometry/step_import.h"

#include <BOPAlgo_CellsBuilder.hxx>
#include <BOPAlgo_Splitter.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <GProp_GProps.hxx>
#include <IntTools_CommonPrt.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <IntTools_EdgeFace.hxx>
#include <IntTools_Range.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_State.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
	using Clock = std::chrono::steady_clock;
	using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;

	struct Options
	{
		std::string step = "Test-Data/Two-Cells.stp";
		std::set<std::size_t> excluded_faces;
		std::size_t contact_offset = 0;
		std::size_t max_contacts = 4;
		double cell_mm = 62.5;
		double deflection_mm = 2.0;
		bool run_global = false;
		bool run_global_native = false;
		bool run_growing_native = false;
		bool run_batch_target_prototype = false;
		bool inspect_overlaps = false;
		bool scan_only = false;
	};

	void usage(const char* exe)
	{
		std::printf("Usage: %s [--step FILE] [--exclude-faces 8,9] "
			"[--contact-offset N] [--max-contacts N] [--cell-mm H] "
			"[--deflection-mm D] [--global] [--global-native] [--growing-native] "
			"[--batch-target-prototype] "
			"[--inspect-overlaps] [--scan-only]\n", exe);
	}

	bool parse_face_list(const std::string& text, std::set<std::size_t>& result)
	{
		std::size_t begin = 0;
		while (begin < text.size())
		{
			const std::size_t end = text.find(',', begin);
			const std::string token = text.substr(begin,
				end == std::string::npos ? std::string::npos : end - begin);
			if (token.empty()) return false;
			char* tail = nullptr;
			const unsigned long long value = std::strtoull(token.c_str(), &tail, 10);
			if (!tail || *tail != '\0') return false;
			result.insert(static_cast<std::size_t>(value));
			if (end == std::string::npos) break;
			begin = end + 1;
		}
		return true;
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
			else if (arg == "--exclude-faces")
			{
				const char* v = value("--exclude-faces");
				if (!v || !parse_face_list(v, options.excluded_faces)) return false;
			}
			else if (arg == "--contact-offset")
			{
				const char* v = value("--contact-offset"); if (!v) return false;
				options.contact_offset = static_cast<std::size_t>(std::strtoull(v, nullptr, 10));
			}
			else if (arg == "--max-contacts")
			{
				const char* v = value("--max-contacts"); if (!v) return false;
				options.max_contacts = static_cast<std::size_t>(std::strtoull(v, nullptr, 10));
			}
			else if (arg == "--cell-mm")
			{
				const char* v = value("--cell-mm"); if (!v) return false;
				options.cell_mm = std::strtod(v, nullptr);
			}
			else if (arg == "--deflection-mm")
			{
				const char* v = value("--deflection-mm"); if (!v) return false;
				options.deflection_mm = std::strtod(v, nullptr);
			}
			else if (arg == "--global") options.run_global = true;
			else if (arg == "--global-native") options.run_global_native = true;
			else if (arg == "--growing-native") options.run_growing_native = true;
			else if (arg == "--batch-target-prototype") options.run_batch_target_prototype = true;
			else if (arg == "--inspect-overlaps") options.inspect_overlaps = true;
			else if (arg == "--scan-only") options.scan_only = true;
			else if (arg == "--help" || arg == "-h") { usage(argv[0]); std::exit(0); }
			else { std::fprintf(stderr, "unknown option: %s\n", arg.c_str()); return false; }
		}
		return options.max_contacts > 0 && options.cell_mm > 0.0
			&& options.deflection_mm > 0.0;
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

	TopoDS_Shape compound(const std::vector<TopoDS_Shape>& shapes)
	{
		TopoDS_Compound result;
		BRep_Builder builder;
		builder.MakeCompound(result);
		for (const TopoDS_Shape& shape : shapes) builder.Add(result, shape);
		return result;
	}

	Bnd_Box bounds(const TopoDS_Shape& shape)
	{
		Bnd_Box box;
		BRepBndLib::AddOptimal(shape, box, false, true);
		return box;
	}

	bool exact_full_edge_on_face(const TopoDS_Edge& edge, const TopoDS_Face& face,
		std::uint32_t& added_sectors)
	{
		BRepAdaptor_Curve curve(edge);
		const double first = curve.FirstParameter(), last = curve.LastParameter();
		if (!std::isfinite(first) || !std::isfinite(last) || !(last > first)) return false;
		IntTools_EdgeFace intersection;
		intersection.SetEdge(edge);
		intersection.SetFace(face);
		intersection.SetRange(first, last);
		intersection.SetFuzzyValue(0.0);
		intersection.Perform();
		if (!intersection.IsDone()) return false;
		std::vector<std::pair<double, double>> ranges;
		for (NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
			intersection.CommonParts()); common.More(); common.Next())
		{
			if (common.Value().Type() != TopAbs_EDGE) continue;
			double lo = 0.0, hi = 0.0;
			common.Value().Range1(lo, hi);
			if (hi < lo) std::swap(lo, hi);
			ranges.emplace_back(std::max(lo, first), std::min(hi, last));
		}
		if (ranges.empty()) return false;
		std::sort(ranges.begin(), ranges.end());
		const double parameter_tolerance = std::max(Precision::PConfusion(),
			64.0 * std::numeric_limits<double>::epsilon()
			* std::max({ 1.0, std::abs(first), std::abs(last) }));
		double covered = first;
		for (const auto& range : ranges)
		{
			if (range.first > covered + parameter_tolerance) return false;
			covered = std::max(covered, range.second);
		}
		if (covered < last - parameter_tolerance) return false;
		const gp_Pnt midpoint = curve.Value(0.5 * (first + last));
		const BRepClass_FaceClassifier classifier(face, midpoint, 0.0, true);
		if (classifier.State() == TopAbs_IN) added_sectors = 2;
		else if (classifier.State() == TopAbs_ON) added_sectors = 1;
		else return false;
		return true;
	}

	bool ranges_cover(std::vector<std::pair<double, double>> ranges, double first,
		double last)
	{
		if (last < first) std::swap(first, last);
		if (ranges.empty()) return false;
		for (auto& range : ranges) if (range.second < range.first)
			std::swap(range.first, range.second);
		std::sort(ranges.begin(), ranges.end());
		const double tolerance = std::max(Precision::PConfusion(),
			64.0 * std::numeric_limits<double>::epsilon()
			* std::max({ 1.0, std::abs(first), std::abs(last) }));
		double covered = first;
		for (const auto& range : ranges)
		{
			if (range.second < first || range.first > last) continue;
			if (range.first > covered + tolerance) return false;
			covered = std::max(covered, std::min(last, range.second));
		}
		return covered >= last - tolerance;
	}

	bool exact_same_edge_extent(const TopoDS_Edge& a, const TopoDS_Edge& b)
	{
		BRepAdaptor_Curve curve_a(a), curve_b(b);
		const double first_a = curve_a.FirstParameter(), last_a = curve_a.LastParameter();
		const double first_b = curve_b.FirstParameter(), last_b = curve_b.LastParameter();
		if (!std::isfinite(first_a) || !std::isfinite(last_a) || !(last_a > first_a)
			|| !std::isfinite(first_b) || !std::isfinite(last_b) || !(last_b > first_b))
			return false;
		IntTools_EdgeEdge intersection(a, b);
		intersection.SetFuzzyValue(0.0);
		intersection.Perform();
		if (!intersection.IsDone()) return false;
		std::vector<std::pair<double, double>> ranges_a, ranges_b;
		for (NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
			intersection.CommonParts()); common.More(); common.Next())
		{
			if (common.Value().Type() != TopAbs_EDGE) continue;
			double lo = 0.0, hi = 0.0;
			common.Value().Range1(lo, hi);
			ranges_a.emplace_back(lo, hi);
			for (NCollection_Sequence<IntTools_Range>::Iterator range(
				common.Value().Ranges2()); range.More(); range.Next())
				ranges_b.emplace_back(range.Value().First(), range.Value().Last());
		}
		return ranges_cover(std::move(ranges_a), first_a, last_a)
			&& ranges_cover(std::move(ranges_b), first_b, last_b);
	}

	struct Contact
	{
		std::size_t edge = 0;
		std::size_t owner_face = 0;
		std::size_t target_face = 0;
		std::uint32_t target_sectors = 0;
		std::uint32_t expected_fan = 0;
	};

	struct ContactGroup
	{
		std::size_t edge = 0;
		std::size_t owner_face = 0;
		std::vector<Contact> targets;
		std::uint32_t expected_fan = 1;
	};

	struct Model
	{
		TopoDS_Shape shape;
		std::vector<TopoDS_Face> faces;
		ShapeMap edges;
		std::vector<std::vector<std::size_t>> edge_owners;
		std::vector<Bnd_Box> face_bounds;
	};

	Model make_model(TopoDS_Shape shape)
	{
		Model model;
		model.shape = std::move(shape);
		model.faces = unique_faces(model.shape);
		TopExp::MapShapes(model.shape, TopAbs_EDGE, model.edges);
		model.edge_owners.resize(static_cast<std::size_t>(model.edges.Extent()));
		model.face_bounds.reserve(model.faces.size());
		for (std::size_t face = 0; face < model.faces.size(); ++face)
		{
			model.face_bounds.push_back(bounds(model.faces[face]));
			for (TopExp_Explorer exp(model.faces[face], TopAbs_EDGE); exp.More(); exp.Next())
			{
				const Standard_Integer edge = model.edges.FindIndex(exp.Current());
				if (edge > 0) model.edge_owners[static_cast<std::size_t>(edge - 1)].push_back(face);
			}
		}
		return model;
	}

	std::vector<ContactGroup> find_contacts(const Model& model, const Options& options)
	{
		std::vector<ContactGroup> result;
		const std::size_t stop_after = options.max_contacts
			> std::numeric_limits<std::size_t>::max() - options.contact_offset
			? std::numeric_limits<std::size_t>::max()
			: options.contact_offset + options.max_contacts;
		for (std::size_t edge_id = 0; edge_id < model.edge_owners.size(); ++edge_id)
		{
			const auto& owners = model.edge_owners[edge_id];
			if (owners.size() != 1 || options.excluded_faces.contains(owners.front())) continue;
			const TopoDS_Edge edge = TopoDS::Edge(model.edges(static_cast<Standard_Integer>(edge_id + 1)));
			const Bnd_Box edge_bounds = bounds(edge);
			ContactGroup group;
			group.edge = edge_id;
			group.owner_face = owners.front();
			for (std::size_t target = 0; target < model.faces.size(); ++target)
			{
				if (target == owners.front() || options.excluded_faces.contains(target)
					|| edge_bounds.IsOut(model.face_bounds[target])) continue;
				std::uint32_t sectors = 0;
				if (!exact_full_edge_on_face(edge, model.faces[target], sectors)) continue;
				group.targets.push_back({ edge_id, owners.front(), target, sectors, 1u + sectors });
				group.expected_fan += sectors;
			}
			if (!group.targets.empty()) result.push_back(std::move(group));
			if (result.size() >= stop_after) return result;
		}
		return result;
	}

	template<class Builder, class ResultShape>
	std::vector<TopoDS_Face> face_images(Builder& builder, const TopoDS_Face& source,
		const ResultShape& result_shape)
	{
		TopTools_IndexedMapOfShape all_faces;
		TopExp::MapShapes(result_shape, TopAbs_FACE, all_faces);
		TopTools_MapOfShape seen;
		std::vector<TopoDS_Face> result;
		const NCollection_List<TopoDS_Shape>& modified = builder.Modified(source);
		for (NCollection_List<TopoDS_Shape>::Iterator it(modified); it.More(); it.Next())
			if (it.Value().ShapeType() == TopAbs_FACE && all_faces.FindIndex(it.Value()) > 0
				&& seen.Add(it.Value())) result.push_back(TopoDS::Face(it.Value()));
		if (all_faces.FindIndex(source) > 0 && seen.Add(source)) result.push_back(source);
		return result;
	}

	template<class Builder, class ResultShape>
	std::vector<TopoDS_Edge> edge_images(Builder& builder, const TopoDS_Edge& source,
		const ResultShape& result_shape)
	{
		TopTools_IndexedMapOfShape all_edges;
		TopExp::MapShapes(result_shape, TopAbs_EDGE, all_edges);
		TopTools_MapOfShape seen;
		std::vector<TopoDS_Edge> result;
		const NCollection_List<TopoDS_Shape>& modified = builder.Modified(source);
		for (NCollection_List<TopoDS_Shape>::Iterator it(modified); it.More(); it.Next())
			if (it.Value().ShapeType() == TopAbs_EDGE && all_edges.FindIndex(it.Value()) > 0
				&& seen.Add(it.Value())) result.push_back(TopoDS::Edge(it.Value()));
		if (all_edges.FindIndex(source) > 0 && seen.Add(source)) result.push_back(source);
		return result;
	}

	std::size_t edge_face_incidence(const TopoDS_Edge& edge,
		const std::vector<TopoDS_Face>& faces)
	{
		std::size_t count = 0;
		for (const TopoDS_Face& face : faces)
		{
			bool found = false;
			for (TopExp_Explorer exp(face, TopAbs_EDGE); exp.More(); exp.Next())
				if (exp.Current().IsSame(edge)) { found = true; break; }
			if (found) ++count;
		}
		return count;
	}

	std::size_t geometric_edge_fan(const TopoDS_Edge& source,
		const TopoDS_Shape& result_shape, const std::vector<TopoDS_Face>& faces,
		std::size_t* matching_edges = nullptr)
	{
		ShapeMap edges;
		TopExp::MapShapes(result_shape, TopAbs_EDGE, edges);
		std::size_t fan = 0, matches = 0;
		for (Standard_Integer i = 1; i <= edges.Extent(); ++i)
		{
			const TopoDS_Edge candidate = TopoDS::Edge(edges(i));
			if (!exact_same_edge_extent(source, candidate)) continue;
			++matches;
			fan = std::max(fan, edge_face_incidence(candidate, faces));
		}
		if (matching_edges) *matching_edges = matches;
		return fan;
	}

	template<class Builder>
	std::size_t reproduced_fan(Builder& builder, const TopoDS_Edge& source,
		const TopoDS_Shape& result_shape, const std::vector<TopoDS_Face>& retained_faces)
	{
		std::size_t fan = 0;
		for (const TopoDS_Edge& image : edge_images(builder, source, result_shape))
			fan = std::max(fan, edge_face_incidence(image, retained_faces));
		return fan;
	}

	struct Provenance
	{
		std::size_t result_faces = 0;
		std::size_t one_parent = 0;
		std::size_t multi_parent = 0;
		std::size_t orphan = 0;
		std::size_t missing_source = 0;
	};

	template<class Builder>
	Provenance provenance(Builder& builder, const std::vector<TopoDS_Face>& sources,
		const TopoDS_Shape& result_shape, const std::vector<TopoDS_Face>* selected = nullptr)
	{
		TopTools_IndexedMapOfShape result_faces;
		TopExp::MapShapes(result_shape, TopAbs_FACE, result_faces);
		std::vector<std::set<std::size_t>> parents(static_cast<std::size_t>(result_faces.Extent()));
		Provenance stats;
		stats.result_faces = selected ? selected->size() : static_cast<std::size_t>(result_faces.Extent());
		for (std::size_t source = 0; source < sources.size(); ++source)
		{
			bool found = false;
			for (const TopoDS_Face& image : face_images(builder, sources[source], result_shape))
			{
				const Standard_Integer index = result_faces.FindIndex(image);
				if (index <= 0) continue;
				parents[static_cast<std::size_t>(index - 1)].insert(source);
				found = true;
			}
			if (!found) ++stats.missing_source;
		}
		auto classify = [&](const TopoDS_Face& face)
		{
			const Standard_Integer index = result_faces.FindIndex(face);
			const std::size_t count = index > 0 ? parents[static_cast<std::size_t>(index - 1)].size() : 0;
			if (count == 0) ++stats.orphan;
			else if (count == 1) ++stats.one_parent;
			else ++stats.multi_parent;
		};
		if (selected) for (const TopoDS_Face& face : *selected) classify(face);
		else for (Standard_Integer i = 1; i <= result_faces.Extent(); ++i)
			classify(TopoDS::Face(result_faces(i)));
		return stats;
	}

	struct Tessellation
	{
		std::size_t triangles = 0;
		std::size_t shared_edges = 0;
		std::size_t conforming_edges = 0;
		std::size_t nonconforming_edges = 0;
		std::size_t missing_chains = 0;
	};

	bool same_chain(const std::vector<gp_Pnt>& a, const std::vector<gp_Pnt>& b,
		double tolerance)
	{
		if (a.size() != b.size()) return false;
		auto matches = [&](bool reverse)
		{
			for (std::size_t i = 0; i < a.size(); ++i)
			{
				const std::size_t j = reverse ? b.size() - 1 - i : i;
				if (a[i].Distance(b[j]) > tolerance) return false;
			}
			return true;
		};
		return matches(false) || matches(true);
	}

	Tessellation tessellate(TopoDS_Shape shape, double deflection_mm)
	{
		BRepTools::Clean(shape);
		BRepMesh_IncrementalMesh mesher(shape, deflection_mm, false, 0.1, false);
		Tessellation stats;
		ShapeMap faces, edges;
		TopExp::MapShapes(shape, TopAbs_FACE, faces);
		TopExp::MapShapes(shape, TopAbs_EDGE, edges);
		std::vector<std::vector<std::vector<gp_Pnt>>> chains(static_cast<std::size_t>(edges.Extent()));
		std::vector<std::size_t> uses(static_cast<std::size_t>(edges.Extent()), 0);
		for (Standard_Integer face_id = 1; face_id <= faces.Extent(); ++face_id)
		{
			const TopoDS_Face face = TopoDS::Face(faces(face_id));
			TopLoc_Location location;
			const Handle(Poly_Triangulation) mesh = BRep_Tool::Triangulation(face, location);
			if (mesh.IsNull()) continue;
			stats.triangles += static_cast<std::size_t>(mesh->NbTriangles());
			for (TopExp_Explorer exp(face, TopAbs_EDGE); exp.More(); exp.Next())
			{
				const TopoDS_Edge edge = TopoDS::Edge(exp.Current());
				const Standard_Integer edge_id = edges.FindIndex(edge);
				if (edge_id <= 0) continue;
				++uses[static_cast<std::size_t>(edge_id - 1)];
				const Handle(Poly_PolygonOnTriangulation) polygon =
					BRep_Tool::PolygonOnTriangulation(edge, mesh, location);
				if (polygon.IsNull()) { ++stats.missing_chains; continue; }
				std::vector<gp_Pnt> chain;
				chain.reserve(static_cast<std::size_t>(polygon->NbNodes()));
				const gp_Trsf transform = location.Transformation();
				for (Standard_Integer node = 1; node <= polygon->NbNodes(); ++node)
				{
					gp_Pnt point = mesh->Node(polygon->Node(node));
					point.Transform(transform);
					chain.push_back(point);
				}
				chains[static_cast<std::size_t>(edge_id - 1)].push_back(std::move(chain));
			}
		}
		for (std::size_t edge = 0; edge < uses.size(); ++edge)
		{
			if (uses[edge] < 2) continue;
			++stats.shared_edges;
			bool conforming = chains[edge].size() == uses[edge] && !chains[edge].empty();
			const double tolerance = std::max(1.0e-8,
				BRep_Tool::Tolerance(TopoDS::Edge(edges(static_cast<Standard_Integer>(edge + 1)))));
			for (std::size_t chain = 1; conforming && chain < chains[edge].size(); ++chain)
				conforming = same_chain(chains[edge][0], chains[edge][chain], tolerance);
			if (conforming) ++stats.conforming_edges;
			else ++stats.nonconforming_edges;
		}
		return stats;
	}

	std::vector<TopoDS_Face> result_faces(const TopoDS_Shape& shape)
	{
		return unique_faces(shape);
	}

	void print_transaction(const char* label, const TopoDS_Shape& shape,
		const Provenance& history, std::size_t fan, std::uint32_t expected,
		double milliseconds, double deflection)
	{
		const Tessellation mesh = tessellate(shape, deflection);
		std::printf("    %-18s %7.2f ms valid=%s faces=%zu provenance 1/multi/orphan/missing=%zu/%zu/%zu/%zu "
			"fan=%zu/%u tess tris=%zu shared conform/non=%zu/%zu\n", label, milliseconds,
			BRepCheck_Analyzer(shape, true).IsValid() ? "yes" : "NO", history.result_faces,
			history.one_parent, history.multi_parent, history.orphan, history.missing_source,
			fan, expected, mesh.triangles, mesh.conforming_edges, mesh.nonconforming_edges);
	}

	void run_boundary_relink(const TopoDS_Face& owner, const TopoDS_Face& target,
		const TopoDS_Edge& owner_edge, const Options& options)
	{
		TopoDS_Edge target_edge;
		for (TopExp_Explorer exp(target, TopAbs_EDGE); exp.More(); exp.Next())
		{
			const TopoDS_Edge candidate = TopoDS::Edge(exp.Current());
			if (!exact_same_edge_extent(owner_edge, candidate)) continue;
			if (!target_edge.IsNull())
			{
				std::printf("    boundary-relink    ambiguous target edge\n");
				return;
			}
			target_edge = candidate;
		}
		if (target_edge.IsNull())
		{
			std::printf("    boundary-relink    not applicable (contact is not one equal target edge)\n");
			return;
		}

		const auto begin = Clock::now();
		BRepBuilderAPI_Copy copy(owner_edge, false, false);
		if (!copy.IsDone() || copy.Shape().IsNull())
		{
			std::printf("    boundary-relink    could not copy canonical edge\n");
			return;
		}
		TopoDS_Edge canonical = TopoDS::Edge(copy.Shape());
		Standard_Real first = 0.0, last = 0.0;
		if (BRep_Tool::CurveOnSurface(target_edge, target, first, last).IsNull())
		{
			std::printf("    boundary-relink    target edge has no pcurve\n");
			return;
		}
		BRep_Builder edge_builder;
		const double tolerance = std::max({ BRep_Tool::Tolerance(owner_edge),
			BRep_Tool::Tolerance(target_edge), Precision::Confusion() });
		// Copy the complete curve-on-surface representation, including its native
		// parameterization, rather than projecting or expanding a tolerance.
		edge_builder.Transfert(target_edge, canonical);
		edge_builder.UpdateEdge(canonical, tolerance);
		BRepLib::SameParameter(canonical, tolerance, true);

		Handle(BRepTools_ReShape) owner_reshape = new BRepTools_ReShape;
		owner_reshape->Replace(owner_edge, canonical);
		const TopoDS_Face new_owner = TopoDS::Face(owner_reshape->Apply(owner));
		Handle(BRepTools_ReShape) target_reshape = new BRepTools_ReShape;
		target_reshape->Replace(target_edge, canonical);
		const TopoDS_Face new_target = TopoDS::Face(target_reshape->Apply(target));
		const TopoDS_Shape joined = compound({ new_owner, new_target });
		const double elapsed = std::chrono::duration<double, std::milli>(
			Clock::now() - begin).count();
		Provenance history;
		history.result_faces = 2;
		history.one_parent = 2;
		print_transaction("boundary-relink", joined, history,
			edge_face_incidence(canonical, { new_owner, new_target }), 2,
			elapsed, options.deflection_mm);
	}

	void run_local_sew(const TopoDS_Face& owner, const TopoDS_Face& target,
		const TopoDS_Edge& owner_edge, const Options& options)
	{
		double tolerance = std::max({ BRep_Tool::Tolerance(owner_edge),
			BRep_Tool::Tolerance(owner), BRep_Tool::Tolerance(target), Precision::Confusion() });
		for (TopExp_Explorer edge(target, TopAbs_EDGE); edge.More(); edge.Next())
			tolerance = std::max(tolerance, BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
		// This is not gap repair: the transaction contains only an already exact-certified
		// pair and never exceeds the source shapes' own tolerance.
		BRepBuilderAPI_Sewing sewing(tolerance, true, true, true, true);
		sewing.SetLocalTolerancesMode(true);
		sewing.Add(owner);
		sewing.Add(target);
		const auto begin = Clock::now(); sewing.Perform();
		const double elapsed = std::chrono::duration<double, std::milli>(
			Clock::now() - begin).count();
		const TopoDS_Shape joined = sewing.SewedShape();
		if (joined.IsNull())
		{
			std::printf("    local-native-sew   FAILED after %.2f ms\n", elapsed);
			return;
		}
		TopoDS_Edge canonical = owner_edge;
		const TopoDS_Shape modified_edge = sewing.ModifiedSubShape(owner_edge);
		if (!modified_edge.IsNull() && modified_edge.ShapeType() == TopAbs_EDGE)
			canonical = TopoDS::Edge(modified_edge);
		const auto faces = result_faces(joined);
		Provenance history;
		history.result_faces = faces.size();
		history.one_parent = faces.size();
		print_transaction("local-native-sew", joined, history,
			edge_face_incidence(canonical, faces), 2, elapsed, options.deflection_mm);
	}

	void run_local_sew_group(const Model& model, const ContactGroup& contact,
		const Options& options)
	{
		std::vector<TopoDS_Face> sources{ model.faces[contact.owner_face] };
		for (const Contact& target : contact.targets)
			sources.push_back(model.faces[target.target_face]);
		double tolerance = Precision::Confusion();
		for (const TopoDS_Face& face : sources)
		{
			tolerance = std::max(tolerance, BRep_Tool::Tolerance(face));
			for (TopExp_Explorer edge(face, TopAbs_EDGE); edge.More(); edge.Next())
				tolerance = std::max(tolerance,
					BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
		}
		BRepBuilderAPI_Sewing sewing(tolerance, true, true, true, true);
		sewing.SetLocalTolerancesMode(true);
		for (const TopoDS_Face& face : sources) sewing.Add(face);
		const auto begin = Clock::now(); sewing.Perform();
		const double elapsed = std::chrono::duration<double, std::milli>(
			Clock::now() - begin).count();
		const TopoDS_Shape joined = sewing.SewedShape();
		if (joined.IsNull())
		{
			std::printf("    contact-group-sew  FAILED after %.2f ms\n", elapsed);
			return;
		}
		const auto faces = result_faces(joined);
		TopoDS_Edge canonical = TopoDS::Edge(model.edges(
			static_cast<Standard_Integer>(contact.edge + 1)));
		const TopoDS_Shape modified_edge = sewing.ModifiedSubShape(canonical);
		if (!modified_edge.IsNull() && modified_edge.ShapeType() == TopAbs_EDGE)
			canonical = TopoDS::Edge(modified_edge);

		TopTools_IndexedMapOfShape result_face_map;
		TopExp::MapShapes(joined, TopAbs_FACE, result_face_map);
		std::vector<std::set<std::size_t>> parents(static_cast<std::size_t>(
			result_face_map.Extent()));
		std::size_t missing = 0;
		for (std::size_t source = 0; source < sources.size(); ++source)
		{
			TopoDS_Shape image = sewing.ModifiedSubShape(sources[source]);
			if (image.IsNull()) image = sources[source];
			bool found = false;
			for (TopExp_Explorer face(image, TopAbs_FACE); face.More(); face.Next())
			{
				const Standard_Integer index = result_face_map.FindIndex(face.Current());
				if (index <= 0) continue;
				parents[static_cast<std::size_t>(index - 1)].insert(source);
				found = true;
			}
			if (!found && result_face_map.FindIndex(sources[source]) > 0)
			{
				parents[static_cast<std::size_t>(result_face_map.FindIndex(sources[source]) - 1)]
					.insert(source);
				found = true;
			}
			if (!found) ++missing;
		}
		Provenance history;
		history.result_faces = faces.size();
		history.missing_source = missing;
		for (const auto& parent : parents)
			if (parent.empty()) ++history.orphan;
			else if (parent.size() == 1) ++history.one_parent;
			else ++history.multi_parent;
		print_transaction("contact-group-sew", joined, history,
			edge_face_incidence(canonical, faces), contact.expected_fan,
			elapsed, options.deflection_mm);
		double output_tolerance = 0.0;
		for (TopExp_Explorer edge(joined, TopAbs_EDGE); edge.More(); edge.Next())
			output_tolerance = std::max(output_tolerance,
				BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
		std::printf("      native input/output max edge tolerance %.9g / %.9g mm; "
			"sewing free/contiguous/multiple=%d/%d/%d\n", tolerance, output_tolerance,
			sewing.NbFreeEdges(), sewing.NbContigousEdges(), sewing.NbMultipleEdges());
	}

	void run_iterative_contact_sew(const Model& model, const ContactGroup& contact,
		const Options& options)
	{
		TopoDS_Shape joined = model.faces[contact.owner_face];
		std::vector<TopoDS_Face> source_images{ model.faces[contact.owner_face] };
		TopoDS_Edge canonical = TopoDS::Edge(model.edges(
			static_cast<Standard_Integer>(contact.edge + 1)));
		double maximum_input_tolerance = 0.0;
		const auto begin = Clock::now();
		for (const Contact& target_contact : contact.targets)
		{
			const TopoDS_Face target = model.faces[target_contact.target_face];
			double tolerance = Precision::Confusion();
			for (TopExp_Explorer edge(joined, TopAbs_EDGE); edge.More(); edge.Next())
				tolerance = std::max(tolerance,
					BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
			for (TopExp_Explorer edge(target, TopAbs_EDGE); edge.More(); edge.Next())
				tolerance = std::max(tolerance,
					BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
			maximum_input_tolerance = std::max(maximum_input_tolerance, tolerance);
			BRepBuilderAPI_Sewing sewing(tolerance, true, true, true, true);
			sewing.SetLocalTolerancesMode(true);
			sewing.Add(joined);
			sewing.Add(target);
			sewing.Perform();
			if (sewing.SewedShape().IsNull())
			{
				std::printf("    iterative-sew      FAILED while adding target %zu\n",
					target_contact.target_face);
				return;
			}
			for (TopoDS_Face& image : source_images)
			{
				const TopoDS_Shape modified = sewing.ModifiedSubShape(image);
				if (!modified.IsNull() && modified.ShapeType() == TopAbs_FACE)
					image = TopoDS::Face(modified);
			}
			TopoDS_Shape target_image = sewing.ModifiedSubShape(target);
			if (target_image.IsNull()) target_image = target;
			if (target_image.ShapeType() != TopAbs_FACE)
			{
				std::printf("    iterative-sew      target %zu lost one-face provenance\n",
					target_contact.target_face);
				return;
			}
			source_images.push_back(TopoDS::Face(target_image));
			const TopoDS_Shape modified_edge = sewing.ModifiedSubShape(canonical);
			if (!modified_edge.IsNull() && modified_edge.ShapeType() == TopAbs_EDGE)
				canonical = TopoDS::Edge(modified_edge);
			joined = sewing.SewedShape();
		}
		const double elapsed = std::chrono::duration<double, std::milli>(
			Clock::now() - begin).count();
		const auto faces = result_faces(joined);
		TopTools_IndexedMapOfShape face_map;
		TopExp::MapShapes(joined, TopAbs_FACE, face_map);
		std::vector<std::set<std::size_t>> parents(static_cast<std::size_t>(face_map.Extent()));
		Provenance history;
		history.result_faces = faces.size();
		for (std::size_t source = 0; source < source_images.size(); ++source)
		{
			const Standard_Integer index = face_map.FindIndex(source_images[source]);
			if (index <= 0) ++history.missing_source;
			else parents[static_cast<std::size_t>(index - 1)].insert(source);
		}
		for (const auto& parent : parents)
			if (parent.empty()) ++history.orphan;
			else if (parent.size() == 1) ++history.one_parent;
			else ++history.multi_parent;
		print_transaction("iterative-sew", joined, history,
			edge_face_incidence(canonical, faces), contact.expected_fan,
			elapsed, options.deflection_mm);
		double output_tolerance = 0.0;
		for (TopExp_Explorer edge(joined, TopAbs_EDGE); edge.More(); edge.Next())
			output_tolerance = std::max(output_tolerance,
				BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
		std::printf("      iterative native max input/output edge tolerance %.9g / %.9g mm\n",
			maximum_input_tolerance, output_tolerance);
	}

	void run_pair(const Model& model, const Contact& contact, const Options& options)
	{
		const TopoDS_Edge edge = TopoDS::Edge(model.edges(
			static_cast<Standard_Integer>(contact.edge + 1)));
		const std::vector<TopoDS_Face> sources{
			model.faces[contact.owner_face], model.faces[contact.target_face] };
		run_local_sew(sources[0], sources[1], edge, options);
		run_boundary_relink(sources[0], sources[1], edge, options);

		NCollection_List<TopoDS_Shape> arguments;
		arguments.Append(sources[0]); arguments.Append(sources[1]);
		BOPAlgo_Splitter pair;
		pair.SetArguments(arguments);
		pair.SetNonDestructive(true);
		pair.SetRunParallel(false);
		pair.SetUseOBB(true);
		pair.SetFuzzyValue(0.0);
		const auto begin = Clock::now(); pair.Perform();
		const double elapsed = std::chrono::duration<double, std::milli>(Clock::now() - begin).count();
		if (pair.HasErrors() || pair.Shape().IsNull())
		{
			std::printf("    pair-fuse          FAILED after %.2f ms\n", elapsed);
			return;
		}
		const auto faces = result_faces(pair.Shape());
		print_transaction("pair-fuse", pair.Shape(), provenance(pair, sources, pair.Shape()),
			reproduced_fan(pair, edge, pair.Shape(), faces), contact.expected_fan,
			elapsed, options.deflection_mm);

		NCollection_List<TopoDS_Shape> target_arguments, tools;
		target_arguments.Append(sources[1]); tools.Append(edge);
		BOPAlgo_Splitter target;
		target.SetArguments(target_arguments);
		target.SetTools(tools);
		target.SetNonDestructive(true);
		target.SetRunParallel(false);
		target.SetUseOBB(true);
		target.SetFuzzyValue(0.0);
		const auto target_begin = Clock::now(); target.Perform();
		const double target_elapsed = std::chrono::duration<double, std::milli>(
			Clock::now() - target_begin).count();
		if (target.HasErrors() || target.Shape().IsNull())
		{
			std::printf("    edge-on-target     FAILED after %.2f ms\n", target_elapsed);
			return;
		}
		std::vector<TopoDS_Shape> combined_shapes{ sources[0] };
		for (const TopoDS_Face& face : result_faces(target.Shape())) combined_shapes.push_back(face);
		const TopoDS_Shape combined = compound(combined_shapes);
		const auto combined_faces = result_faces(combined);
		std::vector<TopoDS_Face> target_source{ sources[1] };
		print_transaction("edge-on-target", combined,
			provenance(target, target_source, target.Shape()),
			reproduced_fan(target, edge, combined, combined_faces), contact.expected_fan,
			target_elapsed, options.deflection_mm);

		// A target-only split creates the required interior edge but intentionally does
		// not make it the same topological edge as the owner.  Join only this already
		// exact-certified local transaction at the source shapes' native tolerance.
		// This is the candidate production primitive for an interior T-junction: split
		// the target, then unify the coincident edges without any fuzzy gap repair.
		const std::vector<TopoDS_Face> target_pieces =
			face_images(target, sources[1], target.Shape());
		double native_tolerance = Precision::Confusion();
		auto accumulate_tolerance = [&](const TopoDS_Shape& shape)
		{
			for (TopExp_Explorer exp(shape, TopAbs_EDGE); exp.More(); exp.Next())
				native_tolerance = std::max(native_tolerance,
					BRep_Tool::Tolerance(TopoDS::Edge(exp.Current())));
		};
		accumulate_tolerance(sources[0]);
		for (const TopoDS_Face& piece : target_pieces) accumulate_tolerance(piece);
		BRepBuilderAPI_Sewing split_sewing(native_tolerance, true, true, true, true);
		split_sewing.SetLocalTolerancesMode(true);
		split_sewing.Add(sources[0]);
		for (const TopoDS_Face& piece : target_pieces) split_sewing.Add(piece);
		const auto split_sew_begin = Clock::now();
		split_sewing.Perform();
		const double split_sew_elapsed = target_elapsed
			+ std::chrono::duration<double, std::milli>(Clock::now() - split_sew_begin).count();
		const TopoDS_Shape split_joined = split_sewing.SewedShape();
		if (split_joined.IsNull())
		{
			std::printf("    split-then-sew     FAILED after %.2f ms\n", split_sew_elapsed);
			return;
		}

		TopTools_IndexedMapOfShape split_face_map;
		TopExp::MapShapes(split_joined, TopAbs_FACE, split_face_map);
		std::vector<std::set<std::size_t>> split_parents(
			static_cast<std::size_t>(split_face_map.Extent()));
		Provenance split_history;
		split_history.result_faces = static_cast<std::size_t>(split_face_map.Extent());
		auto map_source_image = [&](const TopoDS_Face& input, std::size_t source)
		{
			TopoDS_Shape image = split_sewing.ModifiedSubShape(input);
			if (image.IsNull()) image = input;
			bool found = false;
			for (TopExp_Explorer exp(image, TopAbs_FACE); exp.More(); exp.Next())
			{
				const Standard_Integer index = split_face_map.FindIndex(exp.Current());
				if (index <= 0) continue;
				split_parents[static_cast<std::size_t>(index - 1)].insert(source);
				found = true;
			}
			if (!found)
			{
				const Standard_Integer index = split_face_map.FindIndex(input);
				if (index > 0)
				{
					split_parents[static_cast<std::size_t>(index - 1)].insert(source);
					found = true;
				}
			}
			if (!found) ++split_history.missing_source;
		};
		map_source_image(sources[0], 0);
		for (const TopoDS_Face& piece : target_pieces) map_source_image(piece, 1);
		for (const auto& parents : split_parents)
			if (parents.empty()) ++split_history.orphan;
			else if (parents.size() == 1) ++split_history.one_parent;
			else ++split_history.multi_parent;

		const auto split_faces = result_faces(split_joined);
		std::size_t matching_edges = 0;
		const std::size_t split_fan = geometric_edge_fan(edge, split_joined,
			split_faces, &matching_edges);
		print_transaction("split-then-sew", split_joined, split_history,
			split_fan, contact.expected_fan,
			split_sew_elapsed, options.deflection_mm);
		double output_tolerance = 0.0;
		for (TopExp_Explorer exp(split_joined, TopAbs_EDGE); exp.More(); exp.Next())
			output_tolerance = std::max(output_tolerance,
				BRep_Tool::Tolerance(TopoDS::Edge(exp.Current())));
		std::printf("      split/native max input/output edge tolerance %.9g / %.9g mm; "
			"geometric edge images=%zu\n", native_tolerance, output_tolerance, matching_edges);
	}

	void run_cell(const Model& model, const ContactGroup& contact, const Options& options)
	{
		const TopoDS_Edge edge = TopoDS::Edge(model.edges(
			static_cast<Standard_Integer>(contact.edge + 1)));
		BRepAdaptor_Curve curve(edge);
		const gp_Pnt centre = curve.Value(0.5 * (curve.FirstParameter() + curve.LastParameter()));
		const double half = 0.5 * options.cell_mm;
		const gp_Pnt lower(centre.X() - half, centre.Y() - half, centre.Z() - half);
		const gp_Pnt upper(centre.X() + half, centre.Y() + half, centre.Z() + half);
		BRepPrimAPI_MakeBox make_box(lower, upper);
		const TopoDS_Solid box = make_box.Solid();
		const Bnd_Box box_bounds = bounds(box);
		std::vector<TopoDS_Face> sources;
		std::vector<std::size_t> source_ids;
		for (std::size_t face = 0; face < model.faces.size(); ++face)
			if (!options.excluded_faces.contains(face) && !box_bounds.IsOut(model.face_bounds[face]))
			{
				sources.push_back(model.faces[face]);
				source_ids.push_back(face);
			}

		NCollection_List<TopoDS_Shape> arguments;
		arguments.Append(box);
		for (const TopoDS_Face& face : sources) arguments.Append(face);
		BOPAlgo_CellsBuilder builder;
		builder.SetArguments(arguments);
		builder.SetNonDestructive(true);
		builder.SetRunParallel(false);
		builder.SetUseOBB(true);
		builder.SetFuzzyValue(0.0);
		const auto begin = Clock::now(); builder.Perform();
		if (!builder.HasErrors()) builder.AddAllToResult();
		const double elapsed = std::chrono::duration<double, std::milli>(Clock::now() - begin).count();
		if (builder.HasErrors() || builder.Shape().IsNull())
		{
			std::printf("    exact-BRep cell     FAILED after %.2f ms (%zu candidate faces)\n",
				elapsed, sources.size());
			return;
		}

		TopTools_MapOfShape solid_faces;
		std::size_t solids = 0;
		for (TopExp_Explorer solid(builder.Shape(), TopAbs_SOLID); solid.More(); solid.Next())
		{
			++solids;
			for (TopExp_Explorer face(solid.Current(), TopAbs_FACE); face.More(); face.Next())
				solid_faces.Add(face.Current());
		}
		std::vector<TopoDS_Face> retained;
		TopTools_MapOfShape seen;
		for (const TopoDS_Face& source : sources)
			for (const TopoDS_Face& image : face_images(builder, source, builder.Shape()))
				if (solid_faces.Contains(image) && seen.Add(image)) retained.push_back(image);

		const Provenance history = provenance(builder, sources, builder.Shape(), &retained);
		const std::size_t fan = reproduced_fan(builder, edge, builder.Shape(), retained);
		std::vector<TopoDS_Shape> retained_shapes;
		retained_shapes.reserve(retained.size());
		for (const TopoDS_Face& face : retained) retained_shapes.push_back(face);
		const TopoDS_Shape fabric = compound(retained_shapes);
		const Tessellation mesh = tessellate(fabric, options.deflection_mm);
		std::printf("    exact-BRep cell     %7.2f ms valid=%s candidates=%zu solids=%zu retained-fabric=%zu "
			"provenance 1/multi/orphan/missing=%zu/%zu/%zu/%zu fan=%zu/%u "
			"tess tris=%zu shared conform/non=%zu/%zu\n", elapsed,
			BRepCheck_Analyzer(builder.Shape(), true).IsValid() ? "yes" : "NO",
			sources.size(), solids, retained.size(), history.one_parent, history.multi_parent,
			history.orphan, history.missing_source, fan, contact.expected_fan,
			mesh.triangles, mesh.conforming_edges, mesh.nonconforming_edges);
	}

	void run_global(const Model& model, const Options& options)
	{
		std::vector<TopoDS_Face> sources;
		for (std::size_t face = 0; face < model.faces.size(); ++face)
			if (!options.excluded_faces.contains(face)) sources.push_back(model.faces[face]);
		NCollection_List<TopoDS_Shape> arguments;
		for (const TopoDS_Face& face : sources) arguments.Append(face);
		BOPAlgo_Splitter splitter;
		splitter.SetArguments(arguments);
		splitter.SetNonDestructive(true);
		splitter.SetRunParallel(false);
		splitter.SetUseOBB(true);
		splitter.SetFuzzyValue(0.0);
		const auto begin = Clock::now(); splitter.Perform();
		const double elapsed = std::chrono::duration<double, std::milli>(Clock::now() - begin).count();
		if (splitter.HasErrors() || splitter.Shape().IsNull())
		{
			std::printf("\n[global %zu-face split] FAILED after %.2f ms\n", sources.size(), elapsed);
			return;
		}
		const Provenance history = provenance(splitter, sources, splitter.Shape());
		const Tessellation mesh = tessellate(splitter.Shape(), options.deflection_mm);
		std::printf("\n[global %zu-face split] %.2f ms valid=%s result-faces=%zu "
			"provenance 1/multi/orphan/missing=%zu/%zu/%zu/%zu; "
			"tess tris=%zu shared conform/non=%zu/%zu\n", sources.size(), elapsed,
			BRepCheck_Analyzer(splitter.Shape(), true).IsValid() ? "yes" : "NO",
			history.result_faces, history.one_parent, history.multi_parent, history.orphan,
			history.missing_source, mesh.triangles, mesh.conforming_edges,
			mesh.nonconforming_edges);
	}

	void run_global_native_sew(const Model& model, const Options& options)
	{
		Options all_options = options;
		all_options.contact_offset = 0;
		all_options.max_contacts = std::numeric_limits<std::size_t>::max();
		const auto contact_begin = Clock::now();
		const std::vector<ContactGroup> contacts = find_contacts(model, all_options);
		const double contact_seconds = std::chrono::duration<double>(
			Clock::now() - contact_begin).count();
		std::vector<TopoDS_Face> sources;
		std::vector<std::size_t> source_ids;
		double tolerance = Precision::Confusion();
		for (std::size_t face = 0; face < model.faces.size(); ++face)
		{
			if (options.excluded_faces.contains(face)) continue;
			sources.push_back(model.faces[face]);
			source_ids.push_back(face);
			for (TopExp_Explorer edge(model.faces[face], TopAbs_EDGE); edge.More(); edge.Next())
				tolerance = std::max(tolerance,
					BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
		}
		BRepBuilderAPI_Sewing sewing(tolerance, true, true, true, true);
		sewing.SetLocalTolerancesMode(true);
		for (const TopoDS_Face& face : sources) sewing.Add(face);
		const auto begin = Clock::now(); sewing.Perform();
		const double elapsed = std::chrono::duration<double, std::milli>(
			Clock::now() - begin).count();
		const TopoDS_Shape result = sewing.SewedShape();
		if (result.IsNull())
		{
			std::printf("\n[global native sewing] FAILED after %.2f ms\n", elapsed);
			return;
		}
		const auto faces = result_faces(result);
		TopTools_IndexedMapOfShape result_face_map;
		TopExp::MapShapes(result, TopAbs_FACE, result_face_map);
		std::vector<std::set<std::size_t>> parents(static_cast<std::size_t>(
			result_face_map.Extent()));
		Provenance history;
		history.result_faces = faces.size();
		for (std::size_t source = 0; source < sources.size(); ++source)
		{
			TopoDS_Shape image = sewing.ModifiedSubShape(sources[source]);
			if (image.IsNull()) image = sources[source];
			bool found = false;
			for (TopExp_Explorer face(image, TopAbs_FACE); face.More(); face.Next())
			{
				const Standard_Integer index = result_face_map.FindIndex(face.Current());
				if (index <= 0) continue;
				parents[static_cast<std::size_t>(index - 1)].insert(source_ids[source]);
				found = true;
			}
			if (!found && result_face_map.FindIndex(sources[source]) > 0)
			{
				parents[static_cast<std::size_t>(result_face_map.FindIndex(sources[source]) - 1)]
					.insert(source_ids[source]);
				found = true;
			}
			if (!found) ++history.missing_source;
		}
		for (const auto& parent : parents)
			if (parent.empty()) ++history.orphan;
			else if (parent.size() == 1) ++history.one_parent;
			else ++history.multi_parent;

		std::vector<std::uint32_t> expected(model.edge_owners.size(), 0);
		for (const ContactGroup& contact : contacts)
			expected[contact.edge] = contact.expected_fan;
		std::size_t certified_ok = 0, certified_bad = 0, free_ok = 0, false_join = 0;
		std::size_t mapping_failures = 0, printed = 0;
		for (std::size_t edge_id = 0; edge_id < model.edge_owners.size(); ++edge_id)
		{
			const auto& owners = model.edge_owners[edge_id];
			if (owners.size() != 1 || options.excluded_faces.contains(owners.front())) continue;
			const TopoDS_Edge source_edge = TopoDS::Edge(model.edges(
				static_cast<Standard_Integer>(edge_id + 1)));
			TopoDS_Shape image = sewing.ModifiedSubShape(source_edge);
			if (image.IsNull()) image = source_edge;
			if (image.ShapeType() != TopAbs_EDGE)
			{
				++mapping_failures;
				continue;
			}
			const std::size_t fan = edge_face_incidence(TopoDS::Edge(image), faces);
			if (expected[edge_id] > 1)
			{
				if (fan == expected[edge_id]) ++certified_ok;
				else
				{
					++certified_bad;
					if (printed++ < 12) std::printf("  unreproduced certified edge %zu fan=%zu/%u\n",
						edge_id, fan, expected[edge_id]);
				}
			}
			else if (fan == 1) ++free_ok;
			else
			{
				++false_join;
				if (printed++ < 12) std::printf("  uncertified free edge %zu was joined with fan=%zu\n",
					edge_id, fan);
			}
		}
		const Tessellation mesh = tessellate(result, options.deflection_mm);
		double output_tolerance = 0.0;
		for (TopExp_Explorer edge(result, TopAbs_EDGE); edge.More(); edge.Next())
			output_tolerance = std::max(output_tolerance,
				BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
		std::printf("\n[global exact-native sewing audit] contacts=%zu scan=%.3f s sew=%.2f ms "
			"valid=%s faces=%zu provenance 1/multi/orphan/missing=%zu/%zu/%zu/%zu\n",
			contacts.size(), contact_seconds, elapsed,
			BRepCheck_Analyzer(result, true).IsValid() ? "yes" : "NO", faces.size(),
			history.one_parent, history.multi_parent, history.orphan, history.missing_source);
		std::printf("  certified fan ok/bad=%zu/%zu; untouched free ok/false-joined=%zu/%zu; "
			"edge-map failures=%zu\n", certified_ok, certified_bad, free_ok, false_join,
			mapping_failures);
		std::printf("  native input/output max tolerance %.9g / %.9g mm; "
			"tess tris=%zu shared conform/non=%zu/%zu; tool contiguous/multiple=%d/%d\n",
			tolerance, output_tolerance, mesh.triangles, mesh.conforming_edges,
			mesh.nonconforming_edges, sewing.NbContigousEdges(), sewing.NbMultipleEdges());
	}

	void run_growing_native_sew(const Model& model, const Options& options)
	{
		Options all_options = options;
		all_options.contact_offset = 0;
		all_options.max_contacts = std::numeric_limits<std::size_t>::max();
		const auto contact_begin = Clock::now();
		const std::vector<ContactGroup> contacts = find_contacts(model, all_options);
		const double contact_seconds = std::chrono::duration<double>(
			Clock::now() - contact_begin).count();
		std::vector<std::vector<std::size_t>> adjacency(model.faces.size());
		auto join = [&](std::size_t a, std::size_t b)
		{
			if (a == b || options.excluded_faces.contains(a)
				|| options.excluded_faces.contains(b)) return;
			adjacency[a].push_back(b);
			adjacency[b].push_back(a);
		};
		for (const auto& owners : model.edge_owners)
		{
			std::vector<std::size_t> unique = owners;
			std::sort(unique.begin(), unique.end());
			unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
			for (std::size_t i = 1; i < unique.size(); ++i) join(unique[0], unique[i]);
		}
		for (const ContactGroup& contact : contacts)
			for (const Contact& target : contact.targets)
				join(contact.owner_face, target.target_face);
		for (auto& neighbors : adjacency)
		{
			std::sort(neighbors.begin(), neighbors.end());
			neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
		}

		std::vector<TopoDS_Face> face_images = model.faces;
		std::vector<TopoDS_Edge> edge_images;
		edge_images.reserve(static_cast<std::size_t>(model.edges.Extent()));
		for (Standard_Integer edge = 1; edge <= model.edges.Extent(); ++edge)
			edge_images.push_back(TopoDS::Edge(model.edges(edge)));
		std::vector<std::uint8_t> visited(model.faces.size(), 0), active_edge(edge_images.size(), 0);
		std::vector<TopoDS_Shape> components;
		std::size_t transactions = 0, map_failures = 0;
		double maximum_input_tolerance = 0.0;
		const auto begin = Clock::now();
		for (std::size_t seed = 0; seed < model.faces.size(); ++seed)
		{
			if (visited[seed] || options.excluded_faces.contains(seed)) continue;
			std::vector<std::size_t> order, queue{ seed };
			visited[seed] = 1;
			for (std::size_t cursor = 0; cursor < queue.size(); ++cursor)
			{
				const std::size_t face = queue[cursor];
				order.push_back(face);
				for (std::size_t neighbor : adjacency[face])
					if (!visited[neighbor]) { visited[neighbor] = 1; queue.push_back(neighbor); }
			}
			TopoDS_Shape current = face_images[order.front()];
			std::vector<std::size_t> added{ order.front() };
			for (TopExp_Explorer edge(model.faces[order.front()], TopAbs_EDGE);
				edge.More(); edge.Next())
			{
				const Standard_Integer id = model.edges.FindIndex(edge.Current());
				if (id > 0) active_edge[static_cast<std::size_t>(id - 1)] = 1;
			}
			for (std::size_t order_index = 1; order_index < order.size(); ++order_index)
			{
				const std::size_t next = order[order_index];
				double tolerance = Precision::Confusion();
				for (TopExp_Explorer edge(current, TopAbs_EDGE); edge.More(); edge.Next())
					tolerance = std::max(tolerance,
						BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
				for (TopExp_Explorer edge(face_images[next], TopAbs_EDGE); edge.More(); edge.Next())
					tolerance = std::max(tolerance,
						BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
				maximum_input_tolerance = std::max(maximum_input_tolerance, tolerance);
				BRepBuilderAPI_Sewing sewing(tolerance, true, true, true, true);
				sewing.SetLocalTolerancesMode(true);
				sewing.Add(current);
				sewing.Add(face_images[next]);
				sewing.Perform();
				if (sewing.SewedShape().IsNull())
				{
					std::printf("\n[growing native sewing] FAILED adding source face %zu\n", next);
					return;
				}
				for (std::size_t face : added)
				{
					const TopoDS_Shape image = sewing.ModifiedSubShape(face_images[face]);
					if (!image.IsNull() && image.ShapeType() == TopAbs_FACE)
						face_images[face] = TopoDS::Face(image);
				}
				const TopoDS_Shape next_image = sewing.ModifiedSubShape(face_images[next]);
				if (!next_image.IsNull() && next_image.ShapeType() == TopAbs_FACE)
					face_images[next] = TopoDS::Face(next_image);
				for (TopExp_Explorer edge(model.faces[next], TopAbs_EDGE); edge.More(); edge.Next())
				{
					const Standard_Integer id = model.edges.FindIndex(edge.Current());
					if (id > 0) active_edge[static_cast<std::size_t>(id - 1)] = 1;
				}
				for (std::size_t edge = 0; edge < edge_images.size(); ++edge)
				{
					if (!active_edge[edge]) continue;
					const TopoDS_Shape image = sewing.ModifiedSubShape(edge_images[edge]);
					if (image.IsNull()) continue;
					if (image.ShapeType() == TopAbs_EDGE) edge_images[edge] = TopoDS::Edge(image);
					else ++map_failures;
				}
				current = sewing.SewedShape();
				added.push_back(next);
				++transactions;
			}
			components.push_back(current);
			for (std::size_t edge = 0; edge < active_edge.size(); ++edge) active_edge[edge] = 0;
		}
		const TopoDS_Shape result = compound(components);
		const double elapsed = std::chrono::duration<double, std::milli>(
			Clock::now() - begin).count();
		const auto faces = result_faces(result);
		TopTools_IndexedMapOfShape result_face_map;
		TopExp::MapShapes(result, TopAbs_FACE, result_face_map);
		Provenance history;
		history.result_faces = faces.size();
		std::vector<std::set<std::size_t>> parents(static_cast<std::size_t>(
			result_face_map.Extent()));
		for (std::size_t source = 0; source < face_images.size(); ++source)
		{
			if (options.excluded_faces.contains(source)) continue;
			const Standard_Integer index = result_face_map.FindIndex(face_images[source]);
			if (index <= 0) ++history.missing_source;
			else parents[static_cast<std::size_t>(index - 1)].insert(source);
		}
		for (const auto& parent : parents)
			if (parent.empty()) ++history.orphan;
			else if (parent.size() == 1) ++history.one_parent;
			else ++history.multi_parent;

		std::vector<std::uint32_t> expected(edge_images.size(), 0);
		for (const ContactGroup& contact : contacts) expected[contact.edge] = contact.expected_fan;
		std::size_t certified_ok = 0, certified_bad = 0, free_ok = 0, false_join = 0;
		std::size_t printed = 0;
		for (std::size_t edge = 0; edge < edge_images.size(); ++edge)
		{
			const auto& owners = model.edge_owners[edge];
			if (owners.size() != 1 || options.excluded_faces.contains(owners.front())) continue;
			const std::size_t fan = edge_face_incidence(edge_images[edge], faces);
			if (expected[edge] > 1)
			{
				if (fan == expected[edge]) ++certified_ok;
				else
				{
					++certified_bad;
					if (printed++ < 12) std::printf("  growing unresolved edge %zu fan=%zu/%u\n",
						edge, fan, expected[edge]);
				}
			}
			else if (fan == 1) ++free_ok;
			else
			{
				++false_join;
				if (printed++ < 12) std::printf("  growing false join edge %zu fan=%zu\n",
					edge, fan);
			}
		}
		const Tessellation mesh = tessellate(result, options.deflection_mm);
		double output_tolerance = 0.0;
		for (TopExp_Explorer edge(result, TopAbs_EDGE); edge.More(); edge.Next())
			output_tolerance = std::max(output_tolerance,
				BRep_Tool::Tolerance(TopoDS::Edge(edge.Current())));
		std::printf("\n[growing exact-contact sewing audit] components=%zu contacts=%zu scan=%.3f s "
			"transactions=%zu time=%.2f ms valid=%s faces=%zu provenance 1/multi/orphan/missing=%zu/%zu/%zu/%zu\n",
			components.size(), contacts.size(), contact_seconds, transactions, elapsed,
			BRepCheck_Analyzer(result, true).IsValid() ? "yes" : "NO", faces.size(),
			history.one_parent, history.multi_parent, history.orphan, history.missing_source);
		std::printf("  certified fan ok/bad=%zu/%zu; untouched free ok/false-joined=%zu/%zu; "
			"edge-map failures=%zu\n", certified_ok, certified_bad, free_ok, false_join,
			map_failures);
		std::printf("  native max input/output tolerance %.9g / %.9g mm; tess tris=%zu "
			"shared conform/non=%zu/%zu\n", maximum_input_tolerance, output_tolerance,
			mesh.triangles, mesh.conforming_edges, mesh.nonconforming_edges);
	}

	// True when every positive-length portion of `candidate` is an exact OCCT
	// edge/edge common with `support`.  This is intentionally one-way: a batch split
	// may divide one certified contact curve into several result edges.
	bool exact_subedge_of(const TopoDS_Edge& candidate, const TopoDS_Edge& support)
	{
		BRepAdaptor_Curve candidate_curve(candidate), support_curve(support);
		const double candidate_first = candidate_curve.FirstParameter();
		const double candidate_last = candidate_curve.LastParameter();
		const double support_first = support_curve.FirstParameter();
		const double support_last = support_curve.LastParameter();
		if (!std::isfinite(candidate_first) || !std::isfinite(candidate_last)
			|| !(candidate_last > candidate_first) || !std::isfinite(support_first)
			|| !std::isfinite(support_last) || !(support_last > support_first)) return false;
		IntTools_EdgeEdge intersection(candidate, support);
		intersection.SetFuzzyValue(0.0);
		intersection.Perform();
		if (!intersection.IsDone()) return false;
		std::vector<std::pair<double, double>> candidate_ranges;
		for (NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
			intersection.CommonParts()); common.More(); common.Next())
		{
			if (common.Value().Type() != TopAbs_EDGE) continue;
			double lo = 0.0, hi = 0.0;
			common.Value().Range1(lo, hi);
			candidate_ranges.emplace_back(lo, hi);
		}
		return ranges_cover(std::move(candidate_ranges), candidate_first, candidate_last);
	}

	struct BatchContactCurve
	{
		std::size_t source_edge = 0;
		TopoDS_Edge edge;
		// A boundary face contributes one half-sheet sector.  A curve in the
		// interior of a face contributes two.  max() during reduction removes the
		// reciprocal reports made by independently modelled face boundaries.
		std::map<std::size_t, std::uint32_t> face_sectors;
	};

	struct BatchFaceImage
	{
		TopoDS_Face face;
		std::size_t source_face = 0;
	};

	using QuantizedPoint = std::array<std::int64_t, 3>;
	using QuantizedSegment = std::array<std::int64_t, 6>;

	QuantizedPoint quantized_point(const gp_Pnt& point, double quantum)
	{
		return { static_cast<std::int64_t>(std::llround(point.X() / quantum)),
			static_cast<std::int64_t>(std::llround(point.Y() / quantum)),
			static_cast<std::int64_t>(std::llround(point.Z() / quantum)) };
	}

	QuantizedSegment quantized_segment(gp_Pnt a, gp_Pnt b, double quantum)
	{
		QuantizedPoint qa = quantized_point(a, quantum), qb = quantized_point(b, quantum);
		if (qb < qa) std::swap(qa, qb);
		return { qa[0], qa[1], qa[2], qb[0], qb[1], qb[2] };
	}

	std::array<std::int64_t, 6> edge_endpoint_key(const TopoDS_Edge& edge)
	{
		BRepAdaptor_Curve curve(edge);
		gp_Pnt a = curve.Value(curve.FirstParameter());
		gp_Pnt b = curve.Value(curve.LastParameter());
		// Broad-phase only.  exact_same_edge_extent() remains the decision, so this
		// key cannot join nearby fabric or close a represented opening.
		return quantized_segment(a, b, 1.0e-5);
	}

	void run_batch_target_prototype(const Model& model, const Options& options)
	{
		Options all_options = options;
		all_options.contact_offset = 0;
		all_options.max_contacts = std::numeric_limits<std::size_t>::max();
		const auto contact_begin = Clock::now();
		const std::vector<ContactGroup> raw_contacts = find_contacts(model, all_options);
		const double contact_seconds = std::chrono::duration<double>(
			Clock::now() - contact_begin).count();

		std::vector<BatchContactCurve> curves;
		std::map<std::array<std::int64_t, 6>, std::vector<std::size_t>> endpoint_index;
		std::vector<std::set<std::size_t>> target_curves(model.faces.size());
		std::size_t raw_pairs = 0;
		for (const ContactGroup& contact : raw_contacts)
		{
			const TopoDS_Edge edge = TopoDS::Edge(model.edges(
				static_cast<Standard_Integer>(contact.edge + 1)));
			const auto key = edge_endpoint_key(edge);
			std::size_t curve_id = std::numeric_limits<std::size_t>::max();
			for (std::size_t candidate : endpoint_index[key])
				if (exact_same_edge_extent(edge, curves[candidate].edge))
				{
					curve_id = candidate;
					break;
				}
			if (curve_id == std::numeric_limits<std::size_t>::max())
			{
				curve_id = curves.size();
				curves.push_back({ contact.edge, edge, {} });
				endpoint_index[key].push_back(curve_id);
			}
			BatchContactCurve& curve = curves[curve_id];
			curve.face_sectors[contact.owner_face] = std::max<std::uint32_t>(
				curve.face_sectors[contact.owner_face], 1u);
			for (const Contact& target : contact.targets)
			{
				++raw_pairs;
				curve.face_sectors[target.target_face] = std::max(
					curve.face_sectors[target.target_face], target.target_sectors);
				target_curves[target.target_face].insert(curve_id);
			}
		}

		std::printf("\n[batch target-face contact prototype]\n");
		std::printf("  complete exact scan %.3f s: raw groups/pairs=%zu/%zu; "
			"canonical full-extent curves=%zu\n", contact_seconds, raw_contacts.size(),
			raw_pairs, curves.size());

		std::vector<BatchFaceImage> images;
		images.reserve(model.faces.size() + curves.size());
		std::size_t transactions = 0, succeeded = 0, warned = 0, failed = 0;
		std::size_t tools_total = 0, input_faces = 0, output_faces = 0;
		double split_milliseconds = 0.0;
		for (std::size_t face_id = 0; face_id < model.faces.size(); ++face_id)
		{
			if (options.excluded_faces.contains(face_id)) continue;
			++input_faces;
			if (target_curves[face_id].empty())
			{
				images.push_back({ model.faces[face_id], face_id });
				++output_faces;
				continue;
			}

			NCollection_List<TopoDS_Shape> arguments, tools;
			arguments.Append(model.faces[face_id]);
			for (std::size_t curve : target_curves[face_id]) tools.Append(curves[curve].edge);
			tools_total += target_curves[face_id].size();
			BOPAlgo_Splitter splitter;
			splitter.SetArguments(arguments);
			splitter.SetTools(tools);
			splitter.SetNonDestructive(true);
			splitter.SetRunParallel(false);
			splitter.SetUseOBB(true);
			splitter.SetFuzzyValue(0.0);
			const auto begin = Clock::now();
			splitter.Perform();
			const double elapsed = std::chrono::duration<double, std::milli>(
				Clock::now() - begin).count();
			split_milliseconds += elapsed;
			++transactions;
			if (splitter.HasErrors() || splitter.Shape().IsNull())
			{
				++failed;
				std::ostringstream report;
				splitter.DumpErrors(report);
				std::printf("  target face %zu FAILED (%zu tools, %.2f ms): %s\n", face_id,
					target_curves[face_id].size(), elapsed, report.str().c_str());
				// Keep the original only so the final mesh audit can show every missing
				// contact instead of terminating at the first OCCT transaction failure.
				images.push_back({ model.faces[face_id], face_id });
				++output_faces;
				continue;
			}
			++succeeded;
			if (splitter.HasWarnings())
			{
				++warned;
				std::ostringstream report;
				splitter.DumpWarnings(report);
				std::printf("  target face %zu WARNING (%zu tools, %.2f ms): %s\n", face_id,
					target_curves[face_id].size(), elapsed, report.str().c_str());
			}
			std::vector<TopoDS_Face> split_faces = face_images(splitter,
				model.faces[face_id], splitter.Shape());
			if (split_faces.empty())
			{
				++failed;
				--succeeded;
				std::printf("  target face %zu produced no provenance image; retaining source\n",
					face_id);
				images.push_back({ model.faces[face_id], face_id });
				++output_faces;
				continue;
			}
			for (const TopoDS_Face& face : split_faces)
			{
				images.push_back({ face, face_id });
				++output_faces;
			}
		}

		std::vector<TopoDS_Shape> assembled_faces;
		assembled_faces.reserve(images.size());
		for (const BatchFaceImage& image : images) assembled_faces.push_back(image.face);
		TopoDS_Shape assembled = compound(assembled_faces);
		BRepTools::Clean(assembled);
		const auto mesh_begin = Clock::now();
		BRepMesh_IncrementalMesh mesher(assembled, options.deflection_mm, false, 0.1, false);
		const double mesh_seconds = std::chrono::duration<double>(Clock::now() - mesh_begin).count();
		const Tessellation mesh = tessellate(assembled, options.deflection_mm);
		std::printf("  target transactions=%zu success/warning/failure=%zu/%zu/%zu; "
			"tools=%zu; faces %zu -> %zu; split %.2f ms\n", transactions, succeeded,
			warned, failed, tools_total, input_faces, output_faces, split_milliseconds);
		std::printf("  assembled valid=%s; mesh %.3f s; triangles=%zu; "
			"topologically shared conform/non=%zu/%zu\n",
			BRepCheck_Analyzer(assembled, true).IsValid() ? "yes" : "NO", mesh_seconds,
			mesh.triangles, mesh.conforming_edges, mesh.nonconforming_edges);

		std::size_t curves_passed = 0, curves_missing = 0, curves_nonconforming = 0;
		std::size_t boundary_curves = 0, boundary_passed = 0;
		std::size_t interior_curves = 0, interior_passed = 0;
		std::size_t printed = 0;
		for (std::size_t curve_id = 0; curve_id < curves.size(); ++curve_id)
		{
			const BatchContactCurve& curve = curves[curve_id];
			std::uint32_t expected_fan = 0;
			bool has_interior_participant = false;
			for (const auto& [face, sectors] : curve.face_sectors)
				if (!options.excluded_faces.contains(face))
				{
					expected_fan += sectors;
					has_interior_participant = has_interior_participant || sectors > 1;
				}
			if (expected_fan < 2) continue;
			if (has_interior_participant) ++interior_curves;
			else ++boundary_curves;
			const double audit_quantum = std::max(1.0e-6,
				10.0 * std::max(Precision::Confusion(), BRep_Tool::Tolerance(curve.edge)));
			std::map<QuantizedSegment, std::set<std::uint64_t>> segment_occurrences;
			std::size_t edge_uses = 0, missing_polygons = 0;
			for (std::size_t image_id = 0; image_id < images.size(); ++image_id)
			{
				const TopoDS_Face& face = images[image_id].face;
				TopLoc_Location location;
				const Handle(Poly_Triangulation) triangulation =
					BRep_Tool::Triangulation(face, location);
				if (triangulation.IsNull()) continue;
				std::uint32_t edge_use = 0;
				for (TopExp_Explorer exp(face, TopAbs_EDGE); exp.More(); exp.Next(), ++edge_use)
				{
					const TopoDS_Edge edge = TopoDS::Edge(exp.Current());
					if (!exact_subedge_of(edge, curve.edge)) continue;
					++edge_uses;
					const Handle(Poly_PolygonOnTriangulation) polygon =
						BRep_Tool::PolygonOnTriangulation(edge, triangulation, location);
					if (polygon.IsNull() || polygon->NbNodes() < 2)
					{
						++missing_polygons;
						continue;
					}
					const gp_Trsf transform = location.Transformation();
					const std::uint64_t occurrence = (static_cast<std::uint64_t>(image_id) << 32)
						| edge_use;
					for (Standard_Integer node = 1; node < polygon->NbNodes(); ++node)
					{
						gp_Pnt a = triangulation->Node(polygon->Node(node));
						gp_Pnt b = triangulation->Node(polygon->Node(node + 1));
						a.Transform(transform);
						b.Transform(transform);
						if (a.Distance(b) <= audit_quantum) continue;
						segment_occurrences[quantized_segment(a, b, audit_quantum)]
							.insert(occurrence);
					}
				}
			}
			std::size_t below = 0, exact = 0, above = 0;
			std::size_t minimum_fan = std::numeric_limits<std::size_t>::max();
			std::size_t maximum_fan = 0;
			for (const auto& [segment, occurrences] : segment_occurrences)
			{
				const std::size_t fan = occurrences.size();
				minimum_fan = std::min(minimum_fan, fan);
				maximum_fan = std::max(maximum_fan, fan);
				if (fan < expected_fan) ++below;
				else if (fan == expected_fan) ++exact;
				else ++above;
			}
			const bool pass = !segment_occurrences.empty() && missing_polygons == 0
				&& below == 0 && above == 0;
			if (pass)
			{
				++curves_passed;
				if (has_interior_participant) ++interior_passed;
				else ++boundary_passed;
			}
			else if (segment_occurrences.empty()) ++curves_missing;
			else ++curves_nonconforming;
			if (!pass && printed++ < 24)
			{
				std::printf("    curve %zu edge=%zu kind=%s expected-fan=%u faces=%zu edge-uses=%zu "
					"mesh-segments=%zu fan[min,max]=%zu/%zu below/exact/above=%zu/%zu/%zu "
					"missing-chains=%zu\n", curve_id, curve.source_edge,
					has_interior_participant ? "interior" : "boundary", expected_fan,
					curve.face_sectors.size(), edge_uses, segment_occurrences.size(),
					segment_occurrences.empty() ? 0 : minimum_fan, maximum_fan,
					below, exact, above, missing_polygons);
			}
		}
		std::printf("  certified mesh-curve audit pass/missing/nonconforming=%zu/%zu/%zu "
			"(audit coordinate quantum >= 1e-6 mm)\n", curves_passed, curves_missing,
			curves_nonconforming);
		std::printf("  boundary curves passed=%zu/%zu; interior-constraint curves passed=%zu/%zu\n",
			boundary_passed, boundary_curves, interior_passed, interior_curves);
		std::printf("  verdict: %s\n", failed == 0 && curves_missing == 0
			&& curves_nonconforming == 0 ? "BATCH TARGET SPLIT IS SUFFICIENT"
			: "BATCH TARGET SPLIT IS NOT YET A CONFORMING SURFACE MESH");
	}

	double face_area(const TopoDS_Face& face)
	{
		GProp_GProps properties;
		BRepGProp::SurfaceProperties(face, properties);
		return std::abs(properties.Mass());
	}

	std::array<double, 6> numeric_bounds(const TopoDS_Shape& shape)
	{
		std::array<double, 6> result{};
		Bnd_Box box = bounds(shape);
		box.Get(result[0], result[1], result[2], result[3], result[4], result[5]);
		return result;
	}

	bool coplanar(const TopoDS_Face& a, const TopoDS_Face& b)
	{
		BRepAdaptor_Surface surface_a(a, true), surface_b(b, true);
		if (surface_a.GetType() != GeomAbs_Plane || surface_b.GetType() != GeomAbs_Plane)
			return false;
		const gp_Pln plane_a = surface_a.Plane(), plane_b = surface_b.Plane();
		const double alignment = std::abs(plane_a.Axis().Direction().Dot(
			plane_b.Axis().Direction()));
		if (1.0 - alignment > 64.0 * std::numeric_limits<double>::epsilon()) return false;
		const double tolerance = std::max({ BRep_Tool::Tolerance(a), BRep_Tool::Tolerance(b),
			Precision::Confusion() });
		return plane_a.Distance(plane_b.Location()) <= tolerance;
	}

	void print_names(const paracfd::core::StepGeometry& geometry, std::size_t face)
	{
		std::printf(" names=");
		if (face >= geometry.source_faces.size()
			|| geometry.source_faces[face].representation_names.empty())
		{
			std::printf("<none>");
			return;
		}
		for (const std::string& name : geometry.source_faces[face].representation_names)
			std::printf("[%s]", name.c_str());
	}

	void inspect_planar_overlaps(const Model& model, const Options& options)
	{
		std::printf("\n[exact positive-area coplanar overlap audit]\n");
		std::string import_error;
		const paracfd::core::StepGeometry metadata = paracfd::core::load_step_geometry(
			options.step, options.deflection_mm, &import_error);
		if (metadata.mesh.empty())
			std::printf("  metadata import failed: %s\n", import_error.c_str());
		const auto begin = Clock::now();
		std::size_t coplanar_pairs = 0, overlap_pairs = 0;
		double total_overlap = 0.0;
		for (std::size_t a = 0; a < model.faces.size(); ++a)
		{
			if (options.excluded_faces.contains(a)) continue;
			for (std::size_t b = a + 1; b < model.faces.size(); ++b)
			{
				if (options.excluded_faces.contains(b)
					|| model.face_bounds[a].IsOut(model.face_bounds[b])
					|| !coplanar(model.faces[a], model.faces[b])) continue;
				++coplanar_pairs;
				BRepAlgoAPI_Common common(model.faces[a], model.faces[b]);
				common.SetFuzzyValue(0.0);
				common.SetRunParallel(false);
				common.Build();
				if (!common.IsDone() || common.HasErrors() || common.Shape().IsNull()) continue;
				double overlap = 0.0;
				for (TopExp_Explorer face(common.Shape(), TopAbs_FACE); face.More(); face.Next())
					overlap += face_area(TopoDS::Face(face.Current()));
				const double area_tolerance = std::max(Precision::Confusion()
					* std::sqrt(std::max(face_area(model.faces[a]), face_area(model.faces[b]))),
					1024.0 * std::numeric_limits<double>::epsilon()
					* std::max(face_area(model.faces[a]), face_area(model.faces[b])));
				if (!(overlap > area_tolerance)) continue;
				++overlap_pairs;
				total_overlap += overlap;
				const auto box_a = numeric_bounds(model.faces[a]);
				const auto box_b = numeric_bounds(model.faces[b]);
				std::printf("  faces %zu / %zu overlap=%.12g mm^2 (areas %.12g / %.12g)\n",
					a, b, overlap, face_area(model.faces[a]), face_area(model.faces[b]));
				std::printf("    face %zu bbox=[%.6g %.6g %.6g]-[%.6g %.6g %.6g]",
					a, box_a[0], box_a[1], box_a[2], box_a[3], box_a[4], box_a[5]);
				print_names(metadata, a); std::printf("\n");
				std::printf("    face %zu bbox=[%.6g %.6g %.6g]-[%.6g %.6g %.6g]",
					b, box_b[0], box_b[1], box_b[2], box_b[3], box_b[4], box_b[5]);
				print_names(metadata, b); std::printf("\n");
			}
		}
		std::printf("  coplanar candidates=%zu; positive-area overlaps=%zu; summed overlap=%.12g mm^2; %.3f s\n",
			coplanar_pairs, overlap_pairs, total_overlap,
			std::chrono::duration<double>(Clock::now() - begin).count());
	}
}

int main(int argc, char** argv)
{
	Options options;
	if (!parse_options(argc, argv, options)) { usage(argv[0]); return 2; }
	try
	{
		std::printf("OCCT exact CAD-contact imprint diagnostic (no sewing/fuzzy repair)\n");
		std::printf("STEP: %s; excluded faces:", options.step.c_str());
		if (options.excluded_faces.empty()) std::printf(" none");
		for (std::size_t face : options.excluded_faces) std::printf(" %zu", face);
		std::printf("\n");
		const auto read_begin = Clock::now();
		const Model model = make_model(read_step(options.step));
		std::printf("read/model %.3f s; source faces=%zu edges=%d\n",
			std::chrono::duration<double>(Clock::now() - read_begin).count(),
			model.faces.size(), model.edges.Extent());
		const auto contact_begin = Clock::now();
		const std::vector<ContactGroup> contacts = find_contacts(model, options);
		std::printf("contact scan %.3f s; found through requested window=%zu\n",
			std::chrono::duration<double>(Clock::now() - contact_begin).count(), contacts.size());
		if (options.scan_only)
		{
			for (std::size_t index = 0; index < contacts.size(); ++index)
			{
				const ContactGroup& contact = contacts[index];
				std::printf("  group %zu edge=%zu owner=%zu expected=%u targets=", index,
					contact.edge, contact.owner_face, contact.expected_fan);
				for (const Contact& target : contact.targets)
					std::printf(" %zu(+%u)", target.target_face, target.target_sectors);
				std::printf("\n");
			}
			return 0;
		}
		if (contacts.size() <= options.contact_offset)
			throw std::runtime_error("contact scan did not reach requested offset");
		const std::size_t end = std::min(contacts.size(),
			options.contact_offset + options.max_contacts);
		for (std::size_t index = options.contact_offset; index < end; ++index)
		{
			const ContactGroup& contact = contacts[index];
			std::printf("\n[contact %zu] edge=%zu owner-face=%zu targets=", index,
				contact.edge, contact.owner_face);
			for (const Contact& target : contact.targets)
				std::printf(" %zu(+%u)", target.target_face, target.target_sectors);
			std::printf(" aggregate expected fan=%u\n", contact.expected_fan);
			run_local_sew_group(model, contact, options);
			run_iterative_contact_sew(model, contact, options);
			for (const Contact& target : contact.targets) run_pair(model, target, options);
			run_cell(model, contact, options);
		}
		if (options.run_global) run_global(model, options);
		if (options.run_global_native) run_global_native_sew(model, options);
		if (options.run_growing_native) run_growing_native_sew(model, options);
		if (options.run_batch_target_prototype) run_batch_target_prototype(model, options);
		if (options.inspect_overlaps) inspect_planar_overlaps(model, options);
		std::printf("\nInterpretation: acceptance requires valid BRep, one source parent per retained "
			"fabric face, reproduced fan, and zero nonconforming shared edges.\n");
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
