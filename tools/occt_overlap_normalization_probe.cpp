// Fixture driver for the reusable OCCT planar surface-overlap normalizer.

#include "core/geometry/occt_surface_normalization.h"

#include <BRepBuilderAPI_MakeFace.hxx>
#include <Geom_BezierSurface.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <NCollection_Array2.hxx>
#include <NCollection_IndexedMap.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pln.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <gp.hxx>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
	using Clock = std::chrono::steady_clock;
	using ShapeMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
	using paracfd::core::OcctSurfaceNormalizationClusterAudit;
	using paracfd::core::OcctSurfaceNormalizationResult;

	struct Options
	{
		std::string step = "Test-Data/Two-Cells-fixed.stp";
		std::optional<std::size_t> expected_clusters;
		std::optional<std::size_t> expected_pairs;
		bool verbose = false;
		bool manufactured_only = false;
	};

	void usage(const char* executable)
	{
		std::printf("Usage: %s [--step FILE] [--expect-clusters N] "
			"[--expect-overlap-pairs N] [--manufactured-only] [--verbose]\n", executable);
	}

	bool parse_count(const char* text, std::size_t& result)
	{
		if (!text || !*text) return false;
		char* tail = nullptr;
		const unsigned long long value = std::strtoull(text, &tail, 10);
		if (!tail || *tail != '\0') return false;
		result = static_cast<std::size_t>(value);
		return static_cast<unsigned long long>(result) == value;
	}

	bool parse_options(int argc, char** argv, Options& options)
	{
		for (int argument = 1; argument < argc; ++argument)
		{
			const std::string option = argv[argument];
			auto value = [&](const char* name) -> const char*
			{
				if (argument + 1 >= argc)
				{
					std::fprintf(stderr, "%s requires a value\n", name);
					return nullptr;
				}
				return argv[++argument];
			};
			if (option == "--step")
			{
				const char* next = value("--step");
				if (!next) return false;
				options.step = next;
			}
			else if (option == "--expect-clusters")
			{
				const char* next = value("--expect-clusters");
				std::size_t count = 0;
				if (!next || !parse_count(next, count)) return false;
				options.expected_clusters = count;
			}
			else if (option == "--expect-overlap-pairs")
			{
				const char* next = value("--expect-overlap-pairs");
				std::size_t count = 0;
				if (!next || !parse_count(next, count)) return false;
				options.expected_pairs = count;
			}
			else if (option == "--verbose") options.verbose = true;
			else if (option == "--manufactured-only") options.manufactured_only = true;
			else if (option == "--help" || option == "-h")
			{
				usage(argv[0]);
				std::exit(0);
			}
			else
			{
				std::fprintf(stderr, "unknown option: %s\n", option.c_str());
				return false;
			}
		}
		return !options.step.empty();
	}

	TopoDS_Shape read_step(const std::string& path)
	{
		STEPControl_Reader reader;
		if (reader.ReadFile(path.c_str()) != IFSelect_RetDone)
			throw std::runtime_error("STEPControl_Reader::ReadFile failed for '" + path + "'");
		if (reader.TransferRoots() <= 0)
			throw std::runtime_error("STEPControl_Reader::TransferRoots transferred no roots");
		TopoDS_Shape shape = reader.OneShape();
		if (shape.IsNull()) throw std::runtime_error("STEP reader returned a null shape");
		return shape;
	}

	std::vector<TopoDS_Face> unique_faces(const TopoDS_Shape& shape)
	{
		ShapeMap faces;
		TopExp::MapShapes(shape, TopAbs_FACE, faces);
		std::vector<TopoDS_Face> result;
		result.reserve(static_cast<std::size_t>(faces.Extent()));
		for (Standard_Integer face = 1; face <= faces.Extent(); ++face)
			result.push_back(TopoDS::Face(faces(face)));
		return result;
	}

	bool manufactured_curved_overlap_regression()
	{
		NCollection_Array2<gp_Pnt> poles(1, 4, 1, 4);
		for (int u = 1; u <= 4; ++u)
			for (int v = 1; v <= 4; ++v)
				poles(u, v) = gp_Pnt((u - 1) * (10.0 / 3.0),
					(v - 1) * (10.0 / 3.0),
					(u > 1 && u < 4 && v > 1 && v < 4) ? 2.0 : 0.0);
		const occ::handle<Geom_BezierSurface> support_a = new Geom_BezierSurface(poles);
		const occ::handle<Geom_BezierSurface> support_b =
			occ::down_cast<Geom_BezierSurface>(support_a->Copy());
		BRepBuilderAPI_MakeFace make_a(support_a, 0.0, 0.6, 0.0, 1.0,
			Precision::Confusion());
		BRepBuilderAPI_MakeFace make_b(support_b, 0.4, 1.0, 0.0, 1.0,
			Precision::Confusion());
		if (!make_a.IsDone() || !make_b.IsDone())
		{
			std::printf("[manufactured curved preflight] FAIL: Bezier face construction\n");
			return false;
		}
		TopoDS_Face face_a = make_a.Face();
		TopoDS_Face face_b = make_b.Face();
		face_b.Reverse();
		gp_Trsf translation;
		translation.SetTranslation(gp_Vec(17.0, -9.0, 4.0));
		const TopLoc_Location location(translation);
		face_a.Location(location);
		face_b.Location(location);
		const auto begin = Clock::now();
		const OcctSurfaceNormalizationResult generic =
			paracfd::core::normalize_occt_surface_faces(
				{ face_a, face_b });
		const bool has_limitation_warning = std::any_of(generic.warnings.begin(),
			generic.warnings.end(), [](const std::string& warning)
			{ return warning.find("preflight is incomplete") != std::string::npos; });
		const bool generic_visible = generic.valid() && !generic.fully_verified()
			&& generic.normalized_faces.size() == 2
			&& generic.unsupported_overlaps.empty()
			&& generic.stats.nonplanar_unverified_pair_count == 1
			&& generic.stats.nonplanar_exact_common_pair_count == 0
			&& has_limitation_warning;

		const gp_Cylinder cylinder(gp_Ax3(gp::Origin(), gp::DZ()), 10.0);
		BRepBuilderAPI_MakeFace make_cylinder_a(cylinder, 0.0, std::numbers::pi,
			0.0, 20.0);
		BRepBuilderAPI_MakeFace make_cylinder_b(cylinder, 0.25 * std::numbers::pi,
			0.75 * std::numbers::pi, 5.0, 15.0);
		if (!make_cylinder_a.IsDone() || !make_cylinder_b.IsDone())
		{
			std::printf("[manufactured curved preflight] FAIL: cylinder construction\n");
			return false;
		}
		const OcctSurfaceNormalizationResult elementary =
			paracfd::core::normalize_occt_surface_faces(
				{ make_cylinder_a.Face(), make_cylinder_b.Face() });
		const double milliseconds = std::chrono::duration<double, std::milli>(
			Clock::now() - begin).count();
		const double expected_area = 10.0 * (0.5 * std::numbers::pi) * 10.0;
		const bool elementary_caught = !elementary.valid()
			&& elementary.normalized_faces.empty()
			&& elementary.unsupported_overlaps.size() == 1
			&& std::abs(elementary.unsupported_overlaps.front().common_area - expected_area)
				<= 1.0e-10 * expected_area
			&& elementary.stats.nonplanar_exact_common_pair_count == 1
			&& elementary.stats.nonplanar_unverified_pair_count == 0;
		const bool passed = generic_visible && elementary_caught;
		std::printf("[manufactured curved preflight] %s; Bezier unverified=%zu; "
			"elementary common=%.12g mm^2; exact=%zu; %.2f ms\n",
			passed ? "PASS" : "FAIL", generic.stats.nonplanar_unverified_pair_count,
			elementary.unsupported_overlaps.empty() ? 0.0
				: elementary.unsupported_overlaps.front().common_area,
			elementary.stats.nonplanar_exact_common_pair_count, milliseconds);
		if (!passed)
		{
			for (const std::string& warning : generic.warnings)
				std::printf("  generic diagnostic: %s\n", warning.c_str());
			for (const std::string& error : elementary.errors)
				std::printf("  manufactured diagnostic: %s\n", error.c_str());
		}
		return passed;
	}

	bool manufactured_planar_overlap_regression(double scale)
	{
		const gp_Pln plane(gp::Origin(),gp::DZ());
		BRepBuilderAPI_MakeFace make_a(plane,0.0,2.0*scale,0.0,scale);
		BRepBuilderAPI_MakeFace make_b(plane,scale,3.0*scale,0.0,scale);
		if(!make_a.IsDone()||!make_b.IsDone())
		{
			std::printf("[manufactured planar scale %.3g] FAIL: face construction\n",scale);
			return false;
		}
		TopoDS_Face face_b=make_b.Face();face_b.Reverse();
		const OcctSurfaceNormalizationResult result=
			paracfd::core::normalize_occt_surface_faces({make_a.Face(),face_b});
		const double expected_input=4.0*scale*scale;
		const double expected_union=3.0*scale*scale;
		const double expected_duplicate=scale*scale;
		const double tolerance=5.0e-8*expected_input;
		bool saw_shared=false,parity_ok=true;
		for(const auto& output:result.normalized_faces)
		{
			saw_shared=saw_shared||output.contributors.size()==2;
			if(output.contributors.size()==2)
				parity_ok=parity_ok&&output.contributors[0].orientation_parity
					==-output.contributors[1].orientation_parity;
		}
		const bool passed=result.valid()&&result.clusters.size()==1
			&&result.overlap_pairs.size()==1&&result.normalized_faces.size()==3
			&&saw_shared&&parity_ok
			&&std::abs(result.stats.input_area_sum-expected_input)<=tolerance
			&&std::abs(result.stats.once_covered_area-expected_union)<=tolerance
			&&std::abs(result.stats.duplicate_area_removed-expected_duplicate)<=tolerance;
		std::printf("[manufactured planar scale %.3g] %s; input/union/duplicate "
			"%.12g / %.12g / %.12g\n",scale,passed?"PASS":"FAIL",
			result.stats.input_area_sum,result.stats.once_covered_area,
			result.stats.duplicate_area_removed);
		if(!passed)
		{
			for(const std::string& warning:result.warnings)
				std::printf("  warning: %s\n",warning.c_str());
			for(const std::string& error:result.errors)
				std::printf("  diagnostic: %s\n",error.c_str());
		}
		return passed;
	}

	void print_cluster(std::size_t id, const OcctSurfaceNormalizationClusterAudit& cluster,
		const OcctSurfaceNormalizationResult& result, bool verbose)
	{
		std::printf("\n[cluster %zu] sources", id);
		for (std::size_t source : cluster.source_face_ids) std::printf(" %zu", source);
		std::printf("\n  local arrangement %.2f ms; atoms=%zu; valid=%s\n",
			cluster.arrangement_milliseconds, cluster.atomic_face_count,
			cluster.valid() ? "yes" : "NO");
		std::printf("  input area sum / once-covered union / removed duplicate: "
			"%.12g / %.12g / %.12g mm^2\n", cluster.input_area_sum,
			cluster.once_covered_area, cluster.duplicate_area_removed);
		std::printf("  max input/output tolerance: %.12g / %.12g mm\n",
			cluster.input_maximum_tolerance, cluster.output_maximum_tolerance);
		std::printf("  atom provenance single/shared: %zu/%zu\n",
			cluster.single_source_atom_count, cluster.shared_source_atom_count);
		if (verbose)
			for (std::size_t output_id : cluster.normalized_face_indices)
			{
				const auto& face = result.normalized_faces[output_id];
				std::printf("    output %zu area=%.12g contributors", output_id, face.area);
				for (const auto& contributor : face.contributors)
					std::printf(" %zu(%c)", contributor.source_face_id,
						contributor.orientation_parity > 0 ? '+' : '-');
				std::printf("\n");
			}
		for (const std::string& warning : cluster.warnings)
			std::printf("  warning: %s\n", warning.c_str());
		for (const std::string& error : cluster.errors)
			std::printf("  FAIL: %s\n", error.c_str());
	}
}

int main(int argc, char** argv)
{
	Options options;
	if (!parse_options(argc, argv, options))
	{
		usage(argv[0]);
		return 2;
	}
	try
	{
		const bool manufactured_ok = manufactured_curved_overlap_regression()
			&& manufactured_planar_overlap_regression(1.0)
			&& manufactured_planar_overlap_regression(1.0e-3);
		if (options.manufactured_only) return manufactured_ok ? 0 : 1;
		if (!manufactured_ok) return 1;
		std::printf("OCCT exact coplanar-overlap normalization diagnostic\nSTEP: %s\n",
			options.step.c_str());
		const auto read_begin = Clock::now();
		const std::vector<TopoDS_Face> faces = unique_faces(read_step(options.step));
		const double read_seconds = std::chrono::duration<double>(
			Clock::now() - read_begin).count();
		const OcctSurfaceNormalizationResult result =
			paracfd::core::normalize_occt_surface_faces(faces);
		std::printf("read %.3f s; source faces=%zu; planar=%zu\n", read_seconds,
			result.stats.source_face_count, result.stats.planar_source_face_count);
		std::printf("exact detection %.3f s (support screen %.3f s); "
			"broad/coplanar/nonplanar candidates=%zu/%zu/%zu; "
			"support rejected/unverified/exact=%zu/%zu/%zu; "
			"overlap pairs=%zu; clusters=%zu\n", result.stats.detection_milliseconds / 1000.0,
			result.stats.support_screen_milliseconds / 1000.0,
			result.stats.broad_phase_candidate_pair_count,
			result.stats.coplanar_candidate_pair_count,
			result.stats.nonplanar_screened_pair_count,
			result.stats.nonplanar_support_rejected_pair_count,
			result.stats.nonplanar_unverified_pair_count,
			result.stats.nonplanar_exact_common_pair_count, result.overlap_pairs.size(),
			result.clusters.size());
		for (const auto& pair : result.overlap_pairs)
			std::printf("  pair %zu / %zu common=%.12g mm^2 (decision tol %.6g)\n",
				pair.source_a, pair.source_b, pair.common_area,
				pair.positive_area_tolerance);
		for (const auto& overlap : result.unsupported_overlaps)
			std::printf("  UNSUPPORTED pair %zu / %zu common=%.12g mm^2 "
				"(planar=%s/%s; decision tol %.6g)\n", overlap.source_a,
				overlap.source_b, overlap.common_area,
				overlap.source_a_planar ? "yes" : "no",
				overlap.source_b_planar ? "yes" : "no",
				overlap.positive_area_tolerance);

		std::size_t passed = 0;
		for (std::size_t cluster = 0; cluster < result.clusters.size(); ++cluster)
		{
			print_cluster(cluster, result.clusters[cluster], result, options.verbose);
			passed += result.clusters[cluster].valid();
		}
		for (const std::string& warning : result.warnings)
			std::printf("WARNING: %s\n", warning.c_str());
		for (const std::string& error : result.errors)
			std::printf("NORMALIZATION FAIL: %s\n", error.c_str());

		bool success = result.valid();
		if (options.expected_clusters && result.clusters.size() != *options.expected_clusters)
		{
			std::printf("EXPECTATION FAIL: clusters expected/actual=%zu/%zu\n",
				*options.expected_clusters, result.clusters.size());
			success = false;
		}
		if (options.expected_pairs && result.overlap_pairs.size() != *options.expected_pairs)
		{
			std::printf("EXPECTATION FAIL: overlap pairs expected/actual=%zu/%zu\n",
				*options.expected_pairs, result.overlap_pairs.size());
			success = false;
		}
		std::printf("\nverdict: %s (%zu/%zu clusters passed; %zu output faces)\n",
			success ? "PASS" : "FAIL", passed, result.clusters.size(),
			result.normalized_faces.size());
		return success ? 0 : 1;
	}
	catch (const Standard_Failure& failure)
	{
		std::fprintf(stderr, "OpenCascade exception: %s\n", failure.GetMessageString());
	}
	catch (const std::exception& exception)
	{
		std::fprintf(stderr, "error: %s\n", exception.what());
	}
	return 1;
}
