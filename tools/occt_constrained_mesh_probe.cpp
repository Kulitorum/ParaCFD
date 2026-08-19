// Diagnostic only: inject an explicitly noded interior polyline graph into OCCT's
// face-local Delaunay data structure and prove that every requested segment survives
// as an edge of the emitted Poly_Triangulation.
//
// This is deliberately an isolated diagnostic and is not part of STEP import.  Its
// default fixture is a curved bicubic BSpline patch with a six-segment interior
// T-junction.  In real-STEP mode it audits every eligible curved interior-contact
// face unless --face-id explicitly narrows the run.  OCCT does not report failed
// constraint recovery, so both transient BRepMesh links and final triangle edges
// are audited.  A missing segment or failed face makes the 100% gate fail.

#include "core/geometry/step_import.h"

#include <BRepAdaptor_Surface.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepMesh_Context.hxx>
#include <BRepMesh_DataStructureOfDelaun.hxx>
#include <BRepMesh_Delaun.hxx>
#include <BRepMesh_DelaunayBaseMeshAlgo.hxx>
#include <BRepMesh_DelaunayDeflectionControlMeshAlgo.hxx>
#include <BRepMesh_Edge.hxx>
#include <BRepMesh_FaceDiscret.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepMesh_NURBSRangeSplitter.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineSurface.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IMeshData_Face.hxx>
#include <IMeshData_Status.hxx>
#include <IMeshTools_MeshAlgo.hxx>
#include <IMeshTools_MeshAlgoFactory.hxx>
#include <IMeshTools_Parameters.hxx>
#include <Message_ProgressRange.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <STEPControl_Reader.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
	struct Options
	{
		double deflection_mm = 2.0;
		double angular_deflection_rad = 0.25;
		std::string step_path;
		std::optional<std::uint32_t> face_id;
	};

	void usage(const char* executable)
	{
		std::printf("Usage: %s [--deflection-mm D] [--angle-rad A] "
			"[--step FILE [--face-id N]]\n", executable);
	}

	bool parse_positive_double(int& index, int argc, char** argv, double& value)
	{
		if (index + 1 >= argc) return false;
		char* end = nullptr;
		const double parsed = std::strtod(argv[++index], &end);
		if (!end || *end != '\0' || !std::isfinite(parsed) || parsed <= 0.0) return false;
		value = parsed;
		return true;
	}

	bool parse_options(int argc, char** argv, Options& options)
	{
		for (int index = 1; index < argc; ++index)
		{
			const std::string argument = argv[index];
			if (argument == "--deflection-mm")
			{
				if (!parse_positive_double(index, argc, argv, options.deflection_mm)) return false;
			}
			else if (argument == "--angle-rad")
			{
				if (!parse_positive_double(index, argc, argv,
					options.angular_deflection_rad)) return false;
			}
			else if (argument == "--step")
			{
				if (index + 1 >= argc) return false;
				options.step_path = argv[++index];
			}
			else if (argument == "--face-id")
			{
				if (index + 1 >= argc) return false;
				char* end = nullptr;
				const unsigned long parsed = std::strtoul(argv[++index], &end, 10);
				if (!end || *end != '\0'
					|| parsed > std::numeric_limits<std::uint32_t>::max()) return false;
				options.face_id = static_cast<std::uint32_t>(parsed);
			}
			else if (argument == "--help" || argument == "-h")
			{
				usage(argv[0]);
				std::exit(EXIT_SUCCESS);
			}
			else
			{
				std::fprintf(stderr, "unknown option: %s\n", argument.c_str());
				return false;
			}
		}
		return true;
	}

	struct ConstraintVertex
	{
		std::string name;
		gp_Pnt2d uv;
		gp_Pnt point_mm;
		double uv_tolerance = 1.0e-9;
		double position_tolerance_mm = 1.0e-7;
	};

	struct ConstraintSegment
	{
		std::string name;
		std::size_t first = 0;
		std::size_t last = 0;
	};

	struct ConstraintFixture
	{
		std::vector<ConstraintVertex> vertices;
		std::vector<ConstraintSegment> segments;
		std::uint32_t source_face_id = std::numeric_limits<std::uint32_t>::max();
		std::vector<std::uint64_t> source_contact_ids;
		std::size_t maximum_graph_degree = 0;
		std::uint32_t maximum_fan_degree = 0;
		double normal_variation_rad = 0.0;
	};

	struct ConstraintRecoveryReport
	{
		bool factory_called = false;
		bool initialization_called = false;
		bool post_process_called = false;
		bool injection_ok = true;
		bool data_structure_ok = true;
		GeomAbs_SurfaceType requested_surface_type = GeomAbs_OtherSurface;
		std::vector<int> node_ids;
		std::vector<int> link_ids;
		std::vector<bool> link_is_fixed;
		std::vector<int> link_incidences;
		std::vector<std::string> errors;
	};

	void record_error(const std::shared_ptr<ConstraintRecoveryReport>& report,
		const std::string& message)
	{
		report->errors.push_back(message);
	}

	template <class BaseAlgorithm>
	class GuardedConstraintMeshAlgorithm final : public BaseAlgorithm
	{
	public:
		GuardedConstraintMeshAlgorithm(std::shared_ptr<const ConstraintFixture> fixture,
			std::shared_ptr<ConstraintRecoveryReport> report)
			: fixture_(std::move(fixture)), report_(std::move(report))
		{
		}

		~GuardedConstraintMeshAlgorithm() override = default;

	protected:
		bool initDataStructure() override
		{
			report_->initialization_called = true;
			if (!BaseAlgorithm::initDataStructure())
			{
				report_->injection_ok = false;
				record_error(report_, "native BRepMesh data-structure initialization failed");
				return false;
			}

			if (!fixture_ || fixture_->vertices.empty() || fixture_->segments.empty())
			{
				report_->injection_ok = false;
				record_error(report_, "constraint fixture is empty");
				return false;
			}

			report_->node_ids.clear();
			report_->link_ids.clear();
			report_->node_ids.reserve(fixture_->vertices.size());
			report_->link_ids.reserve(fixture_->segments.size());

			for (const ConstraintVertex& vertex : fixture_->vertices)
			{
				// Pass unscaled face UV.  The standard NodeInsertion algorithm routes this
				// through its virtual addNodeToStructure() override, which applies the
				// surface-specific range scaling exactly once.
				const int node_id = this->registerNode(vertex.point_mm, vertex.uv,
					BRepMesh_Fixed, false);
				if (node_id <= 0)
				{
					report_->injection_ok = false;
					record_error(report_, std::string("failed to register constraint vertex ")
						+ vertex.name);
					return false;
				}
				report_->node_ids.push_back(node_id);
			}

			const occ::handle<BRepMesh_DataStructureOfDelaun>& structure =
				this->getStructure();
			for (const ConstraintSegment& segment : fixture_->segments)
			{
				if (segment.first >= report_->node_ids.size()
					|| segment.last >= report_->node_ids.size()
					|| segment.first == segment.last)
				{
					report_->injection_ok = false;
					record_error(report_, std::string("invalid endpoint indices for segment ")
						+ segment.name);
					return false;
				}

				const int signed_link_id = structure->AddLink(BRepMesh_Edge(
					report_->node_ids[segment.first], report_->node_ids[segment.last],
					BRepMesh_Fixed));
				const int link_id = std::abs(signed_link_id);
				if (link_id <= 0 || structure->GetLink(link_id).Movability() != BRepMesh_Fixed)
				{
					report_->injection_ok = false;
					record_error(report_, std::string("failed to register fixed link ")
						+ segment.name);
					return false;
				}
				report_->link_ids.push_back(link_id);
			}

			return true;
		}

		void postProcessMesh(BRepMesh_Delaun& mesher,
			const Message_ProgressRange& range) override
		{
			// Let the ordinary NURBS deflection-control path finish inserting and
			// optimizing surface nodes before checking constraint incidence.
			BaseAlgorithm::postProcessMesh(mesher, range);
			report_->post_process_called = true;
			report_->link_is_fixed.clear();
			report_->link_incidences.clear();
			report_->data_structure_ok = true;

			const occ::handle<BRepMesh_DataStructureOfDelaun>& structure =
				this->getStructure();
			if (report_->link_ids.size() != fixture_->segments.size())
			{
				report_->data_structure_ok = false;
				record_error(report_, "not all requested links reached post-processing");
				return;
			}

			for (std::size_t index = 0; index < report_->link_ids.size(); ++index)
			{
				const int link_id = report_->link_ids[index];
				const BRepMesh_Edge& link = structure->GetLink(link_id);
				const bool is_fixed = link.Movability() == BRepMesh_Fixed;
				const int incidence = structure->ElementsConnectedTo(link_id).Extent();
				report_->link_is_fixed.push_back(is_fixed);
				report_->link_incidences.push_back(incidence);
				if (!is_fixed || incidence != 2)
				{
					report_->data_structure_ok = false;
					record_error(report_, std::string("fixed-link audit failed for ")
						+ fixture_->segments[index].name + " (fixed="
						+ (is_fixed ? "true" : "false") + ", incidence="
						+ std::to_string(incidence) + ")");
				}
			}
		}

	private:
		std::shared_ptr<const ConstraintFixture> fixture_;
		std::shared_ptr<ConstraintRecoveryReport> report_;
	};

	using NativeNurbsMeshAlgorithm = BRepMesh_DelaunayDeflectionControlMeshAlgo<
		BRepMesh_NURBSRangeSplitter, BRepMesh_DelaunayBaseMeshAlgo>;
	using ConstrainedNurbsMeshAlgorithm =
		GuardedConstraintMeshAlgorithm<NativeNurbsMeshAlgorithm>;

	class ConstraintMeshAlgorithmFactory final : public IMeshTools_MeshAlgoFactory
	{
	public:
		ConstraintMeshAlgorithmFactory(std::shared_ptr<const ConstraintFixture> fixture,
			std::shared_ptr<ConstraintRecoveryReport> report)
			: fixture_(std::move(fixture)), report_(std::move(report))
		{
		}

		~ConstraintMeshAlgorithmFactory() override = default;

		occ::handle<IMeshTools_MeshAlgo> GetAlgo(GeomAbs_SurfaceType surface_type,
			const IMeshTools_Parameters&) const override
		{
			report_->factory_called = true;
			report_->requested_surface_type = surface_type;
			if (surface_type != GeomAbs_BSplineSurface)
			{
				report_->injection_ok = false;
				record_error(report_, "probe factory received a non-BSpline face");
				return occ::handle<IMeshTools_MeshAlgo>();
			}
			return new ConstrainedNurbsMeshAlgorithm(fixture_, report_);
		}

	private:
		std::shared_ptr<const ConstraintFixture> fixture_;
		std::shared_ptr<ConstraintRecoveryReport> report_;
	};

	occ::handle<Geom_BSplineSurface> make_manufactured_surface()
	{
		NCollection_Array2<gp_Pnt> poles(1, 4, 1, 4);
		for (int u_index = 1; u_index <= 4; ++u_index)
		{
			const double u = static_cast<double>(u_index - 1) / 3.0;
			for (int v_index = 1; v_index <= 4; ++v_index)
			{
				const double v = static_cast<double>(v_index - 1) / 3.0;
				const double x = 120.0 * u;
				const double y = 90.0 * v;
				const double z = 18.0 * u * (1.0 - u)
					- 11.0 * v * (1.0 - v)
					+ 7.0 * (u - 0.5) * (v - 0.5);
				poles.SetValue(u_index, v_index, gp_Pnt(x, y, z));
			}
		}

		NCollection_Array1<double> u_knots(1, 2);
		NCollection_Array1<double> v_knots(1, 2);
		NCollection_Array1<int> u_multiplicities(1, 2);
		NCollection_Array1<int> v_multiplicities(1, 2);
		u_knots.SetValue(1, 0.0);
		u_knots.SetValue(2, 1.0);
		v_knots.SetValue(1, 0.0);
		v_knots.SetValue(2, 1.0);
		u_multiplicities.SetValue(1, 4);
		u_multiplicities.SetValue(2, 4);
		v_multiplicities.SetValue(1, 4);
		v_multiplicities.SetValue(2, 4);

		return new Geom_BSplineSurface(poles, u_knots, v_knots,
			u_multiplicities, v_multiplicities, 3, 3, false, false);
	}

	std::shared_ptr<ConstraintFixture> make_t_junction_fixture(
		const occ::handle<Geom_BSplineSurface>& surface)
	{
		auto fixture = std::make_shared<ConstraintFixture>();
		const auto add_vertex = [&](const char* name, double u, double v)
		{
			fixture->vertices.push_back({name, gp_Pnt2d(u, v), surface->Value(u, v)});
		};

		add_vertex("left", 0.20, 0.50);
		add_vertex("left-mid", 0.35, 0.50);
		add_vertex("junction", 0.50, 0.50);
		add_vertex("right-mid", 0.65, 0.50);
		add_vertex("right", 0.80, 0.50);
		add_vertex("stem-mid", 0.50, 0.65);
		add_vertex("stem-end", 0.50, 0.80);

		fixture->segments = {
			{"left-0", 0, 1},
			{"left-1", 1, 2},
			{"right-0", 2, 3},
			{"right-1", 3, 4},
			{"stem-0", 2, 5},
			{"stem-1", 5, 6},
		};
		fixture->maximum_graph_degree = 3;
		return fixture;
	}

	std::uint64_t undirected_edge_key(int first, int last)
	{
		if (last < first) std::swap(first, last);
		return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(first)) << 32)
			| static_cast<std::uint32_t>(last);
	}

	struct ContactUseRef
	{
		const paracfd::core::StepContactCurve* curve = nullptr;
		const paracfd::core::StepContactFaceUse* use = nullptr;
	};

	bool read_unmeshed_step_faces(const std::string& path, TopoDS_Shape& shape,
		std::vector<TopoDS_Face>& faces, std::string& error)
	{
		STEPControl_Reader reader;
		if (reader.ReadFile(path.c_str()) != IFSelect_RetDone)
		{
			error = "STEPControl_Reader::ReadFile failed";
			return false;
		}
		if (reader.TransferRoots() <= 0)
		{
			error = "STEPControl_Reader::TransferRoots transferred no roots";
			return false;
		}
		shape = reader.OneShape();
		if (shape.IsNull())
		{
			error = "STEPControl_Reader::OneShape returned a null shape";
			return false;
		}

		faces.clear();
		for (TopExp_Explorer face(shape, TopAbs_FACE); face.More(); face.Next())
			faces.push_back(TopoDS::Face(face.Current()));
		if (faces.empty())
		{
			error = "transferred STEP shape contains no faces";
			return false;
		}
		return true;
	}

	double fixture_normal_variation(const TopoDS_Face& face,
		const ConstraintFixture& fixture)
	{
		BRepAdaptor_Surface surface(face, false);
		gp_Vec reference;
		bool has_reference = false;
		double maximum_angle = 0.0;
		for (const ConstraintVertex& vertex : fixture.vertices)
		{
			try
			{
				gp_Pnt point;
				gp_Vec derivative_u;
				gp_Vec derivative_v;
				surface.D1(vertex.uv.X(), vertex.uv.Y(), point, derivative_u, derivative_v);
				gp_Vec normal = derivative_u.Crossed(derivative_v);
				if (normal.SquareMagnitude() <= gp::Resolution()) continue;
				normal.Normalize();
				if (!has_reference)
				{
					reference = normal;
					has_reference = true;
				}
				else
				{
					maximum_angle = std::max(maximum_angle, reference.Angle(normal));
				}
			}
			catch (const Standard_Failure&)
			{
				// A single singular sample is not useful for curvature scoring.  Meshing
				// and the final segment audit remain authoritative.
			}
		}
		return maximum_angle;
	}

	std::shared_ptr<ConstraintFixture> build_real_contact_fixture(const TopoDS_Face& face,
		std::uint32_t face_id, const std::vector<ContactUseRef>& contact_uses)
	{
		auto fixture = std::make_shared<ConstraintFixture>();
		fixture->source_face_id = face_id;
		BRepAdaptor_Surface surface(face, false);
		std::set<std::uint64_t> source_contacts;
		std::set<std::uint64_t> emitted_segments;

		for (const ContactUseRef& ref : contact_uses)
		{
			if (!ref.curve || !ref.use) continue;
			const std::size_t sample_count = std::min(ref.curve->sample_positions_m.size(),
				ref.use->sample_uv.size());
			if (sample_count < 2) continue;
			source_contacts.insert(ref.curve->id);
			fixture->maximum_fan_degree = std::max(fixture->maximum_fan_degree,
				ref.curve->fan_degree);

			const double position_tolerance_mm = std::max(1.0e-6,
				8.0 * ref.curve->tolerance_m * 1000.0);
			const double u_tolerance = std::max(1.0e-11,
				8.0 * std::abs(surface.UResolution(position_tolerance_mm)));
			const double v_tolerance = std::max(1.0e-11,
				8.0 * std::abs(surface.VResolution(position_tolerance_mm)));
			const double uv_tolerance = std::max(u_tolerance, v_tolerance);
			std::vector<std::size_t> chain;
			chain.reserve(sample_count);

			for (std::size_t sample = 0; sample < sample_count; ++sample)
			{
				const auto& source_uv = ref.use->sample_uv[sample];
				const auto& source_point = ref.curve->sample_positions_m[sample];
				if (!std::isfinite(source_uv[0]) || !std::isfinite(source_uv[1])
					|| !std::isfinite(source_point[0]) || !std::isfinite(source_point[1])
					|| !std::isfinite(source_point[2]))
					continue;

				const gp_Pnt2d uv(source_uv[0], source_uv[1]);
				const gp_Pnt point_mm(source_point[0] * 1000.0,
					source_point[1] * 1000.0, source_point[2] * 1000.0);
				std::size_t vertex_id = fixture->vertices.size();
				for (std::size_t existing = 0; existing < fixture->vertices.size(); ++existing)
				{
					const ConstraintVertex& candidate = fixture->vertices[existing];
					const double accepted_uv = std::max(uv_tolerance, candidate.uv_tolerance);
					const double accepted_position = std::max(position_tolerance_mm,
						candidate.position_tolerance_mm);
					if (std::abs(candidate.uv.X() - uv.X()) <= accepted_uv
						&& std::abs(candidate.uv.Y() - uv.Y()) <= accepted_uv
						&& candidate.point_mm.Distance(point_mm) <= accepted_position)
					{
						vertex_id = existing;
						break;
					}
				}

				if (vertex_id == fixture->vertices.size())
				{
					ConstraintVertex vertex;
					vertex.name = "c" + std::to_string(ref.curve->id) + "-p"
						+ std::to_string(sample);
					vertex.uv = uv;
					vertex.point_mm = point_mm;
					vertex.uv_tolerance = uv_tolerance;
					vertex.position_tolerance_mm = position_tolerance_mm;
					fixture->vertices.push_back(std::move(vertex));
				}
				else
				{
					ConstraintVertex& vertex = fixture->vertices[vertex_id];
					vertex.uv_tolerance = std::max(vertex.uv_tolerance, uv_tolerance);
					vertex.position_tolerance_mm = std::max(vertex.position_tolerance_mm,
						position_tolerance_mm);
				}

				if (chain.empty() || chain.back() != vertex_id) chain.push_back(vertex_id);
			}

			for (std::size_t segment = 1; segment < chain.size(); ++segment)
			{
				if (chain[segment - 1] == chain[segment]) continue;
				const std::uint64_t key = undirected_edge_key(
					static_cast<int>(chain[segment - 1]), static_cast<int>(chain[segment]));
				if (!emitted_segments.insert(key).second) continue;
				fixture->segments.push_back({"c" + std::to_string(ref.curve->id) + "-s"
					+ std::to_string(segment - 1), chain[segment - 1], chain[segment]});
			}
		}

		fixture->source_contact_ids.assign(source_contacts.begin(), source_contacts.end());
		std::vector<std::size_t> degree(fixture->vertices.size(), 0);
		for (const ConstraintSegment& segment : fixture->segments)
		{
			if (segment.first < degree.size()) ++degree[segment.first];
			if (segment.last < degree.size()) ++degree[segment.last];
		}
		for (std::size_t value : degree)
			fixture->maximum_graph_degree = std::max(fixture->maximum_graph_degree, value);
		fixture->normal_variation_rad = fixture_normal_variation(face, *fixture);
		return fixture;
	}

	struct RealFixtureSelection
	{
		TopoDS_Face face;
		std::shared_ptr<ConstraintFixture> fixture;
	};

	std::optional<std::vector<RealFixtureSelection>> select_real_fixtures(
		const Options& options,
		std::string& error)
	{
		paracfd::core::StepGeometry geometry = paracfd::core::load_step_geometry(
			options.step_path, options.deflection_mm, &error);
		if (geometry.mesh.empty())
		{
			if (error.empty()) error = "production STEP pass returned an empty mesh";
			return std::nullopt;
		}

		TopoDS_Shape raw_shape;
		std::vector<TopoDS_Face> faces;
		if (!read_unmeshed_step_faces(options.step_path, raw_shape, faces, error))
			return std::nullopt;
		if (faces.size() != geometry.source_faces.size())
		{
			error = "independent STEP transfer changed deterministic source-face count";
			return std::nullopt;
		}

		std::map<std::uint32_t, std::vector<ContactUseRef>> uses_by_face;
		for (const paracfd::core::StepContactCurve& curve : geometry.contacts.curves)
		{
			for (const paracfd::core::StepContactFaceUse& use : curve.uses)
			{
				if (use.kind != paracfd::core::StepContactUseKind::face_interior) continue;
				uses_by_face[use.source_face_id].push_back({&curve, &use});
			}
		}

		std::vector<RealFixtureSelection> selections;
		for (const auto& [face_id, uses] : uses_by_face)
		{
			if (options.face_id && *options.face_id != face_id) continue;
			if (face_id >= faces.size()) continue;
			BRepAdaptor_Surface surface(faces[face_id], false);
			if (surface.GetType() != GeomAbs_BSplineSurface) continue;

			std::shared_ptr<ConstraintFixture> fixture = build_real_contact_fixture(
				faces[face_id], face_id, uses);
			std::printf("  real candidate face=%u curves=%zu vertices=%zu segments=%zu "
				"max_graph_degree=%zu max_fan_degree=%u normal_variation_rad=%.6g\n", face_id,
				fixture->source_contact_ids.size(), fixture->vertices.size(),
				fixture->segments.size(), fixture->maximum_graph_degree,
				fixture->maximum_fan_degree, fixture->normal_variation_rad);
			if (fixture->segments.empty() || fixture->maximum_fan_degree < 3
				|| fixture->normal_variation_rad <= 1.0e-6)
				continue;

			selections.push_back(RealFixtureSelection{faces[face_id], std::move(fixture)});
		}

		if (selections.empty())
		{
			error = options.face_id
				? "requested face is not a curved BSpline carrying an interior contact"
				: "no curved BSpline face carrying an interior contact was found";
			return std::nullopt;
		}
		std::printf("  real candidate_count=%zu\n", selections.size());
		return selections;
	}

	struct PolyAuditResult
	{
		bool ok = true;
		std::vector<int> node_ids;
		std::vector<double> uv_errors;
		std::vector<double> position_errors_mm;
		std::vector<int> edge_incidences;
		std::vector<std::string> errors;
	};

	PolyAuditResult audit_poly_triangulation(const TopoDS_Face& face,
		const ConstraintFixture& fixture)
	{
		PolyAuditResult result;
		TopLoc_Location location;
		const occ::handle<Poly_Triangulation> triangulation =
			BRep_Tool::Triangulation(face, location);
		if (triangulation.IsNull())
		{
			result.ok = false;
			result.errors.push_back("face has no Poly_Triangulation");
			return result;
		}
		if (!triangulation->HasUVNodes())
		{
			result.ok = false;
			result.errors.push_back("Poly_Triangulation has no UV nodes");
			return result;
		}

		result.node_ids.reserve(fixture.vertices.size());
		result.uv_errors.reserve(fixture.vertices.size());
		result.position_errors_mm.reserve(fixture.vertices.size());

		for (const ConstraintVertex& requested : fixture.vertices)
		{
			int closest_node = 0;
			double closest_uv_error = std::numeric_limits<double>::infinity();
			for (int node = 1; node <= triangulation->NbNodes(); ++node)
			{
				const double error = requested.uv.Distance(triangulation->UVNode(node));
				if (error < closest_uv_error)
				{
					closest_uv_error = error;
					closest_node = node;
				}
			}

			result.node_ids.push_back(closest_node);
			result.uv_errors.push_back(closest_uv_error);
			gp_Pnt emitted_point = triangulation->Node(closest_node);
			emitted_point.Transform(location.Transformation());
			const double position_error = emitted_point.Distance(requested.point_mm);
			result.position_errors_mm.push_back(position_error);
			if (closest_node <= 0 || closest_uv_error > requested.uv_tolerance
				|| position_error > requested.position_tolerance_mm)
			{
				result.ok = false;
				result.errors.push_back(std::string("output-node audit failed for ")
					+ requested.name + " (uv_error=" + std::to_string(closest_uv_error)
					+ ", position_error_mm=" + std::to_string(position_error) + ")");
			}
		}

		std::vector<int> distinct_nodes = result.node_ids;
		std::sort(distinct_nodes.begin(), distinct_nodes.end());
		distinct_nodes.erase(std::unique(distinct_nodes.begin(), distinct_nodes.end()),
			distinct_nodes.end());
		if (distinct_nodes.size() != fixture.vertices.size())
		{
			result.ok = false;
			result.errors.push_back("distinct requested vertices collapsed to one output node");
		}

		std::unordered_map<std::uint64_t, int> edge_incidences;
		for (int triangle = 1; triangle <= triangulation->NbTriangles(); ++triangle)
		{
			int nodes[3] = {0, 0, 0};
			triangulation->Triangle(triangle).Get(nodes[0], nodes[1], nodes[2]);
			if (nodes[0] <= 0 || nodes[1] <= 0 || nodes[2] <= 0
				|| nodes[0] > triangulation->NbNodes()
				|| nodes[1] > triangulation->NbNodes()
				|| nodes[2] > triangulation->NbNodes()
				|| nodes[0] == nodes[1] || nodes[1] == nodes[2] || nodes[2] == nodes[0])
			{
				result.ok = false;
				result.errors.push_back("Poly_Triangulation contains an invalid triangle");
				continue;
			}
			for (int edge = 0; edge < 3; ++edge)
			{
				++edge_incidences[undirected_edge_key(nodes[edge], nodes[(edge + 1) % 3])];
			}
		}

		result.edge_incidences.reserve(fixture.segments.size());
		for (const ConstraintSegment& segment : fixture.segments)
		{
			const int first = result.node_ids[segment.first];
			const int last = result.node_ids[segment.last];
			const auto found = edge_incidences.find(undirected_edge_key(first, last));
			const int incidence = found == edge_incidences.end() ? 0 : found->second;
			result.edge_incidences.push_back(incidence);
			if (incidence != 2)
			{
				result.ok = false;
				result.errors.push_back(std::string("requested segment is not a two-sided "
					"triangle edge: ") + segment.name + " (incidence="
					+ std::to_string(incidence) + ")");
			}
		}

		std::printf("  Poly_Triangulation: nodes=%d triangles=%d\n",
			triangulation->NbNodes(), triangulation->NbTriangles());
		return result;
	}

	void print_report(const ConstraintFixture& fixture,
		const ConstraintRecoveryReport& recovery, const PolyAuditResult& poly)
	{
		std::printf("  constraint vertices=%zu segments=%zu\n",
			fixture.vertices.size(), fixture.segments.size());
		for (std::size_t index = 0; index < fixture.vertices.size(); ++index)
		{
			const ConstraintVertex& vertex = fixture.vertices[index];
			const int structure_node = index < recovery.node_ids.size()
				? recovery.node_ids[index] : 0;
			const int poly_node = index < poly.node_ids.size() ? poly.node_ids[index] : 0;
			const double uv_error = index < poly.uv_errors.size()
				? poly.uv_errors[index] : std::numeric_limits<double>::infinity();
			const double position_error = index < poly.position_errors_mm.size()
				? poly.position_errors_mm[index] : std::numeric_limits<double>::infinity();
			std::printf("  vertex %-10s uv=(%.6f, %.6f) ds=%d poly=%d "
				"uv_err=%.3e xyz_err_mm=%.3e\n", vertex.name.c_str(), vertex.uv.X(),
				vertex.uv.Y(), structure_node, poly_node, uv_error, position_error);
		}

		for (std::size_t index = 0; index < fixture.segments.size(); ++index)
		{
			const ConstraintSegment& segment = fixture.segments[index];
			const int link_id = index < recovery.link_ids.size() ? recovery.link_ids[index] : 0;
			const bool link_is_fixed = index < recovery.link_is_fixed.size()
				&& recovery.link_is_fixed[index];
			const int data_structure_incidence = index < recovery.link_incidences.size()
				? recovery.link_incidences[index] : 0;
			const int poly_incidence = index < poly.edge_incidences.size()
				? poly.edge_incidences[index] : 0;
			std::printf("  segment %-8s link=%d fixed=%s ds_incidence=%d "
				"poly_incidence=%d %s\n", segment.name.c_str(), link_id,
				link_is_fixed ? "true" : "false", data_structure_incidence, poly_incidence,
				(link_is_fixed && data_structure_incidence == 2 && poly_incidence == 2)
					? "PASS" : "FAIL");
		}

		for (const std::string& error : recovery.errors)
			std::fprintf(stderr, "  recovery error: %s\n", error.c_str());
		for (const std::string& error : poly.errors)
			std::fprintf(stderr, "  output error: %s\n", error.c_str());
	}

	struct ConstraintCaseResult
	{
		bool ok = false;
		bool mesher_done = false;
		bool lifecycle_ok = false;
		bool injection_ok = false;
		bool data_structure_ok = false;
		bool poly_ok = false;
		std::size_t requested_segments = 0;
		std::size_t data_structure_passed_segments = 0;
		std::size_t poly_passed_segments = 0;
		std::size_t passed_segments = 0;
	};

	ConstraintCaseResult run_constrained_case(const TopoDS_Face& face,
		const std::shared_ptr<ConstraintFixture>& fixture, const Options& options,
		const std::string& label)
	{
		const std::shared_ptr<ConstraintRecoveryReport> recovery =
			std::make_shared<ConstraintRecoveryReport>();
		const occ::handle<IMeshTools_MeshAlgoFactory> factory =
			new ConstraintMeshAlgorithmFactory(fixture, recovery);
		const occ::handle<BRepMesh_Context> context = new BRepMesh_Context(
			IMeshTools_MeshAlgoType_Watson);
		context->SetFaceDiscret(new BRepMesh_FaceDiscret(factory));

		IMeshTools_Parameters parameters;
		parameters.MeshAlgo = IMeshTools_MeshAlgoType_Watson;
		parameters.Deflection = options.deflection_mm;
		parameters.DeflectionInterior = options.deflection_mm;
		parameters.Angle = options.angular_deflection_rad;
		parameters.AngleInterior = options.angular_deflection_rad;
		parameters.MinSize = std::max(Precision::Confusion(),
			0.05 * options.deflection_mm);
		parameters.InParallel = false;
		parameters.Relative = false;
		parameters.InternalVerticesMode = true;
		parameters.ControlSurfaceDeflection = true;
		parameters.EnableControlSurfaceDeflectionAllSurfaces = true;
		parameters.CleanModel = false;

		BRepMesh_IncrementalMesh mesher;
		mesher.SetShape(face);
		mesher.ChangeParameters() = parameters;
		mesher.Perform(context);

		std::printf("OCCT constrained-mesh %s\n", label.c_str());
		std::printf("  deflection_mm=%.6g angle_rad=%.6g mesher_done=%s status=0x%x\n",
			options.deflection_mm, options.angular_deflection_rad,
			mesher.IsDone() ? "true" : "false", mesher.GetStatusFlags());
		if (fixture->source_face_id != std::numeric_limits<std::uint32_t>::max())
		{
			std::printf("  source_face=%u source_contacts=%zu max_graph_degree=%zu "
				"max_fan_degree=%u normal_variation_rad=%.6g\n", fixture->source_face_id,
				fixture->source_contact_ids.size(), fixture->maximum_graph_degree,
				fixture->maximum_fan_degree, fixture->normal_variation_rad);
		}

		const PolyAuditResult poly = audit_poly_triangulation(face, *fixture);
		print_report(*fixture, *recovery, poly);
		const bool lifecycle_ok = recovery->factory_called
			&& recovery->initialization_called && recovery->post_process_called;

		ConstraintCaseResult result;
		result.mesher_done = mesher.IsDone();
		result.lifecycle_ok = lifecycle_ok;
		result.injection_ok = recovery->injection_ok;
		result.data_structure_ok = recovery->data_structure_ok;
		result.poly_ok = poly.ok;
		result.requested_segments = fixture->segments.size();
		for (std::size_t index = 0; index < fixture->segments.size(); ++index)
		{
			const bool data_structure_passed = index < recovery->link_is_fixed.size()
				&& recovery->link_is_fixed[index]
				&& index < recovery->link_incidences.size()
				&& recovery->link_incidences[index] == 2;
			const bool poly_passed = index < poly.edge_incidences.size()
				&& poly.edge_incidences[index] == 2;
			if (data_structure_passed) ++result.data_structure_passed_segments;
			if (poly_passed) ++result.poly_passed_segments;
			if (data_structure_passed && poly_passed) ++result.passed_segments;
		}
		result.ok = result.mesher_done && result.lifecycle_ok && result.injection_ok
			&& result.data_structure_ok && result.poly_ok
			&& result.passed_segments == result.requested_segments;
		std::printf("RESULT: %s\n", result.ok ? "PASS" : "FAIL");
		return result;
	}
}

int main(int argc, char** argv)
{
	Options options;
	if (!parse_options(argc, argv, options))
	{
		usage(argv[0]);
		return EXIT_FAILURE;
	}

	try
	{
		if (!options.step_path.empty())
		{
			std::string error;
			std::optional<std::vector<RealFixtureSelection>> selections =
				select_real_fixtures(options, error);
			if (!selections)
			{
				std::fprintf(stderr, "real STEP fixture selection failed: %s\n", error.c_str());
				return EXIT_FAILURE;
			}

			std::size_t passed_faces = 0;
			std::size_t requested_segments = 0;
			std::size_t data_structure_passed_segments = 0;
			std::size_t poly_passed_segments = 0;
			std::size_t passed_segments = 0;
			for (const RealFixtureSelection& selection : *selections)
			{
				ConstraintCaseResult result;
				result.requested_segments = selection.fixture->segments.size();
				const std::uint32_t face_id = selection.fixture->source_face_id;
				try
				{
					// Each candidate starts from an independent topology copy so triangulation
					// state on shared STEP edges cannot leak from one audit into the next.
					BRepBuilderAPI_Copy face_copy(selection.face, true, false);
					if (!face_copy.IsDone() || face_copy.Shape().IsNull())
					{
						std::fprintf(stderr, "face %u isolation copy failed\n", face_id);
					}
					else
					{
						const TopoDS_Face isolated_face = TopoDS::Face(face_copy.Shape());
						const std::string constraint_kind =
							selection.fixture->maximum_graph_degree >= 3
								? "branched interior constraint"
								: "unbranched interior contact";
						const std::string label = "real STEP face " + std::to_string(face_id)
							+ " with " + constraint_kind;
						result = run_constrained_case(isolated_face, selection.fixture,
							options, label);
					}
				}
				catch (const Standard_Failure& failure)
				{
					std::fprintf(stderr, "face %u OCCT failure: %s\n", face_id,
						failure.what());
				}
				catch (const std::exception& exception)
				{
					std::fprintf(stderr, "face %u failure: %s\n", face_id,
						exception.what());
				}

				if (result.ok) ++passed_faces;
				requested_segments += result.requested_segments;
				data_structure_passed_segments += result.data_structure_passed_segments;
				poly_passed_segments += result.poly_passed_segments;
				passed_segments += result.passed_segments;
				std::printf("FACE AUDIT face=%u result=%s segments=%zu/%zu "
					"ds_fixed_two_sided=%zu/%zu poly_two_sided=%zu/%zu\n", face_id,
					result.ok ? "PASS" : "FAIL", result.passed_segments,
					result.requested_segments, result.data_structure_passed_segments,
					result.requested_segments, result.poly_passed_segments,
					result.requested_segments);
			}

			const bool all_passed = passed_faces == selections->size()
				&& passed_segments == requested_segments;
			std::printf("REAL STEP AUDIT SUMMARY candidates=%zu faces=%zu/%zu "
				"segments=%zu/%zu ds_fixed_two_sided=%zu/%zu poly_two_sided=%zu/%zu\n",
				selections->size(), passed_faces, selections->size(), passed_segments,
				requested_segments, data_structure_passed_segments, requested_segments,
				poly_passed_segments, requested_segments);
			std::printf("100%% GATE: %s\n", all_passed ? "PASS" : "FAIL");
			return all_passed ? EXIT_SUCCESS : EXIT_FAILURE;
		}
		else
		{
			const occ::handle<Geom_BSplineSurface> surface = make_manufactured_surface();
			BRepBuilderAPI_MakeFace face_builder(surface, 0.0, 1.0, 0.0, 1.0,
				Precision::Confusion());
			if (!face_builder.IsDone())
			{
				std::fprintf(stderr, "failed to build manufactured BSpline face\n");
				return EXIT_FAILURE;
			}
			const TopoDS_Face face = face_builder.Face();
			const BRepAdaptor_Surface adaptor(face, false);
			if (adaptor.GetType() != GeomAbs_BSplineSurface)
			{
				std::fprintf(stderr,
					"manufactured face is not reported as a BSpline surface\n");
				return EXIT_FAILURE;
			}
			const std::shared_ptr<ConstraintFixture> fixture =
				make_t_junction_fixture(surface);
			const ConstraintCaseResult result = run_constrained_case(face, fixture, options,
				"manufactured BSpline T-junction");
			return result.ok ? EXIT_SUCCESS : EXIT_FAILURE;
		}
	}
	catch (const Standard_Failure& failure)
	{
		std::fprintf(stderr, "OCCT failure: %s\n", failure.what());
	}
	catch (const std::exception& exception)
	{
		std::fprintf(stderr, "failure: %s\n", exception.what());
	}
	return EXIT_FAILURE;
}
