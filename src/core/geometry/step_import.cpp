// step_import.cpp — see step_import.h. The ONLY translation unit that touches OpenCascade.
//
// Paraglider CAD pipeline:
//   STEPControl_Reader → ReadFile → TransferRoots → OneShape → BRepMesh_IncrementalMesh →
//   per-face BRep_Tool::Triangulation (apply TopLoc_Location) → 1-based→0-based indices.
//
// Three correctness-critical traps, all handled here:
//   (a) mm→m — OCC coordinates are millimetres; every node is multiplied by 0.001.
//   (b) winding — OCC stores triangles in the surface's natural (u,v) sense. We apply the
//       TopoDS face orientation so every patch has a stable local minus-to-plus normal. No
//       globally outward or watertight shell orientation is assumed by the CFD model.
//   (c) smooth normals — per-vertex normals are the area-weighted average of incident face
//       normals (summing UN-normalised cross products naturally area-weights, then normalise).
//   (d) CAD edges — BRep_Tool::PolygonOnTriangulation maps each TopoDS_Edge to face-local
//       tessellation node pairs. Those pairs are attached to the final, orientation-adjusted
//       triangle half-edges; ordinary interior tessellation edges remain explicitly unlabelled.
//
// load_step_mesh() accumulates every face, including open shells and internal fabric, while
// deliberately ignoring standalone STEP wires/edges.
#include "core/geometry/step_import.h"
#include "core/geometry/step_face_uv_projection.h"

#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBndLib.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <Bnd_Box.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IntTools_CommonPrt.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <IntTools_EdgeFace.hxx>
#include <IntTools_Range.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_PolygonOnTriangulation.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_ErrorHandler.hxx>
#include <Standard_Failure.hxx>
#include <StepData_StepModel.hxx>
#include <StepRepr_RepresentationItem.hxx>
#include <TCollection_HAsciiString.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopAbs_State.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <Transfer_TransientProcess.hxx>
#include <TransferBRep.hxx>
#include <XSControl_TransferReader.hxx>
#include <XSControl_WorkSession.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

namespace paracfd::core
{
	namespace
	{
		constexpr double kMmToM = 0.001;

		using ShapeIndexedMap = NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>;
		using ShapeAncestorMap = NCollection_IndexedDataMap<TopoDS_Shape,
			NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>;

		struct CadEdgeProvenance
		{
			std::uint32_t id = TriMesh::kNoCadEdgeId;
			std::uint32_t incident_face_count = 0;
			double tolerance_m = 0.0;
			bool periodic_seam = false;
			bool ambiguous = false;
			std::uint64_t contact_id = TriMesh::kNoCadContactId;
			std::uint32_t certified_fan_degree = 0;
		};

		struct CadEdgeContactCertificate
		{
			std::uint64_t contact_id = TriMesh::kNoCadContactId;
			std::uint32_t fan_degree = 0;
			// Zero-based IDs into source_faces. These are exact OCCT full edge-on-face
			// contacts, not tessellation/proximity matches.
			std::vector<std::uint32_t> target_face_ids;
			// Parallel to target_face_ids: one sector for a target trim boundary,
			// two when the contact lies in the target face interior.
			std::vector<std::uint8_t> target_sector_counts;
		};

		using FaceNodePairMap = std::unordered_map<std::uint64_t, CadEdgeProvenance>;

		bool same_polygon_chain(const Handle(Poly_PolygonOnTriangulation)& a,
			const Handle(Poly_PolygonOnTriangulation)& b)
		{
			if (a.IsNull() || b.IsNull() || a->NbNodes() != b->NbNodes()) return false;
			bool same = true, reversed = true;
			for (Standard_Integer node = 1; node <= a->NbNodes(); ++node)
			{
				same = same && a->Node(node) == b->Node(node);
				reversed = reversed && a->Node(node) == b->Node(a->NbNodes() + 1 - node);
			}
			return same || reversed;
		}

		std::uint64_t undirected_node_pair(Standard_Integer a, Standard_Integer b)
		{
			const std::uint32_t lo = static_cast<std::uint32_t>(std::min(a, b));
			const std::uint32_t hi = static_cast<std::uint32_t>(std::max(a, b));
			return (static_cast<std::uint64_t>(lo) << 32) | hi;
		}

		void insert_face_node_pair(FaceNodePairMap& pairs, Standard_Integer a, Standard_Integer b,
			std::uint32_t edge_id, std::uint32_t incident_face_count, double tolerance_m,
			bool periodic_seam, const CadEdgeContactCertificate& contact)
		{
			if (a <= 0 || b <= 0 || a == b) return;
			const std::uint64_t key = undirected_node_pair(a, b);
			auto [it, inserted] = pairs.emplace(key,
				CadEdgeProvenance{ edge_id, incident_face_count, tolerance_m, periodic_seam, false,
					contact.contact_id, contact.fan_degree });
			if (!inserted && !it->second.ambiguous
				&& (it->second.id != edge_id || it->second.incident_face_count != incident_face_count))
			{
				// Two different BRep edges claim the same tessellation segment. Do not guess:
				// unavailable provenance is safer than assigning the segment to either edge.
				it->second = CadEdgeProvenance{ TriMesh::kNoCadEdgeId, 0u, 0.0, false, true,
					TriMesh::kNoCadContactId, 0u };
			}
			else if (!inserted && !it->second.ambiguous)
			{
				it->second.tolerance_m = std::max(it->second.tolerance_m, tolerance_m);
				it->second.periodic_seam = it->second.periodic_seam || periodic_seam;
			}
		}

		void insert_unknown_face_node_pair(FaceNodePairMap& pairs, Standard_Integer a,
			Standard_Integer b)
		{
			if (a <= 0 || b <= 0 || a == b) return;
			pairs[undirected_node_pair(a, b)] = CadEdgeProvenance{
				TriMesh::kNoCadEdgeId, 0u, 0.0, false, true, TriMesh::kNoCadContactId, 0u };
		}

		bool is_owner_face(const ShapeAncestorMap& edge_faces,const TopoDS_Edge& edge,
			const TopoDS_Face& face)
		{
			if(!edge_faces.Contains(edge))return false;
			const auto& owners=edge_faces.FindFromKey(edge);
			for(NCollection_List<TopoDS_Shape>::Iterator owner(owners);owner.More();owner.Next())
				if(face.IsSame(owner.Value()))return true;
			return false;
		}

		struct ContactCertificationStats
		{
			std::uint64_t free_edges=0;
			std::uint64_t candidate_pairs=0;
			std::uint64_t inttools_calls=0;
			std::uint64_t midpoint_classifications=0;
			std::uint64_t certified_pairs=0;
			double index_seconds=0.0;
			double query_seconds=0.0;
			double inttools_seconds=0.0;
			double accepted_inttools_seconds=0.0;
			double classifier_seconds=0.0;
		};

		bool full_edge_common_with_face(const TopoDS_Edge& edge,const TopoDS_Face& face,
			double fuzzy_tolerance,std::uint32_t& added_sectors,ContactCertificationStats* stats)
		{
			BRepAdaptor_Curve curve(edge);
			const double first=curve.FirstParameter(),last=curve.LastParameter();
			if(!std::isfinite(first)||!std::isfinite(last)||!(last>first))return false;

			IntTools_EdgeFace intersection;
			intersection.SetEdge(edge);intersection.SetFace(face);intersection.SetRange(first,last);
			intersection.SetFuzzyValue(fuzzy_tolerance);
			const auto inttools_begin=std::chrono::steady_clock::now();intersection.Perform();
			const double inttools_elapsed=std::chrono::duration<double>(
				std::chrono::steady_clock::now()-inttools_begin).count();
			if(stats)
			{
				++stats->inttools_calls;
				stats->inttools_seconds+=inttools_elapsed;
			}
			if(!intersection.IsDone())return false;
			std::vector<std::pair<double,double>> ranges;
			for(NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
				intersection.CommonParts());common.More();common.Next())
			{
				if(common.Value().Type()!=TopAbs_EDGE)continue;
				double lo=0.0,hi=0.0;common.Value().Range1(lo,hi);if(hi<lo)std::swap(lo,hi);
				lo=std::max(lo,first);hi=std::min(hi,last);if(hi>lo)ranges.emplace_back(lo,hi);
			}
			if(ranges.empty())return false;
			std::sort(ranges.begin(),ranges.end());
			const double parameter_tolerance=std::max(Precision::PConfusion(),
				64.0*std::numeric_limits<double>::epsilon()*std::max({1.0,std::abs(first),std::abs(last)}));
			double covered=first;
			for(const auto& range:ranges)
			{
				if(range.first>covered+parameter_tolerance)return false;
				covered=std::max(covered,range.second);
			}
			if(covered<last-parameter_tolerance)return false;

			const gp_Pnt midpoint=curve.Value(0.5*(first+last));
			const auto classifier_begin=std::chrono::steady_clock::now();
			const BRepClass_FaceClassifier classifier(face,midpoint,fuzzy_tolerance,true);
			if(stats)
			{
				++stats->midpoint_classifications;
				stats->classifier_seconds+=std::chrono::duration<double>(
					std::chrono::steady_clock::now()-classifier_begin).count();
			}
			if(classifier.State()==TopAbs_IN)added_sectors=2;
			else if(classifier.State()==TopAbs_ON)added_sectors=1;
			else return false;
			if(stats)stats->accepted_inttools_seconds+=inttools_elapsed;
			return true;
		}

		struct NumericAabb
		{
			double lo[3]={std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),
				std::numeric_limits<double>::max()};
			double hi[3]={-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max(),
				-std::numeric_limits<double>::max()};
		};

		bool numeric_aabb(const Bnd_Box& source,NumericAabb& result)
		{
			if(source.IsVoid()||source.IsWhole())return false;
			source.Get(result.lo[0],result.lo[1],result.lo[2],result.hi[0],result.hi[1],result.hi[2]);
			return true;
		}

		bool aabb_overlaps(const NumericAabb& a,const NumericAabb& b)
		{
			for(int axis=0;axis<3;++axis)if(a.hi[axis]<b.lo[axis]||b.hi[axis]<a.lo[axis])return false;
			return true;
		}

		void add_aabb(NumericAabb& destination,const NumericAabb& source)
		{
			for(int axis=0;axis<3;++axis)
			{
				destination.lo[axis]=std::min(destination.lo[axis],source.lo[axis]);
				destination.hi[axis]=std::max(destination.hi[axis],source.hi[axis]);
			}
		}

		struct FaceTriangleAabb
		{
			NumericAabb bounds;
			std::uint32_t face_id=0;
		};

		struct FaceTriangleAabbNode
		{
			NumericAabb bounds;
			std::uint32_t begin=0;
			std::uint32_t count=0;
			int left=-1;
			int right=-1;
		};

		// A face's single whole-shape AABB is extremely loose for a curved canopy and
		// admitted thousands of unrelated edge/face pairs.  The already-created OCCT
		// tessellation supplies a much tighter spatial broad phase.  Each triangle box
		// is enlarged by its recorded meshing deflection plus the CAD face tolerance;
		// therefore an exact point on the represented face remains a candidate.  This
		// BVH never certifies contact: IntTools_EdgeFace below remains the sole exact
		// decision, so nearby fabric is not joined.
		class FaceTriangleAabbIndex
		{
		public:
			explicit FaceTriangleAabbIndex(const std::vector<TopoDS_Face>& faces)
			{
				for(std::size_t face_id=0;face_id<faces.size();++face_id)
				{
					const TopoDS_Face& face=faces[face_id];
					TopLoc_Location location;
					const Handle(Poly_Triangulation) triangulation=BRep_Tool::Triangulation(face,location);
					if(triangulation.IsNull()||triangulation->NbTriangles()<=0)
					{
						Bnd_Box bounds;BRepBndLib::AddOptimal(face,bounds,false,true);
						NumericAabb numeric;
						if(numeric_aabb(bounds,numeric))fallback_.push_back({numeric,
							static_cast<std::uint32_t>(face_id)});
						else unconditional_faces_.push_back(static_cast<std::uint32_t>(face_id));
						continue;
					}
					const gp_Trsf transform=location.Transformation();
					const double scale=std::abs(transform.ScaleFactor());
					const double deflection=std::isfinite(triangulation->Deflection())
						?std::max(0.0,triangulation->Deflection())*scale:0.0;
					const double padding=deflection+std::max(0.0,BRep_Tool::Tolerance(face))*scale
						+Precision::Confusion();
					for(Standard_Integer triangle_id=1;triangle_id<=triangulation->NbTriangles();++triangle_id)
					{
						Standard_Integer node[3];triangulation->Triangle(triangle_id).Get(node[0],node[1],node[2]);
						FaceTriangleAabb primitive;primitive.face_id=static_cast<std::uint32_t>(face_id);
						for(int corner=0;corner<3;++corner)
						{
							gp_Pnt point=triangulation->Node(node[corner]);point.Transform(transform);
							const double coordinate[3]={point.X(),point.Y(),point.Z()};
							for(int axis=0;axis<3;++axis)
							{
								primitive.bounds.lo[axis]=std::min(primitive.bounds.lo[axis],coordinate[axis]-padding);
								primitive.bounds.hi[axis]=std::max(primitive.bounds.hi[axis],coordinate[axis]+padding);
							}
						}
						triangles_.push_back(std::move(primitive));
					}
				}
				nodes_.reserve(triangles_.empty()?0:2*triangles_.size());
				if(!triangles_.empty())root_=build(0,static_cast<std::uint32_t>(triangles_.size()));
			}

			std::vector<std::uint32_t> query(const Bnd_Box& edge_bounds)const
			{
				NumericAabb edge;
				if(numeric_aabb(edge_bounds,edge))return query(std::vector<NumericAabb>{edge});
				std::vector<std::uint32_t> result=unconditional_faces_;
				for(const auto& triangle:triangles_)result.push_back(triangle.face_id);
				for(const auto& fallback:fallback_)result.push_back(fallback.face_id);
				std::sort(result.begin(),result.end());
				result.erase(std::unique(result.begin(),result.end()),result.end());
				return result;
			}

			std::vector<std::uint32_t> query(const std::vector<NumericAabb>& edge_segments)const
			{
				std::vector<std::uint32_t> result;
				bool first_segment=true;
				for(const NumericAabb& edge:edge_segments)
				{
					std::vector<std::uint32_t> segment_faces=unconditional_faces_;
					if(root_>=0)query_node(root_,edge,segment_faces);
					for(const auto& fallback:fallback_)if(aabb_overlaps(edge,fallback.bounds))
						segment_faces.push_back(fallback.face_id);
					std::sort(segment_faces.begin(),segment_faces.end());
					segment_faces.erase(std::unique(segment_faces.begin(),segment_faces.end()),
						segment_faces.end());
					if(first_segment)
					{
						result=std::move(segment_faces);first_segment=false;
					}
					else
					{
						std::vector<std::uint32_t> common;
						common.reserve(std::min(result.size(),segment_faces.size()));
						std::set_intersection(result.begin(),result.end(),segment_faces.begin(),
							segment_faces.end(),std::back_inserter(common));
						result=std::move(common);
						if(result.empty())break;
					}
				}
				return result;
			}

		private:
			int build(std::uint32_t begin,std::uint32_t end)
			{
				FaceTriangleAabbNode node;node.begin=begin;node.count=end-begin;
				NumericAabb centres;
				for(std::uint32_t index=begin;index<end;++index)
				{
					add_aabb(node.bounds,triangles_[index].bounds);
					for(int axis=0;axis<3;++axis)
					{
						const double centre=0.5*(triangles_[index].bounds.lo[axis]+triangles_[index].bounds.hi[axis]);
						centres.lo[axis]=std::min(centres.lo[axis],centre);
						centres.hi[axis]=std::max(centres.hi[axis],centre);
					}
				}
				const int node_id=static_cast<int>(nodes_.size());nodes_.push_back(node);
				if(end-begin<=8)return node_id;
				int axis=0;for(int candidate=1;candidate<3;++candidate)
					if(centres.hi[candidate]-centres.lo[candidate]>centres.hi[axis]-centres.lo[axis])axis=candidate;
				const std::uint32_t middle=begin+(end-begin)/2;
				std::nth_element(triangles_.begin()+begin,triangles_.begin()+middle,triangles_.begin()+end,
					[axis](const FaceTriangleAabb& a,const FaceTriangleAabb& b)
					{return a.bounds.lo[axis]+a.bounds.hi[axis]<b.bounds.lo[axis]+b.bounds.hi[axis];});
				nodes_[node_id].left=build(begin,middle);nodes_[node_id].right=build(middle,end);
				nodes_[node_id].count=0;return node_id;
			}

			void query_node(int node_id,const NumericAabb& edge,std::vector<std::uint32_t>& result)const
			{
				const FaceTriangleAabbNode& node=nodes_[node_id];
				if(!aabb_overlaps(edge,node.bounds))return;
				if(node.count>0)
				{
					for(std::uint32_t offset=0;offset<node.count;++offset)
					{
						const FaceTriangleAabb& triangle=triangles_[node.begin+offset];
						if(aabb_overlaps(edge,triangle.bounds))result.push_back(triangle.face_id);
					}
					return;
				}
				query_node(node.left,edge,result);query_node(node.right,edge,result);
			}

			int root_=-1;
			std::vector<FaceTriangleAabb> triangles_;
			std::vector<FaceTriangleAabbNode> nodes_;
			std::vector<FaceTriangleAabb> fallback_;
			std::vector<std::uint32_t> unconditional_faces_;
		};

		std::vector<NumericAabb> edge_polyline_aabbs(const TopoDS_Edge& edge,
			const TopoDS_Face& owner_face)
		{
			TopLoc_Location location;
			const Handle(Poly_Triangulation) triangulation=BRep_Tool::Triangulation(owner_face,location);
			if(triangulation.IsNull())return {};
			Handle(Poly_PolygonOnTriangulation) polygon=
				BRep_Tool::PolygonOnTriangulation(edge,triangulation,location);
			if(polygon.IsNull()||polygon->NbNodes()<2)return {};
			const gp_Trsf transform=location.Transformation();
			const double scale=std::abs(transform.ScaleFactor());
			const double polygon_deflection=std::isfinite(polygon->Deflection())
				?std::max(0.0,polygon->Deflection())*scale:0.0;
			const double mesh_deflection=std::isfinite(triangulation->Deflection())
				?std::max(0.0,triangulation->Deflection())*scale:0.0;
			const double padding=std::max(polygon_deflection,mesh_deflection)
				+std::max({0.0,BRep_Tool::Tolerance(edge),BRep_Tool::Tolerance(owner_face)})*scale
				+Precision::Confusion();
			std::vector<NumericAabb> result;
			result.reserve(static_cast<std::size_t>(polygon->NbNodes()-1));
			for(Standard_Integer segment=1;segment<polygon->NbNodes();++segment)
			{
				NumericAabb bounds;
				for(const Standard_Integer polygon_node:{polygon->Node(segment),polygon->Node(segment+1)})
				{
					gp_Pnt point=triangulation->Node(polygon_node);point.Transform(transform);
					const double coordinate[3]={point.X(),point.Y(),point.Z()};
					for(int axis=0;axis<3;++axis)
					{
						bounds.lo[axis]=std::min(bounds.lo[axis],coordinate[axis]-padding);
						bounds.hi[axis]=std::max(bounds.hi[axis],coordinate[axis]+padding);
					}
				}
				result.push_back(bounds);
			}
			return result;
		}

		struct EdgeFaceCandidate
		{
			TopoDS_Edge edge;
			TopoDS_Face face;
			std::uint32_t edge_id=0;
			std::uint32_t face_id=0;
			double fuzzy=0.0;
		};

		struct EdgeFaceCandidateResult
		{
			bool common=false;
			std::uint32_t added_sectors=0;
			ContactCertificationStats stats;
			std::exception_ptr error;
		};

		unsigned contact_worker_count(std::size_t jobs)
		{
			if(jobs==0)return 0;
			unsigned limit=8;
			if(const char* configured=std::getenv("PARACFD_STEP_CONTACT_THREADS"))
			{
				char* end=nullptr;const unsigned long parsed=std::strtoul(configured,&end,10);
				if(end!=configured&&*end=='\0'&&parsed>0)limit=static_cast<unsigned>(
					std::min<unsigned long>(parsed,32));
			}
			const unsigned hardware=std::max(1u,std::thread::hardware_concurrency());
			return static_cast<unsigned>(std::min<std::size_t>(jobs,std::min(limit,hardware)));
		}

		void add_contact_stats(ContactCertificationStats& destination,
			const ContactCertificationStats& source)
		{
			destination.inttools_calls+=source.inttools_calls;
			destination.midpoint_classifications+=source.midpoint_classifications;
			destination.inttools_seconds+=source.inttools_seconds;
			destination.accepted_inttools_seconds+=source.accepted_inttools_seconds;
			destination.classifier_seconds+=source.classifier_seconds;
		}

		std::uint64_t contact_certificate_hash(
			const std::vector<CadEdgeContactCertificate>& certificates)
		{
			std::uint64_t hash=1469598103934665603ull;
			auto append=[&](std::uint64_t value)
			{
				for(unsigned byte=0;byte<8;++byte)
				{
					hash^=(value>>(8*byte))&0xffu;hash*=1099511628211ull;
				}
			};
			append(certificates.size());
			for(std::size_t edge=0;edge<certificates.size();++edge)
			{
				const CadEdgeContactCertificate& certificate=certificates[edge];
				append(edge);append(certificate.contact_id);append(certificate.fan_degree);
				append(certificate.target_face_ids.size());
				for(const std::uint32_t face:certificate.target_face_ids)append(face);
			}
			return hash;
		}

		std::vector<CadEdgeContactCertificate> certify_edge_face_contacts(
			const ShapeIndexedMap& shape_edges,const ShapeAncestorMap& edge_faces,
			const std::vector<TopoDS_Face>& faces)
		{
			std::vector<CadEdgeContactCertificate> result(
				static_cast<std::size_t>(shape_edges.Extent()));
			ContactCertificationStats measured_stats;
			ContactCertificationStats* stats=std::getenv("PARACFD_STEP_CONTACT_STATS")?&measured_stats:nullptr;
			const auto index_begin=std::chrono::steady_clock::now();
			const FaceTriangleAabbIndex face_index(faces);
			std::vector<Bnd_Box> exact_face_bounds(faces.size());
			for(std::size_t face=0;face<faces.size();++face)
				BRepBndLib::AddOptimal(faces[face],exact_face_bounds[face],false,true);
			if(stats)stats->index_seconds=std::chrono::duration<double>(
				std::chrono::steady_clock::now()-index_begin).count();
			std::vector<EdgeFaceCandidate> candidates;
			std::vector<std::uint32_t> fan_degree(static_cast<std::size_t>(shape_edges.Extent()),0u);
			for(Standard_Integer edge_index=1;edge_index<=shape_edges.Extent();++edge_index)
			{
				const TopoDS_Edge edge=TopoDS::Edge(shape_edges(edge_index));
				const std::uint32_t owner_sectors=edge_faces.Contains(edge)
					?static_cast<std::uint32_t>(edge_faces.FindFromKey(edge).Size()):0u;
				// Shared/manifold and periodic edges already carry exact TopoDS incidence.
				// The missing information is a nominally free edge terminating on a face
				// that does not own the same topological edge.
				if(owner_sectors!=1)continue;
				const std::uint32_t edge_id=static_cast<std::uint32_t>(edge_index-1);
				fan_degree[edge_id]=owner_sectors;
				if(stats)++stats->free_edges;
				Bnd_Box edge_bounds;BRepBndLib::AddOptimal(edge,edge_bounds,false,true);
				const TopoDS_Face owner_face=TopoDS::Face(edge_faces.FindFromKey(edge).First());
				const std::vector<NumericAabb> segment_bounds=edge_polyline_aabbs(edge,owner_face);
				const auto query_begin=std::chrono::steady_clock::now();
				const std::vector<std::uint32_t> candidate_faces=segment_bounds.empty()
					?face_index.query(edge_bounds):face_index.query(segment_bounds);
				if(stats)
				{
					stats->query_seconds+=std::chrono::duration<double>(
						std::chrono::steady_clock::now()-query_begin).count();
					stats->candidate_pairs+=candidate_faces.size();
				}
				for(const std::uint32_t face_id:candidate_faces)
				{
					const TopoDS_Face& face=faces[face_id];
					if(is_owner_face(edge_faces,edge,face)
						||edge_bounds.IsOut(exact_face_bounds[face_id]))continue;
					const double fuzzy=std::max({Precision::Confusion(),BRep_Tool::Tolerance(edge),
						BRep_Tool::Tolerance(face)});
					candidates.push_back({edge,face,edge_id,face_id,fuzzy});
				}
			}

			std::vector<EdgeFaceCandidateResult> candidate_results(candidates.size());
			const unsigned worker_count=contact_worker_count(candidates.size());
			const auto exact_begin=std::chrono::steady_clock::now();
			std::atomic<std::size_t> next_candidate{0};
			auto worker=[&]()
			{
				for(;;)
				{
					const std::size_t candidate_id=next_candidate.fetch_add(1,std::memory_order_relaxed);
					if(candidate_id>=candidates.size())return;
					const EdgeFaceCandidate& candidate=candidates[candidate_id];
					EdgeFaceCandidateResult& candidate_result=candidate_results[candidate_id];
					try
					{
						candidate_result.common=full_edge_common_with_face(candidate.edge,candidate.face,
							candidate.fuzzy,candidate_result.added_sectors,
							stats?&candidate_result.stats:nullptr);
					}
					catch(...){candidate_result.error=std::current_exception();}
				}
			};
			std::vector<std::thread> workers;workers.reserve(worker_count);
			for(unsigned thread=0;thread<worker_count;++thread)workers.emplace_back(worker);
			for(std::thread& thread:workers)thread.join();
			const double exact_wall_seconds=std::chrono::duration<double>(
				std::chrono::steady_clock::now()-exact_begin).count();

			// Serial, index-ordered reduction makes contact IDs, fan degrees, exception
			// propagation and target face ordering independent of worker scheduling.
			for(std::size_t candidate_id=0;candidate_id<candidates.size();++candidate_id)
			{
				const EdgeFaceCandidate& candidate=candidates[candidate_id];
				const EdgeFaceCandidateResult& candidate_result=candidate_results[candidate_id];
				if(candidate_result.error)std::rethrow_exception(candidate_result.error);
				if(stats)add_contact_stats(*stats,candidate_result.stats);
				if(!candidate_result.common)continue;
				if(stats)++stats->certified_pairs;
				fan_degree[candidate.edge_id]+=candidate_result.added_sectors;
				result[candidate.edge_id].target_face_ids.push_back(candidate.face_id);
				result[candidate.edge_id].target_sector_counts.push_back(static_cast<std::uint8_t>(
					candidate_result.added_sectors));
			}
			for(std::size_t edge_id=0;edge_id<fan_degree.size();++edge_id)
				if(fan_degree[edge_id]>1)
				{
					auto& certificate=result[edge_id];certificate.contact_id=edge_id;
					certificate.fan_degree=fan_degree[edge_id];
					std::vector<std::pair<std::uint32_t,std::uint8_t>> targets;
					targets.reserve(certificate.target_face_ids.size());
					for(std::size_t target=0;target<certificate.target_face_ids.size();++target)
						targets.emplace_back(certificate.target_face_ids[target],
							certificate.target_sector_counts[target]);
					std::sort(targets.begin(),targets.end());
					certificate.target_face_ids.clear();certificate.target_sector_counts.clear();
					for(const auto& target:targets)
					{
						if(!certificate.target_face_ids.empty()
							&&certificate.target_face_ids.back()==target.first)
						{
							certificate.target_sector_counts.back()=std::max(
								certificate.target_sector_counts.back(),target.second);
							continue;
						}
						certificate.target_face_ids.push_back(target.first);
						certificate.target_sector_counts.push_back(target.second);
					}
				}
			if(stats)std::fprintf(stderr,"[STEP contact] free-edges=%llu candidates=%llu "
				"exact-calls=%llu midpoint-tests=%llu certified=%llu workers=%u; index=%.3fs query=%.3fs "
				"IntTools=%.3fs (accepted %.3fs) classifier=%.3fs\n",
				static_cast<unsigned long long>(stats->free_edges),
				static_cast<unsigned long long>(stats->candidate_pairs),
				static_cast<unsigned long long>(stats->inttools_calls),
				static_cast<unsigned long long>(stats->midpoint_classifications),
				static_cast<unsigned long long>(stats->certified_pairs),worker_count,stats->index_seconds,
				stats->query_seconds,stats->inttools_seconds,stats->accepted_inttools_seconds,
				stats->classifier_seconds);
			if(stats)std::fprintf(stderr,"[STEP contact] exact wall %.3fs; aggregate exact CPU %.3fs; "
				"certificate hash=%016llx\n",exact_wall_seconds,
				stats->inttools_seconds+stats->classifier_seconds,
				static_cast<unsigned long long>(contact_certificate_hash(result)));
			return result;
		}

		int face_index_of(const std::vector<TopoDS_Face>& faces,const TopoDS_Shape& candidate)
		{
			for(std::size_t face=0;face<faces.size();++face)
				if(faces[face].IsSame(candidate))return static_cast<int>(face);
			return -1;
		}

		bool parameter_ranges_cover(std::vector<std::pair<double,double>> ranges,
			double first,double last)
		{
			if(last<first)std::swap(first,last);
			if(ranges.empty())return false;
			for(auto& range:ranges)if(range.second<range.first)std::swap(range.first,range.second);
			std::sort(ranges.begin(),ranges.end());
			const double tolerance=std::max(Precision::PConfusion(),
				64.0*std::numeric_limits<double>::epsilon()
				*std::max({1.0,std::abs(first),std::abs(last)}));
			double covered=first;
			for(const auto& range:ranges)
			{
				if(range.second<first||range.first>last)continue;
				if(range.first>covered+tolerance)return false;
				covered=std::max(covered,std::min(last,range.second));
			}
			return covered>=last-tolerance;
		}

		struct ExactEdgeOverlap
		{
			bool has_positive_length_common = false;
			bool covers_a = false;
			bool covers_b = false;
		};

		// This is the authoritative reciprocal/partial-edge decision. Endpoint and AABB
		// tests below are only a broad phase; no distance-only match can create a physical
		// contact curve.
		ExactEdgeOverlap exact_edge_overlap(const TopoDS_Edge& a,const TopoDS_Edge& b)
		{
			ExactEdgeOverlap result;
			BRepAdaptor_Curve curve_a(a),curve_b(b);
			const double first_a=curve_a.FirstParameter(),last_a=curve_a.LastParameter();
			const double first_b=curve_b.FirstParameter(),last_b=curve_b.LastParameter();
			if(!std::isfinite(first_a)||!std::isfinite(last_a)||!(last_a>first_a)
				||!std::isfinite(first_b)||!std::isfinite(last_b)||!(last_b>first_b))return result;
			IntTools_EdgeEdge intersection(a,b);
			intersection.SetFuzzyValue(0.0);
			intersection.Perform();
			if(!intersection.IsDone())return result;
			std::vector<std::pair<double,double>> ranges_a,ranges_b;
			for(NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
				intersection.CommonParts());common.More();common.Next())
			{
				if(common.Value().Type()!=TopAbs_EDGE)continue;
				result.has_positive_length_common=true;
				double lo=0.0,hi=0.0;common.Value().Range1(lo,hi);
				ranges_a.emplace_back(lo,hi);
				for(NCollection_Sequence<IntTools_Range>::Iterator range(
					common.Value().Ranges2());range.More();range.Next())
					ranges_b.emplace_back(range.Value().First(),range.Value().Last());
			}
			result.covers_a=parameter_ranges_cover(std::move(ranges_a),first_a,last_a);
			result.covers_b=parameter_ranges_cover(std::move(ranges_b),first_b,last_b);
			return result;
		}

		double edge_native_tolerance_mm(const TopoDS_Edge& edge)
		{
			return std::max(0.0,BRep_Tool::Tolerance(edge))
				*std::abs(edge.Location().Transformation().ScaleFactor());
		}

		double vertex_native_tolerance_mm(const TopoDS_Vertex& vertex)
		{
			return vertex.IsNull()?0.0:std::max(0.0,BRep_Tool::Tolerance(vertex))
				*std::abs(vertex.Location().Transformation().ScaleFactor());
		}

		struct CanonicalContactEdge
		{
			std::uint32_t edge_id=TriMesh::kNoCadEdgeId;
			std::uint32_t owner_face_id=0;
			TopoDS_Edge edge;
			std::array<gp_Pnt,2> endpoints;
			std::array<double,2> endpoint_tolerance_mm{{0.0,0.0}};
			double edge_tolerance_mm=0.0;
			NumericAabb bounds;
			bool has_bounds=false;
		};

		bool endpoint_pair_within_native_tolerance(const gp_Pnt& a,double tolerance_a,
			const gp_Pnt& b,double tolerance_b)
		{
			const double tolerance=tolerance_a+tolerance_b;
			return a.SquareDistance(b)<=tolerance*tolerance;
		}

		bool reciprocal_edge_broad_phase(const CanonicalContactEdge& a,
			const CanonicalContactEdge& b)
		{
			const bool forward=endpoint_pair_within_native_tolerance(a.endpoints[0],
				a.endpoint_tolerance_mm[0],b.endpoints[0],b.endpoint_tolerance_mm[0])
				&&endpoint_pair_within_native_tolerance(a.endpoints[1],
					a.endpoint_tolerance_mm[1],b.endpoints[1],b.endpoint_tolerance_mm[1]);
			const bool reversed=endpoint_pair_within_native_tolerance(a.endpoints[0],
				a.endpoint_tolerance_mm[0],b.endpoints[1],b.endpoint_tolerance_mm[1])
				&&endpoint_pair_within_native_tolerance(a.endpoints[1],
					a.endpoint_tolerance_mm[1],b.endpoints[0],b.endpoint_tolerance_mm[0]);
			return forward||reversed;
		}

		bool point_in_native_edge_bounds(const gp_Pnt& point,double point_tolerance_mm,
			const CanonicalContactEdge& edge)
		{
			if(!edge.has_bounds)return true;
			const double tolerance=point_tolerance_mm+edge.edge_tolerance_mm;
			const double coordinate[3]={point.X(),point.Y(),point.Z()};
			for(int axis=0;axis<3;++axis)
				if(coordinate[axis]<edge.bounds.lo[axis]-tolerance
					||coordinate[axis]>edge.bounds.hi[axis]+tolerance)return false;
			return true;
		}

		bool endpoint_junction_broad_phase(const CanonicalContactEdge& a,
			const CanonicalContactEdge& b)
		{
			for(std::size_t endpoint=0;endpoint<2;++endpoint)
				if(point_in_native_edge_bounds(a.endpoints[endpoint],
					a.endpoint_tolerance_mm[endpoint],b)
					||point_in_native_edge_bounds(b.endpoints[endpoint],
						b.endpoint_tolerance_mm[endpoint],a))return true;
			return false;
		}

		bool endpoint_lies_on_edge_within_native_tolerance(const gp_Pnt& point,
			double point_tolerance_mm,const CanonicalContactEdge& edge)
		{
			if(!point_in_native_edge_bounds(point,point_tolerance_mm,edge))return false;
			TopLoc_Location location;Standard_Real first=0.0,last=0.0;
			const Handle(Geom_Curve) curve=BRep_Tool::Curve(edge.edge,location,first,last);
			if(curve.IsNull()||!std::isfinite(first)||!std::isfinite(last)||!(last>first))
				return true; // indeterminate broad phase must retain the exact candidate
			gp_Pnt local=point;local.Transform(location.Transformation().Inverted());
			GeomAPI_ProjectPointOnCurve projection(local,curve,first,last);
			if(projection.NbPoints()<=0)return true;
			gp_Pnt projected=curve->Value(projection.LowerDistanceParameter());
			projected.Transform(location.Transformation());
			const double tolerance=point_tolerance_mm+edge.edge_tolerance_mm;
			return projected.SquareDistance(point)<=tolerance*tolerance;
		}

		// A positive-length overlap which is not full/full must terminate at an endpoint
		// of at least one participating trimmed edge. Requiring that endpoint to lie on the
		// other exact support curve is a conservative cheap prerequisite for IntTools.
		bool partial_overlap_broad_phase(const CanonicalContactEdge& a,
			const CanonicalContactEdge& b)
		{
			for(std::size_t endpoint=0;endpoint<2;++endpoint)
				if(endpoint_lies_on_edge_within_native_tolerance(a.endpoints[endpoint],
					a.endpoint_tolerance_mm[endpoint],b)
					||endpoint_lies_on_edge_within_native_tolerance(b.endpoints[endpoint],
						b.endpoint_tolerance_mm[endpoint],a))return true;
			return false;
		}

		struct ContactDisjointSet
		{
			explicit ContactDisjointSet(std::size_t count):parent(count),rank(count,0)
			{
				std::iota(parent.begin(),parent.end(),std::size_t{0});
			}
			std::size_t find(std::size_t value)
			{
				if(parent[value]!=value)parent[value]=find(parent[value]);
				return parent[value];
			}
			void merge(std::size_t a,std::size_t b)
			{
				a=find(a);b=find(b);if(a==b)return;
				if(rank[a]<rank[b])std::swap(a,b);
				parent[b]=a;if(rank[a]==rank[b])++rank[a];
			}
			std::vector<std::size_t> parent;
			std::vector<std::uint8_t> rank;
		};

		struct ProjectedContactSample
		{
			double parameter=0.0;
			gp_Pnt point;
			double tolerance_mm=0.0;
			std::uint32_t source_edge_id=TriMesh::kNoCadEdgeId;
			std::uint32_t source_node_id=0;
			std::vector<std::size_t> junction_ids;
		};

		struct ExactContactJunction
		{
			std::size_t first_group=0;
			std::size_t second_group=0;
			gp_Pnt point;
		};

		std::vector<gp_Pnt> exact_edge_vertex_junctions(const TopoDS_Edge& a,
			const TopoDS_Edge& b)
		{
			std::vector<gp_Pnt> result;
			IntTools_EdgeEdge intersection(a,b);intersection.SetFuzzyValue(0.0);
			intersection.Perform();if(!intersection.IsDone())return result;
			BRepAdaptor_Curve curve_a(a);
			for(NCollection_Sequence<IntTools_CommonPrt>::Iterator common(
				intersection.CommonParts());common.More();common.Next())
			{
				if(common.Value().Type()!=TopAbs_VERTEX)continue;
				const double parameter_a=common.Value().VertexParameter1();
				const double parameter_b=common.Value().VertexParameter2();
				if(!std::isfinite(parameter_a)||!std::isfinite(parameter_b))continue;
				result.push_back(curve_a.Value(parameter_a));
			}
			return result;
		}

		bool project_to_canonical_edge(const gp_Pnt& world_point,
			const Handle(Geom_Curve)& curve,const TopLoc_Location& location,
			double first,double last,double& parameter,gp_Pnt& projected_world)
		{
			if(curve.IsNull())return false;
			gp_Pnt local_point=world_point;
			local_point.Transform(location.Transformation().Inverted());
			GeomAPI_ProjectPointOnCurve projection(local_point,curve,first,last);
			if(projection.NbPoints()<=0)return false;
			parameter=std::clamp(projection.LowerDistanceParameter(),first,last);
			projected_world=curve->Value(parameter);
			projected_world.Transform(location.Transformation());
			return true;
		}

		std::vector<TopoDS_Edge> certified_boundary_occurrences(
			const CanonicalContactEdge& canonical,
			const std::vector<std::size_t>& members,
			const std::vector<CanonicalContactEdge>& contact_edges,
			std::uint32_t face_id,const std::vector<TopoDS_Face>& faces)
		{
			std::vector<TopoDS_Edge> result;
			if(face_id>=faces.size())return result;
			bool has_owned_member=false;
			for(const std::size_t member_id:members)
				if(contact_edges[member_id].owner_face_id==face_id)
					has_owned_member=true;
			for(TopExp_Explorer occurrence(faces[face_id],TopAbs_EDGE);
				occurrence.More();occurrence.Next())
			{
				const TopoDS_Edge oriented=TopoDS::Edge(occurrence.Current());
				if(has_owned_member)
				{
					for(const std::size_t member_id:members)
						if(contact_edges[member_id].owner_face_id==face_id
							&&oriented.IsSame(contact_edges[member_id].edge))
						{
							result.push_back(oriented);
							break;
						}
					continue;
				}
				// A target trim may be shared and therefore absent from the free-edge
				// physical group. It is usable only when both exact extents agree; a
				// longer/split trim still requires the separately scoped atomization step.
				const ExactEdgeOverlap overlap=exact_edge_overlap(canonical.edge,oriented);
				if(overlap.has_positive_length_common&&overlap.covers_a&&overlap.covers_b)
					result.push_back(oriented);
			}
			return result;
		}

		bool make_contact_graph(const ShapeIndexedMap& shape_edges,
			const ShapeAncestorMap& edge_faces,const std::vector<TopoDS_Face>& faces,
			std::vector<CadEdgeContactCertificate>& certificates,
			StepCadContactGraph& output,std::string& error)
		{
			error.clear();
			std::vector<CanonicalContactEdge> contact_edges;
			for(std::size_t edge_id=0;edge_id<certificates.size();++edge_id)
			{
				const CadEdgeContactCertificate& certificate=certificates[edge_id];
				if(certificate.contact_id==TriMesh::kNoCadContactId||certificate.fan_degree<=1
					||edge_id>=static_cast<std::size_t>(shape_edges.Extent()))continue;
				const TopoDS_Edge edge=TopoDS::Edge(shape_edges(static_cast<Standard_Integer>(edge_id+1)));
				if(!edge_faces.Contains(edge)||edge_faces.FindFromKey(edge).Size()!=1)continue;
				const TopoDS_Face owner=TopoDS::Face(edge_faces.FindFromKey(edge).First());
				const int owner_id=face_index_of(faces,owner);if(owner_id<0)continue;
				BRepAdaptor_Curve curve(edge);
				const double first=curve.FirstParameter(),last=curve.LastParameter();
				if(!std::isfinite(first)||!std::isfinite(last)||!(last>first))continue;
				CanonicalContactEdge record;
				record.edge_id=static_cast<std::uint32_t>(edge_id);
				record.owner_face_id=static_cast<std::uint32_t>(owner_id);
				record.edge=edge;record.endpoints={{curve.Value(first),curve.Value(last)}};
				record.edge_tolerance_mm=edge_native_tolerance_mm(edge);
				TopoDS_Vertex first_vertex,last_vertex;TopExp::Vertices(edge,first_vertex,last_vertex,true);
				record.endpoint_tolerance_mm={{std::max(record.edge_tolerance_mm,
					vertex_native_tolerance_mm(first_vertex)),std::max(record.edge_tolerance_mm,
					vertex_native_tolerance_mm(last_vertex))}};
				Bnd_Box bounds;BRepBndLib::AddOptimal(edge,bounds,false,true);
				record.has_bounds=numeric_aabb(bounds,record.bounds);
				contact_edges.push_back(std::move(record));
			}

			ContactDisjointSet groups(contact_edges.size());
			std::vector<StepCadContactGraph::PartialOverlap> unresolved_partial_overlaps;
			for(std::size_t a=0;a<contact_edges.size();++a)
				for(std::size_t b=a+1;b<contact_edges.size();++b)
				{
					const bool reciprocal=reciprocal_edge_broad_phase(contact_edges[a],contact_edges[b]);
					if(!reciprocal&&!partial_overlap_broad_phase(contact_edges[a],contact_edges[b]))
						continue;
					const ExactEdgeOverlap overlap=exact_edge_overlap(contact_edges[a].edge,
						contact_edges[b].edge);
					if(!overlap.has_positive_length_common)continue;
					if(overlap.covers_a&&overlap.covers_b)
					{
						groups.merge(a,b);
						continue;
					}
					unresolved_partial_overlaps.push_back({contact_edges[a].edge_id,
						contact_edges[b].edge_id,contact_edges[a].owner_face_id,
						contact_edges[b].owner_face_id,overlap.covers_a,overlap.covers_b});
				}
			std::map<std::size_t,std::vector<std::size_t>> members_by_group;
			for(std::size_t edge=0;edge<contact_edges.size();++edge)
				members_by_group[groups.find(edge)].push_back(edge);
			std::vector<std::vector<std::size_t>> physical_groups;
			physical_groups.reserve(members_by_group.size());
			for(auto& [unused_root,members]:members_by_group)
			{
				(void)unused_root;
				std::sort(members.begin(),members.end(),[&](std::size_t a,std::size_t b)
					{return contact_edges[a].edge_id<contact_edges[b].edge_id;});
				physical_groups.push_back(std::move(members));
			}
			std::sort(physical_groups.begin(),physical_groups.end(),[&](const auto& a,const auto& b)
				{return contact_edges[a.front()].edge_id<contact_edges[b.front()].edge_id;});

			// Curves which are not reciprocal copies can still meet at an exact CAD vertex or
			// T-node. Insert that OCCT-certified point into both canonical chains, then join
			// their discrete topology IDs below. The endpoint-vs-edge bounds test is only a
			// broad phase; a zero-fuzzy IntTools vertex common is the decision.
			std::vector<ExactContactJunction> junctions;
			std::vector<std::vector<std::size_t>> junctions_by_group(physical_groups.size());
			for(std::size_t a=0;a<physical_groups.size();++a)
				for(std::size_t b=a+1;b<physical_groups.size();++b)
				{
					const CanonicalContactEdge& edge_a=contact_edges[physical_groups[a].front()];
					const CanonicalContactEdge& edge_b=contact_edges[physical_groups[b].front()];
					if(!endpoint_junction_broad_phase(edge_a,edge_b))continue;
					for(const gp_Pnt& point:exact_edge_vertex_junctions(edge_a.edge,edge_b.edge))
					{
						const std::size_t junction_id=junctions.size();
						junctions.push_back({a,b,point});
						junctions_by_group[a].push_back(junction_id);
						junctions_by_group[b].push_back(junction_id);
					}
				}

			StepCadContactGraph graph;
			graph.unresolved_partial_overlaps=std::move(unresolved_partial_overlaps);
			// Kept parallel to graph.curves/curve.uses until shared junction coordinates
			// have been canonicalized. UVs must round-trip to the final public 3D chain,
			// not to a pre-collapse tolerance-near copy.
			std::vector<std::vector<std::vector<TopoDS_Edge>>> pending_boundary_edges;
			std::vector<std::vector<std::size_t>> topology_locations(junctions.size());
			std::size_t total_sample_count=0;
			for(std::size_t group_id=0;group_id<physical_groups.size();++group_id)
			{
				const std::vector<std::size_t>& members=physical_groups[group_id];
				const CanonicalContactEdge& canonical=contact_edges[members.front()];
				StepContactCurve curve;
				curve.id=canonical.edge_id;curve.source_edge_id=canonical.edge_id;
				std::map<std::uint32_t,std::uint8_t> use_sectors;
				for(const std::size_t member_id:members)
				{
					const CanonicalContactEdge& member=contact_edges[member_id];
					curve.source_edge_ids.push_back(member.edge_id);
					curve.tolerance_m=std::max(curve.tolerance_m,
						std::max({member.edge_tolerance_mm,member.endpoint_tolerance_mm[0],
							member.endpoint_tolerance_mm[1]})*kMmToM);
					use_sectors.emplace(member.owner_face_id,1u);
					const CadEdgeContactCertificate& certificate=certificates[member.edge_id];
					for(std::size_t target=0;target<certificate.target_face_ids.size();++target)
					{
						const std::uint32_t face_id=certificate.target_face_ids[target];
						const std::uint8_t sectors=target<certificate.target_sector_counts.size()
							?certificate.target_sector_counts[target]:1u;
						auto [found,inserted]=use_sectors.emplace(face_id,sectors);
						if(!inserted)found->second=std::max(found->second,sectors);
					}
				}
				for(const auto& [face_id,sectors]:use_sectors)
					if(face_id<faces.size())curve.fan_degree+=sectors;

				TopLoc_Location canonical_location;
				Standard_Real canonical_first=0.0,canonical_last=0.0;
				const Handle(Geom_Curve) canonical_curve=BRep_Tool::Curve(canonical.edge,
					canonical_location,canonical_first,canonical_last);
				std::vector<ProjectedContactSample> samples;
				auto append_projected=[&](const gp_Pnt& point,double tolerance_mm,
					std::uint32_t source_edge_id,std::uint32_t source_node_id,
					std::size_t junction_id=std::numeric_limits<std::size_t>::max())
				{
					double parameter=0.0;gp_Pnt projected;
					if(!project_to_canonical_edge(point,canonical_curve,canonical_location,
						canonical_first,canonical_last,parameter,projected))return;
					ProjectedContactSample sample;
					sample.parameter=parameter;sample.point=projected;
					sample.tolerance_mm=tolerance_mm;sample.source_edge_id=source_edge_id;
					sample.source_node_id=source_node_id;
					if(junction_id!=std::numeric_limits<std::size_t>::max())
						sample.junction_ids.push_back(junction_id);
					samples.push_back(std::move(sample));
				};
				BRepAdaptor_Curve canonical_adaptor(canonical.edge);
				append_projected(canonical_adaptor.Value(canonical_adaptor.FirstParameter()),
					canonical.endpoint_tolerance_mm[0],canonical.edge_id,0u);
				append_projected(canonical_adaptor.Value(canonical_adaptor.LastParameter()),
					canonical.endpoint_tolerance_mm[1],canonical.edge_id,
					std::numeric_limits<std::uint32_t>::max());
				for(const std::size_t member_id:members)
				{
					const CanonicalContactEdge& member=contact_edges[member_id];
					const TopoDS_Face& owner=faces[member.owner_face_id];
					TopLoc_Location location;
					const Handle(Poly_Triangulation) triangulation=BRep_Tool::Triangulation(owner,location);
					if(triangulation.IsNull())continue;
					Handle(Poly_PolygonOnTriangulation) polygon=
						BRep_Tool::PolygonOnTriangulation(member.edge,triangulation,location);
					if(polygon.IsNull())
					{
						TopoDS_Edge reversed=member.edge;reversed.Reverse();
						polygon=BRep_Tool::PolygonOnTriangulation(reversed,triangulation,location);
					}
					if(polygon.IsNull())continue;
					for(Standard_Integer node=1;node<=polygon->NbNodes();++node)
					{
						gp_Pnt point=triangulation->Node(polygon->Node(node));
						point.Transform(location.Transformation());
						double tolerance_mm=member.edge_tolerance_mm;
						if(node==1)tolerance_mm=std::max(tolerance_mm,
							member.endpoint_tolerance_mm[0]);
						if(node==polygon->NbNodes())tolerance_mm=std::max(tolerance_mm,
							member.endpoint_tolerance_mm[1]);
						append_projected(point,tolerance_mm,member.edge_id,
							static_cast<std::uint32_t>(node));
					}
				}
				for(const std::size_t junction_id:junctions_by_group[group_id])
				{
					const ExactContactJunction& junction=junctions[junction_id];
					const std::size_t other_group=junction.first_group==group_id
						?junction.second_group:junction.first_group;
					const CanonicalContactEdge& other=
						contact_edges[physical_groups[other_group].front()];
					append_projected(junction.point,
						std::max(canonical.edge_tolerance_mm,other.edge_tolerance_mm),
						canonical.edge_id,static_cast<std::uint32_t>(junction_id),junction_id);
				}
				std::sort(samples.begin(),samples.end(),[](const ProjectedContactSample& a,
					const ProjectedContactSample& b)
				{
					if(a.parameter!=b.parameter)return a.parameter<b.parameter;
					if(a.source_edge_id!=b.source_edge_id)return a.source_edge_id<b.source_edge_id;
					return a.source_node_id<b.source_node_id;
				});
				std::vector<ProjectedContactSample> unique_samples;
				for(const ProjectedContactSample& sample:samples)
				{
					if(!unique_samples.empty())
					{
						ProjectedContactSample& previous=unique_samples.back();
						const double tolerance=std::max(previous.tolerance_mm,sample.tolerance_mm);
						// Both samples were projected onto one exact canonical curve. Native
						// edge/vertex tolerance is the only permitted node de-duplication radius.
						if(previous.point.SquareDistance(sample.point)<=tolerance*tolerance)
						{
							previous.tolerance_mm=tolerance;
							previous.junction_ids.insert(previous.junction_ids.end(),
								sample.junction_ids.begin(),sample.junction_ids.end());
							std::sort(previous.junction_ids.begin(),previous.junction_ids.end());
							previous.junction_ids.erase(std::unique(previous.junction_ids.begin(),
								previous.junction_ids.end()),previous.junction_ids.end());
							continue;
						}
					}
					unique_samples.push_back(sample);
				}
				std::vector<gp_Pnt> world_points;world_points.reserve(unique_samples.size());
				for(const ProjectedContactSample& sample:unique_samples)
				{
					curve.source_parameters.push_back(sample.parameter);
					world_points.push_back(sample.point);
					curve.sample_positions_m.push_back({{sample.point.X()*kMmToM,
						sample.point.Y()*kMmToM,sample.point.Z()*kMmToM}});
				}
				const auto endpoint_less=[](const std::array<double,3>& a,
					const std::array<double,3>& b)
				{
					if(a[0]!=b[0])return a[0]<b[0];if(a[1]!=b[1])return a[1]<b[1];return a[2]<b[2];
				};
				if(curve.sample_positions_m.size()>=2&&endpoint_less(
					curve.sample_positions_m.back(),curve.sample_positions_m.front()))
				{
					std::reverse(unique_samples.begin(),unique_samples.end());
					std::reverse(world_points.begin(),world_points.end());
					std::reverse(curve.sample_positions_m.begin(),curve.sample_positions_m.end());
					std::reverse(curve.source_parameters.begin(),curve.source_parameters.end());
				}
				std::vector<std::vector<TopoDS_Edge>> curve_boundary_edges;
				for(const auto& [face_id,sectors]:use_sectors)
				{
					if(face_id>=faces.size())continue;
					StepContactFaceUse use;use.source_face_id=face_id;use.sector_count=sectors;
					use.kind=sectors>1?StepContactUseKind::face_interior:
						StepContactUseKind::trim_boundary;
					curve_boundary_edges.push_back(use.kind==StepContactUseKind::trim_boundary
						?certified_boundary_occurrences(canonical,members,contact_edges,
							face_id,faces):std::vector<TopoDS_Edge>{});
					curve.uses.push_back(std::move(use));
				}
				for(const std::size_t member_id:members)
				{
					CadEdgeContactCertificate& certificate=
						certificates[contact_edges[member_id].edge_id];
					certificate.contact_id=curve.id;certificate.fan_degree=curve.fan_degree;
				}
				for(std::size_t sample=0;sample<unique_samples.size();++sample)
					for(const std::size_t junction_id:unique_samples[sample].junction_ids)
						if(junction_id<topology_locations.size())
							topology_locations[junction_id].push_back(total_sample_count+sample);
				total_sample_count+=unique_samples.size();
				graph.curves.push_back(std::move(curve));
				pending_boundary_edges.push_back(std::move(curve_boundary_edges));
			}
			ContactDisjointSet topology_groups(total_sample_count);
			for(const std::vector<std::size_t>& locations:topology_locations)
				for(std::size_t location=1;location<locations.size();++location)
					topology_groups.merge(locations.front(),locations[location]);
			std::map<std::size_t,std::uint32_t> topology_id_by_root;
			std::uint32_t next_topology_id=0;std::size_t global_sample=0;
			for(StepContactCurve& curve:graph.curves)
			{
				curve.sample_topology_ids.reserve(curve.sample_positions_m.size());
				for(std::size_t sample=0;sample<curve.sample_positions_m.size();++sample)
				{
					const std::size_t root=topology_groups.find(global_sample++);
					auto [found,inserted]=topology_id_by_root.emplace(root,next_topology_id);
					if(inserted)++next_topology_id;
					curve.sample_topology_ids.push_back(found->second);
				}
			}
			// A topology ID is a discrete identity, not a proximity hint. Curves meeting at
			// an exact OCCT vertex/T-junction were projected independently onto their own
			// support curves above, so their floating-point coordinates can legitimately
			// differ by a native CAD tolerance. Collapse every joined node to the first
			// deterministic graph coordinate now. The STEP conformer will copy this exact
			// coordinate into every face-local render vertex carrying the ID; BVH/EB can
			// therefore require bit-identical coordinates without performing a tolerance weld.
			std::map<std::uint32_t,std::array<double,3>> canonical_position_by_topology;
			std::map<std::uint32_t,double> tolerance_by_topology;
			for(StepContactCurve& curve:graph.curves)
				for(std::size_t sample=0;sample<curve.sample_positions_m.size();++sample)
				{
					const std::uint32_t topology=curve.sample_topology_ids[sample];
					const auto [position,inserted]=canonical_position_by_topology.emplace(
						topology,curve.sample_positions_m[sample]);
					if(!inserted)curve.sample_positions_m[sample]=position->second;
					auto [tolerance,tolerance_inserted]=tolerance_by_topology.emplace(
						topology,curve.tolerance_m);
					if(!tolerance_inserted)tolerance->second=std::max(
						tolerance->second,curve.tolerance_m);
				}
			for(StepContactCurve& curve:graph.curves)
				for(const std::uint32_t topology:curve.sample_topology_ids)
					curve.tolerance_m=std::max(curve.tolerance_m,
						tolerance_by_topology[topology]);
			for(std::size_t curve_id=0;curve_id<graph.curves.size();++curve_id)
			{
				StepContactCurve& curve=graph.curves[curve_id];
				if(curve_id>=pending_boundary_edges.size()
					||pending_boundary_edges[curve_id].size()!=curve.uses.size())
				{
					error="internal contact UV staging mismatch";
					return false;
				}
				std::vector<gp_Pnt> final_world_points;
				final_world_points.reserve(curve.sample_positions_m.size());
				for(const std::array<double,3>& position:curve.sample_positions_m)
					final_world_points.emplace_back(position[0]/kMmToM,
						position[1]/kMmToM,position[2]/kMmToM);
				for(std::size_t use_id=0;use_id<curve.uses.size();++use_id)
				{
					StepContactFaceUse& use=curve.uses[use_id];
					if(use.source_face_id>=faces.size())
					{
						error="contact UV use references an invalid source face";
						return false;
					}
					const TopoDS_Face& face=faces[use.source_face_id];
					TopLoc_Location face_location;
					(void)BRep_Tool::Surface(face,face_location);
					const double location_scale=std::abs(
						face_location.Transformation().ScaleFactor());
					const double tolerance_mm=std::max({Precision::Confusion(),
						curve.tolerance_m/kMmToM,
						BRep_Tool::Tolerance(face)*location_scale});
					std::string projection_error;
					bool used_exact_boundary_pcurve=false;
					if(!detail::project_trim_valid_face_uv(final_world_points,face,
						pending_boundary_edges[curve_id][use_id],
						use.kind==StepContactUseKind::trim_boundary,tolerance_mm,
						use.sample_uv,projection_error,&used_exact_boundary_pcurve))
					{
						// Keep the STEP surface available for display and diagnostics, but
						// retain the exact failure as a hard conforming-mesh prerequisite.
						// Publishing a partial UV chain would be more dangerous than an
						// explicit unresolved record, so the failed use remains empty.
						use.sample_uv.clear();
						graph.unresolved_uv_projections.push_back({curve.id,
							use.source_face_id,use.kind,projection_error});
						continue;
					}
					if(use.kind==StepContactUseKind::trim_boundary
						&&!used_exact_boundary_pcurve)
						++graph.trim_surface_fallback_count;
				}
			}
			output=std::move(graph);
			return true;
		}

		FaceNodePairMap face_cad_edge_segments(const TopoDS_Face& face,
			const Handle(Poly_Triangulation)& triangulation, const TopLoc_Location& location,
			const ShapeIndexedMap& shape_edges, const ShapeAncestorMap& edge_faces,
			const std::vector<CadEdgeContactCertificate>& edge_contacts,
			bool& unresolved_unlocated_boundary)
		{
			unresolved_unlocated_boundary = false;
			FaceNodePairMap result;
			for (TopExp_Explorer edge_exp(face, TopAbs_EDGE); edge_exp.More(); edge_exp.Next())
			{
				const TopoDS_Edge edge = TopoDS::Edge(edge_exp.Current());
				const Standard_Integer edge_index = shape_edges.FindIndex(edge);
				if (edge_index <= 0) continue;

				std::uint32_t incident_faces = edge_faces.Contains(edge)
					? static_cast<std::uint32_t>(edge_faces.FindFromKey(edge).Size()) : 0u;
				const std::uint32_t edge_id = static_cast<std::uint32_t>(edge_index - 1);
				const CadEdgeContactCertificate contact=edge_id<edge_contacts.size()
					?edge_contacts[edge_id]:CadEdgeContactCertificate{};
				const double location_scale = std::abs(location.Transformation().ScaleFactor());
				const double edge_tolerance_m = std::max(0.0,
					BRep_Tool::Tolerance(edge) * location_scale * kMmToM);
				const bool periodic_seam = BRep_Tool::IsClosed(edge, face);
				if (periodic_seam) incident_faces = std::max(incident_faces, 2u);
				auto polygon_for = [&](const TopoDS_Edge& oriented_edge)
				{
					Handle(Poly_PolygonOnTriangulation) polygon =
						BRep_Tool::PolygonOnTriangulation(oriented_edge, triangulation, location);
					return polygon;
				};
				auto append_polygon = [&](const Handle(Poly_PolygonOnTriangulation)& polygon)
				{
					for (Standard_Integer node = 1; node < polygon->NbNodes(); ++node)
						insert_face_node_pair(result, polygon->Node(node), polygon->Node(node + 1),
							edge_id, incident_faces, edge_tolerance_m, periodic_seam,contact);
				};
				const Handle(Poly_PolygonOnTriangulation) first = polygon_for(edge);
				if (!periodic_seam)
				{
					if (!first.IsNull() && first->NbNodes() >= 2) append_polygon(first);
					continue;
				}

				// A seam on a periodic face owns two distinct node chains on the same
				// triangulation. OCCT selects PolygonOnTriangulation2() through the edge-use
				// orientation. Query the reversed use explicitly rather than depending on
				// whether TopExp_Explorer returns both wire occurrences. If either chain is
				// absent or aliases the other, explicitly mark every available chain as an
				// unresolved CAD boundary.  The unavailable sentinel alone would be unsafe:
				// it is also the affirmative representation of an internal tessellation edge.
				if (BRep_Tool::IsClosed(edge, triangulation, location))
				{
					TopoDS_Edge opposite_use = edge;
					opposite_use.Reverse();
					const Handle(Poly_PolygonOnTriangulation) second = polygon_for(opposite_use);
					if (!first.IsNull() && !second.IsNull() && first->NbNodes() >= 2
						&& second->NbNodes() >= 2 && !same_polygon_chain(first, second))
					{
						append_polygon(first);
						append_polygon(second);
					}
					else
					{
						auto append_unknown = [&](const Handle(Poly_PolygonOnTriangulation)& polygon)
						{
							if (polygon.IsNull()) return;
							for (Standard_Integer node = 1; node < polygon->NbNodes(); ++node)
								insert_unknown_face_node_pair(result, polygon->Node(node),
									polygon->Node(node + 1));
						};
						append_unknown(first);
						append_unknown(second);
						unresolved_unlocated_boundary = unresolved_unlocated_boundary
							|| (first.IsNull() && second.IsNull());
					}
				}
				else
				{
					unresolved_unlocated_boundary = true;
					if (!first.IsNull())
						for (Standard_Integer node = 1; node < first->NbNodes(); ++node)
							insert_unknown_face_node_pair(result, first->Node(node),
								first->Node(node + 1));
				}
			}
			return result;
		}

		struct FaceDisjointSet
		{
			explicit FaceDisjointSet(std::size_t count) : parent(count), rank(count, 0)
			{
				std::iota(parent.begin(), parent.end(), std::size_t{0});
			}

			std::size_t find(std::size_t value)
			{
				std::size_t root = value;
				while (parent[root] != root) root = parent[root];
				while (parent[value] != value)
				{
					const std::size_t next = parent[value];
					parent[value] = root;
					value = next;
				}
				return root;
			}

			void join(std::size_t a, std::size_t b)
			{
				a = find(a); b = find(b); if (a == b) return;
				if (rank[a] < rank[b]) std::swap(a, b);
				parent[b] = a;
				if (rank[a] == rank[b]) ++rank[a];
			}

			std::vector<std::size_t> parent;
			std::vector<std::uint8_t> rank;
		};

		int source_face_index(const std::vector<TopoDS_Face>& source_faces,
			const TopoDS_Shape& candidate)
		{
			for (std::size_t index = 0; index < source_faces.size(); ++index)
				if (source_faces[index].IsSame(candidate)) return static_cast<int>(index);
			return -1;
		}

		std::string trimmed_representation_name(const Handle(TCollection_HAsciiString)& value)
		{
			if (value.IsNull() || value->IsEmpty()) return {};
			std::string name = value->ToCString();
			const auto is_space = [](unsigned char byte) { return std::isspace(byte) != 0; };
			const auto first = std::find_if_not(name.begin(), name.end(), is_space);
			const auto last = std::find_if_not(name.rbegin(), name.rend(), is_space).base();
			return first < last ? std::string(first, last) : std::string{};
		}

		std::vector<StepSourceFaceMetadata> source_face_representation_metadata(
			const STEPControl_Reader& reader, const std::vector<TopoDS_Face>& source_faces)
		{
			std::vector<StepSourceFaceMetadata> result(source_faces.size());
			for (std::size_t face = 0; face < result.size(); ++face)
				result[face].source_face_id = static_cast<std::uint32_t>(face);

			const Handle(StepData_StepModel) model = reader.StepModel();
			const Handle(XSControl_WorkSession) session = reader.WS();
			if (model.IsNull() || session.IsNull() || session->TransferReader().IsNull()) return result;
			const Handle(Transfer_TransientProcess) process =
				session->TransferReader()->TransientProcess();
			if (process.IsNull()) return result;

			for (Standard_Integer entity_index = 1; entity_index <= model->NbEntities(); ++entity_index)
			{
				const Handle(Standard_Transient) entity = model->Value(entity_index);
				const Handle(StepRepr_RepresentationItem) item =
					Handle(StepRepr_RepresentationItem)::DownCast(entity);
				if (item.IsNull()) continue;
				const std::string name = trimmed_representation_name(item->Name());
				if (name.empty()) continue;

				// Use the same transient transfer map that produced reader.OneShape().  This is an
				// exact source-entity -> BRep correspondence, not a parse of STEP text or a
				// geometry-proximity guess.
				const TopoDS_Shape transferred = TransferBRep::ShapeResult(process, entity);
				if (transferred.IsNull()) continue;

				std::vector<TopoDS_Face> transferred_faces;
				if (transferred.ShapeType() == TopAbs_FACE)
					transferred_faces.push_back(TopoDS::Face(transferred));
				else
					for (TopExp_Explorer face(transferred, TopAbs_FACE); face.More(); face.Next())
						transferred_faces.push_back(TopoDS::Face(face.Current()));

				for (const TopoDS_Face& transferred_face : transferred_faces)
					for (std::size_t source_face = 0; source_face < source_faces.size(); ++source_face)
						if (source_faces[source_face].IsPartner(transferred_face))
							result[source_face].representation_names.push_back(name);
			}

			for (StepSourceFaceMetadata& face : result)
			{
				std::sort(face.representation_names.begin(), face.representation_names.end());
				face.representation_names.erase(std::unique(face.representation_names.begin(),
					face.representation_names.end()), face.representation_names.end());
			}
			return result;
		}

		std::vector<std::uint32_t> edge_owner_face_ids(const ShapeAncestorMap& edge_faces,
			const TopoDS_Edge& edge, const std::vector<TopoDS_Face>& source_faces)
		{
			std::vector<std::uint32_t> result;
			if (!edge_faces.Contains(edge)) return result;
			const auto& owners = edge_faces.FindFromKey(edge);
			for (NCollection_List<TopoDS_Shape>::Iterator owner(owners); owner.More(); owner.Next())
			{
				const int face = source_face_index(source_faces, owner.Value());
				if (face >= 0) result.push_back(static_cast<std::uint32_t>(face));
			}
			std::sort(result.begin(), result.end());
			result.erase(std::unique(result.begin(), result.end()), result.end());
			return result;
		}

		std::vector<std::size_t> exact_face_components(const ShapeIndexedMap& shape_edges,
			const ShapeAncestorMap& edge_faces, const std::vector<TopoDS_Face>& source_faces,
			const std::vector<CadEdgeContactCertificate>& edge_contacts)
		{
			FaceDisjointSet components(source_faces.size());
			for (Standard_Integer edge_index = 1; edge_index <= shape_edges.Extent(); ++edge_index)
			{
				const TopoDS_Edge edge = TopoDS::Edge(shape_edges(edge_index));
				const std::vector<std::uint32_t> owners = edge_owner_face_ids(edge_faces, edge,
					source_faces);
				for (std::size_t owner = 1; owner < owners.size(); ++owner)
					components.join(owners.front(), owners[owner]);
				const std::size_t contact_index = static_cast<std::size_t>(edge_index - 1);
				if (owners.empty() || contact_index >= edge_contacts.size()) continue;
				for (std::uint32_t target : edge_contacts[contact_index].target_face_ids)
					if (target < source_faces.size()) components.join(owners.front(), target);
			}

			std::vector<std::size_t> result(source_faces.size());
			for (std::size_t face = 0; face < source_faces.size(); ++face)
				result[face] = components.find(face);
			return result;
		}

		double triangle_area(const TriMesh& mesh, std::uint32_t triangle)
		{
			const auto a = mesh.vertex_position_double(mesh.indices[3 * triangle]);
			const auto b = mesh.vertex_position_double(mesh.indices[3 * triangle + 1]);
			const auto c = mesh.vertex_position_double(mesh.indices[3 * triangle + 2]);
			const double abx = b[0] - a[0], aby = b[1] - a[1], abz = b[2] - a[2];
			const double acx = c[0] - a[0], acy = c[1] - a[1], acz = c[2] - a[2];
			const double nx = aby * acz - abz * acy;
			const double ny = abz * acx - abx * acz;
			const double nz = abx * acy - aby * acx;
			return 0.5 * std::sqrt(nx * nx + ny * ny + nz * nz);
		}

		GeometryIssuePolyline tessellated_edge_polyline(const TopoDS_Edge& edge,
			const TopoDS_Face& owner_face, std::uint32_t source_edge_id)
		{
			GeometryIssuePolyline result;
			result.source_edge_id = source_edge_id;
			TopoDS_Edge face_edge;
			for (TopExp_Explorer it(owner_face, TopAbs_EDGE); it.More(); it.Next())
				if (edge.IsSame(it.Current())) { face_edge = TopoDS::Edge(it.Current()); break; }
			if (face_edge.IsNull()) return result;

			TopLoc_Location location;
			const Handle(Poly_Triangulation) triangulation =
				BRep_Tool::Triangulation(owner_face, location);
			const Handle(Poly_PolygonOnTriangulation) polygon = triangulation.IsNull()
				? Handle(Poly_PolygonOnTriangulation){}
				: BRep_Tool::PolygonOnTriangulation(face_edge, triangulation, location);
			if (!polygon.IsNull() && polygon->NbNodes() >= 2)
			{
				const gp_Trsf transform = location.Transformation();
				result.points.reserve(static_cast<std::size_t>(polygon->NbNodes()));
				for (Standard_Integer node = 1; node <= polygon->NbNodes(); ++node)
				{
					gp_Pnt point = triangulation->Node(polygon->Node(node));
					point.Transform(transform);
					result.points.push_back({static_cast<float>(point.X() * kMmToM),
						static_cast<float>(point.Y() * kMmToM),
						static_cast<float>(point.Z() * kMmToM)});
				}
				return result;
			}

			// A successfully meshed ordinary boundary normally owns a polygon chain. Keep a
			// deterministic exact-curve fallback so the warning remains visible if an OCCT
			// triangulation omits that optional representation; this does not alter connectivity.
			BRepAdaptor_Curve curve(face_edge);
			const double first = curve.FirstParameter(), last = curve.LastParameter();
			if (!std::isfinite(first) || !std::isfinite(last) || !(last > first)) return result;
			constexpr int fallback_samples = 33;
			result.points.reserve(fallback_samples);
			for (int sample = 0; sample < fallback_samples; ++sample)
			{
				const double fraction = static_cast<double>(sample) / (fallback_samples - 1);
				const gp_Pnt point = curve.Value(first + fraction * (last - first));
				result.points.push_back({static_cast<float>(point.X() * kMmToM),
					static_cast<float>(point.Y() * kMmToM),
					static_cast<float>(point.Z() * kMmToM)});
			}
			return result;
		}

		GeometryQualityReport disconnected_component_report(const ShapeIndexedMap& shape_edges,
			const ShapeAncestorMap& edge_faces, const std::vector<TopoDS_Face>& source_faces,
			const std::vector<CadEdgeContactCertificate>& edge_contacts, const TriMesh& mesh)
		{
			GeometryQualityReport report;
			if (source_faces.empty() || mesh.empty() || !mesh.has_face_provenance()) return report;
			const std::vector<std::size_t> face_component = exact_face_components(shape_edges,
				edge_faces, source_faces, edge_contacts);

			struct Component
			{
				std::size_t root = 0;
				double area_m2 = 0.0;
				std::vector<std::uint32_t> faces;
				std::vector<std::uint32_t> triangles;
			};
			std::map<std::size_t, Component> by_root;
			for (std::uint32_t triangle = 0; triangle < mesh.triangle_count(); ++triangle)
			{
				const std::uint32_t face = mesh.source_face_ids[triangle];
				if (face >= face_component.size()) continue;
				Component& component = by_root[face_component[face]];
				component.root = face_component[face];
				component.area_m2 += triangle_area(mesh, triangle);
				component.faces.push_back(face);
				component.triangles.push_back(triangle);
			}
			std::vector<Component> components;
			components.reserve(by_root.size());
			for (auto& [root, component] : by_root)
			{
				(void)root;
				std::sort(component.faces.begin(), component.faces.end());
				component.faces.erase(std::unique(component.faces.begin(), component.faces.end()),
					component.faces.end());
				components.push_back(std::move(component));
			}
			std::sort(components.begin(), components.end(), [](const Component& a, const Component& b)
			{
				if (a.area_m2 != b.area_m2) return a.area_m2 > b.area_m2;
				if (a.triangles.size() != b.triangles.size()) return a.triangles.size() > b.triangles.size();
				return a.faces < b.faces;
			});
			report.exact_connected_component_count = components.size();
			if (components.empty()) return report;
			report.largest_component_area_m2 = components.front().area_m2;

			std::map<std::size_t, std::size_t> issue_by_root;
			for (std::size_t component_index = 1; component_index < components.size(); ++component_index)
			{
				const Component& component = components[component_index];
				GeometryIssue issue;
				issue.stage = GeometryIssueStage::cad_input;
				issue.kind = GeometryIssueKind::disconnected_fabric_component;
				issue.severity = GeometryIssueSeverity::warning_run;
				issue.state = GeometryIssueState::detected;
				issue.surface_area_m2 = component.area_m2;
				issue.source_face_ids = component.faces;
				issue.source_triangle_ids = component.triangles;
				std::uint64_t hash = 1469598103934665603ull;
				for (std::uint32_t face : component.faces)
				{
					hash ^= static_cast<std::uint64_t>(face) + 1;
					hash *= 1099511628211ull;
				}
				issue.id = hash;
				std::ostringstream summary;
				summary << "Disconnected fabric component retained: " << component.faces.size()
					<< " source face(s), " << component.triangles.size() << " triangle(s), "
					<< component.area_m2 << " m^2. It has no exact CAD edge/contact path to the "
					<< "largest fabric component.";
				issue.summary = summary.str();
				issue_by_root.emplace(component.root, report.issues.size());
				report.issues.push_back(std::move(issue));
			}

			// Outline only the true BRep boundary of each disconnected component. Shared edges,
			// periodic seams, and exact full edge-on-face contacts are internal connectivity and
			// are deliberately not painted as gaps. Openings within the main component do not
			// create an issue at all.
			for (Standard_Integer edge_index = 1; edge_index <= shape_edges.Extent(); ++edge_index)
			{
				const TopoDS_Edge edge = TopoDS::Edge(shape_edges(edge_index));
				if (!edge_faces.Contains(edge) || edge_faces.FindFromKey(edge).Size() != 1) continue;
				const std::size_t contact_index = static_cast<std::size_t>(edge_index - 1);
				if (contact_index < edge_contacts.size()
					&& !edge_contacts[contact_index].target_face_ids.empty()) continue;
				const std::vector<std::uint32_t> owners = edge_owner_face_ids(edge_faces, edge,
					source_faces);
				if (owners.size() != 1 || owners.front() >= face_component.size()) continue;
				const auto found = issue_by_root.find(face_component[owners.front()]);
				if (found == issue_by_root.end()) continue;
				GeometryIssuePolyline polyline = tessellated_edge_polyline(edge,
					source_faces[owners.front()], static_cast<std::uint32_t>(edge_index - 1));
				if (polyline.points.size() < 2) continue;
				GeometryIssue& issue = report.issues[found->second];
				for (std::size_t point = 1; point < polyline.points.size(); ++point)
				{
					const auto& a = polyline.points[point - 1];
					const auto& b = polyline.points[point];
					const double dx = static_cast<double>(b[0]) - a[0];
					const double dy = static_cast<double>(b[1]) - a[1];
					const double dz = static_cast<double>(b[2]) - a[2];
					issue.boundary_length_m += std::sqrt(dx * dx + dy * dy + dz * dz);
				}
				issue.boundary_polylines.push_back(std::move(polyline));
			}
			return report;
		}

		void set_error(std::string* error, const std::string& msg)
		{
			if (error) *error = msg;
		}

		// Transfer a loaded reader's roots to a single shape and triangulate it in place (BRepMesh).
		// Shared tail of the file- and stream-based readers below. Returns false with *error set on
		// transfer/null-shape failure. deflection is in OCC's native millimetres.
		bool transfer_and_mesh(STEPControl_Reader& reader, double deflection_mm, TopoDS_Shape& out, std::string* error)
		{
			Standard_Integer nroots = reader.TransferRoots();
			if (nroots <= 0)
			{
				set_error(error, "STEPControl_Reader::TransferRoots transferred no roots");
				return false;
			}

			out = reader.OneShape();
			if (out.IsNull())
			{
				set_error(error, "reader.OneShape() returned a null shape");
				return false;
			}

			// Triangulate (deflection is in OCC's native millimetres).
			if (deflection_mm <= 0.0) deflection_mm = 0.1;
			BRepMesh_IncrementalMesh mesher(out, deflection_mm);
			mesher.Perform();
			return true;
		}

		// Read + transfer a STEP file to its single shape and triangulate it in place (BRepMesh).
		// Returns false with *error set on any read/transfer/null-shape failure. deflection is in
		// OCC's native millimetres. May throw Standard_Failure (OCC) — callers wrap in try/catch.
		bool read_step_shape(const std::string& path, double deflection_mm,
			STEPControl_Reader& reader, TopoDS_Shape& out, std::string* error)
		{
			IFSelect_ReturnStatus status = reader.ReadFile(path.c_str());
			if (status != IFSelect_RetDone)
			{
				set_error(error, "STEPControl_Reader::ReadFile failed for '" + path + "'");
				return false;
			}
			return transfer_and_mesh(reader, deflection_mm, out, error);
		}

		// Accumulate every triangulated FACE of `shape` into one TriMesh (metres, face-local winding,
		// area-weighted per-vertex normals, bbox). `shape` must ALREADY be meshed (BRepMesh run on
		// it or an ancestor). Returns empty() if the shape contributes no triangles.
		TriMesh mesh_from_faces(const TopoDS_Shape& shape,
			const std::vector<TopoDS_Face>& source_faces,
			GeometryQualityReport* quality, StepCadContactGraph* contact_graph,
			std::string* error)
		{
			TriMesh mesh;
			// IDs are topology traversal IDs, deliberately independent of tessellation density.
			// TopTools_ShapeMapHasher ignores orientation but retains TShape + location identity,
			// so the same seam edge encountered with opposite face-use orientation keeps one ID.
			ShapeIndexedMap shape_edges;
			ShapeAncestorMap edge_faces;
			TopExp::MapShapes(shape, TopAbs_EDGE, shape_edges);
			// Deliberately retain face-use multiplicity: a periodic seam is the same
			// TopoDS_Edge used twice by one face and therefore bounds two sheet sectors.
			TopExp::MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_faces);
			std::vector<CadEdgeContactCertificate> edge_contacts=
				certify_edge_face_contacts(shape_edges,edge_faces,source_faces);
			if(contact_graph)
			{
				StepCadContactGraph staged_graph;
				std::string graph_error;
				if(!make_contact_graph(shape_edges,edge_faces,source_faces,edge_contacts,
					staged_graph,graph_error))
				{
					set_error(error,graph_error);
					return TriMesh{};
				}
				*contact_graph=std::move(staged_graph);
			}

			// --- Accumulate faces into one flat mesh (metres) ----------------------
			for (std::size_t source_face_index=0;source_face_index<source_faces.size();++source_face_index)
			{
				const std::uint32_t source_face_id=static_cast<std::uint32_t>(source_face_index);
				TopoDS_Face face = source_faces[source_face_index];
				TopLoc_Location loc;
				Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
				if (tri.IsNull() || tri->NbNodes() <= 0 || tri->NbTriangles() <= 0) continue;

				const bool reversed = (face.Orientation() == TopAbs_REVERSED);
				const std::uint32_t base = static_cast<std::uint32_t>(mesh.positions.size() / 3);
				const gp_Trsf trsf = loc.Transformation();
				bool unresolved_unlocated_boundary = false;
				const FaceNodePairMap cad_segments = face_cad_edge_segments(face, tri, loc,
					shape_edges, edge_faces, edge_contacts, unresolved_unlocated_boundary);
				if (unresolved_unlocated_boundary)
				{
					set_error(error, "periodic CAD boundary on source face "
						+ std::to_string(source_face_id)
						+ " has no distinct triangulation chains; refusing to infer an opening");
					return TriMesh{};
				}
				std::unordered_map<std::uint64_t, std::uint32_t> triangulation_segment_counts;
				for (Standard_Integer triangle = 1; triangle <= tri->NbTriangles(); ++triangle)
				{
					Standard_Integer a, b, c;
					tri->Triangle(triangle).Get(a, b, c);
					const Standard_Integer node[3] = { a, b, c };
					for (int edge = 0; edge < 3; ++edge)
						++triangulation_segment_counts[undirected_node_pair(
							node[edge], node[(edge + 1) % 3])];
				}

				for (Standard_Integer i = 1; i <= tri->NbNodes(); ++i)
				{
					gp_Pnt p = tri->Node(i);
					p.Transform(trsf);
					const double x = p.X() * kMmToM;
					const double y = p.Y() * kMmToM;
					const double z = p.Z() * kMmToM;
					mesh.positions_fp64.insert(mesh.positions_fp64.end(), { x, y, z });
					mesh.positions.push_back(static_cast<float>(x));
					mesh.positions.push_back(static_cast<float>(y));
					mesh.positions.push_back(static_cast<float>(z));
					if (tri->HasUVNodes())
					{
						const gp_Pnt2d uv = tri->UVNode(i);
						mesh.vertex_uv.push_back(static_cast<float>(uv.X()));
						mesh.vertex_uv.push_back(static_cast<float>(uv.Y()));
					}
					else
					{
						const float nan = std::numeric_limits<float>::quiet_NaN();
						mesh.vertex_uv.push_back(nan);
						mesh.vertex_uv.push_back(nan);
					}
				}

				for (Standard_Integer i = 1; i <= tri->NbTriangles(); ++i)
				{
					Standard_Integer n1, n2, n3;
					tri->Triangle(i).Get(n1, n2, n3);
					// Keep the face-local 1-based node numbers until after orientation is applied:
					// the per-half-edge metadata must be parallel to the FINAL index ordering.
					if (reversed) std::swap(n1, n2); // outward winding for a solid's face
					const Standard_Integer local_nodes[3] = { n1, n2, n3 };
					for (int corner = 0; corner < 3; ++corner)
						mesh.indices.push_back(base + static_cast<std::uint32_t>(local_nodes[corner] - 1));
					for (int half_edge = 0; half_edge < 3; ++half_edge)
					{
						const std::uint64_t segment = undirected_node_pair(local_nodes[half_edge],
							local_nodes[(half_edge + 1) % 3]);
						const auto found = cad_segments.find(segment);
						const auto count = triangulation_segment_counts.find(segment);
						const std::uint32_t segment_incidence = count == triangulation_segment_counts.end()
							? 0u : count->second;
						const bool face_boundary = segment_incidence == 1;
						const bool regular_interior = segment_incidence == 2;
						const bool known = face_boundary && found != cad_segments.end()
							&& !found->second.ambiguous
							&& found->second.id != TriMesh::kNoCadEdgeId;
						const bool unknown_boundary = (face_boundary && !known)
							|| (!face_boundary && found != cad_segments.end())
							|| (!face_boundary && !regular_interior);
						const CadEdgeProvenance provenance = known ? found->second : CadEdgeProvenance{};
						mesh.triangle_cad_edge_provenance_states.push_back(static_cast<std::uint8_t>(
							known ? CadEdgeProvenanceState::known
							: (unknown_boundary ? CadEdgeProvenanceState::unknown_boundary
								: CadEdgeProvenanceState::none)));
						mesh.triangle_cad_edge_ids.push_back(known ? provenance.id : TriMesh::kNoCadEdgeId);
						mesh.triangle_cad_edge_incident_face_counts.push_back(
							provenance.id == TriMesh::kNoCadEdgeId ? 0u : provenance.incident_face_count);
						mesh.triangle_cad_edge_tolerances.push_back(
							provenance.id == TriMesh::kNoCadEdgeId ? 0.0 : provenance.tolerance_m);
						mesh.triangle_cad_edge_is_periodic_seam.push_back(
							provenance.id == TriMesh::kNoCadEdgeId || !provenance.periodic_seam ? 0u : 1u);
						mesh.triangle_cad_edge_contact_ids.push_back(
							provenance.id == TriMesh::kNoCadEdgeId ? TriMesh::kNoCadContactId
								: provenance.contact_id);
						mesh.triangle_cad_edge_certified_fan_degrees.push_back(
							provenance.id == TriMesh::kNoCadEdgeId ? 0u
								: provenance.certified_fan_degree);
					}
					mesh.source_face_ids.push_back(source_face_id);
				}
			}

			if (mesh.indices.empty()) return mesh; // empty() — caller reports the error/skips

			// --- Per-vertex normals: area-weighted average of incident face normals ----
			const std::size_t nverts = mesh.vertex_count();
			std::vector<double> acc(nverts * 3, 0.0);
			for (std::size_t t = 0; t + 2 < mesh.indices.size(); t += 3)
			{
				const std::uint32_t i0 = mesh.indices[t + 0];
				const std::uint32_t i1 = mesh.indices[t + 1];
				const std::uint32_t i2 = mesh.indices[t + 2];
				const double* p0 = &mesh.positions_fp64[3 * i0];
				const double* p1 = &mesh.positions_fp64[3 * i1];
				const double* p2 = &mesh.positions_fp64[3 * i2];
				const double e1x = p1[0] - p0[0], e1y = p1[1] - p0[1], e1z = p1[2] - p0[2];
				const double e2x = p2[0] - p0[0], e2y = p2[1] - p0[1], e2z = p2[2] - p0[2];
				// Cross(e1,e2): magnitude == 2*area, so accumulating raw crosses area-weights.
				const double nx = e1y * e2z - e1z * e2y;
				const double ny = e1z * e2x - e1x * e2z;
				const double nz = e1x * e2y - e1y * e2x;
				for (std::uint32_t idx : { i0, i1, i2 })
				{
					acc[3 * idx + 0] += nx;
					acc[3 * idx + 1] += ny;
					acc[3 * idx + 2] += nz;
				}
			}

			mesh.normals.resize(nverts * 3);
			for (std::size_t v = 0; v < nverts; ++v)
			{
				double nx = acc[3 * v + 0], ny = acc[3 * v + 1], nz = acc[3 * v + 2];
				const double len = std::sqrt(nx * nx + ny * ny + nz * nz);
				if (len > 1e-20) { nx /= len; ny /= len; nz /= len; }
				else { nx = 0.0; ny = 0.0; nz = 1.0; } // degenerate fallback: z-up
				mesh.normals[3 * v + 0] = static_cast<float>(nx);
				mesh.normals[3 * v + 1] = static_cast<float>(ny);
				mesh.normals[3 * v + 2] = static_cast<float>(nz);
			}

			// --- Bounding box (metres) --------------------------------------------------
			double lo[3] = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
			double hi[3] = { -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max() };
			for (std::size_t v = 0; v < nverts; ++v)
				for (int c = 0; c < 3; ++c)
				{
					const double x = mesh.positions_fp64[3 * v + c];
					lo[c] = std::min(lo[c], x);
					hi[c] = std::max(hi[c], x);
				}
			mesh.bbox_min = { { static_cast<float>(lo[0]), static_cast<float>(lo[1]), static_cast<float>(lo[2]) } };
			mesh.bbox_max = { { static_cast<float>(hi[0]), static_cast<float>(hi[1]), static_cast<float>(hi[2]) } };
			if (quality)
			{
				*quality = disconnected_component_report(shape_edges, edge_faces, source_faces,
					edge_contacts, mesh);
				// Representation names describe the producer's intent; they are not a geometric
				// validity test. In particular, a correctly attached mini-rib is real fabric.
				// Disconnected artifacts are already reported above from exact CAD connectivity.
			}

			return mesh;
		}
	}

	StepGeometry load_step_geometry(const std::string& path, double deflection_mm,
		std::string* error)
	{
		try
		{
			STEPControl_Reader reader;
			TopoDS_Shape shape;
			if (!read_step_shape(path, deflection_mm, reader, shape, error)) return {};

			std::vector<TopoDS_Face> source_faces;
			for (TopExp_Explorer face(shape, TopAbs_FACE); face.More(); face.Next())
				source_faces.push_back(TopoDS::Face(face.Current()));

			StepGeometry result;
			result.source_faces = source_face_representation_metadata(reader, source_faces);
			result.mesh = mesh_from_faces(shape, source_faces, &result.quality,
				&result.contacts, error);
			if (result.mesh.empty())
			{
				if (!error || error->empty()) set_error(error, "shape produced no triangulable faces");
				return {};
			}
			if (error) error->clear();
			return result;
		}
		catch (const Standard_Failure& f)
		{
			set_error(error, std::string("OpenCascade exception: ") + f.GetMessageString());
			return {};
		}
	}

	TriMesh load_step_mesh(const std::string& path, double deflection_mm, std::string* error)
	{
		return load_step_geometry(path, deflection_mm, error).mesh;
	}

}
