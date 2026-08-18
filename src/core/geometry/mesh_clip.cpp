#include "core/geometry/mesh_clip.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace paracfd::core
{
	namespace
	{
		struct Vertex
		{
			double p[3]{};
			double uv[2]{};
			double bary[3]{};
		};

		Vertex interpolate(const Vertex& a,const Vertex& b,double t)
		{
			Vertex out;for(int q=0;q<3;++q){out.p[q]=a.p[q]+t*(b.p[q]-a.p[q]);out.bary[q]=a.bary[q]+t*(b.bary[q]-a.bary[q]);}for(int q=0;q<2;++q)out.uv[q]=a.uv[q]+t*(b.uv[q]-a.uv[q]);return out;
		}

		std::vector<Vertex> clip_plane(const std::vector<Vertex>& input,int axis,double bound,bool keep_greater)
		{
			std::vector<Vertex> output;if(input.empty())return output;auto inside=[&](const Vertex& v){return keep_greater?v.p[axis]>=bound:v.p[axis]<=bound;};Vertex previous=input.back();bool previous_inside=inside(previous);
			for(const Vertex& current:input)
			{
				const bool current_inside=inside(current);if(current_inside!=previous_inside){const double denominator=current.p[axis]-previous.p[axis];if(std::abs(denominator)>1e-30){const double t=std::clamp((bound-previous.p[axis])/denominator,0.0,1.0);output.push_back(interpolate(previous,current,t));}}
				if(current_inside)output.push_back(current);previous=current;previous_inside=current_inside;
			}
			return output;
		}

		struct InheritedCadEdge
		{
			CadEdgeProvenanceState state=CadEdgeProvenanceState::none;
			std::uint32_t id=TriMesh::kNoCadEdgeId;
			std::uint32_t incident_faces=0;
			double tolerance_m=0.0;
			bool periodic_seam=false;
			std::uint64_t contact_id=TriMesh::kNoCadContactId;
			std::uint32_t certified_fan_degree=0;
		};

		InheritedCadEdge inherited_cad_edge(const TriMesh& mesh,std::size_t triangle,
			const Vertex& a,const Vertex& b)
		{
			if(!mesh.has_cad_edge_provenance())return {};
			constexpr double bary_tolerance=64.0*std::numeric_limits<double>::epsilon();
			int opposite=-1;
			for(int corner=0;corner<3;++corner)
			{
				if(std::abs(a.bary[corner])<=bary_tolerance&&std::abs(b.bary[corner])<=bary_tolerance)
				{
					if(opposite>=0)return {};
					opposite=corner;
				}
			}
			if(opposite<0)return {};
			const unsigned source_half_edge=static_cast<unsigned>((opposite+1)%3);
			const CadEdgeProvenanceState state=
				mesh.cad_edge_provenance_state(triangle,source_half_edge);
			const std::uint32_t id=mesh.cad_edge_id(triangle,source_half_edge);
			return {state,id,id==TriMesh::kNoCadEdgeId?0u:
				mesh.cad_edge_incident_face_count(triangle,source_half_edge),
				id==TriMesh::kNoCadEdgeId?0.0:mesh.cad_edge_tolerance(triangle,source_half_edge),
				id!=TriMesh::kNoCadEdgeId&&mesh.cad_edge_is_periodic_seam(triangle,source_half_edge),
				id!=TriMesh::kNoCadEdgeId?mesh.cad_edge_contact_id(triangle,source_half_edge)
					:TriMesh::kNoCadContactId,
				id!=TriMesh::kNoCadEdgeId?mesh.cad_edge_certified_fan_degree(triangle,source_half_edge):0u};
		}
	}

	TriMesh clip_mesh_to_axis_slab(const TriMesh& mesh,int axis,double lower,double upper,
		bool discard_lower_coplanar,bool discard_upper_coplanar)
	{
		TriMesh out;if(mesh.empty()||axis<0||axis>2||!(upper>lower))return out;const bool have_fp64=mesh.has_fp64_positions();const bool have_uv=mesh.has_uv();const bool have_faces=mesh.has_face_provenance();const bool have_edges=mesh.has_cad_edge_provenance();const bool have_contacts=mesh.has_cad_edge_contact_provenance();const double nan=std::numeric_limits<double>::quiet_NaN();const double boundary_tolerance=1e-7*std::max({1.0,std::abs(lower),std::abs(upper),std::abs(upper-lower)});std::array<double,3> bbox_lo{{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()}},bbox_hi{{-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()}};
		for(std::size_t triangle=0;triangle<mesh.triangle_count();++triangle)
		{
			std::vector<Vertex> polygon;polygon.reserve(5);for(int corner=0;corner<3;++corner){const std::uint32_t index=mesh.indices[3*triangle+corner];if(index>=mesh.vertex_count()){polygon.clear();break;}Vertex v;const std::array<double,3> position=mesh.vertex_position_double(index);for(int q=0;q<3;++q)v.p[q]=position[q];v.bary[corner]=1.0;v.uv[0]=have_uv?mesh.vertex_uv[2*index]:nan;v.uv[1]=have_uv?mesh.vertex_uv[2*index+1]:nan;polygon.push_back(v);}if(polygon.size()<3)continue;const bool on_lower=std::all_of(polygon.begin(),polygon.end(),[&](const Vertex& vertex){return std::abs(vertex.p[axis]-lower)<=boundary_tolerance;});const bool on_upper=std::all_of(polygon.begin(),polygon.end(),[&](const Vertex& vertex){return std::abs(vertex.p[axis]-upper)<=boundary_tolerance;});if((discard_lower_coplanar&&on_lower)||(discard_upper_coplanar&&on_upper))continue;polygon=clip_plane(polygon,axis,lower,true);polygon=clip_plane(polygon,axis,upper,false);if(polygon.size()<3)continue;
			for(std::size_t fan=1;fan+1<polygon.size();++fan)
			{
				const Vertex vertices[3]={polygon[0],polygon[fan],polygon[fan+1]};const double ab[3]={vertices[1].p[0]-vertices[0].p[0],vertices[1].p[1]-vertices[0].p[1],vertices[1].p[2]-vertices[0].p[2]},ac[3]={vertices[2].p[0]-vertices[0].p[0],vertices[2].p[1]-vertices[0].p[1],vertices[2].p[2]-vertices[0].p[2]};double normal[3]={ab[1]*ac[2]-ab[2]*ac[1],ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]};const double length=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);if(!(length>1e-18))continue;for(double& value:normal)value/=length;const std::uint32_t base=static_cast<std::uint32_t>(out.vertex_count());for(const Vertex& vertex:vertices){for(int q=0;q<3;++q){if(have_fp64)out.positions_fp64.push_back(vertex.p[q]);out.positions.push_back(static_cast<float>(vertex.p[q]));out.normals.push_back(static_cast<float>(normal[q]));bbox_lo[q]=std::min(bbox_lo[q],vertex.p[q]);bbox_hi[q]=std::max(bbox_hi[q],vertex.p[q]);}out.vertex_uv.push_back(static_cast<float>(vertex.uv[0]));out.vertex_uv.push_back(static_cast<float>(vertex.uv[1]));}out.indices.insert(out.indices.end(),{base,base+1,base+2});if(have_faces)out.source_face_ids.push_back(mesh.source_face_ids[triangle]);if(have_edges){for(int half_edge=0;half_edge<3;++half_edge){const InheritedCadEdge edge=inherited_cad_edge(mesh,triangle,vertices[half_edge],vertices[(half_edge+1)%3]);out.triangle_cad_edge_provenance_states.push_back(static_cast<std::uint8_t>(edge.state));out.triangle_cad_edge_ids.push_back(edge.id);out.triangle_cad_edge_incident_face_counts.push_back(edge.incident_faces);out.triangle_cad_edge_tolerances.push_back(edge.tolerance_m);out.triangle_cad_edge_is_periodic_seam.push_back(edge.periodic_seam?1u:0u);if(have_contacts){out.triangle_cad_edge_contact_ids.push_back(edge.contact_id);out.triangle_cad_edge_certified_fan_degrees.push_back(edge.certified_fan_degree);}}}
			}
		}
		if(out.empty())return out;bbox_lo[axis]=lower;bbox_hi[axis]=upper;for(int q=0;q<3;++q){out.bbox_min[q]=static_cast<float>(bbox_lo[q]);out.bbox_max[q]=static_cast<float>(bbox_hi[q]);}return out;
	}
}
