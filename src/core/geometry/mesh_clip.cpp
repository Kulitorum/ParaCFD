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
		};

		Vertex interpolate(const Vertex& a,const Vertex& b,double t)
		{
			Vertex out;for(int q=0;q<3;++q)out.p[q]=a.p[q]+t*(b.p[q]-a.p[q]);for(int q=0;q<2;++q)out.uv[q]=a.uv[q]+t*(b.uv[q]-a.uv[q]);return out;
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
	}

	TriMesh clip_mesh_to_axis_slab(const TriMesh& mesh,int axis,double lower,double upper)
	{
		TriMesh out;if(mesh.empty()||axis<0||axis>2||!(upper>lower))return out;const bool have_uv=mesh.has_uv();const bool have_faces=mesh.has_face_provenance();const double nan=std::numeric_limits<double>::quiet_NaN();std::array<double,3> bbox_lo{{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()}},bbox_hi{{-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()}};
		for(std::size_t triangle=0;triangle<mesh.triangle_count();++triangle)
		{
			std::vector<Vertex> polygon;polygon.reserve(5);for(int corner=0;corner<3;++corner){const std::uint32_t index=mesh.indices[3*triangle+corner];if(3*static_cast<std::size_t>(index)+2>=mesh.positions.size()){polygon.clear();break;}Vertex v;for(int q=0;q<3;++q)v.p[q]=mesh.positions[3*index+q];v.uv[0]=have_uv?mesh.vertex_uv[2*index]:nan;v.uv[1]=have_uv?mesh.vertex_uv[2*index+1]:nan;polygon.push_back(v);}polygon=clip_plane(polygon,axis,lower,true);polygon=clip_plane(polygon,axis,upper,false);if(polygon.size()<3)continue;
			for(std::size_t fan=1;fan+1<polygon.size();++fan)
			{
				const Vertex vertices[3]={polygon[0],polygon[fan],polygon[fan+1]};const double ab[3]={vertices[1].p[0]-vertices[0].p[0],vertices[1].p[1]-vertices[0].p[1],vertices[1].p[2]-vertices[0].p[2]},ac[3]={vertices[2].p[0]-vertices[0].p[0],vertices[2].p[1]-vertices[0].p[1],vertices[2].p[2]-vertices[0].p[2]};double normal[3]={ab[1]*ac[2]-ab[2]*ac[1],ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]};const double length=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);if(!(length>1e-18))continue;for(double& value:normal)value/=length;const std::uint32_t base=static_cast<std::uint32_t>(out.vertex_count());for(const Vertex& vertex:vertices){for(int q=0;q<3;++q){out.positions.push_back(static_cast<float>(vertex.p[q]));out.normals.push_back(static_cast<float>(normal[q]));bbox_lo[q]=std::min(bbox_lo[q],vertex.p[q]);bbox_hi[q]=std::max(bbox_hi[q],vertex.p[q]);}out.vertex_uv.push_back(static_cast<float>(vertex.uv[0]));out.vertex_uv.push_back(static_cast<float>(vertex.uv[1]));}out.indices.insert(out.indices.end(),{base,base+1,base+2});out.source_face_ids.push_back(have_faces?mesh.source_face_ids[triangle]:0u);
			}
		}
		if(out.empty())return out;bbox_lo[axis]=lower;bbox_hi[axis]=upper;for(int q=0;q<3;++q){out.bbox_min[q]=static_cast<float>(bbox_lo[q]);out.bbox_max[q]=static_cast<float>(bbox_hi[q]);}return out;
	}
}
