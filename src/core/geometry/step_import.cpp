#include "core/geometry/step_import.h"
#include "core/geometry/occt_closed_solid.h"

#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Poly_Triangle.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_Failure.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Trsf.hxx>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <string>
#include <vector>

namespace paracfd::core
{
	namespace
	{
		constexpr double mm_to_m = 0.001;

		void set_error(std::string* output,const std::string& message)
		{
			if(output)*output=message;
		}

		struct ValidatedEnvelope
		{
			TopoDS_Solid solid;
			std::vector<TopoDS_Face> faces;
			bool reversed=false;
			double volume_m3=0.0;
		};

		bool read_step(const std::string& path,TopoDS_Shape& shape,std::string* error)
		{
			STEPControl_Reader reader;
			const IFSelect_ReturnStatus status=reader.ReadFile(path.c_str());
			if(status!=IFSelect_RetDone)
			{
				set_error(error,"cannot read STEP file '"+path+"'");
				return false;
			}
			if(reader.TransferRoots()<=0)
			{
				set_error(error,"STEP file transferred no BRep roots");
				return false;
			}
			shape=reader.OneShape();
			if(shape.IsNull())
			{
				set_error(error,"STEP transfer produced an empty shape");
				return false;
			}
			return true;
		}

		bool validate_envelope(const TopoDS_Shape& transferred,ValidatedEnvelope& output,
			std::string* error)
		{
			TopTools_IndexedMapOfShape solids;
			TopExp::MapShapes(transferred,TopAbs_SOLID,solids);
			if(solids.Extent()!=1)
			{
				set_error(error,"solid STEP must contain exactly one OCCT solid; found "+
					std::to_string(solids.Extent()));
				return false;
			}

			TopoDS_Solid solid=TopoDS::Solid(solids(1));
			TopTools_IndexedMapOfShape shells;
			TopExp::MapShapes(solid,TopAbs_SHELL,shells);
			if(shells.Extent()!=1)
			{
				set_error(error,"solid STEP must have one exterior shell and no cavities; found "+
					std::to_string(shells.Extent()));
				return false;
			}
			if(!BRep_Tool::IsClosed(TopoDS::Shell(shells(1))))
			{
				set_error(error,"the solid's exterior shell is open");
				return false;
			}

			TopTools_IndexedMapOfShape transferred_faces,solid_faces;
			TopExp::MapShapes(transferred,TopAbs_FACE,transferred_faces);
			TopExp::MapShapes(solid,TopAbs_FACE,solid_faces);
			bool same_faces=transferred_faces.Extent()==solid_faces.Extent();
			for(Standard_Integer face=1;same_faces&&face<=transferred_faces.Extent();++face)
				same_faces=solid_faces.Contains(transferred_faces(face));
			if(!same_faces)
			{
				set_error(error,"solid STEP contains faces outside its sole solid; internal ribs, "
					"baffles, and loose construction surfaces are not aerodynamic geometry");
				return false;
			}
			if(solid_faces.IsEmpty())
			{
				set_error(error,"the sole OCCT solid has no boundary faces");
				return false;
			}

			TopTools_IndexedMapOfShape edges;
			TopTools_IndexedDataMapOfShapeListOfShape edge_faces;
			TopExp::MapShapes(solid,TopAbs_EDGE,edges);
			TopExp::MapShapesAndAncestors(solid,TopAbs_EDGE,TopAbs_FACE,edge_faces);
			for(Standard_Integer index=1;index<=edges.Extent();++index)
			{
				const TopoDS_Edge edge=TopoDS::Edge(edges(index));
				if(BRep_Tool::Degenerated(edge))continue;
				const Standard_Integer incident=edge_faces.Contains(edge)
					?edge_faces.FindFromKey(edge).Size():0;
				if(incident!=2)
				{
					set_error(error,"the solid shell is not a two-manifold boundary: edge "+
						std::to_string(index-1)+" has "+std::to_string(incident)+
						" incident faces");
					return false;
				}
			}

			BRepClass3d_SolidClassifier classifier(solid);
			classifier.PerformInfinitePoint(Precision::Confusion());
			if(classifier.State()!=TopAbs_OUT&&classifier.State()!=TopAbs_IN)
			{
				set_error(error,"OCCT could not classify the infinite point against the solid");
				return false;
			}
			const bool reversed=classifier.State()==TopAbs_IN;
			if(reversed)solid.Reverse();

			GProp_GProps properties;
			BRepGProp::VolumeProperties(solid,properties,true,false,false);
			const double volume_mm3=properties.Mass();
			if(!std::isfinite(volume_mm3)||!(volume_mm3>0.0))
			{
				set_error(error,"the outward-oriented OCCT solid has non-positive volume");
				return false;
			}

			ValidatedEnvelope staged;
			staged.solid=solid;
			staged.reversed=reversed;
			staged.volume_m3=volume_mm3*1.0e-9;
			for(TopExp_Explorer face(solid,TopAbs_FACE);face.More();face.Next())
				staged.faces.push_back(TopoDS::Face(face.Current()));
			if(staged.faces.empty())
			{
				set_error(error,"the oriented solid has no traversable boundary faces");
				return false;
			}
			output=std::move(staged);
			return true;
		}

		TriMesh display_mesh(TopoDS_Solid& solid,const std::vector<TopoDS_Face>& faces,
			double deflection_mm,std::string* error)
		{
			if(!(deflection_mm>0.0)||!std::isfinite(deflection_mm))
			{
				set_error(error,"display tessellation deflection must be positive");
				return {};
			}
			BRepTools::Clean(solid);
			BRepMesh_IncrementalMesh mesh_operation(solid,deflection_mm,false,0.35,true);
			mesh_operation.Perform();
			if(!mesh_operation.IsDone())
			{
				set_error(error,"OCCT display tessellation did not complete");
				return {};
			}

			TriMesh mesh;
			double lo[3]={std::numeric_limits<double>::infinity(),
				std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
			double hi[3]={-std::numeric_limits<double>::infinity(),
				-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()};
			for(std::size_t face_id=0;face_id<faces.size();++face_id)
			{
				const TopoDS_Face& face=faces[face_id];
				TopLoc_Location location;
				const Handle(Poly_Triangulation) triangulation=
					BRep_Tool::Triangulation(face,location);
				if(triangulation.IsNull()||triangulation->NbTriangles()<=0)continue;
				const std::uint32_t base=static_cast<std::uint32_t>(mesh.vertex_count());
				const gp_Trsf transform=location.Transformation();
				for(Standard_Integer node=1;node<=triangulation->NbNodes();++node)
				{
					gp_Pnt point=triangulation->Node(node);point.Transform(transform);
					const double coordinate[3]={point.X()*mm_to_m,point.Y()*mm_to_m,
						point.Z()*mm_to_m};
					for(int axis=0;axis<3;++axis)
					{
						mesh.positions_fp64.push_back(coordinate[axis]);
						mesh.positions.push_back(static_cast<float>(coordinate[axis]));
						mesh.normals.push_back(0.0f);
						lo[axis]=std::min(lo[axis],coordinate[axis]);
						hi[axis]=std::max(hi[axis],coordinate[axis]);
					}
					if(triangulation->HasUVNodes())
					{
						const gp_Pnt2d uv=triangulation->UVNode(node);
						mesh.vertex_uv.push_back(static_cast<float>(uv.X()));
						mesh.vertex_uv.push_back(static_cast<float>(uv.Y()));
					}
					else
					{
						const float nan=std::numeric_limits<float>::quiet_NaN();
						mesh.vertex_uv.push_back(nan);mesh.vertex_uv.push_back(nan);
					}
				}
				for(Standard_Integer triangle=1;triangle<=triangulation->NbTriangles();++triangle)
				{
					Standard_Integer a=0,b=0,c=0;
					triangulation->Triangle(triangle).Get(a,b,c);
					if(face.Orientation()==TopAbs_REVERSED)std::swap(b,c);
					const std::uint32_t indices[3]={base+static_cast<std::uint32_t>(a-1),
						base+static_cast<std::uint32_t>(b-1),base+static_cast<std::uint32_t>(c-1)};
					const auto p=[&](std::uint32_t vertex)
					{
						const auto value=mesh.vertex_position_double(vertex);
						return Vec3d{value[0],value[1],value[2]};
					};
					const Vec3d normal=cross(p(indices[1])-p(indices[0]),p(indices[2])-p(indices[0]));
					if(!(length2(normal)>1.0e-40))continue;
					for(std::uint32_t index:indices)
					{
						mesh.indices.push_back(index);
						mesh.normals[3*index]+=static_cast<float>(normal.x);
						mesh.normals[3*index+1]+=static_cast<float>(normal.y);
						mesh.normals[3*index+2]+=static_cast<float>(normal.z);
					}
					mesh.source_face_ids.push_back(static_cast<std::uint32_t>(face_id));
				}
			}
			if(mesh.empty())
			{
				set_error(error,"solid produced no display triangles");
				return {};
			}
			for(std::size_t vertex=0;vertex<mesh.vertex_count();++vertex)
			{
				const Vec3d value=normalized({mesh.normals[3*vertex],mesh.normals[3*vertex+1],
					mesh.normals[3*vertex+2]});
				mesh.normals[3*vertex]=static_cast<float>(value.x);
				mesh.normals[3*vertex+1]=static_cast<float>(value.y);
				mesh.normals[3*vertex+2]=static_cast<float>(value.z);
			}
			for(int axis=0;axis<3;++axis)
			{
				mesh.bbox_min[axis]=static_cast<float>(lo[axis]);
				mesh.bbox_max[axis]=static_cast<float>(hi[axis]);
			}
			return mesh;
		}
	}

	StepGeometry load_step_solid_geometry(const std::string& path,double deflection_mm,
		std::string* error)
	{
		try
		{
			TopoDS_Shape transferred;
			if(!read_step(path,transferred,error))return {};
			ValidatedEnvelope envelope;
			std::string validation_error;
			if(!validate_envelope(transferred,envelope,&validation_error))
			{
				set_error(error,"solid STEP rejected: "+validation_error);
				return {};
			}

			StepGeometry result;
			result.mesh=display_mesh(envelope.solid,envelope.faces,deflection_mm,error);
			if(result.mesh.empty())return {};
			result.closed_solid=make_validated_closed_solid_source(envelope.solid,result.mesh);
			result.solid_envelope.validated=true;
			result.solid_envelope.source_orientation_reversed=envelope.reversed;
			result.solid_envelope.face_count=static_cast<std::uint32_t>(envelope.faces.size());
			result.solid_envelope.volume_m3=envelope.volume_m3;
			if(error)error->clear();
			return result;
		}
		catch(const Standard_Failure& failure)
		{
			set_error(error,std::string("solid STEP OpenCascade exception: ")+
				(failure.GetMessageString()?failure.GetMessageString():"unknown failure"));
		}
		catch(const std::exception& failure)
		{
			set_error(error,std::string("solid STEP import failed: ")+failure.what());
		}
		return {};
	}
}
