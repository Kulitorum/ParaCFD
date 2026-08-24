#include "core/geometry/occt_closed_solid.h"

#include <BRepClass3d_SolidClassifier.hxx>
#include <Precision.hxx>
#include <TopAbs_State.hxx>
#include <TopoDS_Solid.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace paracfd::core
{
	namespace
	{
		class PlacedClosedSolid final:public ClosedSolidGeometry
		{
		public:
			explicit PlacedClosedSolid(TopoDS_Solid solid,Aabb3d bounds,
				ModelPlacement world_to_source_metres)
				:solid_(std::move(solid)),bounds_(bounds),world_to_source_metres_(world_to_source_metres){}
			const Aabb3d& bounds()const override{return bounds_;}
			void classify_points(const Vec3d* points,std::size_t count,double world_tolerance,
				ClosedSolidPointLocation* locations)const override
			{
				if((count&&(!points||!locations))||world_tolerance<0.0)
					throw std::invalid_argument("invalid closed-solid classification request");
				const double source_tolerance_mm=std::max(Precision::Confusion(),
					1000.0*world_tolerance*world_to_source_metres_.maximum_linear_scale());
				BRepClass3d_SolidClassifier classifier(solid_);
				for(std::size_t index=0;index<count;++index)
				{
					double x=0,y=0,z=0;world_to_source_metres_.apply(points[index].x,
						points[index].y,points[index].z,x,y,z);
					classifier.Perform(gp_Pnt(1000.0*x,1000.0*y,1000.0*z),source_tolerance_mm);
					switch(classifier.State())
					{
					case TopAbs_OUT:locations[index]=ClosedSolidPointLocation::outside;break;
					case TopAbs_IN:locations[index]=ClosedSolidPointLocation::inside;break;
					case TopAbs_ON:locations[index]=ClosedSolidPointLocation::boundary;break;
					default:locations[index]=ClosedSolidPointLocation::unknown;break;
					}
				}
			}
		private:
			TopoDS_Solid solid_;
			Aabb3d bounds_{};
			ModelPlacement world_to_source_metres_{};
		};

		class ValidatedClosedSolidSource final:public ClosedSolidSource
		{
		public:
			explicit ValidatedClosedSolidSource(TopoDS_Solid solid,Aabb3d bounds)
				:solid_(std::move(solid)),bounds_(bounds){}

			ClosedSolidGeometryPtr placed(const ModelPlacement& placement)const override
			{
				ModelPlacement inverse;
				if(!bounds_.valid()||!inverse_placement(placement,inverse))
					throw std::invalid_argument("closed-solid placement is singular");
				Aabb3d output;
				for(int z=0;z<2;++z)for(int y=0;y<2;++y)for(int x=0;x<2;++x)
				{
					const Vec3d source{x?bounds_.hi.x:bounds_.lo.x,
						y?bounds_.hi.y:bounds_.lo.y,z?bounds_.hi.z:bounds_.lo.z};
					Vec3d world;placement.apply(source.x,source.y,source.z,
						world.x,world.y,world.z);
					output.lo.x=std::min(output.lo.x,world.x);
					output.lo.y=std::min(output.lo.y,world.y);
					output.lo.z=std::min(output.lo.z,world.z);
					output.hi.x=std::max(output.hi.x,world.x);
					output.hi.y=std::max(output.hi.y,world.y);
					output.hi.z=std::max(output.hi.z,world.z);
				}
				return std::make_shared<PlacedClosedSolid>(solid_,output,inverse);
			}

		private:
			TopoDS_Solid solid_;
			Aabb3d bounds_{};
		};
	}

	ClosedSolidSourcePtr make_validated_closed_solid_source(
		const TopoDS_Solid& source_solid,const Aabb3d& source_bounds_metres)
	{
		if(source_solid.IsNull()||!source_bounds_metres.valid())
			throw std::invalid_argument("cannot certify empty closed-solid bounds");
		return std::make_shared<ValidatedClosedSolidSource>(source_solid,source_bounds_metres);
	}
}
