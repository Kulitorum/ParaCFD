// closed_solid.h -- OCCT-free interface to a validated aerodynamic solid.
//
// STEP import certifies the closed BRep. CFD volumes and loads use its discrete
// triangle envelope consistently; the BRep remains available for independent
// material checks. The full envelope is retained when the display mesh is cropped.
#pragma once

#include "core/geometry/model_placement.h"
#include "core/geometry/triangle_bvh.h"

#include <cstddef>
#include <memory>

namespace paracfd::core
{
	enum class ClosedSolidPointLocation : unsigned char
	{
		outside,
		inside,
		boundary,
		unknown
	};

	struct SolidFaceTriangle
	{
		Vec3d a{},b{},c{};
	};

	class ClosedSolidGeometry
	{
	public:
		virtual ~ClosedSolidGeometry()=default;
		virtual const Aabb3d& bounds()const=0;
		virtual void classify_points(const Vec3d* points,std::size_t count,
			double world_tolerance,ClosedSolidPointLocation* locations)const=0;
        virtual void classify_discrete_points(const Vec3d* points,std::size_t count,
            double world_tolerance,ClosedSolidPointLocation* locations)const
        {
            classify_points(points,count,world_tolerance,locations);
        }
	};
	using ClosedSolidGeometryPtr=std::shared_ptr<const ClosedSolidGeometry>;

	class ClosedSolidSource
	{
	public:
		virtual ~ClosedSolidSource()=default;
		virtual ClosedSolidGeometryPtr placed(
			const ModelPlacement& source_metres_to_world)const=0;
	};
	using ClosedSolidSourcePtr=std::shared_ptr<const ClosedSolidSource>;
}
