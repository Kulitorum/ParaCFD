// closed_solid.h -- OCCT-free interface to a validated aerodynamic solid.
//
// STEP import certifies that the model is one closed solid. The implementation retains
// that BRep behind this interface so CFD material classification remains authoritative;
// tessellation is an acceleration/load surface, not the definition of inside and outside.
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
