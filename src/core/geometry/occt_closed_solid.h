// Private bridge from the validated STEP import to the OCCT-free solid certificate.
#pragma once

#include "core/geometry/closed_solid.h"

class TopoDS_Solid;

namespace paracfd::core
{
	ClosedSolidSourcePtr make_validated_closed_solid_source(
		const TopoDS_Solid& source_solid,const TriMesh& source_mesh);
}
