// Closed-solid STEP import. OpenCascade stays behind this Qt-free interface.
#pragma once

#include "core/geometry/closed_solid.h"
#include "core/geometry/tri_mesh.h"

#include <cstdint>
#include <string>

namespace paracfd::core
{
	struct StepGeometry
	{
		TriMesh mesh; // oriented CFD shell, BVH geometry, display, and load mapping
		struct SolidEnvelopeCertificate
		{
			bool validated = false;
			bool source_orientation_reversed = false;
			std::uint32_t face_count = 0;
			double volume_m3 = 0.0;
		} solid_envelope;
		ClosedSolidSourcePtr closed_solid; // validated BRep retained behind an OCCT-free interface
	};

	// Read exactly one valid, closed, two-manifold OCCT solid. Compounds with loose
	// faces, construction ribs, baffles, or multiple solids are rejected. OCCT is
	// remains authoritative for material classification; tessellation accelerates cutting.
	StepGeometry load_step_solid_geometry(const std::string& path,
		double display_deflection_mm = 2.0, std::string* error = nullptr);
}
