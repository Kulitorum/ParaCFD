#pragma once

#include "core/fluid/amr_fields.h"
#include "core/geometry/triangle_bvh.h"

#include <cstddef>
#include <vector>

namespace paracfd::core
{
	// First production AMR advection layer. Backtraces locate the finest active brick
	// through the integer-coordinate GPU hash. A static protection mask covers every
	// face centre within `protection_cells * h` of fabric; those values use a local
	// first-order fallback, so no sample can jump to the opposite side of a sheet.
	// Trilinear and diffusion stencils resolve same-level brick crossings through the
	// GPU coordinate hash. Coarse/fine samples are interpolated but are not yet the
	// final conservative face reconstruction. This remains first-order, not MacCormack.
	class DeviceAmrAdvection
	{
	public:
		DeviceAmrAdvection(const AmrHierarchy& hierarchy, const TriangleBvh& fabric,
			double protection_cells = 2.5);
		~DeviceAmrAdvection();
		DeviceAmrAdvection(const DeviceAmrAdvection&) = delete;
		DeviceAmrAdvection& operator=(const DeviceAmrAdvection&) = delete;

		void advect(DeviceAmrFields& fields, Real dt);
		void diffuse_smagorinsky(DeviceAmrFields& fields, Real molecular_nu, Real cs, Real dt);
		std::size_t protected_face_count() const { return protected_faces_; }
		std::size_t active_face_count() const { return active_faces_; }
		std::size_t bytes() const { return bytes_ + locator_.bytes(); }

	private:
		struct Level
		{
			BrickFieldLayout layout;
			int brick_count = 0;
			Real *u = nullptr, *v = nullptr, *w = nullptr;
			unsigned char *protect_u = nullptr, *protect_v = nullptr, *protect_w = nullptr;
		};
		DeviceAmrLocator locator_;
		DeviceAmrFieldLevelView* device_views_ = nullptr;
		std::vector<Level> levels_;
		std::size_t protected_faces_ = 0, active_faces_ = 0, bytes_ = 0;
	};
}
