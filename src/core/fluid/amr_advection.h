#pragma once

#include "core/fluid/amr_fields.h"
#include "core/geometry/triangle_bvh.h"

#include <cstddef>
#include <vector>

namespace paracfd::core
{
	// Production AMR advection layer. Backtraces locate the finest active brick
	// through the integer-coordinate GPU hash. Static preprocessing marks faces within
	// `protection_cells * h` of fabric and stores their six Cartesian same-side links.
	// Those faces use bounded, minmod-limited same-side transport (with a first-order
	// fallback wherever the required links are unavailable) plus link-restricted molecular
	// diffusion; far-field faces retain bounded RK2/MacCormack transport and LES. This
	// permits tangential transport near a sheet without per-step triangle traversal or
	// opposite-side stencil sampling.
	// Trilinear and diffusion stencils resolve same-level brick crossings through the
	// GPU coordinate hash. At 2:1 transitions, staggered face and cell values use a
	// one-sided linear prolongation rather than brick-local clamping; the reverse
	// MacCormack correction falls back to the bounded forward RK2 value when its trace
	// changes lattice. This is linearly consistent, while the semi-Lagrangian momentum
	// update itself is not yet a globally conservative finite-volume/refluxed scheme.
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
			Real *forward_u = nullptr, *forward_v = nullptr, *forward_w = nullptr;
			// Bits 0..5 are -/+ xyz links; bit 6 marks the near-fabric band.
			unsigned char *links_u = nullptr, *links_v = nullptr, *links_w = nullptr;
		};
		DeviceAmrLocator locator_;
		DeviceAmrFieldLevelView *device_views_ = nullptr, *device_forward_views_ = nullptr;
		std::vector<Level> levels_;
		std::size_t protected_faces_ = 0, active_faces_ = 0, bytes_ = 0;
	};
}
