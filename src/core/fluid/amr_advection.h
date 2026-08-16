#pragma once

#include "core/fluid/amr_fields.h"
#include "core/geometry/triangle_bvh.h"

#include <cstddef>
#include <vector>

namespace paracfd::core
{
	// Geometry-independent conservative transport contract used by the composite
	// momentum path. A node is one velocity-component dual control volume (a regular
	// MAC face, a 2:1 tile, or a compact EB aperture state). Every connection is
	// oriented a -> b and contributes one shared upwind momentum flux to both nodes.
	// Fabric is impermeable by topology: no connection is emitted through a patch.
	struct PairwiseMomentumConnection
	{
		int a = -1;
		int b = -1;
		double open_area = 0.0;
		double normal_velocity = 0.0; // signed a -> b
	};

	void conservative_pairwise_momentum_cpu(const std::vector<double>& dual_volume,
		const std::vector<PairwiseMomentumConnection>& connections, double dt,
		std::vector<double>& velocity_x, std::vector<double>& velocity_y,
		std::vector<double>& velocity_z);

	// GPU twin of conservative_pairwise_momentum_cpu. Static volumes and topology are
	// uploaded once; only the per-step signed connection velocities are supplied by the
	// caller. Integrated momentum increments are accumulated pairwise in persistent SoA
	// scratch, so the two endpoints receive exactly opposite transfers.
	class DevicePairwiseMomentumTransport
	{
	public:
		DevicePairwiseMomentumTransport(const std::vector<double>& dual_volume,
			const std::vector<PairwiseMomentumConnection>& connections);
		~DevicePairwiseMomentumTransport();
		DevicePairwiseMomentumTransport(const DevicePairwiseMomentumTransport&) = delete;
		DevicePairwiseMomentumTransport& operator=(const DevicePairwiseMomentumTransport&) = delete;

		void step(Real* velocity_x, Real* velocity_y, Real* velocity_z,
			const Real* connection_normal_velocity, Real dt);
		int node_count() const { return node_count_; }
		int connection_count() const { return connection_count_; }
		std::size_t bytes() const { return bytes_; }

	private:
		int *a_ = nullptr, *b_ = nullptr;
		Real *volume_ = nullptr, *area_ = nullptr;
		Real *delta_x_ = nullptr, *delta_y_ = nullptr, *delta_z_ = nullptr;
		int node_count_ = 0, connection_count_ = 0;
		std::size_t bytes_ = 0;
	};

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
		// Incremental conservative-momentum foundation. This finite-volume MAC update
		// uses one shared upwind flux per same-level lattice link and is currently
		// restricted to one uniform level. Fabric-band links can suppress a flux, but
		// compact EB aperture coupling and 2:1 flux registers are not yet included. It
		// remains separate from production advect() until both pieces are complete.
		void advect_conservative_uniform(DeviceAmrFields& fields, Real dt);
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
