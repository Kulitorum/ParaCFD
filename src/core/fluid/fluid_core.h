// fluid_core.h — the swappable fluid-core interface (PLAN §3). Sediment code may
// consume ONLY {u, v, w, nu_t, tau_b} through this interface; the Stam core is one
// implementation, an LBM core (plan-B, RESEARCH §11.6) must stay drop-in.
//
// M1 provides StamFluidCore (MAC + MacCormack + Smagorinsky + MGPCG). tau_b arrives
// in M3 (fluid/bedshear.cu); for M1 the interface exposes the velocity/eddy-viscosity
// fields on the host for validation.
#pragma once

#include "core/fluid/mac_grid.h"

#include <cstdint>
#include <vector>

namespace windcfd::core
{
	// Host-side snapshot of the MAC fields (for validation/IO). Face-staggered.
	struct FluidSnapshot
	{
		MacGrid grid;
		std::vector<double> u, v, w; // face velocities [m/s]
		std::vector<double> nut;     // cell-centered eddy viscosity [m^2/s]
	};

	class FluidCore
	{
	public:
		virtual ~FluidCore() = default;

		// Advance the fluid state by one adaptive time step; returns the dt taken [s].
		virtual double step() = 0;

		// Copy the current device state to a host snapshot.
		virtual FluidSnapshot snapshot() const = 0;

		virtual const MacGrid& grid() const = 0;
	};
}
