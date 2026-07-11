// inlet_fluct.h — abstract inlet-fluctuation hook so the M2 open-channel core stays
// decoupled from the M3 turbulence generators (SEM / precursor replay). ChannelFluidCore
// holds an optional InletFluct*; when null the core is byte-identical to M2 (no regression).
// advance() convects/updates the generator state and stores this step's fluctuation planes;
// add_u()/set_vw() write them onto the inlet MAC faces around the Dirichlet mean.
#pragma once

#include "core/fluid/mac_grid.h"

namespace windcfd::core
{
	class InletFluct
	{
	public:
		virtual ~InletFluct() = default;
		virtual void advance(double dt) = 0;                 // update state, compute fluct planes
		virtual void add_u(double* u, MacGrid g) = 0;        // u(0,j,k) += u'
		virtual void set_vw(double* v, double* w, MacGrid g) = 0; // v(0,·)=v', w(0,·)=w'
	};
}
