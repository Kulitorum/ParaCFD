#pragma once

#include <string>

namespace paracfd::gui
{
	// Viewer-only Cartesian envelope for the coarse visualization resampling. The
	// production solver remains block-AMR and owns its per-brick origins/spacing.
	struct SimInfo
	{
		int nx=0,ny=0,nz=0;
		// h/coarse_h describe the inexpensive 3-D arrow/tracer snapshot; finest_h
		// controls the independent finest-brick-aware scalar-plane sampling density.
		double h=0.25,coarse_h=0.25,finest_h=0.25;
		double Lx=0,Ly=0,Lz=0;
		double U=10,rho=1.225,nu=1.5e-5;
		std::string name="paraglider";
	};
}
