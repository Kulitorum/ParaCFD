// plane_ops.h — generic inlet-plane apply launchers shared by the SEM and precursor
// inlet-fluctuation generators. up[ny·nz], vp[(ny+1)·nz], wp[ny·(nz+1)] are device
// fluctuation planes indexed k·(dimj)+j to match sem_inlet.cu / the MAC i=0 faces.
#pragma once

#include "core/fluid/mac_grid.h"

namespace paracfd::core
{
	void plane_add_u_gpu(double* u, const double* up, MacGrid g); // u(0,j,k) += up
	void plane_set_vw_gpu(double* v, double* w, const double* vp, const double* wp, MacGrid g); // v/w(0,·)=vp/wp
}
