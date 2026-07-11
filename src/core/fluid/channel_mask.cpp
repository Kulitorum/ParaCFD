// channel_mask.cpp — see channel_mask.h.
#include "core/fluid/channel_mask.h"

#include <cmath>

namespace scour::core
{
	int build_cylinder_mask(MacGrid g, double xc, double yc, double R, std::vector<unsigned char>& solid)
	{
		solid.assign((size_t)g.p_count(), 0);
		int count = 0;
		double R2 = R * R;
		for (int k = 0; k < g.nz; ++k)
			for (int j = 0; j < g.ny; ++j)
				for (int i = 0; i < g.nx; ++i)
				{
					double x = (i + 0.5) * g.h - xc, y = (j + 0.5) * g.h - yc;
					if (x * x + y * y <= R2) { solid[g.pidx(i, j, k)] = 1; ++count; }
				}
		return count;
	}

	void dilate_mask(MacGrid g, const std::vector<unsigned char>& solid, int band, std::vector<unsigned char>& nearsolid)
	{
		nearsolid.assign((size_t)g.p_count(), 0);
		if (band < 0) band = 0;
		for (int k = 0; k < g.nz; ++k)
			for (int j = 0; j < g.ny; ++j)
				for (int i = 0; i < g.nx; ++i)
				{
					bool near = false;
					for (int dk = -band; dk <= band && !near; ++dk)
						for (int dj = -band; dj <= band && !near; ++dj)
							for (int di = -band; di <= band && !near; ++di)
							{
								int ii = i + di, jj = j + dj, kk = k + dk;
								if (ii < 0 || ii >= g.nx || jj < 0 || jj >= g.ny || kk < 0 || kk >= g.nz) continue;
								if (solid[g.pidx(ii, jj, kk)]) near = true;
							}
					nearsolid[g.pidx(i, j, k)] = near ? 1 : 0;
				}
	}
}
