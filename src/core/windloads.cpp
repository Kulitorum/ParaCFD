// windloads.cpp — pressure integration over a voxelized building. See windloads.h.
#include "core/windloads.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace windcfd::core
{
	WindLoads compute_wind_loads(const double* p, const unsigned char* solid, MacGrid g,
		const WindLoadParams& prm, std::vector<float>* out_cell_cp)
	{
		WindLoads L;
		const int nx = g.nx, ny = g.ny, nz = g.nz;
		const double h = g.h, area = h * h;
		if (out_cell_cp) out_cell_cp->assign((std::size_t)g.p_count(), std::numeric_limits<float>::quiet_NaN());
		if (nx <= 0 || ny <= 0 || nz <= 0) return L;

		// --- reference (free-stream) pressure: mean over fluid cells in the upstream inlet slab.
		{
			const int slab = std::min(2, nx);
			double psum = 0.0;
			long long pcnt = 0;
			for (int k = 0; k < nz; ++k)
				for (int j = 0; j < ny; ++j)
					for (int i = 0; i < slab; ++i)
					{
						const int c = g.pidx(i, j, k);
						if (!solid[c]) { psum += p[c]; ++pcnt; }
					}
			L.p_ref = pcnt > 0 ? psum / pcnt : 0.0;
		}

		// --- centroid, z-extent, and projected reference areas (frontal = y-z, plan = x-y).
		std::vector<char> front((std::size_t)ny * nz, 0), plan((std::size_t)nx * ny, 0);
		double cx = 0, cy = 0, cz = 0;
		long long scnt = 0;
		int kmin = nz, kmax = -1;
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
					if (solid[g.pidx(i, j, k)])
					{
						cx += (i + 0.5) * h;
						cy += (j + 0.5) * h;
						cz += (k + 0.5) * h;
						++scnt;
						if (k < kmin) kmin = k;
						if (k > kmax) kmax = k;
						front[(std::size_t)k * ny + j] = 1;
						plan[(std::size_t)j * nx + i] = 1;
					}
		if (scnt == 0) return L;
		cx /= scnt;
		cy /= scnt;
		cz /= scnt;
		L.solid_cells = scnt;
		long long fcnt = 0, plcnt = 0;
		for (char c : front) fcnt += c ? 1 : 0;
		for (char c : plan) plcnt += c ? 1 : 0;
		L.A_frontal = fcnt * area;
		L.A_plan = plcnt * area;
		L.L_ref = (kmax - kmin + 1) * h;

		const double q = 0.5 * prm.rho * prm.u_ref * prm.u_ref; // dynamic pressure [Pa]
		const double qinv = q > 1e-12 ? 1.0 / q : 0.0;
		const int off[6][3] = { { +1, 0, 0 }, { -1, 0, 0 }, { 0, +1, 0 }, { 0, -1, 0 }, { 0, 0, +1 }, { 0, 0, -1 } };

		double Fx = 0, Fy = 0, Fz = 0, Mx = 0, My = 0, Mz = 0;
		double cpmin = 1e300, cpmax = -1e300;
		long long faces = 0;
		for (int k = 0; k < nz; ++k)
			for (int j = 0; j < ny; ++j)
				for (int i = 0; i < nx; ++i)
				{
					const int c = g.pidx(i, j, k);
					if (!solid[c]) continue;
					double cellcp = 0.0;
					int cellfaces = 0;
					for (int f = 0; f < 6; ++f)
					{
						const int ni = i + off[f][0], nj = j + off[f][1], nk = k + off[f][2];
						if (ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz)
							continue; // building shouldn't touch a domain face; skip if it does
						const int nc = g.pidx(ni, nj, nk);
						if (solid[nc]) continue; // interior face
						// Exposed face: outward normal = off[f], surface pressure = fluid neighbour.
						const double pg = p[nc] - L.p_ref; // gauge pressure [Pa]
						const double nX = off[f][0], nY = off[f][1], nZ = off[f][2];
						const double dfx = -pg * nX * area, dfy = -pg * nY * area, dfz = -pg * nZ * area;
						Fx += dfx;
						Fy += dfy;
						Fz += dfz;
						// Moment about the centroid, r x dF, r = face centre - centroid.
						const double rx = (i + 0.5) * h + 0.5 * nX * h - cx;
						const double ry = (j + 0.5) * h + 0.5 * nY * h - cy;
						const double rz = (k + 0.5) * h + 0.5 * nZ * h - cz;
						Mx += ry * dfz - rz * dfy;
						My += rz * dfx - rx * dfz;
						Mz += rx * dfy - ry * dfx;
						const double cp = pg * qinv;
						cellcp += cp;
						++cellfaces;
						++faces;
						if (cp < cpmin) cpmin = cp;
						if (cp > cpmax) cpmax = cp;
					}
					if (out_cell_cp && cellfaces > 0) (*out_cell_cp)[c] = (float)(cellcp / cellfaces);
				}

		L.Fx = Fx;
		L.Fy = Fy;
		L.Fz = Fz;
		L.Mx = Mx;
		L.My = My;
		L.Mz = Mz;
		L.exposed_faces = faces;
		L.cp_min = faces ? cpmin : 0.0;
		L.cp_max = faces ? cpmax : 0.0;

		const double Af = std::max(1e-9, L.A_frontal), Ap = std::max(1e-9, L.A_plan), Lr = std::max(1e-9, L.L_ref);
		L.Cd = Fx * qinv / Af;
		L.Cs = Fy * qinv / Af;
		L.Cl = Fz * qinv / Ap;
		L.CMx = Mx * qinv / (Af * Lr);
		L.CMy = My * qinv / (Af * Lr);
		L.CMz = Mz * qinv / (Ap * Lr);
		return L;
	}
}
