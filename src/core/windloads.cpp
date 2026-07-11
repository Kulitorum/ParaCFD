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

	// --- LoadAverager: instantaneous samples -> converged time-averaged statistics -------------------
	namespace
	{
		// One Welford update: fold `x` into (n already incremented) running mean `m` + M2 accumulator `s`.
		inline void welford(double x, long long n, double& m, double& s)
		{
			const double d = x - m;
			m += d / (double)n;
			s += d * (x - m);
		}
	}

	void LoadAverager::reset(int p_count, double t_start)
	{
		n_ = 0;
		t_start_ = t_last_ = t_start;
		m_cd_ = m_cl_ = m_cs_ = 0;
		s_cd_ = s_cl_ = s_cs_ = 0;
		min_cd_ = max_cd_ = min_cl_ = max_cl_ = min_cs_ = max_cs_ = 0;
		const std::size_t n = p_count > 0 ? (std::size_t)p_count : 0;
		cp_sum_.assign(n, 0.0);
		cp_cnt_.assign(n, 0);
		cd_hist_.clear();
	}

	void LoadAverager::add(const WindLoads& L, const std::vector<float>& cell_cp, double sim_time)
	{
		++n_;
		t_last_ = sim_time;
		welford(L.Cd, n_, m_cd_, s_cd_);
		welford(L.Cl, n_, m_cl_, s_cl_);
		welford(L.Cs, n_, m_cs_, s_cs_);
		if (n_ == 1)
		{
			min_cd_ = max_cd_ = L.Cd;
			min_cl_ = max_cl_ = L.Cl;
			min_cs_ = max_cs_ = L.Cs;
		}
		else
		{
			min_cd_ = std::min(min_cd_, L.Cd); max_cd_ = std::max(max_cd_, L.Cd);
			min_cl_ = std::min(min_cl_, L.Cl); max_cl_ = std::max(max_cl_, L.Cl);
			min_cs_ = std::min(min_cs_, L.Cs); max_cs_ = std::max(max_cs_, L.Cs);
		}
		// Per-cell Cp running sum (surface cells only: finite values; NaN marks off the surface).
		const std::size_t n = std::min(cp_sum_.size(), cell_cp.size());
		for (std::size_t c = 0; c < n; ++c)
		{
			const float cp = cell_cp[c];
			if (cp == cp) // false only for NaN
			{
				cp_sum_[c] += (double)cp;
				++cp_cnt_[c];
			}
		}
		cd_hist_.push_back((float)m_cd_);
	}

	WindLoadStats LoadAverager::result(std::vector<float>* out_mean_cp) const
	{
		WindLoadStats st;
		st.samples = n_;
		st.duration = t_last_ - t_start_;
		st.Cd_mean = m_cd_; st.Cl_mean = m_cl_; st.Cs_mean = m_cs_;
		if (n_ > 0)
		{
			st.Cd_rms = std::sqrt(s_cd_ / (double)n_);
			st.Cl_rms = std::sqrt(s_cl_ / (double)n_);
			st.Cs_rms = std::sqrt(s_cs_ / (double)n_);
		}
		st.Cd_min = min_cd_; st.Cd_max = max_cd_;
		st.Cl_min = min_cl_; st.Cl_max = max_cl_;
		st.Cs_min = min_cs_; st.Cs_max = max_cs_;
		// Convergence hint: cumulative mean of the first half of the window vs the whole window. If the
		// flow is settled these agree (drift -> 0); if it is still trending they differ. cd_hist_[k-1] is
		// the running mean after k samples, so the first-half mean is cd_hist_[n/2 - 1].
		if (n_ >= 20 && (std::size_t)(n_ / 2) <= cd_hist_.size() && (n_ / 2) >= 1)
		{
			const double mean_all = m_cd_;
			const double mean_half = (double)cd_hist_[(std::size_t)(n_ / 2) - 1];
			st.Cd_drift = std::fabs(mean_all - mean_half) / std::max(1e-9, std::fabs(mean_all));
		}
		// Time-averaged per-cell Cp + its range over the surface.
		double cpmin = 1e300, cpmax = -1e300;
		if (out_mean_cp) out_mean_cp->assign(cp_sum_.size(), std::numeric_limits<float>::quiet_NaN());
		for (std::size_t c = 0; c < cp_sum_.size(); ++c)
		{
			if (cp_cnt_[c] > 0)
			{
				const double mean = cp_sum_[c] / (double)cp_cnt_[c];
				if (out_mean_cp) (*out_mean_cp)[c] = (float)mean;
				if (mean < cpmin) cpmin = mean;
				if (mean > cpmax) cpmax = mean;
			}
		}
		if (cpmax >= cpmin)
		{
			st.cp_min = cpmin;
			st.cp_max = cpmax;
		}
		return st;
	}
}
