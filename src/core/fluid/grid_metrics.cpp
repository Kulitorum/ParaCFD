// grid_metrics.cpp — fine-core graded-grid generator + GridMetrics owner (see grid_metrics.h).
#include "core/fluid/grid_metrics.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace windcfd::core
{
	namespace
	{
		// Fill a gap of length `gap` (metres) with cells geometrically growing from `h0` at the core
		// side outward at ratio `g`, summing EXACTLY to `gap` with every adjacent ratio ≤ g. Returned
		// widths are ordered from the CORE outward (w[0] nearest the core ≈ h0). Empty if gap≈0.
		std::vector<double> fill_gap(double gap, double h0, double g)
		{
			std::vector<double> w;
			if (gap <= 1e-12 || h0 <= 0.0) return w;
			if (gap <= h0) { w.push_back(gap); return w; } // gap smaller than one core cell
			if (g <= 1.0 + 1e-9)                           // ~uniform: no growth
			{
				int m = std::max(1, (int)std::lround(gap / h0));
				w.assign(m, gap / m);
				return w;
			}
			// Geometric fill with EVERY adjacent ratio held at exactly g (strict growth bound). Choose the
			// cell count m that best fits, then scale to hit `gap` exactly; scaling preserves the internal
			// ratio g and makes the core↔gap junction ratio = the scale factor s ∈ [g^-½, g^½] ⊂ (1/g, g),
			// so no ratio anywhere exceeds g. The junction cell (h0·s) may be a few % below h0 — fine: it
			// is OUTSIDE the building (the windloads guard only checks the exactly-h_fine core cells).
			double m_real = std::log(1.0 + gap * (g - 1.0) / h0) / std::log(g);
			int m = std::max(1, (int)std::lround(m_real));
			double geom_sum = h0 * (std::pow(g, m) - 1.0) / (g - 1.0);
			double s = gap / geom_sum;
			w.resize(m);
			double cur = h0 * s;
			for (int i = 0; i < m; ++i) { w[i] = cur; cur *= g; }
			return w; // ordered from core outward: w[0] == h0·s, growing by exactly g
		}
	}

	std::vector<double> graded_axis_faces(double L, double a, double b, double h_fine, double growth)
	{
		if (a < 0.0) a = 0.0;
		if (b > L) b = L;
		if (b < a) b = a;
		// Degenerate/absent core ⇒ uniform grid at h_fine across [0,L].
		if (h_fine <= 0.0 || (b - a) < h_fine * 0.5)
		{
			int m = std::max(1, (int)std::lround(L / (h_fine > 0.0 ? h_fine : L)));
			std::vector<double> xf(m + 1);
			for (int i = 0; i <= m; ++i) xf[i] = L * (double)i / m;
			return xf;
		}
		// Core cells are EXACTLY h_fine wide (so the fine core is isotropic h_fine cubes on every axis
		// and g.h==h_fine==dx==dy==dz there — windloads/the corner comparison rely on this). Snap the
		// core box to a whole number of h_fine cells centred on the requested box, then clamp to [0,L].
		int n_fine = std::max(1, (int)std::lround((b - a) / h_fine));
		double half = 0.5 * n_fine * h_fine, centre = 0.5 * (a + b);
		double na = centre - half, nb = centre + half;
		if (na < 0.0) { na = 0.0; nb = std::min(L, n_fine * h_fine); }
		if (nb > L) { nb = L; na = std::max(0.0, L - n_fine * h_fine); }
		std::vector<double> left = fill_gap(na, h_fine, growth);       // core(na) → 0, ordered from core out
		std::vector<double> right = fill_gap(L - nb, h_fine, growth);  // core(nb) → L, ordered from core out

		std::vector<double> widths;
		widths.reserve(left.size() + n_fine + right.size());
		for (int i = (int)left.size() - 1; i >= 0; --i) widths.push_back(left[i]); // 0 → na (coarse → fine)
		for (int i = 0; i < n_fine; ++i) widths.push_back(h_fine);                  // uniform core (exactly h_fine)
		for (double w : right) widths.push_back(w);                                 // nb → L (fine → coarse)

		std::vector<double> xf(widths.size() + 1);
		xf[0] = 0.0;
		for (size_t i = 0; i < widths.size(); ++i) xf[i + 1] = xf[i] + widths[i];
		xf.back() = L; // snap the terminal face (kill accumulated fp drift)
		return xf;
	}

	// --- axis helper: faces → widths + centres ----------------------------------
	static void axis_from_faces(const std::vector<double>& xf, std::vector<double>& dx, std::vector<double>& xc)
	{
		int n = (int)xf.size() - 1;
		dx.resize(n);
		xc.resize(n);
		for (int i = 0; i < n; ++i)
		{
			dx[i] = xf[i + 1] - xf[i];
			xc[i] = 0.5 * (xf[i] + xf[i + 1]);
		}
	}

	void GridMetrics::derive_from_faces()
	{
		axis_from_faces(xf_, dx_, xc_);
		axis_from_faces(yf_, dy_, yc_);
		axis_from_faces(zf_, dz_, zc_);
		nx_ = (int)dx_.size();
		ny_ = (int)dy_.size();
		nz_ = (int)dz_.size();
		hmin_ = 1e300;
		for (double v : dx_) hmin_ = std::min(hmin_, v);
		for (double v : dy_) hmin_ = std::min(hmin_, v);
		for (double v : dz_) hmin_ = std::min(hmin_, v);
	}

	GridMetrics GridMetrics::generate(const FineCoreSpec& s)
	{
		GridMetrics m;
		m.xf_ = graded_axis_faces(s.Lx, s.x0, s.x1, s.h_fine, s.growth);
		m.yf_ = graded_axis_faces(s.Ly, s.y0, s.y1, s.h_fine, s.growth);
		m.zf_ = graded_axis_faces(s.Lz, s.z0, s.z1, s.h_fine, s.growth);
		m.h_ = s.h_fine;
		m.derive_from_faces();
		return m;
	}

	GridMetrics GridMetrics::uniform(int nx, int ny, int nz, double h)
	{
		GridMetrics m;
		m.xf_.resize(nx + 1);
		m.yf_.resize(ny + 1);
		m.zf_.resize(nz + 1);
		for (int i = 0; i <= nx; ++i) m.xf_[i] = i * h;
		for (int j = 0; j <= ny; ++j) m.yf_[j] = j * h;
		for (int k = 0; k <= nz; ++k) m.zf_[k] = k * h;
		m.h_ = h;
		m.derive_from_faces();
		return m;
	}

	MacGrid GridMetrics::host_view() const
	{
		MacGrid g;
		g.nx = nx_; g.ny = ny_; g.nz = nz_; g.h = h_; g.hmin = hmin_;
		g.dxa = dx_.data(); g.dya = dy_.data(); g.dza = dz_.data();
		g.xca = xc_.data(); g.yca = yc_.data(); g.zca = zc_.data();
		g.xfa = xf_.data(); g.yfa = yf_.data(); g.zfa = zf_.data();
		return g;
	}

	MacGrid GridMetrics::device_view()
	{
		upload();
		MacGrid g;
		g.nx = nx_; g.ny = ny_; g.nz = nz_; g.h = h_; g.hmin = hmin_;
		g.dxa = ddx_; g.dya = ddy_; g.dza = ddz_;
		g.xca = dxc_; g.yca = dyc_; g.zca = dzc_;
		g.xfa = dxf_; g.yfa = dyf_; g.zfa = dzf_;
		return g;
	}

	static double* upload_array(const std::vector<double>& h)
	{
		double* d = nullptr;
		if (cudaMalloc(&d, h.size() * sizeof(double)) != cudaSuccess) return nullptr;
		cudaMemcpy(d, h.data(), h.size() * sizeof(double), cudaMemcpyHostToDevice);
		return d;
	}

	void GridMetrics::upload()
	{
		if (ddx_) return; // already uploaded
		ddx_ = upload_array(dx_); ddy_ = upload_array(dy_); ddz_ = upload_array(dz_);
		dxc_ = upload_array(xc_); dyc_ = upload_array(yc_); dzc_ = upload_array(zc_);
		dxf_ = upload_array(xf_); dyf_ = upload_array(yf_); dzf_ = upload_array(zf_);
	}

	void GridMetrics::free_device()
	{
		for (double* p : {ddx_, ddy_, ddz_, dxc_, dyc_, dzc_, dxf_, dyf_, dzf_})
			if (p) cudaFree(p);
		ddx_ = ddy_ = ddz_ = dxc_ = dyc_ = dzc_ = dxf_ = dyf_ = dzf_ = nullptr;
	}

	void GridMetrics::move_from(GridMetrics& o)
	{
		nx_ = o.nx_; ny_ = o.ny_; nz_ = o.nz_; h_ = o.h_; hmin_ = o.hmin_;
		dx_ = std::move(o.dx_); dy_ = std::move(o.dy_); dz_ = std::move(o.dz_);
		xc_ = std::move(o.xc_); yc_ = std::move(o.yc_); zc_ = std::move(o.zc_);
		xf_ = std::move(o.xf_); yf_ = std::move(o.yf_); zf_ = std::move(o.zf_);
		ddx_ = o.ddx_; ddy_ = o.ddy_; ddz_ = o.ddz_;
		dxc_ = o.dxc_; dyc_ = o.dyc_; dzc_ = o.dzc_;
		dxf_ = o.dxf_; dyf_ = o.dyf_; dzf_ = o.dzf_;
		o.ddx_ = o.ddy_ = o.ddz_ = o.dxc_ = o.dyc_ = o.dzc_ = o.dxf_ = o.dyf_ = o.dzf_ = nullptr;
	}

	GridMetrics::~GridMetrics() { free_device(); }
}
