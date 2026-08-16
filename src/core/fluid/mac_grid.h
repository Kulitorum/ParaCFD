// mac_grid.h — MAC staggered-grid geometry, indexing, boundary conditions, and
// BC-aware trilinear sampling. RESEARCH.md §3.1 (staggered MAC), §3.2 (advection),
// §3.3 (projection). Everything here is a pure __host__ __device__ helper so the
// CUDA kernels and their CPU reference implementations share identical arithmetic
// (the GPU-vs-CPU tests then verify launch/index/memory correctness).
//
// Layout (pressure cells nx*ny*nz, cell edge h):
//   p at cell centers          ((i+0.5)h, (j+0.5)h, (k+0.5)h),  i in [0,nx)
//   u on x-faces (nx+1,ny,nz)  ( i   h,   (j+0.5)h, (k+0.5)h),  i in [0,nx]
//   v on y-faces (nx,ny+1,nz)  ((i+0.5)h,  j   h,   (k+0.5)h),  j in [0,ny]
//   w on z-faces (nx,ny,nz+1)  ((i+0.5)h, (j+0.5)h,  k   h ),   k in [0,nz]
// Units: SI (m, s). Velocities m/s, pressure Pa.
#pragma once

#ifdef __CUDACC__
#define PARACFD_HD __host__ __device__
#else
#define PARACFD_HD
#endif

#include <cmath>

namespace paracfd::core
{
	// Wall/face condition. The lid (MOVLID) carries a tangential x-velocity in BC::lid_u.
	enum Wall : int
	{
		WALL_NOSLIP = 0,   // tangential velocity = 0 at the wall (ghost = -interior)
		WALL_FREESLIP = 1, // zero normal gradient of tangential velocity (ghost = +interior)
		WALL_MOVLID = 2    // no-slip lid moving tangentially in +x at speed lid_u
	};

	// Boundary condition for the six domain faces (M1 closed box: cavity).
	struct BC
	{
		int xmin = WALL_NOSLIP, xmax = WALL_NOSLIP;
		int ymin = WALL_FREESLIP, ymax = WALL_FREESLIP;
		int zmin = WALL_NOSLIP, zmax = WALL_MOVLID;
		double lid_u = 1.0; // m/s, tangential (+x) speed of the MOVLID face
	};

	// MAC grid geometry. Spacing is EITHER a single scalar `h` (uniform grid — the metric
	// pointers are null) OR the per-axis "graded" metric arrays below (fine core + coarse far
	// field, graded-structured-grid change). Every accessor falls back to the exact uniform
	// closed form when its array is null, so a null-metric MacGrid is byte-identical to the
	// pre-graded uniform grid and the index topology (pidx/uidx/vidx/widx) is unchanged.
	//
	// The metric pointers are valid ONLY in the memory space of the code holding the struct:
	// a HOST-array copy for the CPU-twin references, a DEVICE-array copy for the GPU kernels.
	// GridMetrics (grid_metrics.h) owns both copies and hands out the correct MacGrid view.
	// Array lengths: cell widths dx[nx]/dy[ny]/dz[nz]; cell centres xc[nx]/yc[ny]/zc[nz];
	// cumulative face coords xf[nx+1]/yf[ny+1]/zf[nz+1] with xf[0]=0, xf[nx]=Lx (SI metres).
	struct MacGrid
	{
		int nx = 0, ny = 0, nz = 0;
		double h = 0.0; // uniform spacing / fine-core h_fine; also the uniform-fallback value
		double hmin = 0.0; // smallest cell dimension across the metric arrays (0 ⇒ uniform ⇒ h)

		const double* dxa = nullptr; // cell widths     [nx] (null ⇒ uniform h)
		const double* dya = nullptr; //                 [ny]
		const double* dza = nullptr; //                 [nz]
		const double* xca = nullptr; // cell centres     [nx]
		const double* yca = nullptr; //                 [ny]
		const double* zca = nullptr; //                 [nz]
		const double* xfa = nullptr; // face coords      [nx+1]
		const double* yfa = nullptr; //                 [ny+1]
		const double* zfa = nullptr; //                 [nz+1]

		PARACFD_HD int p_count() const { return nx * ny * nz; }
		PARACFD_HD int u_count() const { return (nx + 1) * ny * nz; }
		PARACFD_HD int v_count() const { return nx * (ny + 1) * nz; }
		PARACFD_HD int w_count() const { return nx * ny * (nz + 1); }

		PARACFD_HD int pidx(int i, int j, int k) const { return (k * ny + j) * nx + i; }
		PARACFD_HD int uidx(int i, int j, int k) const { return (k * ny + j) * (nx + 1) + i; } // i in [0,nx]
		PARACFD_HD int vidx(int i, int j, int k) const { return (k * (ny + 1) + j) * nx + i; } // j in [0,ny]
		PARACFD_HD int widx(int i, int j, int k) const { return (k * ny + j) * nx + i; }        // k in [0,nz]

		// --- metric accessors (uniform fallback keeps a null-metric grid byte-identical) ---
		// Cell width along each axis.
		PARACFD_HD double dx(int i) const { return dxa ? dxa[i] : h; }
		PARACFD_HD double dy(int j) const { return dya ? dya[j] : h; }
		PARACFD_HD double dz(int k) const { return dza ? dza[k] : h; }
		// Cell-centre world coordinate.
		PARACFD_HD double xc(int i) const { return xca ? xca[i] : (i + 0.5) * h; }
		PARACFD_HD double yc(int j) const { return yca ? yca[j] : (j + 0.5) * h; }
		PARACFD_HD double zc(int k) const { return zca ? zca[k] : (k + 0.5) * h; }
		// Face (staggered-node) world coordinate: i in [0,nx].
		PARACFD_HD double xf(int i) const { return xfa ? xfa[i] : i * h; }
		PARACFD_HD double yf(int j) const { return yfa ? yfa[j] : j * h; }
		PARACFD_HD double zf(int k) const { return zfa ? zfa[k] : k * h; }
		// Domain extent per axis (terminal cumulative face coordinate, NOT n*h on a graded grid).
		PARACFD_HD double Lx() const { return xfa ? xfa[nx] : nx * h; }
		PARACFD_HD double Ly() const { return yfa ? yfa[ny] : ny * h; }
		PARACFD_HD double Lz() const { return zfa ? zfa[nz] : nz * h; }
		// Smallest cell dimension anywhere on the grid — the CFL/diffusive timestep limiter (adaptive_dt)
		// must key off this so the fine core does not go unstable. Uniform (hmin unset) ⇒ h.
		PARACFD_HD double h_min() const { return hmin > 0.0 ? hmin : h; }
		// Centre-to-centre distance between cells (i-1) and i — the distance the pressure gradient
		// / FV-Laplacian uses across the face at u-index i (i in [1,nx-1]). Uniform ⇒ h.
		PARACFD_HD double dxc(int i) const { return xca ? (xca[i] - xca[i - 1]) : h; }
		PARACFD_HD double dyc(int j) const { return yca ? (yca[j] - yca[j - 1]) : h; }
		PARACFD_HD double dzc(int k) const { return zca ? (zca[k] - zca[k - 1]) : h; }
		// Distance between adjacent cell centres i and i+1, OOB-safe: a mirror ghost past the
		// domain edge has the edge cell's width. gapx(i)=xc(i+1)−xc(i) for i in [0,nx−2]. Uniform ⇒ h.
		// Used by the metric diffusion Laplacian (tangential axes) and Smagorinsky cross-derivatives.
		PARACFD_HD double gapx(int i) const { if (i < 0) return dx(0); if (i > nx - 2) return dx(nx - 1); return xc(i + 1) - xc(i); }
		PARACFD_HD double gapy(int j) const { if (j < 0) return dy(0); if (j > ny - 2) return dy(ny - 1); return yc(j + 1) - yc(j); }
		PARACFD_HD double gapz(int k) const { if (k < 0) return dz(0); if (k > nz - 2) return dz(nz - 1); return zc(k + 1) - zc(k); }
	};

	// Non-uniform second derivative along one axis from the three node values (fm, c, fp) at
	// left/right spacings (hL to fm, hR to fp): 2nd-order on a uniform stencil, 1st on stretched.
	// Reduces to (fm + fp − 2c)/h² when hL=hR=h. Used by the metric diffusion Laplacian.
	PARACFD_HD inline double d2_axis(double fm, double c, double fp, double hL, double hR)
	{
		return 2.0 / (hL + hR) * ((fp - c) / hR + (fm - c) / hL);
	}

	// Invertible world↔index map (graded-structured-grid change): fractional index of a physical
	// coordinate x among the monotonically-increasing node positions node[0..N-1] (either the
	// cumulative face coords, N=n+1, or the cell centres, N=n). Returns t with
	// node[floor(t)] <= x < node[floor(t)+1] and t = lo + (x-node[lo])/(node[lo+1]-node[lo]);
	// linearly extrapolates using the end interval outside [node[0], node[N-1]]. Used only on the
	// graded path — the uniform path keeps the exact x/h closed form (see grid_locate_* below).
	PARACFD_HD inline double locate_frac(const double* node, int N, double x)
	{
		int lo = 0, hi = N - 1;
		while (hi - lo > 1)
		{
			int m = (lo + hi) >> 1;
			if (node[m] <= x) lo = m; else hi = m;
		}
		return lo + (x - node[lo]) / (node[lo + 1] - node[lo]); // lo in [0,N-2]
	}
	// Physical x → fractional index in a component's node space (matches the trilerp usage):
	// x-face fields sample at xf[i]; cell-centred axes sample at xc[i] (=(i+0.5)h uniform).
	PARACFD_HD inline double grid_fx(MacGrid g, double x) { return g.xfa ? locate_frac(g.xfa, g.nx + 1, x) : x / g.h; }
	PARACFD_HD inline double grid_fy(MacGrid g, double y) { return g.yfa ? locate_frac(g.yfa, g.ny + 1, y) : y / g.h; }
	PARACFD_HD inline double grid_fz(MacGrid g, double z) { return g.zfa ? locate_frac(g.zfa, g.nz + 1, z) : z / g.h; }
	PARACFD_HD inline double grid_cx(MacGrid g, double x) { return g.xca ? locate_frac(g.xca, g.nx, x) : x / g.h - 0.5; }
	PARACFD_HD inline double grid_cy(MacGrid g, double y) { return g.yca ? locate_frac(g.yca, g.ny, y) : y / g.h - 0.5; }
	PARACFD_HD inline double grid_cz(MacGrid g, double z) { return g.zca ? locate_frac(g.zca, g.nz, z) : z / g.h - 0.5; }

	// --- BC-aware ghost fetch of the three velocity components --------------------
	// Out-of-range indices in a *tangential* direction resolve to ghost values per the
	// wall type; the *normal* index is clamped into the stored face range (the wall
	// faces themselves hold the prescribed normal velocity, 0 for a closed box).

	PARACFD_HD inline double fetch_u(const double* u, MacGrid g, BC bc, int i, int j, int k)
	{
		double factor = 1.0, add = 0.0;
		if (i < 0) i = 0; else if (i > g.nx) i = g.nx; // normal for u
		if (j < 0) { if (bc.ymin == WALL_NOSLIP) factor = -factor; j = 0; }
		else if (j > g.ny - 1) { if (bc.ymax == WALL_NOSLIP) factor = -factor; j = g.ny - 1; }
		if (k < 0) { if (bc.zmin == WALL_NOSLIP) factor = -factor; k = 0; }
		else if (k > g.nz - 1)
		{
			if (bc.zmax == WALL_MOVLID) { add += 2.0 * bc.lid_u; factor = -factor; }
			else if (bc.zmax == WALL_NOSLIP) factor = -factor;
			k = g.nz - 1;
		}
		return factor * u[g.uidx(i, j, k)] + add;
	}

	PARACFD_HD inline double fetch_v(const double* v, MacGrid g, BC bc, int i, int j, int k)
	{
		double factor = 1.0;
		if (j < 0) j = 0; else if (j > g.ny) j = g.ny; // normal for v
		if (i < 0) { if (bc.xmin == WALL_NOSLIP) factor = -factor; i = 0; }
		else if (i > g.nx - 1) { if (bc.xmax == WALL_NOSLIP) factor = -factor; i = g.nx - 1; }
		if (k < 0) { if (bc.zmin == WALL_NOSLIP) factor = -factor; k = 0; }
		else if (k > g.nz - 1)
		{
			// lid tangential v = 0 -> reflect; no-slip likewise
			if (bc.zmax == WALL_MOVLID || bc.zmax == WALL_NOSLIP) factor = -factor;
			k = g.nz - 1;
		}
		return factor * v[g.vidx(i, j, k)];
	}

	PARACFD_HD inline double fetch_w(const double* w, MacGrid g, BC bc, int i, int j, int k)
	{
		double factor = 1.0;
		if (k < 0) k = 0; else if (k > g.nz) k = g.nz; // normal for w
		if (i < 0) { if (bc.xmin == WALL_NOSLIP) factor = -factor; i = 0; }
		else if (i > g.nx - 1) { if (bc.xmax == WALL_NOSLIP) factor = -factor; i = g.nx - 1; }
		if (j < 0) { if (bc.ymin == WALL_NOSLIP) factor = -factor; j = 0; }
		else if (j > g.ny - 1) { if (bc.ymax == WALL_NOSLIP) factor = -factor; j = g.ny - 1; }
		return factor * w[g.widx(i, j, k)];
	}

	// --- Trilinear sampling of a face field at an arbitrary physical point -------
	// grid_x/y/z are the sample point expressed in the component's own node index
	// space (see call sites). Returns the interpolated value; if mn/mx are non-null,
	// also returns the min/max of the 8 corner values (for the MacCormack clamp).

	PARACFD_HD inline double trilerp_u(const double* u, MacGrid g, BC bc, double x, double y, double z, double* mn, double* mx)
	{
		double gx = grid_fx(g, x);    // u node at xf[i]
		double gy = grid_cy(g, y);    // cell centre yc[j]
		double gz = grid_cz(g, z);    // cell centre zc[k]
		int i0 = (int)floor(gx), j0 = (int)floor(gy), k0 = (int)floor(gz);
		double fx = gx - i0, fy = gy - j0, fz = gz - k0;
		double c[2][2][2];
		for (int a = 0; a < 2; ++a)
			for (int b = 0; b < 2; ++b)
				for (int cc = 0; cc < 2; ++cc)
					c[a][b][cc] = fetch_u(u, g, bc, i0 + a, j0 + b, k0 + cc);
		if (mn && mx)
		{
			double lo = c[0][0][0], hi = c[0][0][0];
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int cc = 0; cc < 2; ++cc)
			{ lo = c[a][b][cc] < lo ? c[a][b][cc] : lo; hi = c[a][b][cc] > hi ? c[a][b][cc] : hi; }
			*mn = lo; *mx = hi;
		}
		double c00 = c[0][0][0] * (1 - fx) + c[1][0][0] * fx;
		double c10 = c[0][1][0] * (1 - fx) + c[1][1][0] * fx;
		double c01 = c[0][0][1] * (1 - fx) + c[1][0][1] * fx;
		double c11 = c[0][1][1] * (1 - fx) + c[1][1][1] * fx;
		double c0 = c00 * (1 - fy) + c10 * fy;
		double c1 = c01 * (1 - fy) + c11 * fy;
		return c0 * (1 - fz) + c1 * fz;
	}

	PARACFD_HD inline double trilerp_v(const double* v, MacGrid g, BC bc, double x, double y, double z, double* mn, double* mx)
	{
		double gx = grid_cx(g, x);    // cell centre xc[i]
		double gy = grid_fy(g, y);    // v node at yf[j]
		double gz = grid_cz(g, z);    // cell centre zc[k]
		int i0 = (int)floor(gx), j0 = (int)floor(gy), k0 = (int)floor(gz);
		double fx = gx - i0, fy = gy - j0, fz = gz - k0;
		double c[2][2][2];
		for (int a = 0; a < 2; ++a)
			for (int b = 0; b < 2; ++b)
				for (int cc = 0; cc < 2; ++cc)
					c[a][b][cc] = fetch_v(v, g, bc, i0 + a, j0 + b, k0 + cc);
		if (mn && mx)
		{
			double lo = c[0][0][0], hi = c[0][0][0];
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int cc = 0; cc < 2; ++cc)
			{ lo = c[a][b][cc] < lo ? c[a][b][cc] : lo; hi = c[a][b][cc] > hi ? c[a][b][cc] : hi; }
			*mn = lo; *mx = hi;
		}
		double c00 = c[0][0][0] * (1 - fx) + c[1][0][0] * fx;
		double c10 = c[0][1][0] * (1 - fx) + c[1][1][0] * fx;
		double c01 = c[0][0][1] * (1 - fx) + c[1][0][1] * fx;
		double c11 = c[0][1][1] * (1 - fx) + c[1][1][1] * fx;
		double c0 = c00 * (1 - fy) + c10 * fy;
		double c1 = c01 * (1 - fy) + c11 * fy;
		return c0 * (1 - fz) + c1 * fz;
	}

	PARACFD_HD inline double trilerp_w(const double* w, MacGrid g, BC bc, double x, double y, double z, double* mn, double* mx)
	{
		double gx = grid_cx(g, x);    // cell centre xc[i]
		double gy = grid_cy(g, y);    // cell centre yc[j]
		double gz = grid_fz(g, z);    // w node at zf[k]
		int i0 = (int)floor(gx), j0 = (int)floor(gy), k0 = (int)floor(gz);
		double fx = gx - i0, fy = gy - j0, fz = gz - k0;
		double c[2][2][2];
		for (int a = 0; a < 2; ++a)
			for (int b = 0; b < 2; ++b)
				for (int cc = 0; cc < 2; ++cc)
					c[a][b][cc] = fetch_w(w, g, bc, i0 + a, j0 + b, k0 + cc);
		if (mn && mx)
		{
			double lo = c[0][0][0], hi = c[0][0][0];
			for (int a = 0; a < 2; ++a) for (int b = 0; b < 2; ++b) for (int cc = 0; cc < 2; ++cc)
			{ lo = c[a][b][cc] < lo ? c[a][b][cc] : lo; hi = c[a][b][cc] > hi ? c[a][b][cc] : hi; }
			*mn = lo; *mx = hi;
		}
		double c00 = c[0][0][0] * (1 - fx) + c[1][0][0] * fx;
		double c10 = c[0][1][0] * (1 - fx) + c[1][1][0] * fx;
		double c01 = c[0][0][1] * (1 - fx) + c[1][0][1] * fx;
		double c11 = c[0][1][1] * (1 - fx) + c[1][1][1] * fx;
		double c0 = c00 * (1 - fy) + c10 * fy;
		double c1 = c01 * (1 - fy) + c11 * fy;
		return c0 * (1 - fz) + c1 * fz;
	}

	// Clamp a physical point to the closed domain [0,Lx]x[0,Ly]x[0,Lz] (advection
	// backtraces that would leave the box are clipped to the boundary, RESEARCH §3.2).
	PARACFD_HD inline void clamp_to_domain(MacGrid g, double& x, double& y, double& z)
	{
		double Lx = g.Lx(), Ly = g.Ly(), Lz = g.Lz();
		const double eps = 1e-9;
		if (x < eps) x = eps; else if (x > Lx - eps) x = Lx - eps;
		if (y < eps) y = eps; else if (y > Ly - eps) y = Ly - eps;
		if (z < eps) z = eps; else if (z > Lz - eps) z = Lz - eps;
	}
}
