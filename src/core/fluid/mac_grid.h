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
#define WINDCFD_HD __host__ __device__
#else
#define WINDCFD_HD
#endif

#include <cmath>

namespace windcfd::core
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

	struct MacGrid
	{
		int nx = 0, ny = 0, nz = 0;
		double h = 0.0;

		WINDCFD_HD int p_count() const { return nx * ny * nz; }
		WINDCFD_HD int u_count() const { return (nx + 1) * ny * nz; }
		WINDCFD_HD int v_count() const { return nx * (ny + 1) * nz; }
		WINDCFD_HD int w_count() const { return nx * ny * (nz + 1); }

		WINDCFD_HD int pidx(int i, int j, int k) const { return (k * ny + j) * nx + i; }
		WINDCFD_HD int uidx(int i, int j, int k) const { return (k * ny + j) * (nx + 1) + i; } // i in [0,nx]
		WINDCFD_HD int vidx(int i, int j, int k) const { return (k * (ny + 1) + j) * nx + i; } // j in [0,ny]
		WINDCFD_HD int widx(int i, int j, int k) const { return (k * ny + j) * nx + i; }        // k in [0,nz]
	};

	// --- BC-aware ghost fetch of the three velocity components --------------------
	// Out-of-range indices in a *tangential* direction resolve to ghost values per the
	// wall type; the *normal* index is clamped into the stored face range (the wall
	// faces themselves hold the prescribed normal velocity, 0 for a closed box).

	WINDCFD_HD inline double fetch_u(const double* u, MacGrid g, BC bc, int i, int j, int k)
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

	WINDCFD_HD inline double fetch_v(const double* v, MacGrid g, BC bc, int i, int j, int k)
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

	WINDCFD_HD inline double fetch_w(const double* w, MacGrid g, BC bc, int i, int j, int k)
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

	WINDCFD_HD inline double trilerp_u(const double* u, MacGrid g, BC bc, double x, double y, double z, double* mn, double* mx)
	{
		double gx = x / g.h;          // u node at i*h
		double gy = y / g.h - 0.5;    // (j+0.5)h
		double gz = z / g.h - 0.5;    // (k+0.5)h
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

	WINDCFD_HD inline double trilerp_v(const double* v, MacGrid g, BC bc, double x, double y, double z, double* mn, double* mx)
	{
		double gx = x / g.h - 0.5;
		double gy = y / g.h;
		double gz = z / g.h - 0.5;
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

	WINDCFD_HD inline double trilerp_w(const double* w, MacGrid g, BC bc, double x, double y, double z, double* mn, double* mx)
	{
		double gx = x / g.h - 0.5;
		double gy = y / g.h - 0.5;
		double gz = z / g.h;
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
	WINDCFD_HD inline void clamp_to_domain(MacGrid g, double& x, double& y, double& z)
	{
		double Lx = g.nx * g.h, Ly = g.ny * g.h, Lz = g.nz * g.h;
		const double eps = 1e-9;
		if (x < eps) x = eps; else if (x > Lx - eps) x = Lx - eps;
		if (y < eps) y = eps; else if (y > Ly - eps) y = Ly - eps;
		if (z < eps) z = eps; else if (z > Lz - eps) z = Lz - eps;
	}
}
