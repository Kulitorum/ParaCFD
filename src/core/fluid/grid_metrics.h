// grid_metrics.h — fine-core graded-grid generator and the owner of the per-axis metric arrays
// (graded-structured-grid change, design D2/D5). Decouples feature resolution from domain size:
// a UNIFORM h_fine core around the building with a geometrically GRADED, ratio-bounded coarse far
// field. The 3-D grid is the outer product of three independent 1-D meshes (axis-separable, D1).
//
// GridMetrics owns the host arrays (used by the CPU-twin references and the host-side geometry
// consumers — voxelizer, windloads, tracers) AND a device copy (used by the GPU kernels), and
// hands out the matching MacGrid view for each memory space. A uniform grid is representable as
// either the trivial all-h metrics or, more cheaply, a null-metric MacGrid (mac_grid.h accessors
// fall back to h) — this class is for the graded case and for explicit metric round-trip tests.
#pragma once

#include "core/fluid/mac_grid.h"

#include <vector>

namespace windcfd::core
{
	// Fine-core grid specification (SI metres). A uniform-h_fine core inside the box
	// [x0,x1]×[y0,y1]×[z0,z1] with geometric growth (ratio-bounded) out to the domain
	// [0,Lx]×[0,Ly]×[0,Lz]. The grid dimensions (nx,ny,nz) are DERIVED by the generator.
	struct FineCoreSpec
	{
		double Lx = 0, Ly = 0, Lz = 0;                     // domain extent
		double x0 = 0, x1 = 0, y0 = 0, y1 = 0, z0 = 0, z1 = 0; // fine-core box
		double h_fine = 0;                                 // target fine spacing
		double growth = 1.15;                              // max adjacent-cell size ratio (>1)
	};

	// One axis: cumulative face coordinates for a uniform-[a,b] core inside [0,L] with geometric
	// growth of ratio `growth` outward. Returns xf (length = #cells+1), xf[0]=0, xf.back()=L, with
	// no adjacent-cell ratio exceeding `growth`. Pure host, dependency-free — unit-testable.
	std::vector<double> graded_axis_faces(double L, double a, double b, double h_fine, double growth);

	// --- Corner-resolution workflow (design D7) --------------------------------------------------
	// The rounded-vs-sharp corner signal rides on top of the discretization error, so h_fine is not a
	// free knob. Rule: resolve a corner radius r with ≥ ~10 cells across it ⇒ recommended h_fine ≈ r/10
	// (e.g. 30 mm for a 300 mm radius). Hard FLOOR: at h_fine ≥ r/4 the voxelized rounded corner is
	// indistinguishable from sharp, so any measured difference is a discretization artifact — such a run
	// is NOT comparison-grade and the workflow flags it.
	inline double recommended_h_fine(double corner_radius) { return corner_radius > 0.0 ? corner_radius / 10.0 : 0.0; }
	inline bool below_resolution_floor(double h_fine, double corner_radius) { return corner_radius > 0.0 && h_fine >= 0.25 * corner_radius; }

	class GridMetrics
	{
	public:
		GridMetrics() = default;
		~GridMetrics();
		GridMetrics(const GridMetrics&) = delete;
		GridMetrics& operator=(const GridMetrics&) = delete;
		GridMetrics(GridMetrics&& o) noexcept { move_from(o); }
		GridMetrics& operator=(GridMetrics&& o) noexcept { if (this != &o) { free_device(); move_from(o); } return *this; }

		// Build graded metrics from a fine-core spec (derives nx/ny/nz; h = h_fine).
		static GridMetrics generate(const FineCoreSpec& spec);
		// Build explicit uniform metrics (arrays all = h). For scene compat + collapse tests.
		static GridMetrics uniform(int nx, int ny, int nz, double h);

		int nx() const { return nx_; }
		int ny() const { return ny_; }
		int nz() const { return nz_; }
		double h_fine() const { return h_; }
		double h_min() const { return hmin_; }

		// MacGrid whose metric pointers reference the HOST arrays (CPU twins, voxelizer, windloads).
		MacGrid host_view() const;
		// MacGrid whose metric pointers reference the DEVICE arrays (GPU kernels). Uploads on first use.
		MacGrid device_view();

		const std::vector<double>& xf() const { return xf_; }
		const std::vector<double>& yf() const { return yf_; }
		const std::vector<double>& zf() const { return zf_; }
		const std::vector<double>& dx() const { return dx_; }

	private:
		void derive_from_faces(); // fill dx/xc from xf (per axis) + hmin
		void upload();            // host → device (idempotent)
		void free_device();
		void move_from(GridMetrics& o);

		int nx_ = 0, ny_ = 0, nz_ = 0;
		double h_ = 0.0, hmin_ = 0.0;
		std::vector<double> dx_, dy_, dz_, xc_, yc_, zc_, xf_, yf_, zf_;
		// device copies (owned); null until upload().
		double *ddx_ = nullptr, *ddy_ = nullptr, *ddz_ = nullptr;
		double *dxc_ = nullptr, *dyc_ = nullptr, *dzc_ = nullptr;
		double *dxf_ = nullptr, *dyf_ = nullptr, *dzf_ = nullptr;
	};
}
