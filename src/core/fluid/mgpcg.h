// mgpcg.h — MGPCG pressure solver (McAdams, Sifakis & Teran 2010): conjugate
// gradient preconditioned by one geometric multigrid V-cycle. RESEARCH §3.3.
// Handles the singular all-Neumann closed box by making the RHS zero-mean and
// pinning the pressure gauge via mean-removal (equivalent, for the velocity update,
// to pinning one cell — a global pin composes correctly with the multigrid, PLAN M1).
// A plain damped-Jacobi path (use_mg=false) is kept as the bring-up/debug fallback.
#pragma once

#include "core/fluid/mac_grid.h"

#include <vector>

namespace scour::core
{
	struct SolveResult
	{
		int iters = 0;        // outer CG iterations (or Jacobi sweeps in fallback)
		double relres = 0.0;  // final ||r||_2 / ||b0||_2
		bool converged = false;
	};

	class MgpcgSolver
	{
	public:
		explicit MgpcgSolver(MacGrid finest);
		~MgpcgSolver();
		MgpcgSolver(const MgpcgSolver&) = delete;
		MgpcgSolver& operator=(const MgpcgSolver&) = delete;

		// Solve A x = b (device pointers, cell fields at the finest grid). A is the
		// all-Neumann -Laplacian. If warm_start, x is used as the initial guess,
		// else x is zeroed. RHS is internally shifted to zero mean; x is returned
		// mean-zero. Returns iteration count and final relative residual.
		SolveResult solve(double* x, const double* b, double tol, int max_iter, bool warm_start, bool use_mg = true);

		// V-cycle smoothing knobs (defaults follow RESEARCH §3.3).
		int pre_post_jacobi = 2; // damped-Jacobi sweeps pre and post (symmetric)
		double omega = 2.0 / 3.0;
		int gs_band = 2;         // boundary-band width (cells)
		int coarse_sweeps = 40;  // Jacobi sweeps on the coarsest level

		const std::vector<MacGrid>& levels() const { return grids_; }

	private:
		void vcycle(int level);
		void precondition(const double* r, double* z); // z = M^-1 r, mean-zero
		void apply_finest(const double* p, double* Ap);

		std::vector<MacGrid> grids_;
		std::vector<double*> Lp_, Lrhs_, Lres_, Ltmp_; // per-level work buffers
		int n0_ = 0;                                   // finest cell count
		double *b0_ = nullptr, *r_ = nullptr, *z_ = nullptr, *s_ = nullptr, *As_ = nullptr;
	};
}
