// channel_pressure.h — MGPCG pressure solver for the M2 open channel (McAdams,
// Sifakis & Teran 2010, RESEARCH §3.3), generalised to (a) an interior voxel obstacle
// mask (Neumann ∂p/∂n=0 on solid faces) and (b) a Dirichlet p=0 outlet at xmax (which
// removes the null space — the system is NON-singular, so no mean-removal is done).
// Geometric factor-2 transfer operators are reused from project.cu (mask-agnostic);
// the operator/smoothers are the masked ch_* kernels. Per-level solid masks are built
// by the all-8-solid coarsening rule (keeps fluid connectivity on coarse grids).
#pragma once

#include "core/fluid/channel_bc.h"
#include "core/fluid/mgpcg.h" // SolveResult

#include <vector>

namespace windcfd::core
{
	class ChannelMgpcg
	{
	public:
		explicit ChannelMgpcg(MacGrid finest);
		~ChannelMgpcg();
		ChannelMgpcg(const ChannelMgpcg&) = delete;
		ChannelMgpcg& operator=(const ChannelMgpcg&) = delete;

		// (Re)build the per-level solid masks from the finest-level device mask.
		void build_masks(const unsigned char* finest_solid);

		// Solve A x = b (device cell fields). dir_xmax=1 pins p=0 at the outlet.
		SolveResult solve(double* x, const double* b, double tol, int max_iter, bool warm_start);

		// Tidal reversal (M9): choose which x-face carries the Dirichlet p=0 outlet, from the flow sign.
		// flow_sign≥0 ⇒ xmax (i=nx, the M2/M3 default); flow_sign<0 ⇒ xmin (i=0). The velocity-Dirichlet
		// inlet face is Neumann. Applied to every multigrid level (x coarsens 2× so i=0/i=nx−1 persist).
		void set_outlet_dir(int flow_sign) { dir_xmax_ = flow_sign < 0 ? -1 : 1; }

		int pre_post_jacobi = 2;
		double omega = 2.0 / 3.0;
		int gs_band = 2;
		int coarse_sweeps = 60;

		const std::vector<MacGrid>& levels() const { return grids_; }

	private:
		void vcycle(int level);
		void precondition(const double* r, double* z);
		void apply_finest(const double* p, double* Ap);

		std::vector<MacGrid> grids_;
		std::vector<double*> Lp_, Lrhs_, Ltmp_;
		std::vector<unsigned char*> Lsolid_;
		int n0_ = 0, dir_xmax_ = 1;
		double *r_ = nullptr, *z_ = nullptr, *s_ = nullptr, *As_ = nullptr;
	};
}
