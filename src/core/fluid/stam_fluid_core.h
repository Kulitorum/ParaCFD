// stam_fluid_core.h — the M1 fluid core: MAC staggered grid, MacCormack advection,
// explicit Smagorinsky, MGPCG projection, no-slip/free-slip/lid walls. Implements the
// FluidCore interface (PLAN §3). Per-step pipeline (RESEARCH §3, research/08):
//   MacCormack-advect u,v,w -> (Smagorinsky nu_t) -> explicit diffusion ->
//   MGPCG projection -> subtract (dt/rho) grad p.
#pragma once

#include "core/fluid/fluid_core.h"
#include "core/fluid/mac_grid.h"
#include "core/fluid/mgpcg.h"

#include <memory>

namespace scour::core
{
	struct StamParams
	{
		double rho = 1027.0; // kg/m^3
		double nu = 1.36e-6; // m^2/s molecular kinematic viscosity (treated explicitly)
		double Cs = 0.0;     // Smagorinsky constant (0 = LES off, M1 gate)
		double cfl = 1.0;    // advective CFL target
		double safety = 0.9; // dt safety factor
		double proj_tol = 1e-4;
		int proj_max_iter = 30;
		bool use_mg = true; // false => plain-Jacobi projection fallback (debug)
	};

	// One divergent-field projection probe (cold start) for the pressure gate.
	struct ProjectionProbe
	{
		SolveResult solve;
		double maxdiv_before = 0.0; // 1/s
		double maxdiv_after = 0.0;  // 1/s
	};

	class StamFluidCore final : public FluidCore
	{
	public:
		StamFluidCore(MacGrid grid, BC bc, StamParams params);
		~StamFluidCore() override;
		StamFluidCore(const StamFluidCore&) = delete;
		StamFluidCore& operator=(const StamFluidCore&) = delete;

		double step() override;
		FluidSnapshot snapshot() const override;
		const MacGrid& grid() const override { return g_; }

		int last_solve_iters() const { return last_iters_; }
		double last_solve_relres() const { return last_relres_; }
		double last_dt() const { return last_dt_; }

		// max |u - marked| over all components; re-marks the current state. First call
		// marks and returns a large sentinel. Used for steady-state detection.
		double max_velocity_change_and_remark();

		// Run one advect+diffuse+project on a scratch copy of the current state (does
		// not modify the live field) and report the projection stats. tol/use_mg let
		// the pressure gate exercise both production (1e-4) and validation (1e-6) modes
		// and the Jacobi fallback.
		ProjectionProbe projection_probe(double tol, int max_iter, bool use_mg);

		MgpcgSolver& solver() { return *solver_; }

	private:
		void compute_dt_fields();

		MacGrid g_;
		BC bc_;
		StamParams pr_;
		std::unique_ptr<MgpcgSolver> solver_;

		double *uA_, *vA_, *wA_;   // state (face velocities)
		double *uB_, *vB_, *wB_;   // advection target / scratch
		double *pu_, *pv_, *pw_;   // projection-probe scratch
		double *um_, *vm_, *wm_;   // steady-state mark
		double *nut_;              // cell eddy viscosity
		double *scratch_;          // phi_hat scratch (max face count)
		double *rhs_, *pres_, *divscr_; // cell fields
		int last_iters_ = 0;
		double last_relres_ = 0.0, last_dt_ = 0.0;
		bool marked_ = false;
	};
}
