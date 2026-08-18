#include "core/aero_convergence.h"

#include <cmath>
#include <cstdio>
#include <deque>

namespace
{
	int failures=0;

	void check(bool condition,const char* description)
	{
		std::printf("[aero-convergence] %s: %s\n",condition?"PASS":"FAIL",description);
		if(!condition)++failures;
	}

	bool near(double a,double b,double tolerance=1e-12)
	{
		return std::abs(a-b)<=tolerance;
	}
}

int main()
{
	using namespace paracfd::core;
	AerodynamicRunExitState state;
	check(select_aerodynamic_run_exit(state)==AerodynamicRunExitReason::None,
		"an unconverged run continues");

	state.flow_throughs=state.maximum_flow_throughs;
	check(select_aerodynamic_run_exit(state)==AerodynamicRunExitReason::MaximumFlowThroughs,
		"the configured maximum is an explicit non-converged endpoint");

	state.mean_consecutive=3;
	check(select_aerodynamic_run_exit(state)==AerodynamicRunExitReason::MaximumFlowThroughs,
		"a stale mean counter cannot claim convergence without a ready current mean");
	state.mean_force_ready=true;
	state.mean_force_relative_drift=0.01;
	check(select_aerodynamic_run_exit(state)==AerodynamicRunExitReason::MeanConverged,
		"mean convergence wins over a simultaneous maximum-flow endpoint");

	state.settling_consecutive=3;
	check(select_aerodynamic_run_exit(state)==AerodynamicRunExitReason::Steady,
		"field settling plus a converged mean has highest precedence");

	state.mean_consecutive=0;
	state.flow_throughs=0;
	state.mean_force_relative_drift=state.relative_mean_tolerance;
	check(select_aerodynamic_run_exit(state)==AerodynamicRunExitReason::None,
		"the documented drift threshold is strict");

	std::deque<TimedAerodynamicForce> complete;
	for(int q=0;q<=10;++q)
	{
		const double time=0.1*q;
		complete.push_back({time,{2.0+4.0*time,-3.0,5.0}});
	}
	const BoundedAerodynamicMean complete_mean=bounded_aerodynamic_mean(complete,1.0,1.0);
	check(complete_mean.ready&&complete_mean.complete&&near(complete_mean.coverage,1.0)&&
		near(complete_mean.force.x,4.0)&&near(complete_mean.force.y,-3.0)&&
		near(complete_mean.force.z,5.0),
		"a complete final window is time-weighted exactly");

	const std::deque<TimedAerodynamicForce> partial{{0.75,{3.0,2.0,1.0}},
		{1.0,{5.0,4.0,3.0}}};
	const BoundedAerodynamicMean partial_mean=bounded_aerodynamic_mean(partial,1.0,1.0);
	check(partial_mean.ready&&!partial_mean.complete&&near(partial_mean.coverage,0.25)&&
		near(partial_mean.force.x,4.0)&&near(partial_mean.force.y,3.0)&&
		near(partial_mean.force.z,2.0),
		"a maximum-flow exit retains and labels a partial time-weighted mean");

	const std::deque<TimedAerodynamicForce> insufficient{{1.0,{1.0,2.0,3.0}}};
	const BoundedAerodynamicMean missing_mean=bounded_aerodynamic_mean(insufficient,1.0,1.0);
	check(!missing_mean.ready&&!missing_mean.complete&&near(missing_mean.coverage,0.0),
		"one instantaneous load is not mislabelled as a mean");

	if(failures)
	{
		std::fprintf(stderr,"[aero-convergence] %d deterministic check(s) failed\n",failures);
		return 1;
	}
	std::printf("[aero-convergence] all deterministic checks passed\n");
	return 0;
}
