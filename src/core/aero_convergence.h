#pragma once

#include "core/geometry/triangle_bvh.h"

#include <algorithm>
#include <cmath>
#include <deque>

namespace paracfd::core
{
	struct TimedAerodynamicForce
	{
		double time=0;
		Vec3d force{};
	};

	struct AerodynamicMeanConvergence
	{
		bool ready=false;
		Vec3d previous_mean{},current_mean{};
		double relative_drift=0;
		double current_rms_fraction=0;
	};

	enum class AerodynamicRunExitReason
	{
		None,
		Steady,
		MeanConverged,
		MaximumFlowThroughs
	};

	struct AerodynamicRunExitState
	{
		int settling_consecutive=0;
		int mean_consecutive=0;
		bool mean_force_ready=false;
		double mean_force_relative_drift=0;
		double relative_mean_tolerance=0.02;
		double flow_throughs=0;
		double maximum_flow_throughs=2.5;
	};

	// Single runs and sweep cases use this same ordering. A statistically steady
	// result wins over the mean-only result, and either converged result wins over
	// the explicit maximum-flow endpoint reached on the same sample.
	inline AerodynamicRunExitReason select_aerodynamic_run_exit(
		const AerodynamicRunExitState& state)
	{
		const bool mean_converged=state.mean_force_ready&&
			state.mean_force_relative_drift<state.relative_mean_tolerance;
		if(state.settling_consecutive>=3&&mean_converged)
			return AerodynamicRunExitReason::Steady;
		if(state.mean_consecutive>=3&&mean_converged)
			return AerodynamicRunExitReason::MeanConverged;
		if(state.flow_throughs>=state.maximum_flow_throughs)
			return AerodynamicRunExitReason::MaximumFlowThroughs;
		return AerodynamicRunExitReason::None;
	}

	struct BoundedAerodynamicMean
	{
		bool ready=false;
		bool complete=false;
		Vec3d force{};
		double coverage=0;
	};

	// Time-weighted load over the available tail of the requested final window.
	// This remains usable (and explicitly marked partial) when a bounded run reaches
	// its maximum-flow endpoint before a complete adjacent-window comparison exists.
	inline BoundedAerodynamicMean bounded_aerodynamic_mean(
		const std::deque<TimedAerodynamicForce>& samples,double end_time,double window)
	{
		BoundedAerodynamicMean out;
		if(!(window>0)||samples.size()<2)return out;
		const double oldest=end_time-window;
		Vec3d integral{};
		double covered=0;
		int segments=0;
		for(std::size_t q=0;q+1<samples.size();++q)
		{
			const auto& a=samples[q];
			const auto& b=samples[q+1];
			if(!(b.time>a.time))continue;
			const double left=std::max(oldest,a.time),right=std::min(end_time,b.time);
			if(!(right>left))continue;
			const double inverse=1.0/(b.time-a.time);
			const double ta=(left-a.time)*inverse,tb=(right-a.time)*inverse,dt=right-left;
			const Vec3d fa=a.force+(b.force-a.force)*ta;
			const Vec3d fb=a.force+(b.force-a.force)*tb;
			integral=integral+(fa+fb)*(0.5*dt);
			covered+=dt;
			++segments;
		}
		out.ready=segments>0&&covered>0;
		out.coverage=std::clamp(covered/window,0.0,1.0);
		out.complete=out.ready&&out.coverage>=0.98;
		if(out.ready)out.force=integral/covered;
		return out;
	}

	namespace detail
	{
		struct ForceWindowStatistics
		{
			bool ready=false;
			Vec3d mean{};
			double rms=0;
		};

		inline ForceWindowStatistics aerodynamic_force_window(
			const std::deque<TimedAerodynamicForce>& samples,double lo,double hi)
		{
			ForceWindowStatistics out;
			if(!(hi>lo)||samples.size()<2)return out;
			Vec3d integral{},square_integral{};
			double covered=0;
			int segments=0;
			for(std::size_t q=0;q+1<samples.size();++q)
			{
				const auto& a=samples[q];
				const auto& b=samples[q+1];
				if(!(b.time>a.time))continue;
				const double left=std::max(lo,a.time),right=std::min(hi,b.time);
				if(!(right>left))continue;
				const double inv=1.0/(b.time-a.time);
				const double ta=(left-a.time)*inv,tb=(right-a.time)*inv,dt=right-left;
				const Vec3d fa=a.force+(b.force-a.force)*ta;
				const Vec3d fb=a.force+(b.force-a.force)*tb;
				integral=integral+(fa+fb)*(0.5*dt);
				square_integral.x+=dt*(fa.x*fa.x+fa.x*fb.x+fb.x*fb.x)/3.0;
				square_integral.y+=dt*(fa.y*fa.y+fa.y*fb.y+fb.y*fb.y)/3.0;
				square_integral.z+=dt*(fa.z*fa.z+fa.z*fb.z+fb.z*fb.z)/3.0;
				covered+=dt;
				++segments;
			}
			if(segments<3||covered<0.98*(hi-lo))return out;
			out.ready=true;
			out.mean=integral/covered;
			const double mean_square=(square_integral.x+square_integral.y+square_integral.z)/covered;
			out.rms=std::sqrt(std::max(0.0,mean_square-length2(out.mean)));
			return out;
		}
	}

	// Compare two adjacent time-weighted mean-force windows. Instantaneous RMS is
	// reported but deliberately does not enter relative_drift: a periodic LES wake may
	// be statistically converged without ever becoming instantaneously steady.
	inline AerodynamicMeanConvergence assess_aerodynamic_mean_convergence(
		const std::deque<TimedAerodynamicForce>& samples,double end_time,double total_window)
	{
		AerodynamicMeanConvergence out;
		if(!(total_window>0))return out;
		const double split=end_time-0.5*total_window;
		const auto previous=detail::aerodynamic_force_window(samples,end_time-total_window,split);
		const auto current=detail::aerodynamic_force_window(samples,split,end_time);
		if(!previous.ready||!current.ready)return out;
		out.ready=true;
		out.previous_mean=previous.mean;
		out.current_mean=current.mean;
		const Vec3d average=(previous.mean+current.mean)*0.5;
		const Vec3d difference=current.mean-previous.mean;
		const double vector_scale=std::max(1.0,std::sqrt(length2(average)));
		const double component_floor=0.1*vector_scale;
		out.relative_drift=std::max({
			std::abs(difference.x)/std::max(component_floor,std::abs(average.x)),
			std::abs(difference.y)/std::max(component_floor,std::abs(average.y)),
			std::abs(difference.z)/std::max(component_floor,std::abs(average.z))});
		out.current_rms_fraction=current.rms/vector_scale;
		return out;
	}
}
