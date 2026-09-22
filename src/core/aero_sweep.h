#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace paracfd::core
{
	// Inclusive, deterministic angle sequence for GUI/CLI aerodynamic sweeps. The
	// requested maximum is retained even when the range is not an integer number of
	// steps, so a UI labelled "maximum" does not silently stop short of it.
	inline std::vector<double> inclusive_angle_sweep(double minimum,double maximum,
		double step,std::size_t maximum_cases=101)
	{
		if(!std::isfinite(minimum)||!std::isfinite(maximum)||!std::isfinite(step)||
			step<=0||maximum<minimum||maximum_cases<1)
			throw std::invalid_argument("invalid aerodynamic angle sweep");
		std::vector<double> result;
		const double scale=std::max({1.0,std::abs(minimum),std::abs(maximum)});
		const double tolerance=64.0*std::numeric_limits<double>::epsilon()*scale;
		for(std::size_t index=0;;++index)
		{
			const double angle=minimum+static_cast<double>(index)*step;
			if(angle>maximum+tolerance)break;
			if(result.size()>=maximum_cases)throw std::invalid_argument("aerodynamic angle sweep has too many cases");
			result.push_back(std::min(angle,maximum));
		}
		if(result.empty()||result.back()<maximum-tolerance)
		{
			if(result.size()>=maximum_cases)throw std::invalid_argument("aerodynamic angle sweep has too many cases");
			result.push_back(maximum);
		}
		return result;
	}

	// A horizontal-wind CFD run with the fixed-attitude wing rotated by gamma is the
	// rotated-frame equivalent of leaving the wing at its trim attitude and inclining
	// the relative wind upward by gamma. `drag` and `lift` are the +X/+Z forces reported
	// in that rotated CFD frame.
	struct FixedAttitudeGlideSample
	{
		double flight_path_degrees = 0;
		double drag = 0;
		double lift = 0;
	};

	struct FixedAttitudeGlideTrim
	{
		bool equilibrium = false;
		double flight_path_degrees = 0;
		double glide_ratio = 0;
		double airspeed = 0;
		double horizontal_speed = 0;
		double sink_rate = 0;
		double force_scale = 0;
		double reference_vertical_support = 0;
		double reference_horizontal_residual = 0;
		double lower_trial_degrees = 0, upper_trial_degrees = 0;
	};

	inline FixedAttitudeGlideTrim solve_fixed_attitude_glide(
		const std::vector<FixedAttitudeGlideSample>& source,double reference_speed,
		double all_up_mass,double gravity=9.80665)
	{
		if(!(reference_speed>0)||!(all_up_mass>0)||!(gravity>0)||
			!std::isfinite(reference_speed)||!std::isfinite(all_up_mass)||!std::isfinite(gravity))
			throw std::invalid_argument("invalid fixed-attitude glide input");
		std::vector<FixedAttitudeGlideSample> samples;
		for(const auto& sample:source)
			if(std::isfinite(sample.flight_path_degrees)&&std::isfinite(sample.drag)&&
				std::isfinite(sample.lift)&&sample.flight_path_degrees>=0&&
				sample.flight_path_degrees<89)samples.push_back(sample);
		std::sort(samples.begin(),samples.end(),[](const auto& a,const auto& b)
			{return a.flight_path_degrees<b.flight_path_degrees;});
		FixedAttitudeGlideTrim best;double closest=std::numeric_limits<double>::infinity();
		auto components=[](const FixedAttitudeGlideSample& sample)
		{
			const double radians=sample.flight_path_degrees*3.14159265358979323846/180.0;
			const double cosine=std::cos(radians),sine=std::sin(radians);
			return std::array<double,2>{cosine*sample.drag-sine*sample.lift,
				sine*sample.drag+cosine*sample.lift};
		};
		for(const auto& sample:samples)
		{
			const auto value=components(sample);if(!(value[1]>0))continue;
			if(std::abs(value[0])<closest)
			{
				closest=std::abs(value[0]);best.flight_path_degrees=sample.flight_path_degrees;
				best.reference_horizontal_residual=value[0];best.reference_vertical_support=value[1];
			}
			if(sample.flight_path_degrees>0&&value[0]==0&&
				(!best.equilibrium||sample.flight_path_degrees<best.flight_path_degrees))
			{
				best.equilibrium=true;best.flight_path_degrees=sample.flight_path_degrees;
				best.reference_horizontal_residual=0;best.reference_vertical_support=value[1];
				best.lower_trial_degrees=best.upper_trial_degrees=sample.flight_path_degrees;
			}
		}
		for(std::size_t index=1;index<samples.size();++index)
		{
			const auto lower_value=components(samples[index-1]),upper_value=components(samples[index]);
			if(!(lower_value[1]>0)||!(upper_value[1]>0)||lower_value[0]==0||upper_value[0]==0||
				(lower_value[0]<0)==(upper_value[0]<0))continue;
			const FixedAttitudeGlideSample& lower=samples[index-1];
			const FixedAttitudeGlideSample& upper=samples[index];
			auto interpolate=[&](double fraction)
			{
				FixedAttitudeGlideSample value;
				value.flight_path_degrees=lower.flight_path_degrees+
					fraction*(upper.flight_path_degrees-lower.flight_path_degrees);
				value.drag=lower.drag+fraction*(upper.drag-lower.drag);
				value.lift=lower.lift+fraction*(upper.lift-lower.lift);
				return value;
			};
			// D and L are linearly interpolated between CFD trials, but the rotation
			// into the gravity frame is trigonometric. Solve that interpolant rather
			// than treating the already-rotated residual as linear in gamma.
			double lo=0,hi=1;bool lo_negative=lower_value[0]<0;
			for(int iteration=0;iteration<64;++iteration)
			{
				const double mid=0.5*(lo+hi);
				const auto mid_value=components(interpolate(mid));
				if((mid_value[0]<0)==lo_negative)lo=mid;else hi=mid;
			}
			FixedAttitudeGlideSample candidate=interpolate(0.5*(lo+hi));
			const auto value=components(candidate);if(!(value[1]>0))continue;
			if(best.equilibrium&&candidate.flight_path_degrees>=best.flight_path_degrees)continue;
			best.equilibrium=true;best.flight_path_degrees=candidate.flight_path_degrees;
			best.reference_horizontal_residual=value[0];best.reference_vertical_support=value[1];
			best.lower_trial_degrees=samples[index-1].flight_path_degrees;
			best.upper_trial_degrees=samples[index].flight_path_degrees;
		}
		if(!best.equilibrium)return best;
		const double radians=best.flight_path_degrees*3.14159265358979323846/180.0;
		best.glide_ratio=1.0/std::tan(radians);
		best.force_scale=all_up_mass*gravity/best.reference_vertical_support;
		best.airspeed=reference_speed*std::sqrt(best.force_scale);
		best.horizontal_speed=best.airspeed*std::cos(radians);
		best.sink_rate=best.airspeed*std::sin(radians);
		return best;
	}
}
