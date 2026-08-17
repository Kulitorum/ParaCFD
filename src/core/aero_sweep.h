#pragma once

#include <algorithm>
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
}
