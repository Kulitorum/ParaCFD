#pragma once

namespace paracfd::core
{
#if defined(PARACFD_REAL_DOUBLE)
	using Real = double;
#else
	// Production ParaCFD fields target Ada GPUs, where FP32 throughput and capacity are the
	// appropriate default. Geometry preprocessing and reductions remain explicitly double.
	using Real = float;
#endif
}
