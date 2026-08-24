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
	// Pressure is deliberately mixed precision. The 3-D velocity/advection state stays
	// FP32 in production, while the elliptic solve needs FP64 to resolve the cancellation
	// in A*p on large, highly refined domains.
	using PressureReal = double;
}
