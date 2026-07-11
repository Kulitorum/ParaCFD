// main_stub.cpp — scour-gui fallback entry point when the build is configured WITHOUT
// Qt (-DSCOUR_ENABLE_QT=OFF). The real Qt6 + GL 4.3 slice viewer (main.cpp and the
// gui/ sources) is built when SCOUR_ENABLE_QT=ON (the default from G1 onward). This
// keeps a Qt-less configuration compiling for CI hosts that lack Qt.
#include "core/cuda_probe.h"

#include <cstdio>
#include <string>

int main(int, char**)
{
	std::string err;
	const std::string dev = scour::core::cuda_device_name(&err);
	std::printf("scour-gui built without Qt (SCOUR_ENABLE_QT=OFF).\n");
	std::printf("Reconfigure with -DSCOUR_ENABLE_QT=ON for the G1 slice viewer.\n");
	if (!dev.empty())
		std::printf("CUDA device 0: %s\n", dev.c_str());
	else
		std::printf("CUDA device: unavailable (%s)\n", err.c_str());
	return 0;
}
