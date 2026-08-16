// main_stub.cpp — paracfd-gui fallback entry point when the build is configured WITHOUT
// Qt (-DPARACFD_ENABLE_QT=OFF). The Qt6/OpenGL paraglider application is built when
// PARACFD_ENABLE_QT=ON (the default). This
// keeps a Qt-less configuration compiling for CI hosts that lack Qt.
#include "core/cuda_probe.h"

#include <cstdio>
#include <string>

int main(int, char**)
{
	std::string err;
	const std::string dev = paracfd::core::cuda_device_name(&err);
	std::printf("paracfd-gui built without Qt (PARACFD_ENABLE_QT=OFF).\n");
	std::printf("Reconfigure with -DPARACFD_ENABLE_QT=ON for the paraglider viewer.\n");
	if (!dev.empty())
		std::printf("CUDA device 0: %s\n", dev.c_str());
	else
		std::printf("CUDA device: unavailable (%s)\n", err.c_str());
	return 0;
}
