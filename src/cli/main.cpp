// scour CLI — M0 scaffold entry point.
//
// Usage:
//   scour <config.json>              load + validate config, print a summary
//   scour --emit-test-vti <path>     write the canonical M0 synthetic field .vti
//   scour --device                   print the CUDA device name
#include "core/config.h"
#include "core/cuda_probe.h"
#include "core/vti_synthetic.h"
#include "core/vti_writer.h"

#include <cstdio>
#include <cstring>
#include <exception>
#include <string>

namespace
{
	int usage(const char* argv0)
	{
		std::fprintf(stderr,
			"usage:\n"
			"  %s <config.json>            load + validate config, print summary\n"
			"  %s --emit-test-vti <path>   write the M0 synthetic-field .vti\n"
			"  %s --device                 print CUDA device 0 name\n",
			argv0, argv0, argv0);
		return 2;
	}
}

int main(int argc, char** argv)
{
	if (argc < 2)
		return usage(argv[0]);

	const std::string arg1 = argv[1];

	if (arg1 == "--emit-test-vti")
	{
		if (argc < 3)
			return usage(argv[0]);
		const std::string path = argv[2];
		scour::io::VtiImage img = scour::io::make_synthetic_image();
		std::string err;
		if (!scour::io::write_vti(path, img, &err))
		{
			std::fprintf(stderr, "error: %s\n", err.c_str());
			return 1;
		}
		std::printf("wrote synthetic .vti: %s (%d x %d x %d points)\n", path.c_str(), img.nx, img.ny, img.nz);
		return 0;
	}

	if (arg1 == "--device")
	{
		std::string err;
		const std::string name = scour::core::cuda_device_name(&err);
		if (name.empty())
		{
			std::fprintf(stderr, "CUDA device query failed: %s\n", err.c_str());
			return 1;
		}
		std::printf("CUDA device 0: %s\n", name.c_str());
		return 0;
	}

	if (arg1 == "--help" || arg1 == "-h")
		return usage(argv[0]);

	// Otherwise treat arg1 as a config file path.
	try
	{
		const scour::core::Config c = scour::core::Config::from_file(arg1);
		std::printf("loaded config '%s'\n", c.name.c_str());
		std::printf("  domain      : %g x %g x %g m (h = %g m)\n", c.domain_x, c.domain_y, c.domain_z, c.voxel_h);
		std::printf("  grid cells  : %lld x %lld x %lld = %lld\n",
			static_cast<long long>(c.nx()), static_cast<long long>(c.ny()),
			static_cast<long long>(c.nz()), static_cast<long long>(c.cell_count()));
		std::printf("  fluid       : rho = %g kg/m^3, nu = %g m^2/s, U = %g m/s\n", c.rho, c.nu, c.U);
		std::printf("  sediment    : rho_s = %g kg/m^3, d50 = %g m, p = %g, phi = %g deg\n",
			c.rho_s, c.d50, c.porosity, c.phi_repose_deg);
		return 0;
	}
	catch (const std::exception& e)
	{
		std::fprintf(stderr, "config error: %s\n", e.what());
		return 1;
	}
}
