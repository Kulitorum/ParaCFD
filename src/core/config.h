// config.h — run configuration loaded from JSON (vendored nlohmann single header).
// Milestone M0 scaffold. Canonical defaults come from RESEARCH.md §1 (verbatim);
// every value is overridable from the config file. Units: SI (m, s, kg, Pa).
#pragma once

#include <string>
#include <cstdint>

namespace windcfd::core
{
	// Run configuration. Defaults are the canonical parameters of RESEARCH.md §1.
	// NOTE (RESEARCH §1): rho and nu MUST be config parameters — the physics is
	// temperature dependent. Defaults here are dry air at ~15 C, sea level (ISA).
	struct Config
	{
		std::string name = "default";

		// --- Domain & grid (default 10 x 10 x 5 m domain, h = 5 cm) ---
		double domain_x = 10.0;  // m
		double domain_y = 10.0;  // m
		double domain_z = 5.0;   // m
		double voxel_h = 0.05;   // m (uniform voxel edge)

		// --- Fluid (air; RESEARCH §1) ---
		double rho = 1.225;     // kg/m^3  dry air, 15 C, sea level (ISA)
		double nu = 1.5e-5;     // m^2/s   dry air, ~15 C
		double U = 1.0;          // m/s     free-stream wind speed

		// --- Bed roughness (log-law wall model) ---
		double d50 = 0.35e-3;    // m  grain size -> wall roughness ks = 2.5*d50, z0 = d50/12

		// --- Output ---
		std::string out_dir = "out";

		// Derived: number of voxels per axis (floor of extent / h).
		std::int64_t nx() const;
		std::int64_t ny() const;
		std::int64_t nz() const;
		std::int64_t cell_count() const;

		// Parse from a JSON string. Throws std::runtime_error on malformed JSON,
		// unknown keys, wrong types, or out-of-range (non-positive) values.
		static Config from_json_string(const std::string& text);

		// Parse from a file path. Throws std::runtime_error if the file cannot be
		// read or the contents are invalid.
		static Config from_file(const std::string& path);

		// Serialize back to a pretty JSON string (round-trip friendly).
		std::string to_json_string() const;
	};
}
