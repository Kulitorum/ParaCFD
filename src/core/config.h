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
	// temperature/salinity dependent. Defaults here are 10 C seawater, S=35.
	struct Config
	{
		std::string name = "default";

		// --- Domain & grid (RESEARCH §1: 10 x 10 x 5 m, h = 5 cm ranking) ---
		double domain_x = 10.0;  // m
		double domain_y = 10.0;  // m
		double domain_z = 5.0;   // m
		double voxel_h = 0.05;   // m (uniform voxel edge)

		// --- Fluid (RESEARCH §1) ---
		double rho = 1027.0;     // kg/m^3  seawater S=35, 10 C (ITTC)
		double nu = 1.36e-6;     // m^2/s   10 C seawater
		double U = 1.0;          // m/s     depth-averaged current

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
