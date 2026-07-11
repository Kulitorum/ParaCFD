// config.cpp — JSON config loader implementation (nlohmann single header).
#include "core/config.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <set>

namespace windcfd::core
{
	using nlohmann::json;

	namespace
	{
		// Fetch an optional positive double; throws if present-but-nonpositive or wrong type.
		void get_pos(const json& j, const char* key, double& dst)
		{
			if (!j.contains(key))
				return;
			const auto& v = j.at(key);
			if (!v.is_number())
				throw std::runtime_error(std::string("config key '") + key + "' must be a number");
			double d = v.get<double>();
			if (!(d > 0.0) || !std::isfinite(d))
				throw std::runtime_error(std::string("config key '") + key + "' must be finite and > 0");
			dst = d;
		}

		void get_str(const json& j, const char* key, std::string& dst)
		{
			if (!j.contains(key))
				return;
			const auto& v = j.at(key);
			if (!v.is_string())
				throw std::runtime_error(std::string("config key '") + key + "' must be a string");
			dst = v.get<std::string>();
		}
	}

	std::int64_t Config::nx() const { return static_cast<std::int64_t>(std::floor(domain_x / voxel_h)); }
	std::int64_t Config::ny() const { return static_cast<std::int64_t>(std::floor(domain_y / voxel_h)); }
	std::int64_t Config::nz() const { return static_cast<std::int64_t>(std::floor(domain_z / voxel_h)); }
	std::int64_t Config::cell_count() const { return nx() * ny() * nz(); }

	Config Config::from_json_string(const std::string& text)
	{
		json j = json::parse(text); // throws nlohmann::json::parse_error on malformed input

		if (!j.is_object())
			throw std::runtime_error("config root must be a JSON object");

		// Reject unknown keys so typos surface loudly instead of silently doing nothing.
		static const std::set<std::string> known = {
			"name", "domain_x", "domain_y", "domain_z", "voxel_h",
			"rho", "nu", "U", "d50", "out_dir"};
		for (auto it = j.begin(); it != j.end(); ++it)
		{
			if (known.find(it.key()) == known.end())
				throw std::runtime_error("unknown config key '" + it.key() + "'");
		}

		Config c;
		get_str(j, "name", c.name);
		get_pos(j, "domain_x", c.domain_x);
		get_pos(j, "domain_y", c.domain_y);
		get_pos(j, "domain_z", c.domain_z);
		get_pos(j, "voxel_h", c.voxel_h);
		get_pos(j, "rho", c.rho);
		get_pos(j, "nu", c.nu);
		get_pos(j, "U", c.U);
		get_pos(j, "d50", c.d50);
		get_str(j, "out_dir", c.out_dir);

		if (c.voxel_h > c.domain_x || c.voxel_h > c.domain_y || c.voxel_h > c.domain_z)
			throw std::runtime_error("config 'voxel_h' larger than a domain extent");

		return c;
	}

	Config Config::from_file(const std::string& path)
	{
		std::ifstream in(path, std::ios::binary);
		if (!in)
			throw std::runtime_error("cannot open config file: " + path);
		std::ostringstream ss;
		ss << in.rdbuf();
		return from_json_string(ss.str());
	}

	std::string Config::to_json_string() const
	{
		json j;
		j["name"] = name;
		j["domain_x"] = domain_x;
		j["domain_y"] = domain_y;
		j["domain_z"] = domain_z;
		j["voxel_h"] = voxel_h;
		j["rho"] = rho;
		j["nu"] = nu;
		j["U"] = U;
		j["d50"] = d50;
		j["out_dir"] = out_dir;
		return j.dump(2);
	}
}
