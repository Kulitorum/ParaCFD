// scene_io.cpp — see scene_io.h. Binary .scn container: an 8-byte magic, a version, a JSON
// metadata block, then TLV field blobs. Big arrays are stored float32.
#include "gui/scene_io.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>

namespace windcfd::gui
{
	using namespace windcfd::core;
	namespace fs = std::filesystem;

	namespace
	{
		constexpr char kMagic[8] = { 'S', 'C', 'O', 'U', 'R', 'S', 'C', 'N' };
		constexpr std::uint32_t kVersion = 1;
		enum DType : std::uint8_t { DT_F32 = 0, DT_F64 = 1, DT_U8 = 2, DT_U32 = 3 };

		// --- little-endian scalar IO (x64 workstation) --------------------------------
		void put_u8(std::ostream& o, std::uint8_t v) { o.write(reinterpret_cast<const char*>(&v), 1); }
		void put_u32(std::ostream& o, std::uint32_t v) { o.write(reinterpret_cast<const char*>(&v), 4); }
		void put_u64(std::ostream& o, std::uint64_t v) { o.write(reinterpret_cast<const char*>(&v), 8); }
		bool get_u8(std::istream& i, std::uint8_t& v) { return (bool)i.read(reinterpret_cast<char*>(&v), 1); }
		bool get_u32(std::istream& i, std::uint32_t& v) { return (bool)i.read(reinterpret_cast<char*>(&v), 4); }
		bool get_u64(std::istream& i, std::uint64_t& v) { return (bool)i.read(reinterpret_cast<char*>(&v), 8); }

		void blob_head(std::ostream& o, const std::string& key, std::uint8_t dtype, std::uint64_t count)
		{
			put_u32(o, (std::uint32_t)key.size());
			o.write(key.data(), (std::streamsize)key.size());
			put_u8(o, dtype);
			put_u64(o, count);
		}
		void write_f32(std::ostream& o, const std::string& key, const std::vector<double>& a)
		{
			if (a.empty()) return;
			blob_head(o, key, DT_F32, a.size());
			std::vector<float> f(a.size());
			for (std::size_t n = 0; n < a.size(); ++n) f[n] = (float)a[n];
			o.write(reinterpret_cast<const char*>(f.data()), (std::streamsize)(f.size() * sizeof(float)));
		}
		void write_f32raw(std::ostream& o, const std::string& key, const std::vector<float>& a)
		{
			if (a.empty()) return;
			blob_head(o, key, DT_F32, a.size());
			o.write(reinterpret_cast<const char*>(a.data()), (std::streamsize)(a.size() * sizeof(float)));
		}
		void write_f64(std::ostream& o, const std::string& key, const std::vector<double>& a)
		{
			if (a.empty()) return;
			blob_head(o, key, DT_F64, a.size());
			o.write(reinterpret_cast<const char*>(a.data()), (std::streamsize)(a.size() * sizeof(double)));
		}
		void write_u8(std::ostream& o, const std::string& key, const std::vector<unsigned char>& a)
		{
			if (a.empty()) return;
			blob_head(o, key, DT_U8, a.size());
			o.write(reinterpret_cast<const char*>(a.data()), (std::streamsize)a.size());
		}
		void write_u32(std::ostream& o, const std::string& key, const std::vector<std::uint32_t>& a)
		{
			if (a.empty()) return;
			blob_head(o, key, DT_U32, a.size());
			o.write(reinterpret_cast<const char*>(a.data()), (std::streamsize)(a.size() * sizeof(std::uint32_t)));
		}

		struct Blob { std::uint8_t dtype = DT_U8; std::uint64_t count = 0; std::vector<char> raw; };
		std::size_t dtype_size(std::uint8_t d) { return d == DT_F32 ? 4 : d == DT_F64 ? 8 : d == DT_U32 ? 4 : 1; }

		// Read all TLV blobs until EOF into a keyed map.
		bool read_blobs(std::istream& in, std::map<std::string, Blob>& out)
		{
			while (true)
			{
				std::uint32_t klen = 0;
				if (!get_u32(in, klen)) break; // clean EOF
				if (klen == 0 || klen > (1u << 20)) return false;
				std::string key(klen, '\0');
				if (!in.read(&key[0], klen)) return false;
				std::uint8_t dtype = 0; std::uint64_t count = 0;
				if (!get_u8(in, dtype) || !get_u64(in, count)) return false;
				Blob b; b.dtype = dtype; b.count = count;
				const std::size_t bytes = (std::size_t)count * dtype_size(dtype);
				b.raw.resize(bytes);
				if (bytes && !in.read(b.raw.data(), (std::streamsize)bytes)) return false;
				out[key] = std::move(b);
			}
			return true;
		}

		// Return a blob as f64 (float32 blobs are up-cast). Empty if absent.
		std::vector<double> as_f64(const std::map<std::string, Blob>& m, const std::string& key)
		{
			auto it = m.find(key);
			if (it == m.end()) return {};
			const Blob& b = it->second;
			std::vector<double> out((std::size_t)b.count);
			if (b.dtype == DT_F64)
				std::memcpy(out.data(), b.raw.data(), out.size() * sizeof(double));
			else if (b.dtype == DT_F32)
			{
				const float* f = reinterpret_cast<const float*>(b.raw.data());
				for (std::size_t n = 0; n < out.size(); ++n) out[n] = (double)f[n];
			}
			return out;
		}
		std::vector<float> as_f32(const std::map<std::string, Blob>& m, const std::string& key)
		{
			auto it = m.find(key);
			if (it == m.end() || it->second.dtype != DT_F32) return {};
			std::vector<float> out((std::size_t)it->second.count);
			std::memcpy(out.data(), it->second.raw.data(), out.size() * sizeof(float));
			return out;
		}
		std::vector<unsigned char> as_u8(const std::map<std::string, Blob>& m, const std::string& key)
		{
			auto it = m.find(key);
			if (it == m.end() || it->second.dtype != DT_U8) return {};
			std::vector<unsigned char> out((std::size_t)it->second.count);
			if (!out.empty()) std::memcpy(out.data(), it->second.raw.data(), out.size());
			return out;
		}
		std::vector<std::uint32_t> as_u32(const std::map<std::string, Blob>& m, const std::string& key)
		{
			auto it = m.find(key);
			if (it == m.end() || it->second.dtype != DT_U32) return {};
			std::vector<std::uint32_t> out((std::size_t)it->second.count);
			if (!out.empty()) std::memcpy(out.data(), it->second.raw.data(), out.size() * sizeof(std::uint32_t));
			return out;
		}

		// --- struct <-> json -----------------------------------------------------------
		nlohmann::json bc_json(const ChannelBC& b)
		{
			return { {"ymin", b.ymin}, {"ymax", b.ymax}, {"zmin", b.zmin}, {"zmax", b.zmax},
				{"inlet_mode", b.inlet_mode}, {"U_inlet", b.U_inlet}, {"ustar", b.ustar}, {"z0", b.z0},
				{"kappa", b.kappa}, {"bed_datum", b.bed_datum}, {"Uc", b.Uc}, {"solid_mode", b.solid_mode} };
		}
		ChannelBC bc_from(const nlohmann::json& j)
		{
			ChannelBC b;
			if (j.is_null()) return b;
			b.ymin = j.value("ymin", b.ymin); b.ymax = j.value("ymax", b.ymax);
			b.zmin = j.value("zmin", b.zmin); b.zmax = j.value("zmax", b.zmax);
			b.inlet_mode = j.value("inlet_mode", b.inlet_mode); b.U_inlet = j.value("U_inlet", b.U_inlet);
			b.ustar = j.value("ustar", b.ustar); b.z0 = j.value("z0", b.z0); b.kappa = j.value("kappa", b.kappa);
			b.bed_datum = j.value("bed_datum", b.bed_datum); b.Uc = j.value("Uc", b.Uc);
			b.solid_mode = j.value("solid_mode", b.solid_mode);
			return b;
		}
		nlohmann::json pr_json(const ChannelParams& p)
		{
			return { {"rho", p.rho}, {"nu", p.nu}, {"Cs", p.Cs}, {"cfl", p.cfl}, {"safety", p.safety},
				{"proj_tol", p.proj_tol}, {"proj_max_iter", p.proj_max_iter}, {"fixed_dt", p.fixed_dt},
				{"advect_band", p.advect_band} };
		}
		ChannelParams pr_from(const nlohmann::json& j)
		{
			ChannelParams p;
			if (j.is_null()) return p;
			p.rho = j.value("rho", p.rho); p.nu = j.value("nu", p.nu); p.Cs = j.value("Cs", p.Cs);
			p.cfl = j.value("cfl", p.cfl); p.safety = j.value("safety", p.safety);
			p.proj_tol = j.value("proj_tol", p.proj_tol); p.proj_max_iter = j.value("proj_max_iter", p.proj_max_iter);
			p.fixed_dt = j.value("fixed_dt", p.fixed_dt); p.advect_band = j.value("advect_band", p.advect_band);
			return p;
		}
		nlohmann::json info_json(const SimInfo& s)
		{
			return { {"nx", s.nx}, {"ny", s.ny}, {"nz", s.nz}, {"h", s.h}, {"Lx", s.Lx}, {"Ly", s.Ly},
				{"Lz", s.Lz}, {"U", s.U}, {"rho", s.rho}, {"nu", s.nu}, {"cylinder", s.cylinder}, {"name", s.name} };
		}
		SimInfo info_from(const nlohmann::json& j)
		{
			SimInfo s;
			if (j.is_null()) return s;
			s.nx = j.value("nx", s.nx); s.ny = j.value("ny", s.ny); s.nz = j.value("nz", s.nz);
			s.h = j.value("h", s.h); s.Lx = j.value("Lx", s.Lx); s.Ly = j.value("Ly", s.Ly); s.Lz = j.value("Lz", s.Lz);
			s.U = j.value("U", s.U); s.rho = j.value("rho", s.rho); s.nu = j.value("nu", s.nu);
			s.cylinder = j.value("cylinder", s.cylinder); s.name = j.value("name", s.name);
			return s;
		}
		// --- path helpers --------------------------------------------------------------
		bool ends_with_ci(const std::string& s, const std::string& suf)
		{
			if (s.size() < suf.size()) return false;
			for (std::size_t n = 0; n < suf.size(); ++n)
				if (std::tolower((unsigned char)s[s.size() - suf.size() + n]) != std::tolower((unsigned char)suf[n])) return false;
			return true;
		}
		// Strip ".scn" then a trailing ".<digits>" → the scene base stem (path without extension/step).
		std::string scene_stem(const std::string& path)
		{
			std::string p = path;
			if (ends_with_ci(p, ".scn")) p = p.substr(0, p.size() - 4);
			std::size_t dot = p.find_last_of('.');
			if (dot != std::string::npos && dot + 1 < p.size())
			{
				bool all_digit = true;
				for (std::size_t n = dot + 1; n < p.size(); ++n) if (!std::isdigit((unsigned char)p[n])) { all_digit = false; break; }
				if (all_digit) p = p.substr(0, dot);
			}
			return p;
		}
	}

	bool write_scene(const std::string& path, const SceneFile& sf, std::string& warn)
	{
		std::ofstream out(path, std::ios::binary);
		if (!out) { warn = "cannot open '" + path + "' for writing"; return false; }

		const SceneDefinition& d = sf.def;
		const CheckpointState& s = sf.state;

		nlohmann::json j;
		j["version"] = kVersion;
		j["name"] = d.recipe.info.name;
		j["steps"] = s.steps;
		j["sim_time"] = s.sim_time;
		j["grid"] = { {"nx", s.grid.nx}, {"ny", s.grid.ny}, {"nz", s.grid.nz}, {"h", s.grid.h} };
		j["recipe"] = {
			{"base_solid_mode", d.recipe.base_solid_mode}, {"init_u", d.recipe.init_u},
			{"init_v_blip", d.recipe.init_v_blip}, {"source_config", d.recipe.source_config},
			{"bc", bc_json(d.recipe.bc)},
			{"pr", pr_json(d.recipe.pr)}, {"info", info_json(d.recipe.info)} };
		j["state_bc"] = bc_json(s.bc);
		j["solid_mode"] = s.solid_mode;
		j["bed_inlet_mask"] = s.bed_inlet_mask;
		j["has_mesh"] = d.has_mesh;
		if (d.has_mesh)
		{
			j["mesh"] = {
				{"bbox_min", { d.mesh.bbox_min[0], d.mesh.bbox_min[1], d.mesh.bbox_min[2] }},
				{"bbox_max", { d.mesh.bbox_max[0], d.mesh.bbox_max[1], d.mesh.bbox_max[2] }},
				{"vertices", (std::uint64_t)d.mesh.vertex_count()}, {"triangles", (std::uint64_t)d.mesh.triangle_count()} };
			j["place"] = { {"tx", d.place.tx}, {"ty", d.place.ty}, {"tz", d.place.tz},
				{"m", { d.place.m[0], d.place.m[1], d.place.m[2], d.place.m[3], d.place.m[4],
					d.place.m[5], d.place.m[6], d.place.m[7], d.place.m[8] }} };
		}

		const std::string js = j.dump();
		out.write(kMagic, 8);
		put_u32(out, kVersion);
		put_u64(out, (std::uint64_t)js.size());
		out.write(js.data(), (std::streamsize)js.size());

		// Field blobs. Empty vectors are skipped (write_* is a no-op), so absence == not present.
		write_u8(out, "base_solid", d.recipe.base_solid);
		write_f32(out, "u", s.u); write_f32(out, "v", s.v); write_f32(out, "w", s.w); write_f32(out, "p", s.p);
		write_u8(out, "solid", s.solid);
		if (d.has_mesh)
		{
			write_f32raw(out, "mesh_pos", d.mesh.positions);
			write_f32raw(out, "mesh_norm", d.mesh.normals);
			write_u32(out, "mesh_idx", d.mesh.indices);
		}
		if (!out) { warn = "write error on '" + path + "' (disk full?)"; return false; }
		return true;
	}

	namespace
	{
		// Read magic + version + JSON block. Leaves `in` positioned at the first blob.
		bool read_prefix(std::istream& in, nlohmann::json& j, std::string& warn)
		{
			char magic[8] = {};
			if (!in.read(magic, 8) || std::memcmp(magic, kMagic, 8) != 0) { warn = "not a .scn file (bad magic)"; return false; }
			std::uint32_t ver = 0;
			if (!get_u32(in, ver)) { warn = "truncated header"; return false; }
			if (ver == 0 || ver > kVersion) { warn = "unsupported .scn version " + std::to_string(ver); return false; }
			std::uint64_t jlen = 0;
			if (!get_u64(in, jlen) || jlen == 0 || jlen > (1ull << 30)) { warn = "bad metadata length"; return false; }
			std::string js((std::size_t)jlen, '\0');
			if (!in.read(&js[0], (std::streamsize)jlen)) { warn = "truncated metadata"; return false; }
			try { j = nlohmann::json::parse(js); }
			catch (const std::exception& e) { warn = std::string("metadata parse failed: ") + e.what(); return false; }
			return true;
		}
	}

	bool read_scene_header(const std::string& path, SceneHeader& hdr, std::string& warn)
	{
		std::ifstream in(path, std::ios::binary);
		if (!in) { warn = "cannot open '" + path + "'"; return false; }
		nlohmann::json j;
		if (!read_prefix(in, j, warn)) return false;
		hdr.version = j.value("version", 0);
		hdr.name = j.value("name", std::string());
		hdr.steps = j.value("steps", (long long)0);
		hdr.sim_time = j.value("sim_time", 0.0);
		hdr.has_mesh = j.value("has_mesh", false);
		if (j.contains("grid"))
		{
			const auto& g = j.at("grid");
			hdr.nx = g.value("nx", 0); hdr.ny = g.value("ny", 0); hdr.nz = g.value("nz", 0); hdr.h = g.value("h", 0.0);
		}
		return true;
	}

	bool read_scene(const std::string& path, SceneFile& sf, std::string& warn)
	{
		std::ifstream in(path, std::ios::binary);
		if (!in) { warn = "cannot open '" + path + "'"; return false; }
		nlohmann::json j;
		if (!read_prefix(in, j, warn)) return false;

		std::map<std::string, Blob> blobs;
		if (!read_blobs(in, blobs)) { warn = "corrupt blob section in '" + path + "'"; return false; }

		SceneDefinition& d = sf.def;
		CheckpointState& s = sf.state;

		// grid
		if (j.contains("grid"))
		{
			const auto& g = j.at("grid");
			s.grid.nx = g.value("nx", 0); s.grid.ny = g.value("ny", 0); s.grid.nz = g.value("nz", 0); s.grid.h = g.value("h", 0.05);
		}
		d.recipe.grid = s.grid;

		// recipe
		if (j.contains("recipe"))
		{
			const auto& r = j.at("recipe");
			d.recipe.base_solid_mode = r.value("base_solid_mode", (int)SOLID_NOSLIP);
			d.recipe.init_u = r.value("init_u", 1.0);
			d.recipe.init_v_blip = r.value("init_v_blip", 0.0);
			d.recipe.source_config = r.value("source_config", std::string());
			d.recipe.bc = bc_from(r.contains("bc") ? r.at("bc") : nlohmann::json());
			d.recipe.pr = pr_from(r.contains("pr") ? r.at("pr") : nlohmann::json());
			d.recipe.info = info_from(r.contains("info") ? r.at("info") : nlohmann::json());
		}
		d.recipe.base_solid = as_u8(blobs, "base_solid");

		// dynamic state
		s.steps = j.value("steps", (long long)0);
		s.sim_time = j.value("sim_time", 0.0);
		s.bc = bc_from(j.contains("state_bc") ? j.at("state_bc") : nlohmann::json());
		s.solid_mode = j.value("solid_mode", (int)SOLID_NOSLIP);
		s.bed_inlet_mask = j.value("bed_inlet_mask", false);
		s.u = as_f64(blobs, "u"); s.v = as_f64(blobs, "v"); s.w = as_f64(blobs, "w"); s.p = as_f64(blobs, "p");
		s.solid = as_u8(blobs, "solid");

		// display mesh
		d.has_mesh = j.value("has_mesh", false);
		if (d.has_mesh)
		{
			d.mesh.positions = as_f32(blobs, "mesh_pos");
			d.mesh.normals = as_f32(blobs, "mesh_norm");
			d.mesh.indices = as_u32(blobs, "mesh_idx");
			if (j.contains("mesh"))
			{
				const auto& m = j.at("mesh");
				auto lo = m.value("bbox_min", std::vector<float>{ 0, 0, 0 });
				auto hi = m.value("bbox_max", std::vector<float>{ 0, 0, 0 });
				if (lo.size() == 3) d.mesh.bbox_min = { lo[0], lo[1], lo[2] };
				if (hi.size() == 3) d.mesh.bbox_max = { hi[0], hi[1], hi[2] };
			}
			if (j.contains("place"))
			{
				const auto& p = j.at("place");
				d.place.tx = p.value("tx", 0.0); d.place.ty = p.value("ty", 0.0); d.place.tz = p.value("tz", 0.0);
				// Linear part (rotation·scale); absent in pre-gizmo scenes ⇒ keep the identity default.
				if (p.contains("m") && p.at("m").is_array() && p.at("m").size() == 9)
					for (int i = 0; i < 9; ++i) d.place.m[i] = p.at("m")[i].get<double>();
			}
			if (d.mesh.indices.empty()) d.has_mesh = false; // mesh flagged but blobs absent ⇒ treat as none
		}
		return true;
	}

	std::string checkpoint_path_for(const std::string& scene_path, long long step)
	{
		return scene_stem(scene_path) + "." + std::to_string(step) + ".scn";
	}

	bool is_numbered_checkpoint(const std::string& path)
	{
		std::string name = fs::path(path).filename().string();
		static const std::regex re(R"(.*\.[0-9]+\.scn$)", std::regex::icase);
		return std::regex_match(name, re);
	}

	std::vector<CheckpointRef> list_checkpoints(const std::string& scene_path)
	{
		std::vector<CheckpointRef> out;
		fs::path sp(scene_path);
		fs::path dir = sp.has_parent_path() ? sp.parent_path() : fs::path(".");
		const std::string stem = fs::path(scene_stem(scene_path)).filename().string();
		std::error_code ec;
		if (!fs::exists(dir, ec)) return out;

		// A sibling of `stem` is either the base <stem>.scn or a numbered <stem>.<digits>.scn.
		const std::string base_name = stem + ".scn";
		const std::string prefix = stem + ".";
		for (const auto& e : fs::directory_iterator(dir, ec))
		{
			if (!e.is_regular_file()) continue;
			const std::string fn = e.path().filename().string();
			if (!ends_with_ci(fn, ".scn")) continue;
			bool match = ends_with_ci(fn, base_name) && fn.size() == base_name.size(); // exact base
			if (!match && fn.size() > prefix.size() + 4 &&
				ends_with_ci(fn.substr(0, prefix.size()), prefix)) // <stem>.<digits>.scn
			{
				const std::string mid = fn.substr(prefix.size(), fn.size() - prefix.size() - 4); // between "." and ".scn"
				match = !mid.empty() && std::all_of(mid.begin(), mid.end(), [](char c) { return std::isdigit((unsigned char)c) != 0; });
			}
			if (!match) continue;
			SceneHeader h; std::string w;
			CheckpointRef r; r.path = e.path().string();
			if (read_scene_header(r.path, h, w)) r.step = h.steps;
			out.push_back(std::move(r));
		}
		std::sort(out.begin(), out.end(), [](const CheckpointRef& a, const CheckpointRef& b) { return a.step < b.step; });
		return out;
	}
}
