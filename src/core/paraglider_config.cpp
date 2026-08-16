#include "core/paraglider_config.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <sstream>

namespace paracfd::core
{
	namespace
	{
		using json = nlohmann::json;
		double number(const json& j, const char* key, double fallback) { auto it = j.find(key); return it != j.end() && it->is_number() ? it->get<double>() : fallback; }
		int integer(const json& j, const char* key, int fallback) { auto it = j.find(key); return it != j.end() && it->is_number_integer() ? it->get<int>() : fallback; }
	}

	bool load_paraglider_config(const std::string& path, ParagliderConfig& out, std::string* error)
	{
		try
		{
			std::ifstream in(path);
			if (!in) { if (error) *error = "cannot open paraglider config '" + path + "'"; return false; }
			nlohmann::json root; in >> root;
			ParagliderConfig c;
			if (auto it = root.find("step"); it != root.end() && it->is_object())
			{
				c.step_path = it->value("path", std::string{});
				c.tessellation_deflection_mm = number(*it, "tessellation_deflection_mm", c.tessellation_deflection_mm);
			}
			if (auto it = root.find("placement"); it != root.end() && it->is_object())
			{
				c.placement.tx = number(*it, "tx", 0); c.placement.ty = number(*it, "ty", 0); c.placement.tz = number(*it, "tz", 0);
				if (auto m = it->find("matrix"); m != it->end() && m->is_array() && m->size() == 9)
					for (int i = 0; i < 9; ++i) c.placement.m[i] = (*m)[i].get<double>();
			}
			if (auto it = root.find("freestream"); it != root.end() && it->is_object())
			{
				c.freestream.speed = number(*it, "speed", c.freestream.speed); c.freestream.rho = number(*it, "rho", c.freestream.rho); c.freestream.nu = number(*it, "nu", c.freestream.nu);
			}
			if (auto it = root.find("domain"); it != root.end() && it->is_object())
			{
				c.domain.upstream_margin = number(*it, "upstream_margin", c.domain.upstream_margin);
				c.domain.downstream_margin = number(*it, "downstream_margin", c.domain.downstream_margin);
				c.domain.lateral_margin = number(*it, "lateral_margin", c.domain.lateral_margin);
				c.domain.vertical_margin = number(*it, "vertical_margin", c.domain.vertical_margin);
			}
			if (auto it = root.find("amr"); it != root.end() && it->is_object())
			{
				c.amr.base_cell_size = number(*it, "base_cell_size", c.amr.base_cell_size);
				c.amr.max_levels = integer(*it, "max_levels", c.amr.max_levels); c.amr.brick_size = integer(*it, "brick_size", c.amr.brick_size); c.amr.ghost_cells = integer(*it, "ghost_cells", c.amr.ghost_cells);
				c.amr.wing_refinement_distance = number(*it, "wing_refinement_distance", c.amr.wing_refinement_distance);
				c.amr.surface_refinement_distance = number(*it, "surface_refinement_distance", c.amr.surface_refinement_distance);
				c.amr.wake_length = number(*it, "wake_length", c.amr.wake_length); c.amr.wake_radius = number(*it, "wake_radius", c.amr.wake_radius);
				c.amr.complex_subdivisions = integer(*it, "complex_subdivisions", c.amr.complex_subdivisions);
				c.amr.min_volume_fraction = number(*it, "min_volume_fraction", c.amr.min_volume_fraction);
			}
			if (auto it = root.find("solver"); it != root.end() && it->is_object())
			{
				c.solver.cfl = number(*it, "cfl", c.solver.cfl); c.solver.smagorinsky_cs = number(*it, "smagorinsky_cs", c.solver.smagorinsky_cs);
				c.solver.projection_tolerance = number(*it, "projection_tolerance", c.solver.projection_tolerance); c.solver.projection_max_iterations = integer(*it, "projection_max_iterations", c.solver.projection_max_iterations);
			}
			if (auto it = root.find("reference"); it != root.end() && it->is_object())
			{
				c.reference.area = number(*it, "area", 0); c.reference.length = number(*it, "length", 0);
				if (auto o = it->find("moment_origin"); o != it->end() && o->is_array() && o->size() == 3)
					c.reference.moment_origin = {(*o)[0].get<double>(), (*o)[1].get<double>(), (*o)[2].get<double>()};
			}
			if (!(c.tessellation_deflection_mm > 0 && c.freestream.speed >= 0 && c.freestream.rho > 0 && c.freestream.nu > 0 &&
				c.domain.upstream_margin >= 0 && c.domain.downstream_margin >= 0 && c.domain.lateral_margin >= 0 && c.domain.vertical_margin >= 0 &&
				c.amr.base_cell_size > 0 && c.amr.max_levels >= 1 && c.amr.max_levels <= 10 && c.amr.brick_size >= 4 && c.amr.brick_size <= 128 && c.amr.ghost_cells >= 1 && c.amr.ghost_cells <= 4 &&
				c.amr.wing_refinement_distance >= 0 && c.amr.surface_refinement_distance >= 0 && c.amr.wake_length >= 0 && c.amr.wake_radius >= 0 && c.amr.complex_subdivisions >= 0 && c.amr.complex_subdivisions <= 16 && c.amr.min_volume_fraction > 0 && c.amr.min_volume_fraction < 0.5 &&
				c.solver.cfl > 0 && c.solver.smagorinsky_cs >= 0 && c.solver.projection_tolerance > 0 && c.solver.projection_max_iterations > 0 && c.reference.area >= 0 && c.reference.length >= 0))
			{
				if (error) *error = "invalid non-positive or out-of-range paraglider configuration value"; return false;
			}
			out = c; if (error) error->clear(); return true;
		}
		catch (const std::exception& e) { if (error) *error = std::string("paraglider config parse failed: ") + e.what(); return false; }
	}

	bool save_paraglider_config(const std::string& path,const ParagliderConfig& c,std::string* error)
	{
		try
		{
			nlohmann::json root;
			root["step"]={{"path",c.step_path},{"tessellation_deflection_mm",c.tessellation_deflection_mm}};
			root["placement"]={{"tx",c.placement.tx},{"ty",c.placement.ty},{"tz",c.placement.tz},
				{"matrix",{c.placement.m[0],c.placement.m[1],c.placement.m[2],c.placement.m[3],c.placement.m[4],c.placement.m[5],c.placement.m[6],c.placement.m[7],c.placement.m[8]}}};
			root["freestream"]={{"speed",c.freestream.speed},{"rho",c.freestream.rho},{"nu",c.freestream.nu}};
			root["domain"]={{"upstream_margin",c.domain.upstream_margin},{"downstream_margin",c.domain.downstream_margin},{"lateral_margin",c.domain.lateral_margin},{"vertical_margin",c.domain.vertical_margin}};
			root["amr"]={{"base_cell_size",c.amr.base_cell_size},{"max_levels",c.amr.max_levels},{"brick_size",c.amr.brick_size},{"ghost_cells",c.amr.ghost_cells},{"wing_refinement_distance",c.amr.wing_refinement_distance},{"surface_refinement_distance",c.amr.surface_refinement_distance},{"wake_length",c.amr.wake_length},{"wake_radius",c.amr.wake_radius},{"complex_subdivisions",c.amr.complex_subdivisions},{"min_volume_fraction",c.amr.min_volume_fraction}};
			root["solver"]={{"cfl",c.solver.cfl},{"smagorinsky_cs",c.solver.smagorinsky_cs},{"projection_tolerance",c.solver.projection_tolerance},{"projection_max_iterations",c.solver.projection_max_iterations}};
			root["reference"]={{"area",c.reference.area},{"length",c.reference.length},{"moment_origin",{c.reference.moment_origin.x,c.reference.moment_origin.y,c.reference.moment_origin.z}}};
			std::ofstream out(path);if(!out){if(error)*error="cannot open paraglider config '"+path+"' for writing";return false;}out<<root.dump(2)<<'\n';if(!out){if(error)*error="failed while writing paraglider config '"+path+"'";return false;}if(error)error->clear();return true;
		}
		catch(const std::exception& e){if(error)*error=std::string("paraglider config write failed: ")+e.what();return false;}
	}

	Aabb3d automatic_flow_domain(const TriMesh& m, const DomainConfig& d)
	{
		Aabb3d b;
		if (m.empty()) return b;
		b.lo = {m.bbox_min[0] - d.upstream_margin, m.bbox_min[1] - d.lateral_margin, m.bbox_min[2] - d.vertical_margin};
		b.hi = {m.bbox_max[0] + d.downstream_margin, m.bbox_max[1] + d.lateral_margin, m.bbox_max[2] + d.vertical_margin};
		return b;
	}
}
