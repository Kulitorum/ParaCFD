// sim_setup.cpp — see sim_setup.h. Parses the viewer JSON directly (not via the strict
// core Config loader, so viewer-only keys are allowed) and constructs a ChannelFluidCore
// using the channel_core fluid solver.
#include "gui/sim_setup.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace paracfd::gui
{
	using namespace paracfd::core;

	namespace
	{
		double jd(const nlohmann::json& j, const char* k, double d) { return j.contains(k) ? j.at(k).get<double>() : d; }
		int ji(const nlohmann::json& j, const char* k, int d) { return j.contains(k) ? j.at(k).get<int>() : d; }

		// If an override is active, replace the domain + voxel size (ONLY the grid) in place.
		void apply_override(const GridOverride* ov, double& Lx, double& Ly, double& Lz, double& h)
		{
			if (ov && ov->active)
			{
				Lx = ov->Lx; Ly = ov->Ly; Lz = ov->Lz; h = ov->h;
			}
		}
	}

	void grid_dims_for(double Lx, double Ly, double Lz, double h, int& nx, int& ny, int& nz)
	{
		// Same derivation the sim builder uses — floor to whole cells with the minimum clamps — so a
		// GUI readout computed here matches the grid Apply actually builds.
		if (h <= 0.0) h = 0.05;
		nx = std::max(4, (int)std::floor(Lx / h));
		ny = std::max(4, (int)std::floor(Ly / h));
		nz = std::max(3, (int)std::floor(Lz / h));
	}

	std::unique_ptr<ChannelFluidCore> make_core(const SimRecipe& r, const std::vector<unsigned char>& solid, int solid_mode)
	{
		ChannelBC bc = r.bc;
		bc.solid_mode = solid_mode;
		auto core = std::make_unique<ChannelFluidCore>(r.grid, bc, r.pr, solid);
		core->init_uniform(r.init_u, r.init_v_blip);
		return core;
	}

	std::unique_ptr<ChannelFluidCore> build_sim(const std::string& config_path, SimRecipe& recipe, std::string& warn, const GridOverride* ov)
	{
		nlohmann::json j = nlohmann::json::object();
		if (!config_path.empty())
		{
			try
			{
				std::ifstream in(config_path);
				if (!in) throw std::runtime_error("cannot open");
				std::ostringstream ss; ss << in.rdbuf();
				j = nlohmann::json::parse(ss.str());
			}
			catch (const std::exception& e)
			{
				warn = std::string("config '") + config_path + "' unreadable (" + e.what() + "); using defaults";
				j = nlohmann::json::object();
			}
		}

		// --- Domain & grid --------------------------------------------------------
		double Lx = jd(j, "domain_x", 10.0), Ly = jd(j, "domain_y", 10.0), Lz = jd(j, "domain_z", 5.0);
		double h = jd(j, "voxel_h", 0.05);
		apply_override(ov, Lx, Ly, Lz, h); // GUI Apply: change ONLY the grid, keep the physics below
		double U = jd(j, "U", 1.0);
		if (ov && ov->active && ov->U > 0.0) U = ov->U; // GUI "Input speed": test at a different current
		double rho = jd(j, "rho", 1.0);
		double Cs = jd(j, "Cs", 0.0);
		double v_blip = jd(j, "v_blip", 0.05);

		// --- Fine-core graded grid (graded-structured-grid change): a uniform-h_fine core around the
		// building with a geometrically graded coarse far field, so a rounded corner is resolvable without
		// a fine grid over the whole domain. Opt-in via a "fine_core" object; a missing/disabled/degenerate
		// spec (or h_fine >= the coarse voxel_h) keeps the uniform grid — BYTE-IDENTICAL to before. The box
		// is absolute domain metres and must enclose the placed building (windloads asserts bbox ⊆ core).
		FineCoreSpec fc; fc.Lx = Lx; fc.Ly = Ly; fc.Lz = Lz; fc.h_fine = h; fc.growth = 1.15;
		bool graded = false;
		if (ov && ov->active && ov->set_fine_core)
		{
			// GUI dock is authoritative for the fine core (the config's fine_core is ignored on Apply).
			graded = ov->graded;
			fc = ov->fine_core; fc.Lx = Lx; fc.Ly = Ly; fc.Lz = Lz;
			if (!(fc.growth > 1.0)) fc.growth = 1.15;
			if (!(fc.h_fine > 0.0) || fc.h_fine >= h) graded = false;
		}
		else if (j.contains("fine_core") && j.at("fine_core").is_object())
		{
			const auto& o = j.at("fine_core");
			graded = o.value("enabled", true);
			fc.x0 = jd(o, "x0", 0.0); fc.x1 = jd(o, "x1", 0.0);
			fc.y0 = jd(o, "y0", 0.0); fc.y1 = jd(o, "y1", 0.0);
			fc.z0 = jd(o, "z0", 0.0); fc.z1 = jd(o, "z1", 0.0);
			fc.h_fine = jd(o, "h_fine", h);
			fc.growth = jd(o, "growth", 1.15);
			if (!(fc.growth > 1.0)) fc.growth = 1.15;
			if (!(fc.h_fine > 0.0) || fc.h_fine >= h) graded = false; // no core finer than the coarse grid ⇒ uniform
		}

		MacGrid g;
		std::shared_ptr<GridMetrics> metrics;
		if (graded)
		{
			metrics = std::make_shared<GridMetrics>(GridMetrics::generate(fc));
			g = metrics->device_view(); // device metric pointers for the GPU core (one lazy H2D upload)
		}
		else
		{
			g.h = h;
			grid_dims_for(Lx, Ly, Lz, h, g.nx, g.ny, g.nz);
		}

		// --- No config obstacle: the flow starts as an empty channel. A building is injected
		// later via the centerline Build workflow (or a loaded STEP model). ------------------
		std::vector<unsigned char> solid(g.p_count(), 0);

		// --- Viscosity from the config (moderate default for a visible developing wake) ------
		double nu = jd(j, "nu", 1.0e-3);
		if (nu <= 0.0) nu = 1.0e-3;

		ChannelBC bc;
		bc.inlet_mode = INLET_UNIFORM; bc.U_inlet = U; bc.Uc = U;
		bc.solid_mode = SOLID_NOSLIP;
		bc.ymin = bc.ymax = WALL_FREESLIP;
		bc.zmin = bc.zmax = WALL_FREESLIP;

		ChannelParams pr;
		pr.rho = rho; pr.nu = nu; pr.Cs = Cs; pr.cfl = 1.0; pr.safety = 0.9;
		pr.proj_tol = jd(j, "proj_tol", 1e-3);
		pr.proj_max_iter = ji(j, "proj_max_iter", 40);
		pr.fixed_dt = 0.0; // adaptive dt (viewer needn't lock the sample rate)
		pr.advect_band = 1;

		// --- Fill the recipe (immutable ingredients for re-injection rebuilds) -----
		recipe.grid = g;
		recipe.metrics = metrics; // null ⇒ uniform; shared so every recipe copy keeps the metric arrays alive
		recipe.fine_core = fc;
		recipe.graded = graded;
		recipe.bc = bc;
		recipe.pr = pr;
		recipe.base_solid = solid;
		recipe.base_solid_mode = SOLID_NOSLIP;
		recipe.init_u = U;
		recipe.init_v_blip = v_blip * U;
		recipe.source_config = config_path; // so the GUI can rebuild this exact sim at another grid

		SimInfo& info = recipe.info;
		info.nx = g.nx; info.ny = g.ny; info.nz = g.nz;
		info.h = graded ? fc.h_fine : h; // finest spacing (h_fine on a graded grid)
		info.coarse_h = h;               // the "Voxel/cell size" control value (nominal coarse voxel; = h_fine's parent)
		if (graded) { MacGrid hv = metrics->host_view(); info.Lx = hv.Lx(); info.Ly = hv.Ly(); info.Lz = hv.Lz(); }
		else { info.Lx = g.nx * h; info.Ly = g.ny * h; info.Lz = g.nz * h; }
		info.U = U; info.rho = rho; info.nu = nu;
		if (j.contains("name")) info.name = j.at("name").get<std::string>();

		return make_core(recipe, recipe.base_solid, recipe.base_solid_mode);
	}
}
