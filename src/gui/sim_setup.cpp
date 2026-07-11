// sim_setup.cpp — see sim_setup.h. Parses the viewer JSON directly (not via the strict
// core Config loader, so viewer-only keys are allowed) and constructs a ChannelFluidCore
// mirroring the M2 cylinder benchmark plumbing (channel_core + channel_mask).
#include "gui/sim_setup.h"

#include "core/fluid/channel_mask.h"
#include "core/geometry/step_import.h"
#include "core/geometry/voxelize.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace scour::gui
{
	using namespace scour::core;

	namespace
	{
		double jd(const nlohmann::json& j, const char* k, double d) { return j.contains(k) ? j.at(k).get<double>() : d; }
		int ji(const nlohmann::json& j, const char* k, int d) { return j.contains(k) ? j.at(k).get<int>() : d; }
		bool jb(const nlohmann::json& j, const char* k, bool d) { return j.contains(k) ? j.at(k).get<bool>() : d; }

		// If an override is active, replace the domain + voxel size (ONLY the grid) in place.
		void apply_override(const GridOverride* ov, double& Lx, double& Ly, double& Lz, double& h)
		{
			if (ov && ov->active)
			{
				Lx = ov->Lx; Ly = ov->Ly; Lz = ov->Lz; h = ov->h;
			}
		}

		// Assemble the morphodynamic engine + flow core + recipe/scen from a ready SeabedParams and a
		// rigid structure mask (shared by build_seabed_sim and build_seabed_from_structure). The engine
		// fills sand where z < sand_depth OR structure; the core is masked with that + a bed-inlet mask.
		// Caller sets info.name, scen.spinup, scen.structure_mesh/place and any log line.
		std::unique_ptr<ChannelFluidCore> assemble_seabed(const MacGrid& g, const SeabedParams& sp_in, double U,
			double nu_fluid, double proj_tol, int proj_max_iter, double cfl, bool inlet_loglaw,
			const std::vector<unsigned char>& structure, const std::string& source_config,
			SimRecipe& recipe, SeabedScenario& scen, bool build_core = true)
		{
			SeabedParams sp = sp_in; // local copy so the engine + recipe carry the inlet current …
			sp.U_inlet = U;          // … from which SeabedMorpho derives the open-sea equilibrium inflow u*
			bool have_structure = false;
			for (unsigned char c : structure) if (c) { have_structure = true; break; }

			scen.engine = std::make_unique<SeabedMorpho>(g, sp, have_structure ? structure : std::vector<unsigned char>{});
			scen.bed_z0 = sp.sand_depth;
			scen.params = sp;
			std::vector<unsigned char> init_solid = scen.engine->initial_flow_solid();

			ChannelBC bc;
			bc.inlet_mode = INLET_UNIFORM; bc.U_inlet = U; bc.Uc = U;
			bc.solid_mode = SOLID_NOSLIP;
			bc.ymin = bc.ymax = WALL_FREESLIP; bc.zmin = bc.zmax = WALL_FREESLIP;
			// Log-law boundary-layer inlet (default for sediment runs): the profile rises from ~0 at the
			// sand top (bed_datum = sand_depth) so the leading-edge bed shear is realistic, not the top-hat
			// full-U spike. u* is flux-matched to U over the fluid depth. z0 = grain roughness d50/12.
			if (inlet_loglaw)
			{
				bc.inlet_mode = INLET_LOGLAW;
				bc.z0 = (sp.d50 > 0.0) ? sp.d50 / 12.0 : bc.z0;
				bc.bed_datum = sp.sand_depth;
				bc.kappa = 0.40;
				double H = g.nz * g.h - bc.bed_datum;
				bc.ustar = loglaw_ustar_for_U(U, H, bc.z0, bc.kappa);
				double zp = 1.5 * g.h; // representative near-bed probe height
				double u_nb = (bc.ustar / bc.kappa) * std::log(std::max(zp / bc.z0, 2.718281828459045));
				std::fprintf(stderr, "[seabed] inlet = LOG-LAW boundary layer (leading-edge fix): u*=%.4f m/s, "
					"z0=%.2e m, bed_datum=%.2f m; near-bed u(%.3f m)=%.3f m/s vs top-hat %.2f m/s\n",
					bc.ustar, bc.z0, bc.bed_datum, zp, u_nb, U);
			}

			ChannelParams pr;
			pr.rho = sp.rho; pr.nu = nu_fluid; pr.Cs = sp.Cs; pr.cfl = cfl; pr.safety = 0.9;
			pr.proj_tol = proj_tol; pr.proj_max_iter = proj_max_iter;
			pr.fixed_dt = 0.0; pr.advect_band = 1;

			// Skip the flow-core construction (build_core=false) when the caller reuses an EXISTING,
			// already-developed core and only needs the engine + masks + BCs (the GUI's live "Add sand",
			// which preserves the running flow via ChannelFluidCore::update_solid instead of resetting).
			std::unique_ptr<ChannelFluidCore> core;
			if (build_core)
			{
				core = std::make_unique<ChannelFluidCore>(g, bc, pr, init_solid);
				core->set_bed_inlet_mask(true); // the erodible bed reaches the inlet plane
				core->init_uniform(U, 0.0);
			}

			recipe.grid = g; recipe.bc = bc; recipe.pr = pr;
			recipe.base_solid = init_solid; recipe.base_solid_mode = SOLID_NOSLIP;
			recipe.init_u = U; recipe.init_v_blip = 0.0;
			recipe.source_config = source_config;
			recipe.is_scenario = true;
			SimInfo& info = recipe.info;
			info.nx = g.nx; info.ny = g.ny; info.nz = g.nz; info.h = g.h;
			info.Lx = g.nx * g.h; info.Ly = g.ny * g.h; info.Lz = g.nz * g.h;
			info.U = U; info.rho = sp.rho; info.nu = nu_fluid; info.cylinder = false;

			scen.active = true;
			scen.has_structure = have_structure;
			return core;
		}
	}

	void grid_dims_for(double Lx, double Ly, double Lz, double h, int& nx, int& ny, int& nz)
	{
		// Same derivation both sim builders use — floor to whole cells with the minimum clamps — so a
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
		double Re = jd(j, "Re", 150.0);
		double v_blip = jd(j, "v_blip", 0.05);

		MacGrid g; g.h = h;
		grid_dims_for(Lx, Ly, Lz, h, g.nx, g.ny, g.nz);

		// --- Obstacle (procedural cylinder, on by default for visible dynamics) ----
		nlohmann::json jc = j.contains("cylinder") ? j.at("cylinder") : nlohmann::json::object();
		bool cyl = jb(jc, "enabled", true);
		double D = jd(jc, "diameter_m", 0.10 * Ly); // default: 10% of the span
		double xc = jd(jc, "xc_m", 0.25 * Lx);       // 1/4 downstream
		double yc = jd(jc, "yc_m", 0.5 * Ly);
		double R = 0.5 * D;

		std::vector<unsigned char> solid(g.p_count(), 0);
		if (cyl)
			build_cylinder_mask(g, xc, yc, R, solid);

		// --- Viscosity: artificial nu = U*D/Re (shedding wake) ---------------------
		double nu = (cyl && Re > 0.0) ? (U * D / Re) : jd(j, "nu", 1.0e-3);
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
		recipe.bc = bc;
		recipe.pr = pr;
		recipe.base_solid = solid;
		recipe.base_solid_mode = SOLID_NOSLIP;
		recipe.init_u = U;
		recipe.init_v_blip = v_blip * U;
		recipe.source_config = config_path; // so the GUI can rebuild this exact sim at another grid
		recipe.is_scenario = false;

		SimInfo& info = recipe.info;
		info.nx = g.nx; info.ny = g.ny; info.nz = g.nz; info.h = h;
		info.Lx = g.nx * h; info.Ly = g.ny * h; info.Lz = g.nz * h;
		info.U = U; info.rho = rho; info.nu = nu; info.cylinder = cyl;
		if (j.contains("name")) info.name = j.at("name").get<std::string>();

		return make_core(recipe, recipe.base_solid, recipe.base_solid_mode);
	}

	std::unique_ptr<ChannelFluidCore> build_seabed_sim(const std::string& config_path,
		SimRecipe& recipe, SeabedScenario& scen, std::string& warn, const GridOverride* ov,
		const ModelPlacement* structure_place)
	{
		nlohmann::json j = nlohmann::json::object();
		try
		{
			std::ifstream in(config_path);
			if (!in) throw std::runtime_error("cannot open");
			std::ostringstream ss; ss << in.rdbuf();
			j = nlohmann::json::parse(ss.str());
		}
		catch (const std::exception& e)
		{
			warn = std::string("seabed config '") + config_path + "' unreadable (" + e.what() + ")";
			return build_sim(config_path, recipe, warn, ov); // falls back to a fluid-only viewer
		}

		std::fprintf(stderr, "[seabed] PRE-CALIBRATION DEMO — constants nominal, results qualitative\n");

		// --- Domain & grid --------------------------------------------------------
		double Lx = jd(j, "domain_x", 10.0), Ly = jd(j, "domain_y", 10.0), Lz = jd(j, "domain_z", 5.0);
		double h = jd(j, "voxel_h", 0.05);
		apply_override(ov, Lx, Ly, Lz, h); // GUI Apply: change ONLY the grid; structure re-voxelizes at new h
		double U = jd(j, "U", 0.5);
		if (ov && ov->active && ov->U > 0.0) U = ov->U; // GUI "Input speed": test at a different current
		double rho = jd(j, "rho", 1027.0);

		MacGrid g; g.h = h;
		grid_dims_for(Lx, Ly, Lz, h, g.nx, g.ny, g.nz);

		// --- Seabed parameters (NOMINAL / uncalibrated) ---------------------------
		SeabedParams sp;
		sp.sand_depth = jd(j, "sand_depth", 1.0);
		sp.d50 = jd(j, "d50", 0.2e-3);
		sp.rho = rho; sp.rho_s = jd(j, "rho_s", 2650.0); sp.nu = jd(j, "nu_sed", 1.36e-6);
		sp.Cs = jd(j, "Cs", 0.11);
		sp.alpha = jd(j, "alpha", 0.00033);
		sp.morfac = jd(j, "morfac", 5.0);
		sp.bedload_formula = ji(j, "bedload_formula", 1);
		sp.diffusion_on = jb(j, "suspended_diffusion", true) ? 1 : 0;
		// Suspended-sediment x-boundary (RESEARCH §6): open-sea equilibrium inflow by DEFAULT ("open"),
		// "recycle" = flux-matched recirculating flume, "closed" = the legacy zero-flux box. Key `sediment_bc`.
		{
			std::string sbc = j.contains("sediment_bc") ? j.at("sediment_bc").get<std::string>() : std::string("open");
			sp.sed_bc = (sbc == "closed") ? SED_BC_CLOSED : (sbc == "recycle") ? SED_BC_RECYCLE : SED_BC_OPEN;
		}
		scen.spinup = ji(j, "spinup_steps", 200);
		scen.bed_z0 = sp.sand_depth;

		// --- Structure: voxelize the STEP model, seated on the bed ----------------
		std::vector<unsigned char> structure((size_t)g.p_count(), 0);
		std::string step_path = j.contains("structure_step") ? j.at("structure_step").get<std::string>() : std::string();
		bool have_structure = false;
		ModelPlacement place;
		if (!step_path.empty() && jb(j, "structure_enabled", true))
		{
			std::string err;
			TriMesh mesh = load_step_mesh(step_path, 0.1, &err);
			if (mesh.empty())
				warn = std::string("structure STEP '") + step_path + "' load failed (" + err + "); bare seabed";
			else
			{
				// Gizmo placement override (GUI Apply) already includes the sand-depth seating; otherwise
				// use the default centre-on-bed placement and seat it on top of the sand reservoir.
				if (structure_place)
					place = *structure_place;
				else
				{
					place = place_model_on_bed(mesh, Lx, Ly);
					place.tz += sp.sand_depth; // seat the structure ON TOP of the sand reservoir
				}
				double mv = 0, vv = 0; int thin = 0;
				structure = voxelize_mesh(mesh, g, place, &mv, &vv, nullptr, &thin);
				long long ns = 0; for (unsigned char c : structure) ns += (c != 0);
				std::fprintf(stderr, "[seabed] structure '%s': %lld solid cells on the bed (z=%.3f m)\n",
					step_path.c_str(), ns, place.tz);
				have_structure = ns > 0;
				scen.structure_mesh = std::move(mesh);
				scen.structure_place = place;
			}
		}

		// --- Morphodynamic engine + flow core (shared assembly) -------------------
		double nu_fluid = jd(j, "nu_fluid", 5.0e-3); // artificial FLOW viscosity (coarse-grid stability)
		double proj_tol = jd(j, "proj_tol", 1e-4);
		int proj_max_iter = ji(j, "proj_max_iter", 60);
		double cfl = jd(j, "cfl", 0.7);
		// Inlet profile: "loglaw" (boundary-layer, DEFAULT for sediment — avoids the top-hat leading-edge
		// scour) or "uniform" (top-hat). Config key `inlet_profile`.
		std::string ip = j.contains("inlet_profile") ? j.at("inlet_profile").get<std::string>() : std::string("loglaw");
		bool inlet_loglaw = (ip != "uniform");
		auto core = assemble_seabed(g, sp, U, nu_fluid, proj_tol, proj_max_iter, cfl, inlet_loglaw,
			have_structure ? structure : std::vector<unsigned char>{}, config_path, recipe, scen);
		recipe.info.name = j.contains("name") ? j.at("name").get<std::string>() : std::string("seabed");
		std::fprintf(stderr, "[seabed] %dx%dx%d cells, h=%.3f m, U=%.2f m/s, sand=%.2f m, MORFAC=%.1f, "
			"alpha=%.5f, Cs=%.2f, nu_fluid=%.1e, spinup=%d — UNCALIBRATED\n",
			g.nx, g.ny, g.nz, h, U, sp.sand_depth, sp.morfac, sp.alpha, sp.Cs, nu_fluid, scen.spinup);
		return core;
	}

	std::unique_ptr<ChannelFluidCore> build_seabed_from_structure(
		const MacGrid& g, double U, double sand_depth,
		const std::vector<unsigned char>& structure, const SeabedParams* base,
		const std::string& source_config, int spinup, double nu_fluid,
		SimRecipe& recipe, SeabedScenario& scen, std::string& warn, bool build_core)
	{
		std::fprintf(stderr, "[seabed] PRE-CALIBRATION DEMO — constants nominal, results qualitative\n");
		if ((int)structure.size() != g.p_count())
		{
			warn = "add-sand: structure mask size mismatch";
			return nullptr;
		}
		SeabedParams sp = base ? *base : SeabedParams{}; // struct defaults are the nominal demo values
		if (!base) sp.sed_bc = SED_BC_OPEN; // Add-sand on a plain fluid viewer ⇒ open-sea by default
		sp.sand_depth = sand_depth < 0.0 ? 0.0 : sand_depth;

		scen.spinup = spinup;
		auto core = assemble_seabed(g, sp, U, nu_fluid, /*proj_tol*/1e-4, /*proj_max_iter*/60, /*cfl*/0.7,
			/*inlet_loglaw*/true, structure, source_config, recipe, scen, build_core);
		recipe.info.name = "seabed (added sand)";

		long long ns = 0; for (unsigned char c : structure) ns += (c != 0);
		std::fprintf(stderr, "[seabed] added %.3f m sand on %dx%dx%d @ h=%.3f m, U=%.2f m/s, %lld rigid cells, "
			"MORFAC=%.1f, alpha=%.5f, Cs=%.2f — UNCALIBRATED\n",
			sp.sand_depth, g.nx, g.ny, g.nz, g.h, U, ns, sp.morfac, sp.alpha, sp.Cs);
		return core;
	}
}
