#include "core/fluid/external_aero_core.h"
#include "core/geometry/step_import.h"
#include "core/geometry/triangle_bvh.h"
#include "core/paraglider_config.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace paracfd::core;

namespace
{
	double elapsed_ms(std::chrono::steady_clock::time_point begin,std::chrono::steady_clock::time_point end)
	{
		return std::chrono::duration<double,std::milli>(end-begin).count();
	}

	void report(int step,const ExternalAeroCore& core,const ExternalAeroStepStats& stats,double rho)
	{
		const AerodynamicLoads loads=core.pressure_loads();const ExternalAeroConservationStats conservation=core.conservation_stats();
		std::printf("[paraglider-case] step=%d t=%.6f s dt=%.3e max|u_i|=%.5g [regular=%.5g special=%.5g CFL-source=%.5g] CFL=%.3f Fp=[%.8g %.8g %.8g] N div[max=%.3e rmsV=%.3e] flux-error[max=%.3e abs=%.3e net=%.3e] PCG=%d/%.3e step=%.2f ms projection=%.2f ms EB=%.2f ms\n",step,core.physical_time(),stats.dt,stats.max_abs_velocity,stats.max_abs_regular_velocity,stats.max_abs_special_velocity,stats.cfl_velocity,stats.effective_cfl,loads.pressure_force.x,loads.pressure_force.y,loads.pressure_force.z,conservation.max_abs_divergence,conservation.volume_weighted_rms_divergence,conservation.max_integrated_flux_error,conservation.absolute_integrated_flux_error,conservation.net_integrated_flux_error,stats.pressure.iterations,stats.pressure.relative_residual,stats.gpu_step_ms,stats.projection_ms,stats.embedded_transport_ms);
		CompositeAmrFluxes fluxes;core.download_special_fluxes(fluxes);double coarse_fine_max=0,embedded_max=0;std::size_t embedded_edge=0;for(double value:fluxes.coarse_fine_velocity)coarse_fine_max=std::max(coarse_fine_max,std::abs(value));for(std::size_t edge=0;edge<fluxes.embedded_velocity.size();++edge)if(std::abs(fluxes.embedded_velocity[edge])>embedded_max){embedded_max=std::abs(fluxes.embedded_velocity[edge]);embedded_edge=edge;}const CompositeAmrPressureSystem& system=core.pressure_system();if(embedded_edge<system.embedded.size()){std::vector<double> pressure;core.download_pressure(pressure);const CoarseFinePressureConnection& connection=system.embedded[embedded_edge];const int lower=connection.direction>0?connection.coarse_dof:connection.fine_dof,upper=connection.direction>0?connection.fine_dof:connection.coarse_dof;const double delta_p=pressure[upper]-pressure[lower],correction=-(stats.dt/rho)*delta_p/connection.centre_distance;std::printf("[paraglider-case]   special maxima: coarse/fine=%.5g EB=%.5g; EB edge=%zu axis=%d A=%.3e d=%.3e Va=%.3e Vb=%.3e dp=%.5g Pa du_p=%.5g m/s centroid=[%.5g %.5g %.5g]\n",coarse_fine_max,embedded_max,embedded_edge,static_cast<int>(connection.axis),connection.open_area,connection.centre_distance,system.volume[connection.coarse_dof],system.volume[connection.fine_dof],delta_p,correction,connection.face_centroid.x,connection.face_centroid.y,connection.face_centroid.z);}
	}
}

int main(int argc,char** argv)
{
	std::string config_path="configs/planb_parakite.json";int steps=50,sample_every=5;double projection_tolerance=0,min_volume_fraction=0;
	for(int i=1;i<argc;++i)
	{
		const std::string argument=argv[i];
		if(argument=="--config"&&i+1<argc)config_path=argv[++i];
		else if(argument=="--steps"&&i+1<argc)steps=std::atoi(argv[++i]);
		else if(argument=="--sample-every"&&i+1<argc)sample_every=std::atoi(argv[++i]);
		else if(argument=="--projection-tolerance"&&i+1<argc)projection_tolerance=std::atof(argv[++i]);
		else if(argument=="--min-volume-fraction"&&i+1<argc)min_volume_fraction=std::atof(argv[++i]);
		else{std::fprintf(stderr,"usage: paraglider_case_probe [--config file] [--steps N] [--sample-every N] [--projection-tolerance value] [--min-volume-fraction value]\n");return 2;}
	}
	if(steps<0||sample_every<1){std::fprintf(stderr,"steps must be nonnegative and sample-every must be positive\n");return 2;}
	ParagliderConfig config;std::string error;if(!load_paraglider_config(config_path,config,&error)){std::fprintf(stderr,"[paraglider-case] %s\n",error.c_str());return 2;}if(projection_tolerance>0)config.solver.projection_tolerance=projection_tolerance;if(min_volume_fraction>0)config.amr.min_volume_fraction=min_volume_fraction;
	const auto load_begin=std::chrono::steady_clock::now();TriMesh source=load_step_mesh(config.step_path,config.tessellation_deflection_mm,&error);if(source.empty()){std::fprintf(stderr,"[paraglider-case] STEP load failed: %s\n",error.c_str());return 2;}TriMesh wing=placed_mesh(source,config.placement);TriangleBvh bvh(wing);const auto load_end=std::chrono::steady_clock::now();
	try
	{
		const auto build_begin=std::chrono::steady_clock::now();ExternalAeroCore core(wing,bvh,config);const auto build_end=std::chrono::steady_clock::now();ExternalAeroStepStats stats=core.initialize();
		std::printf("[paraglider-case] triangles=%zu bricks=%zu DOFs=%d EB-edges=%zu patches=%zu load=%.2f ms preprocess=%.2f ms GPU=%.2f MiB tolerance=%.3e\n",wing.triangle_count(),core.hierarchy().active_brick_count(),core.pressure_system().storage_size,core.pressure_system().embedded.size(),core.pressure_system().surface_patches.size(),elapsed_ms(load_begin,load_end),elapsed_ms(build_begin,build_end),core.gpu_bytes()/(1024.0*1024.0),config.solver.projection_tolerance);
		if(!stats.pressure.converged){std::fprintf(stderr,"[paraglider-case] initialization did not converge: iterations=%d residual=%.3e\n",stats.pressure.iterations,stats.pressure.relative_residual);return 3;}report(0,core,stats,config.freestream.rho);
		for(int step=1;step<=steps;++step)
		{
			stats=core.step();if(!stats.pressure.converged){std::fprintf(stderr,"[paraglider-case] step %d did not converge: iterations=%d residual=%.3e\n",step,stats.pressure.iterations,stats.pressure.relative_residual);return 4;}
			if(step==1||step==steps||step%sample_every==0)report(step,core,stats,config.freestream.rho);
		}
	}
	catch(const std::exception& exception){std::fprintf(stderr,"[paraglider-case] failed: %s\n",exception.what());return 5;}
	return 0;
}
