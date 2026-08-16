#include "core/fluid/external_aero_core.h"

#include <cmath>
#include <cstdio>
#include <string>

using namespace paracfd::core;

namespace
{
	int failures=0;
	void check(bool condition,const char* name){std::printf("[paraglider-flow] %-58s %s\n",name,condition?"PASS":"FAIL");if(!condition)++failures;}
	TriMesh quad(Vec3d a,Vec3d b,Vec3d c,Vec3d d)
	{
		TriMesh mesh;for(Vec3d p:{a,b,c,d}){mesh.positions.push_back(static_cast<float>(p.x));mesh.positions.push_back(static_cast<float>(p.y));mesh.positions.push_back(static_cast<float>(p.z));}mesh.indices={0,1,2,0,2,3};mesh.source_face_ids={0,0};return mesh;
	}
	TriMesh rotate_y(const TriMesh& source,double angle)
	{
		TriMesh mesh=source;const double c=std::cos(angle),s=std::sin(angle);for(std::size_t vertex=0;vertex<mesh.vertex_count();++vertex){const double x=source.positions[3*vertex],z=source.positions[3*vertex+2];mesh.positions[3*vertex]=static_cast<float>(c*x+s*z);mesh.positions[3*vertex+2]=static_cast<float>(-s*x+c*z);}return mesh;
	}
	void append_quad(TriMesh& mesh,Vec3d a,Vec3d b,Vec3d c,Vec3d d,std::uint32_t face)
	{
		const std::uint32_t first=static_cast<std::uint32_t>(mesh.vertex_count());for(Vec3d p:{a,b,c,d}){mesh.positions.push_back(static_cast<float>(p.x));mesh.positions.push_back(static_cast<float>(p.y));mesh.positions.push_back(static_cast<float>(p.z));}mesh.indices.insert(mesh.indices.end(),{first,first+1,first+2,first,first+2,first+3});mesh.source_face_ids.insert(mesh.source_face_ids.end(),{face,face});
	}
	TriMesh cavity_box(bool upstream_open)
	{
		TriMesh mesh;std::uint32_t face=0;if(!upstream_open)append_quad(mesh,{-0.5,-0.5,-0.5},{-0.5,-0.5,0.5},{-0.5,0.5,0.5},{-0.5,0.5,-0.5},face++);append_quad(mesh,{0.5,-0.5,-0.5},{0.5,0.5,-0.5},{0.5,0.5,0.5},{0.5,-0.5,0.5},face++);append_quad(mesh,{-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{0.5,-0.5,0.5},{-0.5,-0.5,0.5},face++);append_quad(mesh,{-0.5,0.5,-0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5},{0.5,0.5,-0.5},face++);append_quad(mesh,{-0.5,-0.5,-0.5},{-0.5,0.5,-0.5},{0.5,0.5,-0.5},{0.5,-0.5,-0.5},face++);append_quad(mesh,{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{0.5,0.5,0.5},{-0.5,0.5,0.5},face++);return mesh;
	}
	ParagliderConfig flow_config()
	{
		ParagliderConfig config;config.domain={1.5,2.5,1.5,1.5};config.amr.base_cell_size=0.25;config.amr.max_levels=1;config.amr.brick_size=8;config.amr.ghost_cells=1;config.amr.complex_subdivisions=8;config.solver.cfl=0.25;config.solver.smagorinsky_cs=0.1;config.solver.projection_tolerance=sizeof(Real)==4?1e-5:1e-10;config.solver.projection_max_iterations=1000;config.freestream.speed=1;config.freestream.rho=1;config.freestream.nu=1.5e-5;config.reference.area=1;return config;
	}
	std::size_t pressure_gauges(const TriMesh& mesh)
	{
		TriangleBvh bvh(mesh);const ParagliderConfig config=flow_config();const AmrHierarchy hierarchy=AmrHierarchy::build_static(automatic_flow_domain(mesh,config.domain),mesh,bvh,config.amr);EmbeddedBoundaryBuildOptions options;options.complex_subdivisions=config.amr.complex_subdivisions;options.min_volume_fraction=config.amr.min_volume_fraction;const AmrEmbeddedBoundaryAtlas atlas=build_amr_embedded_boundary_atlas(hierarchy,mesh,bvh,options);if(!atlas.ready_for_flow())return static_cast<std::size_t>(-1);return build_composite_amr_pressure_system(hierarchy,atlas,true).gauges.size();
	}
	struct CaseResult{bool converged=false;AerodynamicLoads loads;double max_divergence=0,last_step_ms=0;int iterations=0;};
	CaseResult run_case(const char* name,const TriMesh& mesh,int steps=8,ParagliderConfig config=flow_config())
	{
		TriangleBvh bvh(mesh);ExternalAeroCore core(mesh,bvh,config);ExternalAeroStepStats stats=core.initialize();bool converged=stats.pressure.converged;for(int step=0;step<steps&&converged;++step){stats=core.step();converged=stats.pressure.converged;}CaseResult result{converged,core.pressure_loads(),core.max_abs_divergence(),stats.gpu_step_ms,stats.pressure.iterations};std::printf("[paraglider-flow] %s: converged=%d F=[%.6g %.6g %.6g] Cd_p=%.6g Cl_p=%.6g maxDiv=%.3e step=%.3f ms it=%d\n",name,converged?1:0,result.loads.pressure_force.x,result.loads.pressure_force.y,result.loads.pressure_force.z,result.loads.cd_pressure,result.loads.cl_pressure,result.max_divergence,result.last_step_ms,result.iterations);return result;
	}
}

int main()
{
	const TriMesh normal=quad({0,-0.5,-0.5},{0,0.5,-0.5},{0,0.5,0.5},{0,-0.5,0.5});
	const TriMesh parallel=quad({-0.5,-0.5,0},{0.5,-0.5,0},{0.5,0.5,0},{-0.5,0.5,0});
	const TriMesh inclined=rotate_y(parallel,15.0*3.14159265358979323846/180.0);
	const CaseResult normal_result=run_case("normal plate",normal),parallel_result=run_case("parallel plate",parallel),inclined_result=run_case("inclined plate",inclined);
	ParagliderConfig uniform_fine=flow_config();uniform_fine.domain={3,4,3,3};uniform_fine.amr.base_cell_size=0.125;ParagliderConfig refined=uniform_fine;refined.amr.base_cell_size=0.25;refined.amr.max_levels=2;refined.amr.wing_refinement_distance=0.1;refined.amr.surface_refinement_distance=0.1;refined.amr.wake_length=0;refined.amr.wake_radius=0;const CaseResult uniform_fine_result=run_case("uniform-fine inclined plate",inclined,8,uniform_fine),refined_result=run_case("2:1 AMR inclined plate",inclined,8,refined);
	check(normal_result.converged&&parallel_result.converged&&inclined_result.converged,"all manufactured plate timesteps converge");
	check(std::isfinite(normal_result.loads.pressure_force.x)&&std::isfinite(parallel_result.loads.pressure_force.x)&&std::isfinite(inclined_result.loads.pressure_force.z),"plate pressure loads are finite");
	check(normal_result.loads.pressure_force.x>0&&normal_result.max_divergence<2e-4,"normal plate has downstream pressure force and projected flow");
	check(std::abs(parallel_result.loads.pressure_force.x)<0.05*normal_result.loads.pressure_force.x,"parallel plate avoids large stair-step pressure drag");
	check(inclined_result.loads.pressure_force.x>=0&&inclined_result.loads.pressure_force.z>0,"positive-angle plate has expected pressure-drag/lift signs");
	const double lift_scale=std::max(1e-12,std::abs(uniform_fine_result.loads.pressure_force.z));check(uniform_fine_result.converged&&refined_result.converged&&std::abs(refined_result.loads.pressure_force.z-uniform_fine_result.loads.pressure_force.z)/lift_scale<0.25,"AMR inclined-plate lift agrees with equal-finest uniform grid");
	const std::size_t closed_gauges=pressure_gauges(cavity_box(false)),open_gauges=pressure_gauges(cavity_box(true));check(closed_gauges>=1&&open_gauges==0,"removing cavity inlet fabric reconnects internal and external air");
	std::printf("[paraglider-flow] dynamic flow probe: %s (%d failures)\n",failures?"FAIL":"PASS",failures);return failures?1:0;
}
