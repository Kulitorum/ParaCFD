#include "core/fluid/external_aero_core.h"

#include <algorithm>
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
	struct OpeningFluxResult{bool converged=false;double positive_flux=0,absolute_flux=0;};
	CaseResult run_case(const char* name,const TriMesh& mesh,int steps=8,ParagliderConfig config=flow_config())
	{
		TriangleBvh bvh(mesh);ExternalAeroCore core(mesh,bvh,config);ExternalAeroStepStats stats=core.initialize();bool converged=stats.pressure.converged;for(int step=0;step<steps&&converged;++step){stats=core.step();converged=stats.pressure.converged;}CaseResult result{converged,core.pressure_loads(),core.max_abs_divergence(),stats.gpu_step_ms,stats.pressure.iterations};std::printf("[paraglider-flow] %s: converged=%d F=[%.6g %.6g %.6g] Cd_p=%.6g Cl_p=%.6g maxDiv=%.3e step=%.3f ms it=%d\n",name,converged?1:0,result.loads.pressure_force.x,result.loads.pressure_force.y,result.loads.pressure_force.z,result.loads.cd_pressure,result.loads.cl_pressure,result.max_divergence,result.last_step_ms,result.iterations);return result;
	}
	OpeningFluxResult run_opening_flux_case(bool upstream_open,int steps=12)
	{
		const TriMesh mesh=cavity_box(upstream_open);TriangleBvh bvh(mesh);const ParagliderConfig config=flow_config();ExternalAeroCore core(mesh,bvh,config);ExternalAeroStepStats stats=core.initialize();bool converged=stats.pressure.converged;for(int step=0;step<steps&&converged;++step){stats=core.step();converged=stats.pressure.converged;}AmrHostFields fields(core.hierarchy());core.download_fields(fields);CompositeAmrFluxes special;core.download_special_fluxes(special);const CompositeAmrPressureSystem& system=core.pressure_system();const AmrHierarchy& hierarchy=core.hierarchy();const int bs=hierarchy.brick_size();OpeningFluxResult result;result.converged=converged;
		for(int level=0;level<(int)hierarchy.levels().size();++level){const AmrLevel& metadata=hierarchy.levels()[level];const auto& values=fields.levels()[level];for(int brick=0;brick<(int)metadata.bricks.size();++brick){const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const double x=record.origin.x+(i+1)*metadata.h,y=record.origin.y+(j+0.5)*metadata.h,z=record.origin.z+(k+0.5)*metadata.h;if(std::abs(x+0.5)>1e-8||y<=-0.5||y>=0.5||z<=-0.5||z>=0.5)continue;const int a=system.dof(level,brick,i,j,k);if(!system.active[a]||(system.cut_face_mask[a]&1u))continue;int b=-1;if(i+1<bs)b=system.dof(level,brick,i+1,j,k);else{const int neighbour=record.same_level_neighbor[1];if(neighbour>=0&&metadata.bricks[neighbour].active())b=system.dof(level,neighbour,0,j,k);}if(b<0||!system.active[b])continue;const double velocity=values.u[values.layout.u_index(brick,i+1,j,k)],area=metadata.h*metadata.h;result.positive_flux+=std::max(0.0,velocity)*area;result.absolute_flux+=std::abs(velocity)*area;}}}
		for(std::size_t edge=0;edge<system.embedded.size();++edge){const auto& connection=system.embedded[edge];if(connection.axis!=0||std::abs(connection.face_centroid.x+0.5)>1e-8||connection.face_centroid.y<=-0.5||connection.face_centroid.y>=0.5||connection.face_centroid.z<=-0.5||connection.face_centroid.z>=0.5)continue;const double velocity=special.embedded_velocity[edge];result.positive_flux+=std::max(0.0,velocity)*connection.open_area;result.absolute_flux+=std::abs(velocity)*connection.open_area;}
		std::printf("[paraglider-flow] %s cavity opening: converged=%d positive-flux=%.6g abs-flux=%.6g\n",upstream_open?"open":"closed",converged?1:0,result.positive_flux,result.absolute_flux);return result;
	}
}

int main()
{
	const TriMesh normal=quad({0,-0.5,-0.5},{0,0.5,-0.5},{0,0.5,0.5},{0,-0.5,0.5});
	const TriMesh parallel=quad({-0.5,-0.5,0},{0.5,-0.5,0},{0.5,0.5,0},{-0.5,0.5,0});
	const TriMesh inclined=rotate_y(parallel,15.0*3.14159265358979323846/180.0);
	// The impulsively initialized normal plate needs a longer pressure adjustment than
	// the tangential cases before its force reaches the downstream-sign regime.
	const CaseResult normal_result=run_case("normal plate",normal,16),parallel_result=run_case("parallel plate",parallel),inclined_result=run_case("inclined plate",inclined);
	ParagliderConfig uniform_fine=flow_config();uniform_fine.domain={3,4,3,3};uniform_fine.amr.base_cell_size=0.125;ParagliderConfig refined=uniform_fine;refined.amr.base_cell_size=0.25;refined.amr.max_levels=2;refined.amr.wing_refinement_distance=0.1;refined.amr.surface_refinement_distance=0.1;refined.amr.wake_length=0;refined.amr.wake_radius=0;const CaseResult uniform_fine_result=run_case("uniform-fine inclined plate",inclined,8,uniform_fine),refined_result=run_case("2:1 AMR inclined plate",inclined,8,refined);
	check(normal_result.converged&&parallel_result.converged&&inclined_result.converged,"all manufactured plate timesteps converge");
	check(std::isfinite(normal_result.loads.pressure_force.x)&&std::isfinite(parallel_result.loads.pressure_force.x)&&std::isfinite(inclined_result.loads.pressure_force.z),"plate pressure loads are finite");
	check(normal_result.loads.pressure_force.x>0&&normal_result.max_divergence<2e-4,"normal plate has downstream pressure force and projected flow");
	check(std::abs(parallel_result.loads.pressure_force.x)<0.05*normal_result.loads.pressure_force.x,"parallel plate avoids large stair-step pressure drag");
	check(inclined_result.loads.pressure_force.x>=0&&inclined_result.loads.pressure_force.z>0,"positive-angle plate has expected pressure-drag/lift signs");
	const double lift_scale=std::max(1e-12,std::abs(uniform_fine_result.loads.pressure_force.z));check(uniform_fine_result.converged&&refined_result.converged&&std::abs(refined_result.loads.pressure_force.z-uniform_fine_result.loads.pressure_force.z)/lift_scale<0.25,"AMR inclined-plate lift agrees with equal-finest uniform grid");
	const std::size_t closed_gauges=pressure_gauges(cavity_box(false)),open_gauges=pressure_gauges(cavity_box(true));check(closed_gauges>=1&&open_gauges==0,"removing cavity inlet fabric reconnects internal and external air");
	const OpeningFluxResult closed_flux=run_opening_flux_case(false),open_flux=run_opening_flux_case(true);check(closed_flux.converged&&open_flux.converged&&closed_flux.absolute_flux<1e-10&&open_flux.positive_flux>1e-3,"developed flow enters the cavity only through the real opening");
	std::printf("[paraglider-flow] dynamic flow probe: %s (%d failures)\n",failures?"FAIL":"PASS",failures);return failures?1:0;
}
