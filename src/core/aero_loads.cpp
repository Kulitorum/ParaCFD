#include "core/aero_loads.h"
#include "core/fluid/amr_pressure.h"

#include <cmath>
#include <stdexcept>

namespace paracfd::core
{
	double EbPressureState::pressure(FragmentRef r) const
	{
		if (fragment_is_regular(r))
		{
			const int i=regular_fragment_cell(r);if(i<0||i>=static_cast<int>(regular.size()))throw std::out_of_range("regular EB pressure reference");return regular[i];
		}
		const int i=irregular_fragment_index(r);if(i<0||i>=static_cast<int>(irregular.size()))throw std::out_of_range("irregular EB pressure reference");return irregular[i];
	}

	AerodynamicLoads compute_pressure_loads(const EmbeddedBoundary& eb,const EbPressureState& p,std::size_t ntri,const FreestreamConfig& fs,const AeroReferenceConfig& ref,double pref)
	{
		AerodynamicLoads out;out.triangles.resize(ntri);std::vector<double> plus_sum(ntri,0),minus_sum(ntri,0),area_sum(ntri,0);const double q=0.5*fs.rho*fs.speed*fs.speed;
		for(const SurfacePatch& patch:eb.patches)
		{
			if(patch.source_triangle_id>=ntri||patch.plus_fragment==invalid_fragment||patch.minus_fragment==invalid_fragment)continue;const double pp=p.pressure(patch.plus_fragment),pm=p.pressure(patch.minus_fragment),dp=pm-pp;const Vec3d force=patch.normal*(dp*patch.area);const Vec3d arm=patch.centroid-ref.moment_origin;
			out.pressure_force=out.pressure_force+force;out.pressure_moment=out.pressure_moment+cross(arm,force);auto& tr=out.triangles[patch.source_triangle_id];tr.pressure_force=tr.pressure_force+force;tr.represented_area+=patch.area;plus_sum[patch.source_triangle_id]+=pp*patch.area;minus_sum[patch.source_triangle_id]+=pm*patch.area;area_sum[patch.source_triangle_id]+=patch.area;
		}
		for(std::size_t i=0;i<ntri;++i)if(area_sum[i]>0)
		{
			auto& tr=out.triangles[i];tr.p_plus=plus_sum[i]/area_sum[i];tr.p_minus=minus_sum[i]/area_sum[i];tr.delta_p=tr.p_minus-tr.p_plus;if(q>0){tr.cp_plus=(tr.p_plus-pref)/q;tr.cp_minus=(tr.p_minus-pref)/q;tr.delta_cp=tr.delta_p/q;}
		}
		out.pressure_drag=out.pressure_force.x;out.pressure_side=out.pressure_force.y;out.pressure_lift=out.pressure_force.z;
		if(q>0&&ref.area>0){out.force_coefficients_valid=true;out.cd_pressure=out.pressure_drag/(q*ref.area);out.cs_pressure=out.pressure_side/(q*ref.area);out.cl_pressure=out.pressure_lift/(q*ref.area);}
		return out;
	}

	AerodynamicLoads compute_pressure_loads(const CompositeAmrPressureSystem& system,const std::vector<double>& p,std::size_t ntri,const FreestreamConfig& fs,const AeroReferenceConfig& ref,double pref)
	{
		if(p.size()<static_cast<std::size_t>(system.storage_size))throw std::invalid_argument("composite aerodynamic pressure vector is too small");AerodynamicLoads out;out.triangles.resize(ntri);std::vector<double> plus_sum(ntri,0),minus_sum(ntri,0),area_sum(ntri,0);const double q=0.5*fs.rho*fs.speed*fs.speed;
		for(const CompositeSurfacePressurePatch& patch:system.surface_patches)
		{
			if(patch.source_triangle_id>=ntri||patch.plus_dof<0||patch.minus_dof<0||patch.plus_dof>=system.storage_size||patch.minus_dof>=system.storage_size)continue;const double pp=p[patch.plus_dof],pm=p[patch.minus_dof],dp=pm-pp;const Vec3d force=patch.normal*(dp*patch.area);const Vec3d arm=patch.centroid-ref.moment_origin;out.pressure_force=out.pressure_force+force;out.pressure_moment=out.pressure_moment+cross(arm,force);auto& tr=out.triangles[patch.source_triangle_id];tr.pressure_force=tr.pressure_force+force;tr.represented_area+=patch.area;plus_sum[patch.source_triangle_id]+=pp*patch.area;minus_sum[patch.source_triangle_id]+=pm*patch.area;area_sum[patch.source_triangle_id]+=patch.area;
		}
		for(std::size_t i=0;i<ntri;++i)if(area_sum[i]>0){auto& tr=out.triangles[i];tr.p_plus=plus_sum[i]/area_sum[i];tr.p_minus=minus_sum[i]/area_sum[i];tr.delta_p=tr.p_minus-tr.p_plus;if(q>0){tr.cp_plus=(tr.p_plus-pref)/q;tr.cp_minus=(tr.p_minus-pref)/q;tr.delta_cp=tr.delta_p/q;}}
		out.pressure_drag=out.pressure_force.x;out.pressure_side=out.pressure_force.y;out.pressure_lift=out.pressure_force.z;if(q>0&&ref.area>0){out.force_coefficients_valid=true;out.cd_pressure=out.pressure_drag/(q*ref.area);out.cs_pressure=out.pressure_side/(q*ref.area);out.cl_pressure=out.pressure_lift/(q*ref.area);}return out;
	}
}
