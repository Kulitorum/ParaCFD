#include "gui/slice_field.h"
#include "core/fluid/amr_sampling.h"
#include "core/aero_sweep.h"

#include <cmath>
#include <cstdio>
#include <vector>

namespace
{
	bool near(float actual,float expected,float tolerance,const char* label)
	{
		if(std::abs(actual-expected)<=tolerance)return true;
		std::fprintf(stderr,"visualization probe: %s expected %.9g, got %.9g\n",label,expected,actual);return false;
	}

	bool composite_amr_sampling()
	{
		using namespace paracfd::core;
		TriMesh marker;marker.positions={0.25f,0.2f,0.2f,0.25f,0.8f,0.2f,0.25f,0.5f,0.8f};
		marker.indices={0,1,2};marker.bbox_min={0.25f,0.2f,0.2f};marker.bbox_max={0.25f,0.8f,0.8f};
		TriangleBvh bvh(marker);AmrConfig config;config.base_cell_size=0.25;config.max_levels=2;
		config.brick_size=4;config.ghost_cells=1;config.wing_refinement_distance=0;
		config.surface_refinement_distance=0;config.wake_length=0;config.wake_radius=0;
		const AmrHierarchy hierarchy=AmrHierarchy::build_static({{0,0,0},{2,1,1}},marker,bvh,config);
		if(hierarchy.levels().size()!=2||hierarchy.levels()[1].active_bricks==0)
		{
			std::fprintf(stderr,"visualization probe: synthetic AMR transition was not built\n");return false;
		}
		AmrHostFields fields(hierarchy);
		for(int level_index=0;level_index<static_cast<int>(hierarchy.levels().size());++level_index)
		{
			const AmrLevel& metadata=hierarchy.levels()[level_index];
			AmrHostLevelFields& values=fields.levels()[level_index];const int bs=hierarchy.brick_size();
			for(int brick=0;brick<static_cast<int>(metadata.bricks.size());++brick)
			{
				const BrickMetadata& record=metadata.bricks[brick];
				for(int component=0;component<3;++component)
				{
					const int ni=component==0?bs+1:bs,nj=component==1?bs+1:bs,nk=component==2?bs+1:bs;
					for(int k=0;k<nk;++k)for(int j=0;j<nj;++j)for(int i=0;i<ni;++i)
					{
						const Vec3d point{record.origin.x+(i+(component==0?0.0:0.5))*record.h,
							record.origin.y+(j+(component==1?0.0:0.5))*record.h,
							record.origin.z+(k+(component==2?0.0:0.5))*record.h};
						double value=component==0?1+2*point.x-0.5*point.y+0.25*point.z:
							(component==1?-2+0.3*point.x+point.y-0.4*point.z:0.5-0.2*point.x+0.7*point.y+point.z);
						const std::size_t index=component==0?values.layout.u_index(brick,i,j,k):
							(component==1?values.layout.v_index(brick,i,j,k):values.layout.w_index(brick,i,j,k));
						if(component==0)values.u[index]=static_cast<Real>(value);
						else if(component==1)values.v[index]=static_cast<Real>(value);else values.w[index]=static_cast<Real>(value);
					}
				}
			}
		}
		bool ok=true;const Vec3d point{1.17,0.37,0.63};const Vec3d sampled=sample_amr_velocity(hierarchy,fields,point);
		ok&=near(static_cast<float>(sampled.x),static_cast<float>(1+2*point.x-0.5*point.y+0.25*point.z),2e-5f,"AMR affine u");
		ok&=near(static_cast<float>(sampled.y),static_cast<float>(-2+0.3*point.x+point.y-0.4*point.z),2e-5f,"AMR affine v");
		ok&=near(static_cast<float>(sampled.z),static_cast<float>(0.5-0.2*point.x+0.7*point.y+point.z),2e-5f,"AMR affine w");

        const double circulation=rectangular_circulation_xz(hierarchy,fields,
            0.213,1.577,0.37,0.173,0.863);
        const double analytic_circulation=(-0.2-0.25)*(1.577-0.213)*(0.863-0.173);
        ok&=near(static_cast<float>(circulation),static_cast<float>(analytic_circulation),
            2e-5f,"AMR off-grid contour circulation");

        // Change d(w)/dx to d(u)/dz, giving a curl-free affine field.
        for(int level_index=0;level_index<static_cast<int>(hierarchy.levels().size());++level_index)
        {
            const auto& metadata=hierarchy.levels()[level_index];auto& values=fields.levels()[level_index];
            const int bs=hierarchy.brick_size();
            for(int brick=0;brick<static_cast<int>(metadata.bricks.size());++brick)
                for(int k=0;k<=bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i)
                {
                    const auto& record=metadata.bricks[brick];
                    const double x=record.origin.x+(i+0.5)*record.h;
                    values.w[values.layout.w_index(brick,i,j,k)]+=static_cast<Real>(0.45*x);
                }
        }
        ok&=near(static_cast<float>(rectangular_circulation_xz(hierarchy,fields,
            0.213,1.577,0.37,0.173,0.863)),0.0f,2e-5f,"AMR curl-free circulation");

		// Give fine normal faces a resolved tangential mode whose coarse conservative face
		// stores only its mean. The visualization limits must still meet at the interface.
		for(int level_index=1;level_index<static_cast<int>(hierarchy.levels().size());++level_index)
		{
			const AmrLevel& metadata=hierarchy.levels()[level_index];AmrHostLevelFields& values=fields.levels()[level_index];
			for(int brick=0;brick<static_cast<int>(metadata.bricks.size());++brick)
			{
				const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;
				for(int k=0;k<hierarchy.brick_size();++k)for(int j=0;j<hierarchy.brick_size();++j)
					for(int i=0;i<=hierarchy.brick_size();++i)
					{
						const double y=record.origin.y+(j+0.5)*record.h;
						values.u[values.layout.u_index(brick,i,j,k)]=static_cast<Real>(10+0.2*(y-0.5));
					}
			}
		}
		const double epsilon=1e-8;const Vec3d fine=sample_amr_velocity(hierarchy,fields,{1-epsilon,0.37,0.63});
		const Vec3d coarse=sample_amr_velocity(hierarchy,fields,{1+epsilon,0.37,0.63});
		ok&=near(static_cast<float>(coarse.x),static_cast<float>(fine.x),2e-5f,"AMR coarse/fine interface u");
		return ok;
	}

	bool fixed_attitude_glide_trim()
	{
		using namespace paracfd::core;
		const std::vector<FixedAttitudeGlideSample> samples={{0,100,1000},{10,100,1000}};
		const FixedAttitudeGlideTrim trim=solve_fixed_attitude_glide(samples,10,100);
		bool ok=trim.equilibrium;
		ok&=near(static_cast<float>(trim.flight_path_degrees),5.71059f,0.01f,"trim flight-path angle");
		ok&=near(static_cast<float>(trim.glide_ratio),10.0f,0.02f,"trim glide ratio");
		ok&=near(static_cast<float>(trim.reference_vertical_support),1004.9875f,0.2f,
			"trim reference support");
		ok&=near(static_cast<float>(trim.sink_rate),0.983f,0.01f,"trim sink rate");
		const FixedAttitudeGlideTrim missing=solve_fixed_attitude_glide({{0,100,100},{10,100,100}},10,100);
		if(missing.equilibrium){std::fprintf(stderr,"visualization probe: false glide equilibrium\n");ok=false;}
		return ok;
	}
}

int main()
{
	using paracfd::core::MacGrid;using paracfd::gui::Field;
	MacGrid grid;grid.nx=grid.ny=grid.nz=4;grid.h=0.5;
	std::vector<double> u(grid.u_count()),v(grid.v_count()),w(grid.w_count()),p(grid.p_count());
	// Rigid rotation u=(-y,x,0) has |curl u|=2 and Q=1 everywhere. Pressure is a
	// separate affine field whose gradient magnitude is sqrt(1+4+9)=sqrt(14).
	for(int k=0;k<grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<=grid.nx;++i)
		u[grid.uidx(i,j,k)]=-(j+0.5)*grid.h;
	for(int k=0;k<grid.nz;++k)for(int j=0;j<=grid.ny;++j)for(int i=0;i<grid.nx;++i)
		v[grid.vidx(i,j,k)]=(i+0.5)*grid.h;
	for(int k=0;k<=grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<grid.nx;++i)
		w[grid.widx(i,j,k)]=0.0;
	for(int k=0;k<grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<grid.nx;++i)
	{
		const double x=(i+0.5)*grid.h,y=(j+0.5)*grid.h,z=(k+0.5)*grid.h;
		p[grid.pidx(i,j,k)]=x+2*y+3*z;
	}
	float range[3]{};bool ok=composite_amr_sampling()&&fixed_attitude_glide_trim();
	paracfd::gui::slice_reduce_cpu(u.data(),v.data(),w.data(),p.data(),nullptr,grid,Field::VorticityMagnitude,range);
	ok&=near(range[0],2.0f,1e-5f,"vorticity min")&&near(range[1],2.0f,1e-5f,"vorticity max");
	paracfd::gui::slice_reduce_cpu(u.data(),v.data(),w.data(),p.data(),nullptr,grid,Field::QCriterion,range);
	ok&=near(range[0],1.0f,1e-5f,"Q min")&&near(range[1],1.0f,1e-5f,"Q max");
	paracfd::gui::slice_reduce_cpu(u.data(),v.data(),w.data(),p.data(),nullptr,grid,Field::PressureGradient,range);
	ok&=near(range[0],std::sqrt(14.0f),1e-5f,"pressure-gradient min")&&near(range[1],std::sqrt(14.0f),1e-5f,"pressure-gradient max");
	paracfd::gui::slice_reduce_cpu(u.data(),v.data(),w.data(),p.data(),nullptr,grid,Field::PressureCoefficient,range,2.0f,4.0f);
	ok&=near(range[0],1.5f/16.0f,1e-6f,"Cp min")&&near(range[1],10.5f/16.0f,1e-6f,"Cp max");
	if(ok)std::printf("visualization probe: derived scalar fields and composite AMR sampling pass\n");return ok?0:1;
}
