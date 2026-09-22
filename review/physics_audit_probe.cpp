// Review-only diagnostics against the production library; no solver changes.
#include "core/fluid/external_aero_core.h"
#include "core/geometry/step_import.h"
#include "core/geometry/mesh_clip.h"
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

using namespace paracfd::core;

double difference(const std::vector<Real>& a,const std::vector<Real>& b)
{
    double result=0;
    for(std::size_t q=0;q<a.size();++q)
        result=std::max(result,std::abs(double(a[q])-double(b[q])));
    return result;
}

int main(int argc,char** argv)
{
    try
    {
        std::setvbuf(stdout,nullptr,_IONBF,0);
        ParagliderConfig config;std::string error;
        if(!load_paraglider_config(argc>1?argv[1]:"configs/planb_parakite.json",config,&error))
            throw std::runtime_error(error);
        const auto imported=load_step_solid_geometry(config.step_path,config.tessellation_deflection_mm,&error);
        if(imported.mesh.empty())throw std::runtime_error(error);
        TriMesh source=imported.mesh;
        const bool slab=argc>2&&std::string(argv[2])=="--naca-slab";
        if(slab)
        {
            source=clip_mesh_to_axis_slab(source,1,-0.0625,0.0625);
            config.domain.lateral_margin=0;
            config.amr.surface_refinement_distance=0.15;
            config.domain.upstream_margin+=0.5*config.amr.base_cell_size/(1<<(config.amr.max_levels-1));
        }
        config.placement=frame_wing_for_external_domain(source,config.placement,
            config.domain.upstream_margin,config.domain.lateral_margin,
            config.domain.vertical_margin,config.amr.base_cell_size*config.amr.brick_size);
        const auto wing=placed_mesh(source,config.placement);
        TriangleBvh bvh(wing);
        ExternalAeroExecutionOptions options;
        options.closed_solid=imported.closed_solid->placed(config.placement);
        ExternalAeroCore core(wing,bvh,config,options);
        const auto& system=core.pressure_system();
        std::vector<Vec3d> sealed_points;
        std::size_t sealed_count=0;double sealed_volume=0;
        for(int q=0;q<system.storage_size;++q)if(system.active[q]&&!system.freestream_connected[q])
        {
            ++sealed_count;sealed_volume+=system.volume[q];
            if(sealed_points.size()<32)sealed_points.push_back(system.centroid[q]);
        }
        std::vector<ClosedSolidPointLocation> locations(sealed_points.size());
        if(!sealed_points.empty())options.closed_solid->classify_points(sealed_points.data(),
            sealed_points.size(),1e-9,locations.data());
        std::size_t inside=0;for(const auto location:locations)inside+=location==ClosedSolidPointLocation::inside;
        std::printf("[audit-solid-mask] sealed-active=%zu volume=%.9g sampled=%zu classified-inside-solid=%zu\n",
            sealed_count,sealed_volume,sealed_points.size(),inside);

        // For p=x/y/z, compare the production compact two-point face derivative
        // with the exact Cartesian normal derivative. No PDE solve is needed.
        for(int component=0;component<3;++component)
        {
            for(int kind=0;kind<2;++kind)
            {
                const auto& edges=kind?system.embedded:system.coarse_fine;
                double maximum=0,square=0,total_area=0;std::size_t bad=0;
                for(const auto& edge:edges)
                {
                    const int lo=edge.direction>0?edge.coarse_dof:edge.fine_dof;
                    const int hi=edge.direction>0?edge.fine_dof:edge.coarse_dof;
                    const double value=pressure_gradient_factor(edge)*
                        (system.centroid[hi][component]-system.centroid[lo][component]);
                    const double defect=value-(component==edge.axis?1.0:0.0);
                    maximum=std::max(maximum,std::abs(defect));
                    square+=edge.open_area*defect*defect;total_area+=edge.open_area;
                    bad+=std::abs(defect)>1e-6;
                }
                std::printf("[audit-affine-pressure] p=%c faces=%s count=%zu bad=%zu max=%.9g area-rms=%.9g\n",
                    'x'+component,kind?"EB":"AMR",edges.size(),bad,maximum,
                    total_area>0?std::sqrt(square/total_area):0);
            }
        }
        const auto closure=composite_cell_pressure_closure_cpu(system);
        double max_closure=0,sum_closure=0;std::size_t bad_rows=0;
        for(int q=0;q<system.storage_size;++q)if(system.active[q])
        {
            const double value=std::sqrt(length2(closure[q]));
            max_closure=std::max(max_closure,value);sum_closure+=value;bad_rows+=value>1e-10;
        }
        std::vector<double> constant(system.storage_size,1.0);
        const auto constant_load=compute_pressure_loads(system,constant,wing.triangle_count(),
            config.freestream,config.reference,0.0);
        std::printf("[audit-geometry] closure-rows=%zu max-area-defect=%.9g sum-area-defect=%.9g unit-pressure-force=[%.9g %.9g %.9g]\n",
            bad_rows,max_closure,sum_closure,constant_load.pressure_force.x,
            constant_load.pressure_force.y,constant_load.pressure_force.z);
        if(slab)return 0;

        const auto initial=core.initialize();
        if(!initial.pressure.converged)throw std::runtime_error("initialization failed");
        AmrHostFields before(core.hierarchy());core.download_fields(before);
        CompositeAmrFluxes before_flux;core.download_special_fluxes(before_flux);
        DeviceAmrFields fields(core.hierarchy());
        DeviceCompositeAmrProjection projection(system,fields);
        for(const int evolved_steps:{0,100})
        {
        if(evolved_steps)
        {
            for(int step=0;step<evolved_steps;++step)
                if(!core.step().pressure.converged)throw std::runtime_error("evolved projection failed");
            core.download_fields(before);core.download_special_fluxes(before_flux);
        }
        std::printf("[audit-state] evolved-steps=%d physical-time=%.9g\n",evolved_steps,core.physical_time());
        std::vector<double> first_trial;
        for(const double dt:{1e-6,1e-9,1e-12,1e-12})
        {
            fields.upload(before);projection.upload_special_fluxes(before_flux);
            projection.transport_embedded_apertures(Real(dt),Real(0),Real(0),false);
            AmrHostFields after(core.hierarchy());fields.download(after);
            CompositeAmrFluxes after_flux;projection.download_special_fluxes(after_flux);
            double regular_change=0;
            for(std::size_t l=0;l<before.levels().size();++l)
            {
                const auto& a=before.levels()[l];const auto& b=after.levels()[l];
                regular_change=std::max({regular_change,difference(a.u,b.u),difference(a.v,b.v),difference(a.w,b.w)});
            }
            double special_change=0,repeat_change=0;
            for(std::size_t q=0;q<before_flux.embedded_velocity.size();++q)
            {
                special_change=std::max(special_change,std::abs(double(after_flux.embedded_velocity[q])-double(before_flux.embedded_velocity[q])));
                if(dt==1e-12&&!first_trial.empty())repeat_change=std::max(repeat_change,
                    std::abs(double(after_flux.embedded_velocity[q])-double(first_trial[q])));
            }
            std::printf("[audit-zero-dt] dt=%.3g nu=0 Cs=0 regular-change=%.9g EB-change=%.9g repeated-EB-difference=%.9g m/s\n",
                dt,regular_change,special_change,repeat_change);
            if(dt==1e-12)first_trial.assign(after_flux.embedded_velocity.begin(),after_flux.embedded_velocity.end());
        }
        }

        // A compactly supported manufactured pressure is zero on every outer
        // boundary. Its fluid impulse and CAD reaction should therefore cancel.
        CompositeCellMomentumState zero;
        zero.x.assign(system.storage_size,0);zero.y=zero.z=zero.x;
        std::vector<Real> p(system.storage_size,0);std::vector<double> p64(system.storage_size,0);
        Vec3d centre{},radius{};
        for(int axis=0;axis<3;++axis)
        {
            centre[axis]=0.5*(wing.bbox_min[axis]+wing.bbox_max[axis]);
            radius[axis]=0.5*(wing.bbox_max[axis]-wing.bbox_min[axis])+0.5;
        }
        for(int q=0;q<system.storage_size;++q)if(system.active[q])
        {
            const Vec3d d=system.centroid[q]-centre;double value=1;
            for(int axis=0;axis<3;++axis)value*=std::max(0.0,1-d[axis]*d[axis]/(radius[axis]*radius[axis]));
            p[q]=Real(100*value*value*(1+0.5*d.z/radius.z));p64[q]=p[q];
        }
        Real* device_pressure=nullptr;
        if(cudaMalloc(&device_pressure,p.size()*sizeof(Real))!=cudaSuccess)throw std::runtime_error("pressure allocation");
        if(cudaMemcpy(device_pressure,p.data(),p.size()*sizeof(Real),cudaMemcpyHostToDevice)!=cudaSuccess)throw std::runtime_error("pressure upload");
        DeviceCompositeCellMomentumTransport cell(system,fields,config.amr.min_volume_fraction);
        cell.upload_state(zero);const Real dt=Real(0.001),rho=Real(config.freestream.rho);
        cell.apply_projected_pressure_gradient(device_pressure,dt,rho);
        const auto momentum=cell.momentum();cudaFree(device_pressure);
        const auto load=compute_pressure_loads(system,p64,wing.triangle_count(),config.freestream,config.reference,0);
        const Vec3d reaction{-double(rho)*momentum[0]/dt,-double(rho)*momentum[1]/dt,-double(rho)*momentum[2]/dt};
        const Vec3d mismatch=reaction-load.pressure_force;
        std::printf("[audit-pressure-reaction] CAD=[%.9g %.9g %.9g] negative-fluid-impulse/dt=[%.9g %.9g %.9g] mismatch=[%.9g %.9g %.9g] N\n",
            load.pressure_force.x,load.pressure_force.y,load.pressure_force.z,
            reaction.x,reaction.y,reaction.z,mismatch.x,mismatch.y,mismatch.z);
        return 0;
    }
    catch(const std::exception& e){std::fprintf(stderr,"[audit-error] %s\n",e.what());return 1;}
}
