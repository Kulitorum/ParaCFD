// Manufactured tests of the production solid/AMR/projection/momentum operators.
#include "core/fluid/external_aero_core.h"
#include "core/geometry/step_import.h"
#include "core/geometry/mesh_clip.h"
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

using namespace paracfd::core;

namespace
{
    void require(bool condition,const char* message)
    {
        if(!condition)throw std::runtime_error(message);
    }
    TriMesh box_mesh(Vec3d half)
    {
        TriMesh mesh;
        for(int vertex=0;vertex<8;++vertex)
            for(int axis=0;axis<3;++axis)
            {
                const double value=(vertex&(1<<axis))?half[axis]:-half[axis];
                mesh.positions_fp64.push_back(value);mesh.positions.push_back(float(value));
                mesh.normals.push_back(0);
            }
        const int faces[6][4]={{0,4,6,2},{1,3,7,5},{0,1,5,4},
            {2,6,7,3},{0,2,3,1},{4,5,7,6}};
        for(int face=0;face<6;++face)for(int triangle=0;triangle<2;++triangle)
        {
            for(int corner:{0,triangle+1,triangle+2})mesh.indices.push_back(faces[face][corner]);
            mesh.source_face_ids.push_back(face);
        }
        for(int axis=0;axis<3;++axis){mesh.bbox_min[axis]=float(-half[axis]);mesh.bbox_max[axis]=float(half[axis]);}
        return mesh;
    }
    class BoxGeometry final:public ClosedSolidGeometry
    {
        Vec3d half_;ModelPlacement inverse_;Aabb3d bounds_;
    public:
        BoxGeometry(Vec3d half,const ModelPlacement& placement):half_(half)
        {
            require(inverse_placement(placement,inverse_),"invalid test placement");
            for(int vertex=0;vertex<8;++vertex)
            {
                Vec3d p;placement.apply((vertex&1)?half.x:-half.x,
                    (vertex&2)?half.y:-half.y,(vertex&4)?half.z:-half.z,p.x,p.y,p.z);
                for(int axis=0;axis<3;++axis)
                    {bounds_.lo[axis]=std::min(bounds_.lo[axis],p[axis]);bounds_.hi[axis]=std::max(bounds_.hi[axis],p[axis]);}
            }
        }
        const Aabb3d& bounds()const override{return bounds_;}
        void classify_points(const Vec3d* points,std::size_t count,double tolerance,
            ClosedSolidPointLocation* locations)const override
        {
            for(std::size_t q=0;q<count;++q)
            {
                Vec3d p;inverse_.apply(points[q].x,points[q].y,points[q].z,p.x,p.y,p.z);
                double distance=-1e30;
                for(int axis=0;axis<3;++axis)distance=std::max(distance,std::abs(p[axis])-half_[axis]);
                locations[q]=distance>tolerance?ClosedSolidPointLocation::outside:
                    (distance<-tolerance?ClosedSolidPointLocation::inside:ClosedSolidPointLocation::boundary);
            }
        }
    };
    double norm(Vec3d value){return std::sqrt(length2(value));}
    double field_difference(const AmrHostFields& a,const AmrHostFields& b)
    {
        double result=0;
        auto compare=[&](const auto& x,const auto& y)
        {
            for(std::size_t q=0;q<x.size();++q)
            {
                require(std::isfinite(x[q])&&std::isfinite(y[q]),"nonfinite velocity");
                result=std::max(result,std::abs(double(x[q])-double(y[q])));
            }
        };
        for(std::size_t l=0;l<a.levels().size();++l)
        {
            compare(a.levels()[l].u,b.levels()[l].u);
            compare(a.levels()[l].v,b.levels()[l].v);
            compare(a.levels()[l].w,b.levels()[l].w);
        }
        return result;
    }
}

int main(int argc,char** argv)
{
    try
    {
        std::setvbuf(stdout,nullptr,_IONBF,0);
        bool planb=false,cropped=false,geometry_only=false,interior_brick=false,local=false,section=false,cell_momentum=false,export_mesh=false,conditioning=false;
        Vec3d local_origin{};
        for(int q=1;q<argc;++q)
        {
            const std::string arg=argv[q];
            if(arg=="--local-planb"&&q+3<argc){planb=local=true;local_origin.x=std::stod(argv[++q]);local_origin.y=std::stod(argv[++q]);local_origin.z=std::stod(argv[++q]);}
            else if(arg=="--cell-momentum")cell_momentum=true;
            else if(arg=="--conditioning")conditioning=true;
            else if(arg=="--section-planb")planb=section=true;
            else if(arg=="--export-planb")planb=export_mesh=true;
            else if(arg=="--planb")planb=true;else if(arg=="--cropped")cropped=true;
            else if(arg=="--geometry-only")geometry_only=true;
            else if(arg=="--interior-brick")interior_brick=true;
            else throw std::runtime_error("unknown physics contract option");
        }
        ParagliderConfig config;TriMesh source;ClosedSolidSourcePtr solid_source;
        const Vec3d half=interior_brick?Vec3d{1.4,1.4,1.4}:Vec3d{0.6,0.8,0.18};
        if(planb)
        {
            std::string error;require(load_paraglider_config("configs/planb_parakite.json",config,&error),error.c_str());
            const auto imported=load_step_solid_geometry(config.step_path,config.tessellation_deflection_mm,&error);
            require(!imported.mesh.empty(),error.c_str());source=imported.mesh;solid_source=imported.closed_solid;
        }
        else
        {
            source=box_mesh(half);config.freestream.speed=2;
            config.domain={1.25,1.25,1.25,1.25,false};
            config.amr.base_cell_size=0.125;config.amr.brick_size=8;config.amr.max_levels=2;
            config.amr.wing_refinement_distance=0.2;config.amr.surface_refinement_distance=0.2;
            config.amr.wake_length=0.5;config.amr.wake_radius=0.5;config.amr.min_volume_fraction=0.005;
            config.solver.projection_tolerance=1e-5;config.solver.projection_max_iterations=600;
            const double pitch=0.31,yaw=cropped?0:0.17,c=std::cos(pitch),s=std::sin(pitch),a=std::cos(yaw),b=std::sin(yaw);
            const double rotation[9]={a*c,-b,a*s,b*c,a,b*s,-s,0,c};
            std::copy(rotation,rotation+9,config.placement.m);
            if(cropped)
            {
                source=clip_mesh_to_axis_slab(source,1,-0.25,0.25);
                config.domain.lateral_margin=0;config.amr.base_cell_size=0.0625;config.amr.max_levels=1;
            }
        }
        config.placement=frame_wing_for_external_domain(source,config.placement,
            config.domain.upstream_margin,config.domain.lateral_margin,
            config.domain.vertical_margin,config.amr.base_cell_size*config.amr.brick_size);
        const auto wing=placed_mesh(source,config.placement);TriangleBvh bvh(wing);
        if(export_mesh)
        {
            FILE* output=std::fopen("review/final-planb-placed.obj","w");
            require(output,"could not write placed geometry diagnostic");
            std::fprintf(output,"# Supplied PlanB placement; no extra rotation\n");
            for(std::size_t vertex=0;vertex<wing.positions.size()/3;++vertex)
            {
                const auto p=wing.vertex_position_double(vertex);
                std::fprintf(output,"v %.17g %.17g %.17g\n",p[0],p[1],p[2]);
            }
            for(std::size_t triangle=0;triangle<wing.triangle_count();++triangle)
                std::fprintf(output,"f %u %u %u\n",wing.indices[3*triangle]+1,
                    wing.indices[3*triangle+1]+1,wing.indices[3*triangle+2]+1);
            std::fclose(output);std::puts("[contracts] exported unchanged PlanB geometry");return 0;
        }
        if(section)
        {
            const double y=0.5*(wing.bbox_min[1]+wing.bbox_max[1])+0.001;
            FILE* output=std::fopen("review/planb-neutral-section.csv","w");require(output,"could not write section diagnostic");
            std::fprintf(output,"x,z_lower,z_upper\n");
            for(int q=0;q<=4000;++q)
            {
                const double x=wing.bbox_min[0]+(wing.bbox_max[0]-wing.bbox_min[0])*q/4000.;
                const Vec3d lo{x,y,wing.bbox_min[2]-1},hi{x,y,wing.bbox_max[2]+1};
                const auto lower=bvh.intersect_segment(lo,hi),upper=bvh.intersect_segment(hi,lo);
                if(lower.hit&&upper.hit)std::fprintf(output,"%.17g,%.17g,%.17g\n",x,lower.position.z,upper.position.z);
            }
            std::fclose(output);
            // Preserve overhangs/intake details for an independent section solver.
            // Vertical lower/upper envelopes above are only a plotting diagnostic.
            output=std::fopen("review/planb-neutral-section-segments.csv","w");
            require(output,"could not write section segments");
            std::fprintf(output,"x0,z0,x1,z1,face\n");
            for(std::size_t triangle=0;triangle<wing.triangle_count();++triangle)
            {
                Vec3d vertex[3],hit[3];int count=0;
                for(int c=0;c<3;++c)
                {
                    const auto p=wing.vertex_position_double(wing.indices[3*triangle+c]);
                    vertex[c]={p[0],p[1],p[2]};
                }
                for(int c=0;c<3;++c)
                {
                    const Vec3d a=vertex[c],b=vertex[(c+1)%3];
                    if((a.y<y&&b.y>=y)||(b.y<y&&a.y>=y))
                        hit[count++]=a+(b-a)*((y-a.y)/(b.y-a.y));
                }
                if(count==2)std::fprintf(output,"%.17g,%.17g,%.17g,%.17g,%u\n",
                    hit[0].x,hit[0].z,hit[1].x,hit[1].z,wing.source_face_ids[triangle]);
            }
            std::fclose(output);std::puts("[contracts] wrote neutral section envelopes and exact segments");return 0;
        }
        ExternalAeroExecutionOptions options;options.conservative_cell_momentum=cell_momentum;
        options.closed_solid=planb?solid_source->placed(config.placement):
            std::make_shared<BoxGeometry>(half,config.placement);
        if(local)
        {
            UniformEbGrid grid{local_origin,3,3,3,0.125};EmbeddedBoundaryBuildOptions eb_options;
            const auto eb=build_closed_solid_embedded_boundary(wing,bvh,*options.closed_solid,grid,eb_options);
            require(eb.unresolved.empty(),"local unresolved geometry");
            std::vector<Vec3d> closure(eb.fragments.size());
            auto add=[&](FragmentRef ref,Vec3d value){if(ref&&!fragment_is_regular(ref))closure[irregular_fragment_index(ref)]=closure[irregular_fragment_index(ref)]+value;};
            for(const auto& face:eb.apertures){Vec3d area{};area[face.axis]=face.area;add(face.fragment_a,area);add(face.fragment_b,area*-1);}
            for(const auto& face:eb.boundary_apertures){Vec3d area{};area[face.axis]=face.direction*face.area;add(face.fragment,area);}
            for(const auto& patch:eb.patches){add(patch.plus_fragment,patch.normal*(-patch.area));add(patch.minus_fragment,patch.normal*patch.area);}
            double maximum=0;for(std::size_t q=0;q<closure.size();++q)
            {
                const auto c=eb.fragments[q].centroid,d=closure[q];maximum=std::max(maximum,norm(d));
                if(norm(d)>1e-9)std::printf("[contracts-local] cell=%d centroid=[%.12g %.12g %.12g] defect=[%.12g %.12g %.12g]\n",eb.fragments[q].parent_cell,c.x,c.y,c.z,d.x,d.y,d.z);
            }
            std::printf("[contracts-local] max defect=%.12g m2\n",maximum);
            require(maximum<1e-7*grid.h*grid.h,"local control-volume closure failed");return 0;
        }
        ExternalAeroCore core(wing,bvh,config,options);
        const auto& system=core.pressure_system();const double h=core.hierarchy().finest_cell_size();
        const auto closure=composite_cell_pressure_closure_cpu(system);double maximum_closure=0;
        std::size_t sealed=0;
        for(int q=0;q<system.storage_size;++q)if(system.active[q])
        {
            require(std::isfinite(norm(closure[q])),"nonfinite control-volume closure");
            maximum_closure=std::max(maximum_closure,norm(closure[q]));
            sealed+=!system.freestream_connected[q];
        }
        std::printf("[contracts] geometry maximum area-vector defect=%.12g m2 (%.6g h2), sealed=%zu\n",maximum_closure,maximum_closure/(h*h),sealed);
        if(maximum_closure>=1e-7*h*h)
        {
            std::vector<int> worst;
            for(int q=0;q<system.storage_size;++q)if(system.active[q]&&norm(closure[q])>1e-7*h*h)worst.push_back(q);
            std::sort(worst.begin(),worst.end(),[&](int a,int b){return norm(closure[a])>norm(closure[b]);});
            std::printf("[contracts] defective volumes=%zu\n",worst.size());
            for(std::size_t n=0;n<std::min<std::size_t>(worst.size(),10);++n)
            {
                const int q=worst[n];const auto c=system.centroid[q],d=closure[q];
                std::printf("[contracts] dof=%d centroid=[%.12g %.12g %.12g] volume=%.12g defect=[%.12g %.12g %.12g]\n",q,c.x,c.y,c.z,system.volume[q],d.x,d.y,d.z);
            }
        }
        require(sealed==0,"solid material retained as sealed fluid");
        require(maximum_closure<1e-7*h*h,"fluid control volumes do not close");
        if(geometry_only){std::puts("[contracts] geometry PASS");return 0;}

        DeviceAmrFields fields(core.hierarchy());DeviceCompositeAmrProjection projection(system,fields);
        AmrHostFields zero_fields(core.hierarchy());CompositeAmrFluxes zero_flux;
        zero_flux.coarse_fine_velocity.resize(system.coarse_fine.size(),0);
        zero_flux.embedded_velocity.resize(system.embedded.size(),0);
        // Exercise the actual GPU flux correction, including both regular AMR
        // tiles and apertures, rather than testing only a copied coefficient formula.
        for(int component=-1;component<3;++component)
        {
            std::vector<double> pressure(system.storage_size,3.0);
            for(int q=0;q<system.storage_size;++q)if(component>=0)pressure[q]+=system.centroid[q][component];
            fields.upload(zero_fields);projection.upload_special_fluxes(zero_flux);projection.upload_pressure(pressure);
            projection.correct_fluxes(1,1);CompositeAmrFluxes corrected;projection.download_special_fluxes(corrected);
            double error=0;
            auto check=[&](const auto& edges,const auto& velocity)
            {
                for(std::size_t q=0;q<edges.size();++q)
                {
                    require(std::isfinite(velocity[q]),"nonfinite manufactured pressure correction");
                    error=std::max(error,std::abs(double(velocity[q])+(component==edges[q].axis?1:0)));
                }
            };
            check(system.coarse_fine,corrected.coarse_fine_velocity);check(system.embedded,corrected.embedded_velocity);
            std::printf("[contracts] GPU pressure field %d maximum derivative error=%.9g\n",component,error);
            require(error<1e-4,"production pressure correction is not affine consistent");
        }

        std::vector<double> pressure(system.storage_size,0);
        Vec3d centre{},radius{};
        for(int axis=0;axis<3;++axis)
        {
            centre[axis]=0.5*(wing.bbox_min[axis]+wing.bbox_max[axis]);
            radius[axis]=0.5*(wing.bbox_max[axis]-wing.bbox_min[axis])+0.25;
            require(centre[axis]-radius[axis]>core.hierarchy().domain().lo[axis]&&
                centre[axis]+radius[axis]<core.hierarchy().domain().hi[axis],"pressure fixture reaches outer boundary");
        }
        for(int q=0;q<system.storage_size;++q)if(system.active[q])
        {
            const Vec3d d=system.centroid[q]-centre;double bump=1;
            for(int axis=0;axis<3;++axis)bump*=std::max(0.0,1-d[axis]*d[axis]/(radius[axis]*radius[axis]));
            pressure[q]=Real(100*bump*bump*(1+0.5*d.z/radius.z));
        }
        projection.upload_pressure(pressure);CompositeCellMomentumState cell_zero;
        cell_zero.x.resize(system.storage_size,0);cell_zero.y=cell_zero.z=cell_zero.x;
        DeviceCompositeCellMomentumTransport cells(system,fields,
            conditioning?0.25:config.amr.min_volume_fraction);
        // The pressure solve's face derivative and the cell momentum impulse
        // must both reproduce a linear pressure field. Testing only the former
        // leaves the wall/face traction quadrature unchecked.
        for(int component=0;component<3;++component)
        {
            std::vector<double> linear(system.storage_size,0);
            for(int q=0;q<system.storage_size;++q)linear[q]=3+system.centroid[q][component];
            projection.upload_pressure(linear);cells.upload_state(cell_zero);
            cells.apply_pressure_impulse(projection.pressure(),Real(1),Real(1),true);
            CompositeCellMomentumState after;cells.download_state(after);
            double error=0;int worst=-1;
            const auto& domain=core.hierarchy().domain();
            for(int q=0;q<system.storage_size;++q)if(system.active[q])
            {
                const auto p=system.centroid[q];bool interior=true;
                for(int axis=0;axis<3;++axis)
                    interior=interior&&p[axis]>domain.lo[axis]+2*h&&p[axis]<domain.hi[axis]-2*h;
                if(!interior)continue;
                const double value=std::max({std::abs(double(after.x[q])+(component==0)),
                    std::abs(double(after.y[q])+(component==1)),
                    std::abs(double(after.z[q])+(component==2))});
                if(value>error){error=value;worst=q;}
            }
            std::printf("[contracts] cell linear pressure axis=%d error=%.9g worst=%d\n",component,error,worst);
            require(error<2e-3,"cell pressure traction is not linear exact");
            if(!planb)
            {
                const auto result=compute_pressure_loads(system,linear,wing.triangle_count(),config.freestream,config.reference,0);
                Vec3d expected{};expected[component]=-8*half.x*half.y*half.z;
                const double body_error=norm(result.pressure_force-expected);
                std::printf("[contracts] linear pressure body force axis=%d error=%.9g N\n",component,body_error);
                require(body_error<1e-6,"linear surface pressure does not integrate to displaced-volume force");
            }
        }
        projection.upload_pressure(pressure);
        cells.upload_state(cell_zero);const Real dt=Real(0.001),rho=Real(config.freestream.rho);
        cells.apply_projected_pressure_gradient(projection.pressure(),dt,rho);
        const auto momentum=cells.momentum();
        const auto loads=compute_pressure_loads(system,pressure,wing.triangle_count(),config.freestream,config.reference,0);
        const Vec3d reaction{-rho*momentum[0]/dt,-rho*momentum[1]/dt,-rho*momentum[2]/dt};
        const double reaction_error=norm(reaction-loads.pressure_force)/std::max(1.0,norm(loads.pressure_force));
        std::printf("[contracts] pressure reaction relative error=%.9g\n",reaction_error);
        require(reaction_error<2e-5,"pressure fluid impulse differs from surface reaction");

        // Exercise the cell-momentum path itself. The original tiny-dt test
        // below only exercised MAC aperture transport, even in --cell-momentum.
        // With no flux and no pressure, any resolved velocity field must survive
        // unchanged; small-cell stabilization must not erase its gradients.
        CompositeCellMomentumState affine=cell_zero;
        for(int q=0;q<system.storage_size;++q)if(system.active[q])
        {
            const Vec3d p=system.centroid[q]-centre;
            affine.x[q]=Real(2+0.3*p.y-0.2*p.z);
            affine.y[q]=Real(-0.4*p.x+0.1*p.z);
            affine.z[q]=Real(0.2*p.x+0.3*p.y);
        }
        fields.upload(zero_fields);
        projection.upload_special_fluxes(zero_flux);
        projection.upload_pressure(std::vector<double>(system.storage_size,0));
        auto cell_change=[&]()
        {
            CompositeCellMomentumState after;cells.download_state(after);double error=0;
            for(int q=0;q<system.storage_size;++q)if(system.active[q])
                error=std::max({error,std::abs(double(after.x[q])-affine.x[q]),
                    std::abs(double(after.y[q])-affine.y[q]),std::abs(double(after.z[q])-affine.z[q])});
            return error;
        };
        double no_op_error=0;
        for(const char* operation:{"zero-flux transport","zero-pressure impulse","tiny-dt diffusion","tiny-dt wall"})
        {
            cells.upload_state(affine);
            if(std::string(operation)=="zero-flux transport")
                cells.step(projection.coarse_fine_velocity_device(),projection.embedded_velocity_device(),Real(1e-12),false);
            else if(std::string(operation)=="zero-pressure impulse")
                cells.apply_projected_pressure_gradient(projection.pressure(),Real(1e-12),rho);
            else if(std::string(operation)=="tiny-dt diffusion")cells.diffuse(Real(1.5e-5),Real(1e-12));
            else cells.apply_smooth_fabric_wall_model(Real(1e-12),Real(1.5e-5));
            const double change=cell_change();no_op_error=std::max(no_op_error,change);
            std::printf("[contracts] cell %s affine change=%.9g m/s\n",operation,change);
        }
        require(no_op_error<1e-5,"cell stabilization erases velocity gradients without a physical impulse");

        auto expected_velocity=[&](Vec3d p)
        {
            p=p-centre;return Vec3d{2+0.3*p.y-0.2*p.z,-0.4*p.x+0.1*p.z,0.2*p.x+0.3*p.y};
        };
        cells.upload_state(affine);
        cells.reconstruct_fluxes(projection.coarse_fine_velocity_device(),projection.embedded_velocity_device());
        CompositeAmrFluxes reconstructed;projection.download_special_fluxes(reconstructed);
        auto flux_error=[&](const CompositeAmrFluxes& flux)
        {
            double error=0;
            auto check_edges=[&](const auto& edges,const auto& values)
            {
                for(std::size_t edge=0;edge<edges.size();++edge)
                {
                    const double exact=expected_velocity(edges[edge].face_centroid)[edges[edge].axis];
                    require(std::isfinite(values[edge]),"nonfinite affine reconstructed flux");
                    error=std::max(error,std::abs(values[edge]-exact));
                }
            };
            check_edges(system.coarse_fine,flux.coarse_fine_velocity);
            check_edges(system.embedded,flux.embedded_velocity);return error;
        };
        const double gpu_flux_error=flux_error(reconstructed);
        reconstruct_composite_cell_fluxes_cpu(system,affine,zero_fields,reconstructed);
        const double cpu_flux_error=flux_error(reconstructed);
        std::printf("[contracts] affine velocity reconstruction error GPU=%.9g CPU=%.9g m/s\n",gpu_flux_error,cpu_flux_error);
        require(std::max(gpu_flux_error,cpu_flux_error)<1e-5,"cell-to-face reconstruction destroys affine velocity");

        require(core.initialize().pressure.converged,"production initialization did not converge");
        const int evolved_steps=cell_momentum?64:8;
        for(int stage=0;stage<2;++stage)
        {
            if(stage)for(int step=0;step<evolved_steps;++step)require(core.step().pressure.converged,"production timestep did not converge");
            AmrHostFields before(core.hierarchy());core.download_fields(before);
            CompositeAmrFluxes initial;core.download_special_fluxes(initial);
            for(double tiny_dt:{1e-6,1e-9,1e-12})
            {
                fields.upload(before);projection.upload_special_fluxes(initial);
                projection.transport_embedded_apertures(Real(tiny_dt),0,0,false);
                AmrHostFields after(core.hierarchy());fields.download(after);
                CompositeAmrFluxes after_flux;projection.download_special_fluxes(after_flux);
                double change=field_difference(before,after);
                for(std::size_t q=0;q<initial.embedded_velocity.size();++q)
                {
                    require(std::isfinite(after_flux.embedded_velocity[q]),"nonfinite zero-time transport");
                    change=std::max(change,std::abs(double(initial.embedded_velocity[q])-double(after_flux.embedded_velocity[q])));
                }
                std::printf("[contracts] evolved=%d dt=%.3g velocity change=%.9g m/s\n",stage*evolved_steps,tiny_dt,change);
                if(tiny_dt==1e-12)require(change<1e-5,"transport applies a finite zero-time filter");
            }
        }
        std::puts("[contracts] production physics contracts PASS");return 0;
    }
    catch(const std::exception& error)
    {
        std::fprintf(stderr,"[contracts] FAIL: %s\n",error.what());return 1;
    }
}
