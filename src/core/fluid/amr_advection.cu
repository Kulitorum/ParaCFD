#include "core/fluid/amr_advection.h"
#include "core/fluid/amr_pressure.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace paracfd::core
{
	namespace
	{
		void check(cudaError_t error,const char* operation){if(error!=cudaSuccess)throw std::runtime_error(std::string(operation)+": "+cudaGetErrorString(error));}
		template<class T>T* allocate(std::size_t count,const char* operation){if(!count)return nullptr;T* pointer=nullptr;check(cudaMalloc(&pointer,count*sizeof(T)),operation);return pointer;}
		template<class T>T* upload(const std::vector<T>& host,const char* operation){T* pointer=allocate<T>(host.size(),operation);if(pointer)check(cudaMemcpy(pointer,host.data(),host.size()*sizeof(T),cudaMemcpyHostToDevice),operation);return pointer;}

		__global__ void accumulate_pairwise_momentum_kernel(const int* a,const int* b,
			const Real* area,const Real* normal_velocity,const Real* velocity_x,
			const Real* velocity_y,const Real* velocity_z,Real dt,Real* delta_x,
			Real* delta_y,Real* delta_z,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;
			const int first=a[edge],second=b[edge];const Real transport=dt*area[edge]*normal_velocity[edge];
			const int donor=transport>=Real(0)?first:second;
			const Real mx=transport*velocity_x[donor],my=transport*velocity_y[donor],mz=transport*velocity_z[donor];
			atomicAdd(delta_x+first,-mx);atomicAdd(delta_x+second,mx);
			atomicAdd(delta_y+first,-my);atomicAdd(delta_y+second,my);
			atomicAdd(delta_z+first,-mz);atomicAdd(delta_z+second,mz);
		}
		__global__ void apply_pairwise_momentum_kernel(const Real* volume,const Real* delta_x,
			const Real* delta_y,const Real* delta_z,Real* velocity_x,Real* velocity_y,
			Real* velocity_z,int count)
		{
			const int node=blockIdx.x*blockDim.x+threadIdx.x;if(node>=count)return;const Real inverse=Real(1)/volume[node];
			velocity_x[node]+=delta_x[node]*inverse;velocity_y[node]+=delta_y[node]*inverse;velocity_z[node]+=delta_z[node]*inverse;
		}
		__global__ void accumulate_pairwise_scalar_kernel(const int* a,const int* b,
			const Real* area,const Real* normal_velocity,const Real* velocity,Real dt,
			Real* delta,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int first=a[edge],second=b[edge];const Real transport=dt*area[edge]*normal_velocity[edge];const Real momentum=transport*velocity[transport>=Real(0)?first:second];atomicAdd(delta+first,-momentum);atomicAdd(delta+second,momentum);
		}
		__global__ void apply_pairwise_scalar_kernel(const Real* volume,const Real* delta,
			Real* velocity,int count)
		{
			const int node=blockIdx.x*blockDim.x+threadIdx.x;if(node<count)velocity[node]+=delta[node]/volume[node];
		}
		__device__ void accumulate_cell_momentum_flux(int a,int b,Real swept,
			const Real* x,const Real* y,const Real* z,Real* dx,Real* dy,Real* dz)
		{
			const int donor=swept>=Real(0)?a:b;const Real mx=swept*x[donor],my=swept*y[donor],mz=swept*z[donor];atomicAdd(dx+a,-mx);atomicAdd(dx+b,mx);atomicAdd(dy+a,-my);atomicAdd(dy+b,my);atomicAdd(dz+a,-mz);atomicAdd(dz+b,mz);
		}
		__global__ void transport_regular_cell_momentum_kernel(DeviceAmrFieldLevelView view,
			Real h,int level_offset,const unsigned char* active,const unsigned char* cut_face_mask,
			const unsigned char* compact_plus_mask,const Real* x,const Real* y,const Real* z,
			Real dt,bool external_aero,Real freestream_speed,Real* dx,Real* dy,Real* dz)
		{
			const int bs=view.layout.brick_size,cells=bs*bs*bs,q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*cells;if(q>=total)return;const int brick=q/cells,local=q%cells;if(view.flags[brick]&BRICK_COVERED)return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs),a=level_offset+q;if(!active[a])return;
			for(int axis=0;axis<3;++axis)
			{
				if((cut_face_mask[a]&(1u<<axis))||(compact_plus_mask[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=view.neighbors[brick*6+2*axis+1];if(other_brick<0||(view.flags[other_brick]&BRICK_COVERED))continue;c[axis]=0;}const int b=level_offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(!active[b])continue;Real normal=0;if(axis==0)normal=view.u[view.layout.u_index(brick,i+1,j,k)];else if(axis==1)normal=view.v[view.layout.v_index(brick,i,j+1,k)];else normal=view.w[view.layout.w_index(brick,i,j,k+1)];accumulate_cell_momentum_flux(a,b,dt*h*h*normal,x,y,z,dx,dy,dz);
			}
			if(!external_aero)return;const Real area=h*h;if(i==0&&(view.flags[brick]&BRICK_XMIN)){const Real normal=view.u[view.layout.u_index(brick,0,j,k)],swept=dt*area*normal;const Real ex=freestream_speed,ey=Real(0),ez=Real(0);if(swept>=Real(0)){atomicAdd(dx+a,swept*ex);atomicAdd(dy+a,swept*ey);atomicAdd(dz+a,swept*ez);}else{atomicAdd(dx+a,swept*x[a]);atomicAdd(dy+a,swept*y[a]);atomicAdd(dz+a,swept*z[a]);}}if(i==bs-1&&(view.flags[brick]&BRICK_XMAX)){const Real normal=view.u[view.layout.u_index(brick,bs,j,k)],swept=dt*area*normal;if(swept>=Real(0)){atomicAdd(dx+a,-swept*x[a]);atomicAdd(dy+a,-swept*y[a]);atomicAdd(dz+a,-swept*z[a]);}else{atomicAdd(dx+a,-swept*freestream_speed);}}
		}
		__global__ void transport_compact_cell_momentum_kernel(const int* a,const int* b,
			const Real* area,const Real* normal_velocity,const Real* x,const Real* y,const Real* z,
			Real dt,Real* dx,Real* dy,Real* dz,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count)accumulate_cell_momentum_flux(a[edge],b[edge],dt*area[edge]*normal_velocity[edge],x,y,z,dx,dy,dz);
		}
		__global__ void apply_cell_momentum_kernel(const unsigned char* active,const Real* volume,
			const Real* dx,const Real* dy,const Real* dz,Real* x,Real* y,Real* z,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<count&&active[q]){const Real inverse=Real(1)/volume[q];x[q]+=dx[q]*inverse;y[q]+=dy[q]*inverse;z[q]+=dz[q]*inverse;}
		}
		__device__ void accumulate_cell_momentum_diffusion(int a,int b,Real scale,
			const Real* x,const Real* y,const Real* z,Real* dx,Real* dy,Real* dz)
		{
			const Real mx=scale*(x[b]-x[a]),my=scale*(y[b]-y[a]),mz=scale*(z[b]-z[a]);atomicAdd(dx+a,mx);atomicAdd(dx+b,-mx);atomicAdd(dy+a,my);atomicAdd(dy+b,-my);atomicAdd(dz+a,mz);atomicAdd(dz+b,-mz);
		}
		__global__ void diffuse_regular_cell_momentum_kernel(DeviceAmrFieldLevelView view,
			Real h,int level_offset,const unsigned char* active,const unsigned char* cut_face_mask,
			const unsigned char* compact_plus_mask,const Real* x,const Real* y,const Real* z,
			Real viscosity_dt,Real* dx,Real* dy,Real* dz)
		{
			const int bs=view.layout.brick_size,cells=bs*bs*bs,q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*cells;if(q>=total)return;const int brick=q/cells,local=q%cells;if(view.flags[brick]&BRICK_COVERED)return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs),a=level_offset+q;if(!active[a])return;for(int axis=0;axis<3;++axis){if((cut_face_mask[a]&(1u<<axis))||(compact_plus_mask[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=view.neighbors[brick*6+2*axis+1];if(other_brick<0||(view.flags[other_brick]&BRICK_COVERED))continue;c[axis]=0;}const int b=level_offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(active[b])accumulate_cell_momentum_diffusion(a,b,viscosity_dt*h,x,y,z,dx,dy,dz);}
		}
		__global__ void diffuse_compact_cell_momentum_kernel(const int* a,const int* b,
			const Real* conductance,const Real* x,const Real* y,const Real* z,Real viscosity_dt,
			Real* dx,Real* dy,Real* dz,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count)accumulate_cell_momentum_diffusion(a[edge],b[edge],viscosity_dt*conductance[edge],x,y,z,dx,dy,dz);
		}
		__device__ Real cell_component(const Real* x,const Real* y,const Real* z,int dof,int axis);
		__device__ int structured_cell_neighbor(const DeviceAmrFieldLevelView& view,
			int level_offset,const unsigned char* active,const unsigned char* cut_face_mask,
			const unsigned char* compact_plus_mask,int brick,int i,int j,int k,int axis,int sign)
		{
			const int bs=view.layout.brick_size,cells=bs*bs*bs;int c[3]={i,j,k},other_brick=brick;c[axis]+=sign;if(c[axis]<0||c[axis]>=bs){other_brick=view.neighbors[brick*6+2*axis+(sign>0)];if(other_brick<0||(view.flags[other_brick]&BRICK_COVERED))return -1;c[axis]=sign>0?0:bs-1;}const int node=level_offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(!active[node])return -1;const int lower=sign>0?level_offset+brick*cells+i+bs*(j+bs*k):node;if((cut_face_mask[lower]&(1u<<axis))||(compact_plus_mask[lower]&(1u<<axis)))return -1;return node;
		}
		__global__ void compute_regular_cell_viscosity_kernel(DeviceAmrFieldLevelView view,
			Real h,int level_offset,const unsigned char* active,const unsigned char* cut_face_mask,
			const unsigned char* compact_plus_mask,const unsigned char* special_mask,
			const Real* volume,const Real* x,const Real* y,const Real* z,Real molecular_viscosity,
			Real smagorinsky_cs,Real* viscosity)
		{
			const int bs=view.layout.brick_size,cells=bs*bs*bs,q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*cells;if(q>=total)return;const int brick=q/cells,local=q%cells;if(view.flags[brick]&BRICK_COVERED)return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs),a=level_offset+q;if(!active[a])return;viscosity[a]=molecular_viscosity;if(special_mask[a]||smagorinsky_cs==Real(0))return;Real gradient[9];for(int component=0;component<3;++component)for(int axis=0;axis<3;++axis){const int lower=structured_cell_neighbor(view,level_offset,active,cut_face_mask,compact_plus_mask,brick,i,j,k,axis,-1),upper=structured_cell_neighbor(view,level_offset,active,cut_face_mask,compact_plus_mask,brick,i,j,k,axis,1);const Real current=cell_component(x,y,z,a,component);if(lower>=0&&upper>=0)gradient[component*3+axis]=(cell_component(x,y,z,upper,component)-cell_component(x,y,z,lower,component))/(Real(2)*h);else if(upper>=0)gradient[component*3+axis]=(cell_component(x,y,z,upper,component)-current)/h;else if(lower>=0)gradient[component*3+axis]=(current-cell_component(x,y,z,lower,component))/h;else gradient[component*3+axis]=Real(0);}const Real strain=sqrt(max(Real(0),detail::compact_strain_magnitude_squared(gradient))),length=smagorinsky_cs*cbrt(volume[a]);viscosity[a]+=length*length*strain;
		}
		__global__ void accumulate_special_cell_gradient_kernel(const int* special_index,
			const int* special_dof,const int* neighbor,const Real* weighted_displacement,const Real* x,const Real* y,
			const Real* z,Real* rhs,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int special=special_index[edge],node=special_dof[special],other=neighbor[edge];for(int component=0;component<3;++component){const Real delta=cell_component(x,y,z,other,component)-cell_component(x,y,z,node,component);for(int axis=0;axis<3;++axis)atomicAdd(rhs+(special*9+component*3+axis),weighted_displacement[3*edge+axis]*delta);}
		}
		__global__ void solve_special_cell_viscosity_kernel(const int* dof,
			const Real* inverse,const Real* rhs,const Real* volume,Real molecular_viscosity,
			Real smagorinsky_cs,Real* viscosity,int count)
		{
			const int special=blockIdx.x*blockDim.x+threadIdx.x;if(special>=count)return;Real gradient[9];for(int component=0;component<3;++component)for(int row=0;row<3;++row){Real value=Real(0);for(int column=0;column<3;++column)value+=inverse[special*9+row*3+column]*rhs[special*9+component*3+column];gradient[component*3+row]=value;}const int node=dof[special];const Real strain=sqrt(max(Real(0),detail::compact_strain_magnitude_squared(gradient))),length=smagorinsky_cs*cbrt(volume[node]);viscosity[node]=molecular_viscosity+length*length*strain;
		}
		__global__ void diffuse_regular_cell_momentum_variable_kernel(DeviceAmrFieldLevelView view,
			Real h,int level_offset,const unsigned char* active,const unsigned char* cut_face_mask,
			const unsigned char* compact_plus_mask,const Real* viscosity,const Real* x,const Real* y,
			const Real* z,Real dt,Real* dx,Real* dy,Real* dz)
		{
			const int bs=view.layout.brick_size,cells=bs*bs*bs,q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*cells;if(q>=total)return;const int brick=q/cells,local=q%cells;if(view.flags[brick]&BRICK_COVERED)return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs),a=level_offset+q;if(!active[a])return;for(int axis=0;axis<3;++axis){if((cut_face_mask[a]&(1u<<axis))||(compact_plus_mask[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=view.neighbors[brick*6+2*axis+1];if(other_brick<0||(view.flags[other_brick]&BRICK_COVERED))continue;c[axis]=0;}const int b=level_offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(active[b])accumulate_cell_momentum_diffusion(a,b,dt*h*Real(0.5)*(viscosity[a]+viscosity[b]),x,y,z,dx,dy,dz);}
		}
		__global__ void diffuse_compact_cell_momentum_variable_kernel(const int* a,const int* b,
			const Real* conductance,const Real* viscosity,const Real* x,const Real* y,const Real* z,
			Real dt,Real* dx,Real* dy,Real* dz,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count)accumulate_cell_momentum_diffusion(a[edge],b[edge],dt*conductance[edge]*Real(0.5)*(viscosity[a[edge]]+viscosity[b[edge]]),x,y,z,dx,dy,dz);
		}
		__device__ Real cell_component(const Real* x,const Real* y,const Real* z,int dof,int axis)
		{
			return axis==0?x[dof]:(axis==1?y[dof]:z[dof]);
		}
		__device__ void set_positive_cell_face(const DeviceAmrFieldLevelView& view,int brick,
			int axis,int i,int j,int k,Real value)
		{
			if(axis==0)view.u[view.layout.u_index(brick,i+1,j,k)]=value;else if(axis==1)view.v[view.layout.v_index(brick,i,j+1,k)]=value;else view.w[view.layout.w_index(brick,i,j,k+1)]=value;
		}
		__global__ void reconstruct_regular_cell_flux_kernel(DeviceAmrFieldLevelView view,
			int level_offset,const unsigned char* active,const unsigned char* cut_face_mask,
			const unsigned char* compact_plus_mask,const Real* x,const Real* y,const Real* z)
		{
			const int bs=view.layout.brick_size,cells=bs*bs*bs,q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*cells;if(q>=total)return;const int brick=q/cells,local=q%cells;if(view.flags[brick]&BRICK_COVERED)return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs),a=level_offset+q;if(!active[a])return;for(int axis=0;axis<3;++axis){if((cut_face_mask[a]&(1u<<axis))||(compact_plus_mask[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=view.neighbors[brick*6+2*axis+1];if(other_brick<0||(view.flags[other_brick]&BRICK_COVERED))continue;c[axis]=0;}const int b=level_offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(!active[b])continue;const Real value=Real(0.5)*(cell_component(x,y,z,a,axis)+cell_component(x,y,z,b,axis));set_positive_cell_face(view,brick,axis,i,j,k,value);if(other_brick!=brick){int lower[3]={c[0],c[1],c[2]};lower[axis]=-1;set_positive_cell_face(view,other_brick,axis,lower[0],lower[1],lower[2],value);}}
		}
		__global__ void reconstruct_compact_cell_flux_kernel(const int* a,const int* b,
			const std::int8_t* axis,const Real* upper_weight,const Real* x,const Real* y,
			const Real* z,Real* velocity,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int component=axis[edge];const Real weight=upper_weight[edge];velocity[edge]=(Real(1)-weight)*cell_component(x,y,z,a[edge],component)+weight*cell_component(x,y,z,b[edge],component);
		}
		__device__ void accumulate_axis_pressure_impulse(int a,int b,int axis,Real impulse,
			Real* dx,Real* dy,Real* dz)
		{
			Real* delta=axis==0?dx:(axis==1?dy:dz);atomicAdd(delta+a,-impulse);atomicAdd(delta+b,impulse);
		}
		__global__ void pressure_regular_cell_impulse_kernel(DeviceAmrFieldLevelView view,
			Real h,int level_offset,const unsigned char* active,const unsigned char* cut_face_mask,
			const unsigned char* compact_plus_mask,const Real* pressure,Real scale,bool outlet_xmax,
			bool include_physical_boundaries,Real* dx,Real* dy,Real* dz)
		{
			const int bs=view.layout.brick_size,cells=bs*bs*bs,q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*cells;if(q>=total)return;const int brick=q/cells,local=q%cells;if(view.flags[brick]&BRICK_COVERED)return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs),a=level_offset+q;if(!active[a])return;const Real area=h*h;
			for(int axis=0;axis<3;++axis){if((cut_face_mask[a]&(1u<<axis))||(compact_plus_mask[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=view.neighbors[brick*6+2*axis+1];if(other_brick<0||(view.flags[other_brick]&BRICK_COVERED))continue;c[axis]=0;}const int b=level_offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(active[b])accumulate_axis_pressure_impulse(a,b,axis,scale*area*Real(0.5)*(pressure[a]+pressure[b]),dx,dy,dz);}
			if(!include_physical_boundaries)return;const Real boundary=scale*area*pressure[a];if(i==0&&(view.flags[brick]&BRICK_XMIN))atomicAdd(dx+a,boundary);if(i==bs-1&&(view.flags[brick]&BRICK_XMAX)&&!outlet_xmax)atomicAdd(dx+a,-boundary);if(j==0&&(view.flags[brick]&BRICK_YMIN))atomicAdd(dy+a,boundary);if(j==bs-1&&(view.flags[brick]&BRICK_YMAX))atomicAdd(dy+a,-boundary);if(k==0&&(view.flags[brick]&BRICK_ZMIN))atomicAdd(dz+a,boundary);if(k==bs-1&&(view.flags[brick]&BRICK_ZMAX))atomicAdd(dz+a,-boundary);
		}
		__global__ void pressure_compact_cell_impulse_kernel(const int* a,const int* b,
			const std::int8_t* axis,const Real* area,const Real* upper_weight,const Real* pressure,
			Real scale,Real* dx,Real* dy,Real* dz,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const Real weight=upper_weight[edge],face_pressure=(Real(1)-weight)*pressure[a[edge]]+weight*pressure[b[edge]];accumulate_axis_pressure_impulse(a[edge],b[edge],axis[edge],scale*area[edge]*face_pressure,dx,dy,dz);
		}
		__global__ void pressure_surface_cell_impulse_kernel(const int* dof,
			const Real* coefficient,const Real* pressure,Real scale,Real* dx,Real* dy,Real* dz,
			int count)
		{
			const int patch=blockIdx.x*blockDim.x+threadIdx.x;if(patch>=count)return;const int node=dof[patch];const Real impulse=scale*pressure[node];atomicAdd(dx+node,impulse*coefficient[3*patch]);atomicAdd(dy+node,impulse*coefficient[3*patch+1]);atomicAdd(dz+node,impulse*coefficient[3*patch+2]);
		}
		__global__ void gather_mac_face_map_kernel(const DeviceAmrFieldLevelView* levels,
			const int* level,const std::uint64_t* index,const std::int8_t* component,
			Real* compact,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;const DeviceAmrFieldLevelView view=levels[level[q]];compact[q]=component[q]==0?view.u[index[q]]:(component[q]==1?view.v[index[q]]:view.w[index[q]]);
		}
		__global__ void scatter_mac_face_map_kernel(const DeviceAmrFieldLevelView* levels,
			const int* level,const std::uint64_t* index,const std::int8_t* component,
			const Real* compact,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;const DeviceAmrFieldLevelView view=levels[level[q]];if(component[q]==0)view.u[index[q]]=compact[q];else if(component[q]==1)view.v[index[q]]=compact[q];else view.w[index[q]]=compact[q];
		}
		__global__ void accumulate_momentum_alias_kernel(const int* node,const int* group,
			const Real* area,const Real* state,Real* sum,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<count)atomicAdd(sum+group[q],area[q]*state[node[q]]);
		}
		__global__ void normalize_momentum_alias_kernel(const Real* sum,const Real* area,
			Real* value,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<count)value[q]=area[q]>Real(0)?sum[q]/area[q]:Real(0);
		}
		__global__ void initialize_momentum_connection_velocity_kernel(const int* lower,const int* upper,
			const std::int8_t* axis,const std::int8_t* component,const Real* node_velocity,
			Real* connection_velocity,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count)connection_velocity[edge]=component[edge]==axis[edge]?Real(0.5)*(node_velocity[lower[edge]]+node_velocity[upper[edge]]):Real(0);
		}
		__global__ void accumulate_momentum_mass_flux_kernel(const int* connection,const int* source,
			const Real* area,const Real* coarse_fine_velocity,Real* connection_velocity,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<count)atomicAdd(connection_velocity+connection[q],area[q]*coarse_fine_velocity[source[q]]);
		}
		__global__ void normalize_tangential_momentum_velocity_kernel(const std::int8_t* axis,
			const std::int8_t* component,const Real* area,Real* velocity,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count&&component[edge]!=axis[edge])velocity[edge]/=area[edge];
		}
		__global__ void accumulate_eb_momentum_carrier_state_kernel(const int* node,const int* carrier,
			const std::int8_t* component,const Real* half_mass,const Real* carrier_state,
			Real* sum,Real* weight,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;const int slot=node[q]*3+component[q],source=carrier[q];const Real mass=half_mass[q];atomicAdd(sum+slot,mass*carrier_state[source]);atomicAdd(weight+slot,mass);
		}
		__global__ void accumulate_eb_momentum_aperture_state_kernel(const int* a,const int* b,
			const std::int8_t* axis,const Real* mass,const Real* velocity,Real* sum,Real* weight,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int component=axis[edge],sa=a[edge]*3+component,sb=b[edge]*3+component;const Real half=Real(0.5)*mass[edge],momentum=half*velocity[edge];atomicAdd(sum+sa,momentum);atomicAdd(sum+sb,momentum);atomicAdd(weight+sa,half);atomicAdd(weight+sb,half);
		}
		__device__ void accumulate_eb_momentum_transfer(int a,int b,Real swept,const Real* sum,
			const Real* weight,Real* delta)
		{
			const int donor=swept>=Real(0)?a:b;for(int component=0;component<3;++component){const int source=donor*3+component,sa=a*3+component,sb=b*3+component;if(!(weight[source]>Real(0)&&weight[sa]>Real(0)&&weight[sb]>Real(0)))continue;const Real transfer=swept*sum[source]/weight[source];atomicAdd(delta+sa,-transfer);atomicAdd(delta+sb,transfer);}
		}
		__global__ void transport_eb_momentum_aperture_kernel(const int* a,const int* b,
			const Real* area,const Real* normal_velocity,const Real* sum,const Real* weight,
			Real dt,Real* delta,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count)accumulate_eb_momentum_transfer(a[edge],b[edge],dt*area[edge]*normal_velocity[edge],sum,weight,delta);
		}
		__global__ void transport_eb_momentum_regular_kernel(const int* a,const int* b,
			const int* normal_carrier,const Real* area,const Real* carrier_state,const Real* sum,
			const Real* weight,Real dt,Real* delta,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count)accumulate_eb_momentum_transfer(a[edge],b[edge],dt*area[edge]*carrier_state[normal_carrier[edge]],sum,weight,delta);
		}
		__global__ void scatter_eb_momentum_carrier_delta_kernel(const int* node,const int* carrier,
			const std::int8_t* component,const Real* weight,const Real* node_delta,
			Real* carrier_delta,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;const int slot=node[q]*3+component[q];if(weight[slot]>Real(0))atomicAdd(carrier_delta+carrier[q],Real(0.5)*node_delta[slot]/weight[slot]);
		}
		__global__ void apply_eb_momentum_carrier_delta_kernel(Real* state,const Real* delta,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<count)state[q]+=delta[q];
		}
		__global__ void gather_eb_momentum_alias_kernel(const int* source,const Real* state,Real* alias,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<count)alias[q]=state[source[q]];
		}
		__global__ void apply_eb_momentum_aperture_delta_kernel(const int* a,const int* b,
			const std::int8_t* axis,const Real* weight,const Real* node_delta,Real* velocity,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int component=axis[edge],sa=a[edge]*3+component,sb=b[edge]*3+component;const Real da=weight[sa]>Real(0)?node_delta[sa]/weight[sa]:Real(0),db=weight[sb]>Real(0)?node_delta[sb]/weight[sb]:Real(0);velocity[edge]+=Real(0.5)*(da+db);
		}

		__host__ __device__ std::uint64_t coordinate_hash(int x,int y,int z)
		{
			std::uint64_t key=(static_cast<std::uint64_t>(static_cast<std::uint32_t>(x))<<42)^(static_cast<std::uint64_t>(static_cast<std::uint32_t>(y))<<21)^static_cast<std::uint32_t>(z);std::uint64_t q=key+0x9e3779b97f4a7c15ull;q=(q^(q>>30))*0xbf58476d1ce4e5b9ull;q=(q^(q>>27))*0x94d049bb133111ebull;return q^(q>>31);
		}
		__device__ GpuBrickLocation locate_finest(GpuAmrHierarchyView hierarchy,GpuAmrPoint point)
		{
			GpuBrickLocation out;if(point.x<hierarchy.domain_lo.x||point.y<hierarchy.domain_lo.y||point.z<hierarchy.domain_lo.z||point.x>=hierarchy.domain_hi.x||point.y>=hierarchy.domain_hi.y||point.z>=hierarchy.domain_hi.z)return out;
			for(int level=hierarchy.level_count-1;level>=0;--level){const GpuAmrLevelView lev=hierarchy.levels[level];const float width=hierarchy.brick_size*lev.h;const int bx=static_cast<int>(floorf((point.x-hierarchy.domain_lo.x)/width)),by=static_cast<int>(floorf((point.y-hierarchy.domain_lo.y)/width)),bz=static_cast<int>(floorf((point.z-hierarchy.domain_lo.z)/width));std::uint32_t slot=static_cast<std::uint32_t>(coordinate_hash(bx,by,bz))&lev.lookup_mask;int brick=-1;for(int probe=0;probe<lev.lookup_size;++probe){const GpuBrickLookupEntry entry=lev.lookup[slot];if(entry.brick<0)break;if(entry.x==bx&&entry.y==by&&entry.z==bz){brick=entry.brick;break;}slot=(slot+1)&lev.lookup_mask;}if(brick<0||brick>=lev.brick_count||(lev.bricks[brick].flags&BRICK_COVERED))continue;const GpuBrickRecord record=lev.bricks[brick];out.level=level;out.brick=brick;out.i=max(0,min(hierarchy.brick_size-1,static_cast<int>(floorf((point.x-record.origin_x)/record.h))));out.j=max(0,min(hierarchy.brick_size-1,static_cast<int>(floorf((point.y-record.origin_y)/record.h))));out.k=max(0,min(hierarchy.brick_size-1,static_cast<int>(floorf((point.z-record.origin_z)/record.h))));return out;}return out;
		}

		__device__ bool valid_face_index(int component,int bs,int i,int j,int k)
		{
			return i>=0&&j>=0&&k>=0&&i<(component==0?bs+1:bs)&&j<(component==1?bs+1:bs)&&k<(component==2?bs+1:bs);
		}
		__device__ Real value_at_clamped(const DeviceAmrFieldLevelView& view,int component,int brick,int i,int j,int k)
		{
			const int bs=view.layout.brick_size;
			if(component==0){i=max(0,min(bs,i));j=max(0,min(bs-1,j));k=max(0,min(bs-1,k));return view.u[view.layout.u_index(brick,i,j,k)];}
			if(component==1){i=max(0,min(bs-1,i));j=max(0,min(bs,j));k=max(0,min(bs-1,k));return view.v[view.layout.v_index(brick,i,j,k)];}
			i=max(0,min(bs-1,i));j=max(0,min(bs-1,j));k=max(0,min(bs,k));return view.w[view.layout.w_index(brick,i,j,k)];
		}
		__device__ Real value_at_extrapolated(const DeviceAmrFieldLevelView& view,int component,int brick,int i,int j,int k)
		{
			const int bs=view.layout.brick_size,n[3]={component==0?bs+1:bs,component==1?bs+1:bs,component==2?bs+1:bs},requested[3]={i,j,k};int index[3][2],count[3];Real weight[3][2];
			for(int axis=0;axis<3;++axis)
			{
				if(requested[axis]<0){index[axis][0]=0;index[axis][1]=1;weight[axis][0]=Real(2);weight[axis][1]=Real(-1);count[axis]=2;}
				else if(requested[axis]>=n[axis]){index[axis][0]=n[axis]-1;index[axis][1]=n[axis]-2;weight[axis][0]=Real(2);weight[axis][1]=Real(-1);count[axis]=2;}
				else{index[axis][0]=requested[axis];weight[axis][0]=Real(1);count[axis]=1;}
			}
			Real result=0;for(int az=0;az<count[2];++az)for(int ay=0;ay<count[1];++ay)for(int ax=0;ax<count[0];++ax)result+=weight[0][ax]*weight[1][ay]*weight[2][az]*value_at_clamped(view,component,brick,index[0][ax],index[1][ay],index[2][az]);return result;
		}
		// Prolong a staggered value from another AMR level. Tangential MAC lattices stop
		// half a cell before a brick boundary, so a plain brick-local clamp is only
		// zeroth order. The one-sided extension above preserves globally linear fields.
		__device__ Real sample_local_extrapolated(const DeviceAmrFieldLevelView& view,const GpuBrickRecord& record,int component,int brick,GpuAmrPoint point)
		{
			Real x=(Real(point.x)-Real(record.origin_x))/Real(record.h)-(component==0?Real(0):Real(0.5));
			Real y=(Real(point.y)-Real(record.origin_y))/Real(record.h)-(component==1?Real(0):Real(0.5));
			Real z=(Real(point.z)-Real(record.origin_z))/Real(record.h)-(component==2?Real(0):Real(0.5));
			const int i0=static_cast<int>(floor(x)),j0=static_cast<int>(floor(y)),k0=static_cast<int>(floor(z));x-=i0;y-=j0;z-=k0;Real result=0;
			for(int dz=0;dz<2;++dz)for(int dy=0;dy<2;++dy)for(int dx=0;dx<2;++dx){const Real wx=dx?x:Real(1)-x,wy=dy?y:Real(1)-y,wz=dz?z:Real(1)-z;result+=wx*wy*wz*value_at_extrapolated(view,component,brick,i0+dx,j0+dy,k0+dz);}return result;
		}
		__device__ GpuAmrPoint face_node_point(const GpuBrickRecord& record,int component,int i,int j,int k)
		{
			return {record.origin_x+(i+(component==0?0.0f:0.5f))*record.h,record.origin_y+(j+(component==1?0.0f:0.5f))*record.h,record.origin_z+(k+(component==2?0.0f:0.5f))*record.h};
		}
		__device__ int nearest_lattice_index(float coordinate,float origin,float h,float offset)
		{
			return static_cast<int>(floorf((coordinate-origin)/h-offset+0.5f));
		}
		__device__ __noinline__ Real value_at_node(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int source_level,const GpuBrickRecord& source_record,int component,int source_brick,int i,int j,int k)
		{
			const DeviceAmrFieldLevelView source_view=views[source_level];const int bs=source_view.layout.brick_size;
			if(valid_face_index(component,bs,i,j,k))return value_at_clamped(source_view,component,source_brick,i,j,k);
			// The overwhelmingly common boundary case is one step into an active
			// same-level neighbour. Resolve it from the six compact neighbour IDs;
			// reserve the hierarchy hash for corners and 2:1 transitions.
			const int nx=component==0?bs+1:bs,ny=component==1?bs+1:bs,nz=component==2?bs+1:bs;int face=-1,outside=0;if(i<0){face=0;++outside;}else if(i>=nx){face=1;++outside;}if(j<0){face=2;++outside;}else if(j>=ny){face=3;++outside;}if(k<0){face=4;++outside;}else if(k>=nz){face=5;++outside;}
			if(outside==1){const int neighbor=source_view.neighbors[source_brick*6+face];if(neighbor>=0&&!(source_view.flags[neighbor]&BRICK_COVERED)){if(face==0)i+=bs;else if(face==1)i-=bs;else if(face==2)j+=bs;else if(face==3)j-=bs;else if(face==4)k+=bs;else k-=bs;return value_at_clamped(source_view,component,neighbor,i,j,k);}}
			const GpuAmrPoint node=face_node_point(source_record,component,i,j,k);
			if(node.x<hierarchy.domain_lo.x||node.y<hierarchy.domain_lo.y||node.z<hierarchy.domain_lo.z||node.x>hierarchy.domain_hi.x||node.y>hierarchy.domain_hi.y||node.z>hierarchy.domain_hi.z)return value_at_clamped(source_view,component,source_brick,i,j,k);
			GpuAmrPoint probe=node;const float epsilon=1e-5f*source_record.h;probe.x=fminf(hierarchy.domain_hi.x-epsilon,fmaxf(hierarchy.domain_lo.x+epsilon,probe.x));probe.y=fminf(hierarchy.domain_hi.y-epsilon,fmaxf(hierarchy.domain_lo.y+epsilon,probe.y));probe.z=fminf(hierarchy.domain_hi.z-epsilon,fmaxf(hierarchy.domain_lo.z+epsilon,probe.z));const GpuBrickLocation location=locate_finest(hierarchy,probe);if(location.brick<0)return value_at_clamped(source_view,component,source_brick,i,j,k);const GpuBrickRecord record=hierarchy.levels[location.level].bricks[location.brick];const DeviceAmrFieldLevelView view=views[location.level];if(location.level!=source_level)return sample_local_extrapolated(view,record,component,location.brick,node);const int ni=nearest_lattice_index(node.x,record.origin_x,record.h,component==0?0.0f:0.5f),nj=nearest_lattice_index(node.y,record.origin_y,record.h,component==1?0.0f:0.5f),nk=nearest_lattice_index(node.z,record.origin_z,record.h,component==2?0.0f:0.5f);return value_at_clamped(view,component,location.brick,ni,nj,nk);
		}
		__device__ Real sample_local_bounded(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,const GpuBrickRecord& record,int component,int brick,GpuAmrPoint point,Real* lower,Real* upper)
		{
			Real x=(Real(point.x)-Real(record.origin_x))/Real(record.h)-(component==0?Real(0):Real(0.5));
			Real y=(Real(point.y)-Real(record.origin_y))/Real(record.h)-(component==1?Real(0):Real(0.5));
			Real z=(Real(point.z)-Real(record.origin_z))/Real(record.h)-(component==2?Real(0):Real(0.5));
			const int i0=static_cast<int>(floor(x)),j0=static_cast<int>(floor(y)),k0=static_cast<int>(floor(z));x-=i0;y-=j0;z-=k0;Real result=0,mn=Real(0),mx=Real(0);bool first=true;const int bs=views[level].layout.brick_size;
			for(int dz=0;dz<2;++dz)for(int dy=0;dy<2;++dy)for(int dx=0;dx<2;++dx){const int ni=i0+dx,nj=j0+dy,nk=k0+dz;const Real wx=dx?x:Real(1)-x,wy=dy?y:Real(1)-y,wz=dz?z:Real(1)-z;const Real value=valid_face_index(component,bs,ni,nj,nk)?value_at_clamped(views[level],component,brick,ni,nj,nk):value_at_node(hierarchy,views,level,record,component,brick,ni,nj,nk);result+=wx*wy*wz*value;if(first){mn=mx=value;first=false;}else{mn=min(mn,value);mx=max(mx,value);}}
			if(lower)*lower=mn;if(upper)*upper=mx;return result;
		}
		__device__ Real sample_local(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,const GpuBrickRecord& record,int component,int brick,GpuAmrPoint point)
		{
			return sample_local_bounded(hierarchy,views,level,record,component,brick,point,nullptr,nullptr);
		}
		__device__ GpuAmrPoint clamp_point(GpuAmrHierarchyView hierarchy,GpuAmrPoint point)
		{
			point.x=fminf(hierarchy.domain_hi.x-1e-5f,fmaxf(hierarchy.domain_lo.x+1e-5f,point.x));point.y=fminf(hierarchy.domain_hi.y-1e-5f,fmaxf(hierarchy.domain_lo.y+1e-5f,point.y));point.z=fminf(hierarchy.domain_hi.z-1e-5f,fmaxf(hierarchy.domain_lo.z+1e-5f,point.z));return point;
		}
		__device__ Real sample_finest_bounded(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int component,GpuAmrPoint point,Real fallback,Real* lower,Real* upper)
		{
			point=clamp_point(hierarchy,point);const GpuBrickLocation location=locate_finest(hierarchy,point);if(location.brick<0){if(lower)*lower=fallback;if(upper)*upper=fallback;return fallback;}return sample_local_bounded(hierarchy,views,location.level,hierarchy.levels[location.level].bricks[location.brick],component,location.brick,point,lower,upper);
		}
		__device__ Real sample_finest(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int component,GpuAmrPoint point,Real fallback)
		{
			return sample_finest_bounded(hierarchy,views,component,point,fallback,nullptr,nullptr);
		}
		__device__ GpuAmrPoint rk2_trace(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,const GpuBrickRecord& record,int brick,GpuAmrPoint point,Real dt)
		{
			const Real ux=sample_local(hierarchy,views,level,record,0,brick,point),uy=sample_local(hierarchy,views,level,record,1,brick,point),uz=sample_local(hierarchy,views,level,record,2,brick,point);GpuAmrPoint midpoint=clamp_point(hierarchy,{point.x-static_cast<float>(Real(0.5)*dt*ux),point.y-static_cast<float>(Real(0.5)*dt*uy),point.z-static_cast<float>(Real(0.5)*dt*uz)});const Real mx=sample_finest(hierarchy,views,0,midpoint,ux),my=sample_finest(hierarchy,views,1,midpoint,uy),mz=sample_finest(hierarchy,views,2,midpoint,uz);return clamp_point(hierarchy,{point.x-static_cast<float>(dt*mx),point.y-static_cast<float>(dt*my),point.z-static_cast<float>(dt*mz)});
		}
		__device__ void decode_face_work(const DeviceAmrFieldLevelView& view,int q,int& component,int& brick,int& i,int& j,int& k)
		{
			const int bs=view.layout.brick_size,per_component=bs*bs*(bs+1);int r=q,local=r%per_component;r/=per_component;component=r%3;brick=r/3;i=j=k=0;if(component==0){i=local%(bs+1);local/=bs+1;j=local%bs;k=local/bs;}else if(component==1){i=local%bs;local/=bs;j=local%(bs+1);k=local/(bs+1);}else{i=local%bs;local/=bs;j=local%bs;k=local/bs;}
		}
		__device__ Real minmod(Real a,Real b)
		{
			return a*b<=Real(0)?Real(0):(abs(a)<abs(b)?a:b);
		}
		__device__ bool step_same_level_face_node(const DeviceAmrFieldLevelView& view,int component,int axis,int direction,int& brick,int& i,int& j,int& k)
		{
			int* coordinate=axis==0?&i:(axis==1?&j:&k);*coordinate+=direction;const int bs=view.layout.brick_size;if(valid_face_index(component,bs,i,j,k))return true;
			const int nx=component==0?bs+1:bs,ny=component==1?bs+1:bs,nz=component==2?bs+1:bs;int face=-1,outside=0;if(i<0){face=0;++outside;}else if(i>=nx){face=1;++outside;}if(j<0){face=2;++outside;}else if(j>=ny){face=3;++outside;}if(k<0){face=4;++outside;}else if(k>=nz){face=5;++outside;}if(outside!=1)return false;const int neighbour=view.neighbors[brick*6+face];if(neighbour<0||(view.flags[neighbour]&BRICK_COVERED))return false;if(face==0)i+=bs;else if(face==1)i-=bs;else if(face==2)j+=bs;else if(face==3)j-=bs;else if(face==4)k+=bs;else k-=bs;brick=neighbour;return valid_face_index(component,bs,i,j,k);
		}
		__device__ std::size_t face_index(const DeviceAmrFieldLevelView& view,int component,int brick,int i,int j,int k)
		{
			return component==0?view.layout.u_index(brick,i,j,k):(component==1?view.layout.v_index(brick,i,j,k):view.layout.w_index(brick,i,j,k));
		}
		__device__ Real side_safe_muscl(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,const GpuBrickRecord& record,int component,int brick,int i,int j,int k,Real current,unsigned char links,const unsigned char* component_links,Real dt)
		{
			const DeviceAmrFieldLevelView view=views[level];const GpuAmrPoint point=face_node_point(record,component,i,j,k);const Real velocity[3]={sample_local(hierarchy,views,level,record,0,brick,point),sample_local(hierarchy,views,level,record,1,brick,point),sample_local(hierarchy,views,level,record,2,brick,point)};Real weighted_delta=0,weight_sum=0,lower=current,upper=current;
			for(int axis=0;axis<3;++axis)
			{
				const int upstream_direction=velocity[axis]>=Real(0)?-1:1,upstream_bit=2*axis+(upstream_direction>0),downstream_bit=2*axis+(upstream_direction<0);if(!(links&(1u<<upstream_bit)))continue;
				int upstream_brick=brick,ui=i,uj=j,uk=k;if(!step_same_level_face_node(view,component,axis,upstream_direction,upstream_brick,ui,uj,uk))continue;const Real upstream=value_at_clamped(view,component,upstream_brick,ui,uj,uk);lower=min(lower,upstream);upper=max(upper,upstream);const Real weight=dt*abs(velocity[axis])/Real(record.h);Real delta=weight*(upstream-current);
				// A limited linear reconstruction is legal only when all three additional
				// segments are represented by same-side links. Otherwise retain the robust
				// first-order donor update locally at a fabric edge, physical boundary, or
				// coarse/fine transition.
				const unsigned char upstream_links=component_links[face_index(view,component,upstream_brick,ui,uj,uk)];int upstream2_brick=upstream_brick,u2i=ui,u2j=uj,u2k=uk,downstream_brick=brick,di=i,dj=j,dk=k;const bool second_order=(links&(1u<<downstream_bit))&&(upstream_links&(1u<<upstream_bit))&&step_same_level_face_node(view,component,axis,upstream_direction,upstream2_brick,u2i,u2j,u2k)&&step_same_level_face_node(view,component,axis,-upstream_direction,downstream_brick,di,dj,dk);
				if(second_order){const Real upstream2=value_at_clamped(view,component,upstream2_brick,u2i,u2j,u2k),downstream=value_at_clamped(view,component,downstream_brick,di,dj,dk);lower=min(lower,downstream);upper=max(upper,downstream);const Real upstream_slope=minmod(upstream-upstream2,current-upstream),current_slope=minmod(current-upstream,downstream-current);const Real incoming=upstream+Real(0.5)*upstream_slope,outgoing=current+Real(0.5)*current_slope;delta=-weight*(outgoing-incoming);}
				weighted_delta+=delta;weight_sum+=weight;
			}
			const Real scale=weight_sum>Real(1)?Real(1)/weight_sum:Real(1);return max(lower,min(upper,current+scale*weighted_delta));
		}
		__global__ void advect_forward_level_kernel(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,Real dt,const unsigned char* links_u,const unsigned char* links_v,const unsigned char* links_w,Real* u_out,Real* v_out,Real* w_out)
		{
			const DeviceAmrFieldLevelView view=views[level];const int bs=view.layout.brick_size,per_component=bs*bs*(bs+1),q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*3*per_component;if(q>=total)return;int component,brick,i,j,k;decode_face_work(view,q,component,brick,i,j,k);if(view.flags[brick]&BRICK_COVERED)return;const GpuBrickRecord record=hierarchy.levels[level].bricks[brick];const GpuAmrPoint point=face_node_point(record,component,i,j,k);const std::size_t index=component==0?view.layout.u_index(brick,i,j,k):(component==1?view.layout.v_index(brick,i,j,k):view.layout.w_index(brick,i,j,k));const Real current=component==0?view.u[index]:(component==1?view.v[index]:view.w[index]);const unsigned char* component_links=component==0?links_u:(component==1?links_v:links_w);const unsigned char links=component_links[index];const Real result=(links&0x40)?side_safe_muscl(hierarchy,views,level,record,component,brick,i,j,k,current,links,component_links,dt):sample_finest(hierarchy,views,component,rk2_trace(hierarchy,views,level,record,brick,point,dt),current);if(component==0)u_out[index]=result;else if(component==1)v_out[index]=result;else w_out[index]=result;
		}
		__global__ void advect_correct_level_kernel(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,const DeviceAmrFieldLevelView* forward_views,int level,Real dt,const unsigned char* links_u,const unsigned char* links_v,const unsigned char* links_w,Real* u_out,Real* v_out,Real* w_out)
		{
			const DeviceAmrFieldLevelView view=views[level],forward=forward_views[level];const int bs=view.layout.brick_size,per_component=bs*bs*(bs+1),q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*3*per_component;if(q>=total)return;int component,brick,i,j,k;decode_face_work(view,q,component,brick,i,j,k);if(view.flags[brick]&BRICK_COVERED)return;const GpuBrickRecord record=hierarchy.levels[level].bricks[brick];const GpuAmrPoint point=face_node_point(record,component,i,j,k);const std::size_t index=component==0?view.layout.u_index(brick,i,j,k):(component==1?view.layout.v_index(brick,i,j,k):view.layout.w_index(brick,i,j,k));const Real current=component==0?view.u[index]:(component==1?view.v[index]:view.w[index]),phf=component==0?forward.u[index]:(component==1?forward.v[index]:forward.w[index]);const unsigned char links=component==0?links_u[index]:(component==1?links_v[index]:links_w[index]);Real result=phf;
			if(!(links&0x40))
			{
				Real mn,mx;const GpuAmrPoint back=rk2_trace(hierarchy,views,level,record,brick,point,dt);sample_finest_bounded(hierarchy,views,component,back,current,&mn,&mx);const GpuAmrPoint forward_point=rk2_trace(hierarchy,views,level,record,brick,point,-dt);
				// The reverse MacCormack correction is only consistent when its old and
				// forward samples use the same uniform lattice as this face. Across a 2:1
				// transition the forward hierarchy contains independently evolved values;
				// mixing those lattices in the reverse trace creates an interface extremum
				// even for a globally linear characteristic. The forward RK2/trilinear value
				// is already bounded and linearly exact, so retain it locally at transitions.
				const GpuBrickLocation owner=locate_finest(hierarchy,clamp_point(hierarchy,point)),back_owner=locate_finest(hierarchy,back),forward_owner=locate_finest(hierarchy,forward_point);const bool one_lattice=owner.level==level&&back_owner.level==level&&forward_owner.level==level;
				if(one_lattice){const Real phn=sample_finest(hierarchy,forward_views,component,forward_point,phf);result=max(mn,min(mx,phf+Real(0.5)*(current-phn)));}
			}
			if(component==0)u_out[index]=result;else if(component==1)v_out[index]=result;else w_out[index]=result;
		}
		__device__ Real conservative_link_flux(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,const GpuBrickRecord& record,int component,int brick,int i,int j,int k,Real current,int transport_axis,int direction,const unsigned char* component_links,bool external_aero,Real freestream_speed)
		{
			const DeviceAmrFieldLevelView view=views[level];int neighbour_brick=brick,ni=i,nj=j,nk=k;const int bit=2*transport_axis+(direction>0),opposite_bit=2*transport_axis+(direction<0);const std::size_t index=face_index(view,component,brick,i,j,k);const unsigned char links=component_links[index];const bool has_neighbour=step_same_level_face_node(view,component,transport_axis,direction,neighbour_brick,ni,nj,nk);const int physical_face=2*transport_axis+(direction>0);const bool physical_boundary=external_aero&&!has_neighbour&&(record.flags&(1u<<physical_face));Real neighbour=current;if(has_neighbour){const std::size_t neighbour_index=face_index(view,component,neighbour_brick,ni,nj,nk);const unsigned char neighbour_links=component_links[neighbour_index];if(!(links&(1u<<bit))||!(neighbour_links&(1u<<opposite_bit)))return Real(0);neighbour=value_at_clamped(view,component,neighbour_brick,ni,nj,nk);}else if(!physical_boundary&&!(links&(1u<<bit)))return Real(0);
			const GpuAmrPoint point=face_node_point(record,component,i,j,k);Real advector=sample_local(hierarchy,views,level,record,transport_axis,brick,point);if(has_neighbour){const GpuBrickRecord neighbour_record=hierarchy.levels[level].bricks[neighbour_brick];const GpuAmrPoint neighbour_point=face_node_point(neighbour_record,component,ni,nj,nk);advector=Real(0.5)*(advector+sample_local(hierarchy,views,level,neighbour_record,transport_axis,neighbour_brick,neighbour_point));}else if(physical_boundary){if(transport_axis!=0)return Real(0);if(direction<0){advector=freestream_speed;const Real exterior=component==0?freestream_speed:Real(0);return advector>=Real(0)?advector*exterior:advector*current;}return advector*current;}const Real lower=direction>0?current:neighbour,upper=direction>0?neighbour:current;return advector*(advector>=Real(0)?lower:upper);
		}
		__global__ void advect_conservative_uniform_kernel(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,Real dt,const unsigned char* links_u,const unsigned char* links_v,const unsigned char* links_w,Real* u_out,Real* v_out,Real* w_out,bool external_aero,Real freestream_speed)
		{
			const int level=0;const DeviceAmrFieldLevelView view=views[level];const int bs=view.layout.brick_size,per_component=bs*bs*(bs+1),q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*3*per_component;if(q>=total)return;int component,brick,i,j,k;decode_face_work(view,q,component,brick,i,j,k);if(view.flags[brick]&BRICK_COVERED)return;const GpuBrickRecord record=hierarchy.levels[level].bricks[brick];const std::size_t index=face_index(view,component,brick,i,j,k);const Real current=component==0?view.u[index]:(component==1?view.v[index]:view.w[index]);const unsigned char* component_links=component==0?links_u:(component==1?links_v:links_w);const int coordinate[3]={i,j,k};const bool component_lower=coordinate[component]==0&&(record.flags&(1u<<(2*component))),component_upper=coordinate[component]==bs&&(record.flags&(1u<<(2*component+1)));const Real component_width=(external_aero&&(component_lower||component_upper))?Real(0.5):Real(1),h=Real(record.h),volume=component_width*h*h*h;Real integrated=0;for(int axis=0;axis<3;++axis){const Real positive=conservative_link_flux(hierarchy,views,level,record,component,brick,i,j,k,current,axis,1,component_links,external_aero,freestream_speed),negative=conservative_link_flux(hierarchy,views,level,record,component,brick,i,j,k,current,axis,-1,component_links,external_aero,freestream_speed),area=(axis==component?Real(1):component_width)*h*h;integrated+=area*(positive-negative);}const Real result=current-dt*integrated/volume;if(component==0)u_out[index]=result;else if(component==1)v_out[index]=result;else w_out[index]=result;
		}
		__device__ Real cell_value_clamped(const DeviceAmrFieldLevelView& view,int brick,int i,int j,int k)
		{
			const int bs=view.layout.brick_size;i=max(0,min(bs-1,i));j=max(0,min(bs-1,j));k=max(0,min(bs-1,k));return view.nut[view.layout.cell_index(brick,i,j,k)];
		}
		__device__ Real cell_value_extrapolated(const DeviceAmrFieldLevelView& view,int brick,int i,int j,int k)
		{
			const int bs=view.layout.brick_size,requested[3]={i,j,k};int index[3][2],count[3];Real weight[3][2];for(int axis=0;axis<3;++axis){if(requested[axis]<0){index[axis][0]=0;index[axis][1]=1;weight[axis][0]=Real(2);weight[axis][1]=Real(-1);count[axis]=2;}else if(requested[axis]>=bs){index[axis][0]=bs-1;index[axis][1]=bs-2;weight[axis][0]=Real(2);weight[axis][1]=Real(-1);count[axis]=2;}else{index[axis][0]=requested[axis];weight[axis][0]=Real(1);count[axis]=1;}}Real result=0;for(int az=0;az<count[2];++az)for(int ay=0;ay<count[1];++ay)for(int ax=0;ax<count[0];++ax)result+=weight[0][ax]*weight[1][ay]*weight[2][az]*cell_value_clamped(view,brick,index[0][ax],index[1][ay],index[2][az]);return result;
		}
		__device__ Real sample_cell_local_extrapolated(const DeviceAmrFieldLevelView& view,const GpuBrickRecord& record,int brick,GpuAmrPoint point)
		{
			Real x=(Real(point.x)-Real(record.origin_x))/Real(record.h)-Real(0.5),y=(Real(point.y)-Real(record.origin_y))/Real(record.h)-Real(0.5),z=(Real(point.z)-Real(record.origin_z))/Real(record.h)-Real(0.5);const int i0=static_cast<int>(floor(x)),j0=static_cast<int>(floor(y)),k0=static_cast<int>(floor(z));x-=i0;y-=j0;z-=k0;Real result=0;for(int dz=0;dz<2;++dz)for(int dy=0;dy<2;++dy)for(int dx=0;dx<2;++dx){const Real wx=dx?x:Real(1)-x,wy=dy?y:Real(1)-y,wz=dz?z:Real(1)-z;result+=wx*wy*wz*cell_value_extrapolated(view,brick,i0+dx,j0+dy,k0+dz);}return result;
		}
		__device__ __noinline__ Real cell_value_at_node(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int source_level,const GpuBrickRecord& source_record,int source_brick,int i,int j,int k)
		{
			const DeviceAmrFieldLevelView source_view=views[source_level];const int bs=source_view.layout.brick_size;if(i>=0&&j>=0&&k>=0&&i<bs&&j<bs&&k<bs)return cell_value_clamped(source_view,source_brick,i,j,k);int face=-1,outside=0;if(i<0){face=0;++outside;}else if(i>=bs){face=1;++outside;}if(j<0){face=2;++outside;}else if(j>=bs){face=3;++outside;}if(k<0){face=4;++outside;}else if(k>=bs){face=5;++outside;}if(outside==1){const int neighbor=source_view.neighbors[source_brick*6+face];if(neighbor>=0&&!(source_view.flags[neighbor]&BRICK_COVERED)){if(face==0)i+=bs;else if(face==1)i-=bs;else if(face==2)j+=bs;else if(face==3)j-=bs;else if(face==4)k+=bs;else k-=bs;return cell_value_clamped(source_view,neighbor,i,j,k);}}const GpuAmrPoint node{source_record.origin_x+(i+0.5f)*source_record.h,source_record.origin_y+(j+0.5f)*source_record.h,source_record.origin_z+(k+0.5f)*source_record.h};if(node.x<hierarchy.domain_lo.x||node.y<hierarchy.domain_lo.y||node.z<hierarchy.domain_lo.z||node.x>=hierarchy.domain_hi.x||node.y>=hierarchy.domain_hi.y||node.z>=hierarchy.domain_hi.z)return cell_value_clamped(source_view,source_brick,i,j,k);const GpuBrickLocation location=locate_finest(hierarchy,node);if(location.brick<0)return cell_value_clamped(source_view,source_brick,i,j,k);const GpuBrickRecord record=hierarchy.levels[location.level].bricks[location.brick];const DeviceAmrFieldLevelView view=views[location.level];if(location.level!=source_level)return sample_cell_local_extrapolated(view,record,location.brick,node);const int ni=nearest_lattice_index(node.x,record.origin_x,record.h,0.5f),nj=nearest_lattice_index(node.y,record.origin_y,record.h,0.5f),nk=nearest_lattice_index(node.z,record.origin_z,record.h,0.5f);return cell_value_clamped(view,location.brick,ni,nj,nk);
		}
		__device__ Real stencil_face_value(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,const GpuBrickRecord& record,int component,int brick,int i,int j,int k)
		{
			return valid_face_index(component,views[level].layout.brick_size,i,j,k)?value_at_clamped(views[level],component,brick,i,j,k):value_at_node(hierarchy,views,level,record,component,brick,i,j,k);
		}
		__device__ Real stencil_cell_value(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,const GpuBrickRecord& record,int brick,int i,int j,int k)
		{
			const int bs=views[level].layout.brick_size;return i>=0&&j>=0&&k>=0&&i<bs&&j<bs&&k<bs?cell_value_clamped(views[level],brick,i,j,k):cell_value_at_node(hierarchy,views,level,record,brick,i,j,k);
		}
		__global__ void smagorinsky_kernel(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,Real cs)
		{
			const DeviceAmrFieldLevelView view=views[level];const int bs=view.layout.brick_size,q=blockIdx.x*blockDim.x+threadIdx.x,cells=bs*bs*bs,total=view.brick_count*cells;if(q>=total)return;const int brick=q/cells,local=q%cells;if(view.flags[brick]&BRICK_COVERED)return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs);const GpuBrickRecord record=hierarchy.levels[level].bricks[brick];const float h=record.h;const GpuAmrPoint centre{record.origin_x+(i+0.5f)*h,record.origin_y+(j+0.5f)*h,record.origin_z+(k+0.5f)*h};Real gradient[3][3];for(int component=0;component<3;++component)for(int axis=0;axis<3;++axis){GpuAmrPoint lower=centre,upper=centre;if(axis==0){lower.x-=h;upper.x+=h;}else if(axis==1){lower.y-=h;upper.y+=h;}else{lower.z-=h;upper.z+=h;}gradient[component][axis]=(sample_local(hierarchy,views,level,record,component,brick,upper)-sample_local(hierarchy,views,level,record,component,brick,lower))/(Real(2)*h);}const Real sxx=gradient[0][0],syy=gradient[1][1],szz=gradient[2][2],sxy=Real(0.5)*(gradient[0][1]+gradient[1][0]),sxz=Real(0.5)*(gradient[0][2]+gradient[2][0]),syz=Real(0.5)*(gradient[1][2]+gradient[2][1]);const Real magnitude=sqrt(Real(2)*(sxx*sxx+syy*syy+szz*szz+Real(2)*(sxy*sxy+sxz*sxz+syz*syz)));const Real length=cs*h;view.nut[view.layout.cell_index(brick,i,j,k)]=length*length*magnitude;
		}
		__global__ void diffuse_level_kernel(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,Real molecular_nu,Real dt,const unsigned char* links_u,const unsigned char* links_v,const unsigned char* links_w,Real* u_out,Real* v_out,Real* w_out)
		{
			const DeviceAmrFieldLevelView view=views[level];const int bs=view.layout.brick_size,per_component=bs*bs*(bs+1),q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*3*per_component;if(q>=total)return;int r=q,local=r%per_component;r/=per_component;const int component=r%3,brick=r/3;if(view.flags[brick]&BRICK_COVERED)return;int i=0,j=0,k=0;if(component==0){i=local%(bs+1);local/=bs+1;j=local%bs;k=local/bs;}else if(component==1){i=local%bs;local/=bs;j=local%(bs+1);k=local/(bs+1);}else{i=local%bs;local/=bs;j=local%bs;k=local/bs;}const GpuBrickRecord record=hierarchy.levels[level].bricks[brick];const std::size_t index=component==0?view.layout.u_index(brick,i,j,k):(component==1?view.layout.v_index(brick,i,j,k):view.layout.w_index(brick,i,j,k));const Real current=component==0?view.u[index]:(component==1?view.v[index]:view.w[index]);const unsigned char links=component==0?links_u[index]:(component==1?links_v[index]:links_w[index]);Real laplacian=0,viscosity=molecular_nu;
			if(!(links&0x40)){laplacian=stencil_face_value(hierarchy,views,level,record,component,brick,i-1,j,k)+stencil_face_value(hierarchy,views,level,record,component,brick,i+1,j,k)+stencil_face_value(hierarchy,views,level,record,component,brick,i,j-1,k)+stencil_face_value(hierarchy,views,level,record,component,brick,i,j+1,k)+stencil_face_value(hierarchy,views,level,record,component,brick,i,j,k-1)+stencil_face_value(hierarchy,views,level,record,component,brick,i,j,k+1)-Real(6)*current;int ci=i,cj=j,ck=k;if(component==0)--ci;else if(component==1)--cj;else --ck;const int ci2=ci+(component==0),cj2=cj+(component==1),ck2=ck+(component==2);viscosity+=Real(0.5)*(stencil_cell_value(hierarchy,views,level,record,brick,ci,cj,ck)+stencil_cell_value(hierarchy,views,level,record,brick,ci2,cj2,ck2));}
			else for(int axis=0;axis<3;++axis)for(int direction=-1;direction<=1;direction+=2){const int bit=2*axis+(direction>0);if(!(links&(1u<<bit)))continue;int c[3]={i,j,k};c[axis]+=direction;laplacian+=(valid_face_index(component,bs,c[0],c[1],c[2])?value_at_clamped(view,component,brick,c[0],c[1],c[2]):value_at_node(hierarchy,views,level,record,component,brick,c[0],c[1],c[2]))-current;}const Real h=hierarchy.levels[level].h,result=current+dt*viscosity*laplacian/(h*h);
			if(component==0)u_out[index]=result;else if(component==1)v_out[index]=result;else w_out[index]=result;
		}

		Vec3d face_centre(const BrickMetadata& brick,double h,int component,int i,int j,int k)
		{
			return brick.origin+Vec3d{(i+(component==0?0.0:0.5))*h,(j+(component==1?0.0:0.5))*h,(k+(component==2?0.0:0.5))*h};
		}
	}

	void conservative_pairwise_momentum_cpu(const std::vector<double>& dual_volume,
		const std::vector<PairwiseMomentumConnection>& connections,double dt,
		std::vector<double>& velocity_x,std::vector<double>& velocity_y,
		std::vector<double>& velocity_z)
	{
		if(!(dt>0))throw std::invalid_argument("pairwise momentum timestep must be positive");const std::size_t n=dual_volume.size();
		if(velocity_x.size()!=n||velocity_y.size()!=n||velocity_z.size()!=n)throw std::invalid_argument("pairwise momentum state size mismatch");
		for(double volume:dual_volume)if(!(volume>0)||!std::isfinite(volume))throw std::invalid_argument("pairwise momentum dual volume must be finite and positive");
		std::vector<double> dx(n,0),dy(n,0),dz(n,0);
		for(const PairwiseMomentumConnection& edge:connections)
		{
			if(edge.a<0||edge.b<0||edge.a==edge.b||edge.a>=static_cast<int>(n)||edge.b>=static_cast<int>(n)||!(edge.open_area>=0)||!std::isfinite(edge.open_area)||!std::isfinite(edge.normal_velocity))throw std::invalid_argument("invalid pairwise momentum connection");
			const double transport=dt*edge.open_area*edge.normal_velocity;const int donor=transport>=0?edge.a:edge.b;
			const double mx=transport*velocity_x[donor],my=transport*velocity_y[donor],mz=transport*velocity_z[donor];
			dx[edge.a]-=mx;dx[edge.b]+=mx;dy[edge.a]-=my;dy[edge.b]+=my;dz[edge.a]-=mz;dz[edge.b]+=mz;
		}
		for(std::size_t node=0;node<n;++node){velocity_x[node]+=dx[node]/dual_volume[node];velocity_y[node]+=dy[node]/dual_volume[node];velocity_z[node]+=dz[node]/dual_volume[node];}
	}

	double regular_mac_dual_volume(const AmrHierarchy& hierarchy,const AmrMacFaceAddress& address)
	{
		if(address.level<0||address.level>=static_cast<int>(hierarchy.levels().size())||address.component<0||address.component>2)throw std::invalid_argument("regular MAC dual-volume address outside hierarchy");const AmrLevel& level=hierarchy.levels()[address.level];if(address.brick<0||address.brick>=static_cast<int>(level.bricks.size())||!level.bricks[address.brick].active())throw std::invalid_argument("regular MAC dual volume requires an active brick");const int bs=hierarchy.brick_size(),coordinate[3]={address.i,address.j,address.k},extent[3]={address.component==0?bs+1:bs,address.component==1?bs+1:bs,address.component==2?bs+1:bs};for(int axis=0;axis<3;++axis)if(coordinate[axis]<0||coordinate[axis]>=extent[axis])throw std::invalid_argument("regular MAC dual-volume local face index outside brick");const double half_volume=0.5*level.h*level.h*level.h;double volume=0;
		for(int side=-1;side<=1;side+=2){int cell[3]={address.i,address.j,address.k};cell[address.component]+=side<0?-1:0;int brick=address.brick;if(cell[address.component]<0||cell[address.component]>=bs){const int face=2*address.component+(side>0),neighbour=level.bricks[brick].same_level_neighbor[face];if(neighbour<0||!level.bricks[neighbour].active())continue;if(cell[address.component]<0)cell[address.component]+=bs;else cell[address.component]-=bs;brick=neighbour;}if(cell[0]>=0&&cell[1]>=0&&cell[2]>=0&&cell[0]<bs&&cell[1]<bs&&cell[2]<bs&&level.bricks[brick].active())volume+=half_volume;}
		return volume;
	}

	AmrMacFaceAddress canonical_amr_mac_face_address(const AmrHierarchy& hierarchy,AmrMacFaceAddress address)
	{
		if(address.level<0||address.level>=static_cast<int>(hierarchy.levels().size())||address.component<0||address.component>2)throw std::invalid_argument("canonical MAC address outside hierarchy");const AmrLevel& level=hierarchy.levels()[address.level];if(address.brick<0||address.brick>=static_cast<int>(level.bricks.size())||!level.bricks[address.brick].active())throw std::invalid_argument("canonical MAC address requires an active brick");const int bs=hierarchy.brick_size(),extent[3]={address.component==0?bs+1:bs,address.component==1?bs+1:bs,address.component==2?bs+1:bs},coordinate[3]={address.i,address.j,address.k};for(int axis=0;axis<3;++axis)if(coordinate[axis]<0||coordinate[axis]>=extent[axis])throw std::invalid_argument("canonical MAC address outside brick lattice");int* normal_coordinate=address.component==0?&address.i:(address.component==1?&address.j:&address.k);if(*normal_coordinate==0){const int neighbour=level.bricks[address.brick].same_level_neighbor[2*address.component];if(neighbour>=0&&level.bricks[neighbour].active()){address.brick=neighbour;*normal_coordinate=bs;}}return address;
	}

	void conservative_pairwise_scalar_cpu(const std::vector<double>& dual_volume,
		const std::vector<PairwiseMomentumConnection>& connections,double dt,
		std::vector<double>& velocity)
	{
		if(!(dt>0))throw std::invalid_argument("pairwise scalar timestep must be positive");const std::size_t n=dual_volume.size();if(velocity.size()!=n)throw std::invalid_argument("pairwise scalar state size mismatch");for(double volume:dual_volume)if(!(volume>0)||!std::isfinite(volume))throw std::invalid_argument("pairwise scalar dual volume must be finite and positive");std::vector<double> delta(n,0);
		for(const PairwiseMomentumConnection& edge:connections){if(edge.a<0||edge.b<0||edge.a==edge.b||edge.a>=static_cast<int>(n)||edge.b>=static_cast<int>(n)||!(edge.open_area>=0)||!std::isfinite(edge.open_area)||!std::isfinite(edge.normal_velocity))throw std::invalid_argument("invalid pairwise scalar connection");const double transport=dt*edge.open_area*edge.normal_velocity,momentum=transport*velocity[transport>=0?edge.a:edge.b];delta[edge.a]-=momentum;delta[edge.b]+=momentum;}for(std::size_t node=0;node<n;++node)velocity[node]+=delta[node]/dual_volume[node];
	}

	std::vector<CompositeEbMomentumRegularConnection>
		build_composite_eb_momentum_regular_connections(const CompositeAmrPressureSystem& system)
	{
		if(!system.hierarchy||system.brick_size<=0)throw std::invalid_argument("EB momentum perimeter requires a composite hierarchy");const AmrHierarchy& hierarchy=*system.hierarchy;const int bs=system.brick_size,cells=bs*bs*bs;
		std::unordered_set<int> embedded_nodes;for(const CoarseFinePressureConnection& connection:system.embedded){if(connection.coarse_dof>=0)embedded_nodes.insert(connection.coarse_dof);if(connection.fine_dof>=0)embedded_nodes.insert(connection.fine_dof);}
		struct CellAddress{int level=-1,brick=-1,i=-1,j=-1,k=-1;};auto cell_address=[&](int dof){CellAddress out;for(int level=0;level<static_cast<int>(hierarchy.levels().size());++level){const int begin=system.level_offset[level],count=static_cast<int>(hierarchy.levels()[level].bricks.size())*cells;if(dof<begin||dof>=begin+count)continue;const int work=dof-begin,local=work%cells;out.level=level;out.brick=work/cells;out.i=local%bs;out.j=(local/bs)%bs;out.k=local/(bs*bs);break;}return out;};
		struct EdgeKey{int lower,upper,axis;bool operator==(const EdgeKey& other)const{return lower==other.lower&&upper==other.upper&&axis==other.axis;}};struct EdgeHash{std::size_t operator()(const EdgeKey& key)const{return (static_cast<std::size_t>(key.lower)*1315423911u^static_cast<std::size_t>(key.upper))*31u+static_cast<unsigned>(key.axis);}};std::unordered_set<EdgeKey,EdgeHash> emitted;std::vector<CompositeEbMomentumRegularConnection> output;
		for(int dof:embedded_nodes)
		{
			const CellAddress cell=cell_address(dof);if(cell.level<0||!system.active[dof])continue;const AmrLevel& level=hierarchy.levels()[cell.level];const BrickMetadata& metadata=level.bricks[cell.brick];if(!metadata.active())continue;
			for(int axis=0;axis<3;++axis)for(int sign=-1;sign<=1;sign+=2)
			{
				int coordinate[3]={cell.i,cell.j,cell.k},other_brick=cell.brick;coordinate[axis]+=sign;if(coordinate[axis]<0||coordinate[axis]>=bs){other_brick=metadata.same_level_neighbor[2*axis+(sign>0)];if(other_brick<0||!level.bricks[other_brick].active())continue;coordinate[axis]=sign>0?0:bs-1;}const int other=system.dof(cell.level,other_brick,coordinate[0],coordinate[1],coordinate[2]);if(other<0||other>=system.storage_size||!system.active[other])continue;const int lower=sign>0?dof:other,upper=sign>0?other:dof;if(system.cut_face_mask[lower]&(1u<<axis))continue;const EdgeKey key{lower,upper,axis};if(!emitted.insert(key).second)continue;
				AmrMacFaceAddress face{cell.level,cell.brick,axis,cell.i,cell.j,cell.k};if(axis==0)face.i+=sign>0;else if(axis==1)face.j+=sign>0;else face.k+=sign>0;face=canonical_amr_mac_face_address(hierarchy,face);const double area=level.h*level.h,volume=composite_mac_carrier_volume(system,lower,upper);if(!(area>0&&volume>0))throw std::runtime_error("EB momentum regular connection has non-positive geometry");output.push_back({lower,upper,face,area,volume,embedded_nodes.count(lower)!=0,embedded_nodes.count(upper)!=0});
			}
		}
		return output;
	}

	struct CompositeCellGradientTopology
	{
		std::vector<unsigned char> special_mask,compact_plus;
		std::vector<int> special_dof,incidence_special,incidence_neighbor;
		std::vector<double> incidence_weighted_displacement,inverse;
	};

	int cell_symmetric_pseudoinverse_3x3(const double* input,double* inverse)
	{
		double a[9],vectors[9]={1,0,0,0,1,0,0,0,1};for(int q=0;q<9;++q)a[q]=input[q];for(int sweep=0;sweep<16;++sweep){bool changed=false;for(int pair=0;pair<3;++pair){const int p=pair==0?0:(pair==1?0:1),q=pair==0?1:(pair==1?2:2);const double apq=a[p*3+q],scale=std::max({std::abs(a[p*3+p]),std::abs(a[q*3+q]),1e-300});if(std::abs(apq)<=1e-14*scale)continue;changed=true;const double phi=0.5*std::atan2(2*apq,a[q*3+q]-a[p*3+p]),c=std::cos(phi),s=std::sin(phi),app=a[p*3+p],aqq=a[q*3+q];for(int k=0;k<3;++k)if(k!=p&&k!=q){const double akp=a[k*3+p],akq=a[k*3+q];a[k*3+p]=a[p*3+k]=c*akp-s*akq;a[k*3+q]=a[q*3+k]=s*akp+c*akq;}a[p*3+p]=c*c*app-2*s*c*apq+s*s*aqq;a[q*3+q]=s*s*app+2*s*c*apq+c*c*aqq;a[p*3+q]=a[q*3+p]=0;for(int k=0;k<3;++k){const double vkp=vectors[k*3+p],vkq=vectors[k*3+q];vectors[k*3+p]=c*vkp-s*vkq;vectors[k*3+q]=s*vkp+c*vkq;}}if(!changed)break;}for(int q=0;q<9;++q)inverse[q]=0;const double maximum=std::max({std::abs(a[0]),std::abs(a[4]),std::abs(a[8])}),threshold=maximum*(sizeof(Real)==4?1e-4:1e-10);int rank=0;for(int mode=0;mode<3;++mode){const double eigenvalue=a[mode*3+mode];if(!(eigenvalue>threshold))continue;++rank;for(int row=0;row<3;++row)for(int column=0;column<3;++column)inverse[row*3+column]+=vectors[row*3+mode]*vectors[column*3+mode]/eigenvalue;}return rank;
	}

	CompositeCellGradientTopology build_composite_cell_gradient_topology(
		const CompositeAmrPressureSystem& system,
		const std::vector<CompositeEbMomentumRegularConnection>& compact_regular)
	{
		if(!system.hierarchy)throw std::invalid_argument("cell gradient topology requires a hierarchy");CompositeCellGradientTopology out;out.special_mask.assign(system.storage_size,0);out.compact_plus.assign(system.storage_size,0);for(int q=0;q<system.storage_size;++q)if(system.active[q]&&system.cut_face_mask[q])out.special_mask[q]=1;auto mark=[&](int q){if(q>=0&&q<system.storage_size&&system.active[q])out.special_mask[q]=1;};for(const auto& connection:compact_regular){mark(connection.lower_dof);mark(connection.upper_dof);out.compact_plus[connection.lower_dof]|=static_cast<unsigned char>(1u<<connection.face.component);}for(const auto& connection:system.coarse_fine){mark(connection.coarse_dof);mark(connection.fine_dof);}for(const auto& connection:system.embedded){mark(connection.coarse_dof);mark(connection.fine_dof);}for(const auto& patch:system.surface_patches){mark(patch.plus_dof);mark(patch.minus_dof);}std::vector<int> special_index(system.storage_size,-1);for(int q=0;q<system.storage_size;++q)if(out.special_mask[q]){special_index[q]=static_cast<int>(out.special_dof.size());out.special_dof.push_back(q);}std::vector<double> normal_matrix(out.special_dof.size()*9,0);auto add_incidence=[&](int node,int neighbor,double area){const int special=special_index[node];if(special<0)return;const Vec3d displacement=system.centroid[neighbor]-system.centroid[node];const double distance2=length2(displacement);if(!(area>0)||!(distance2>1e-20))throw std::runtime_error("invalid cell gradient incidence geometry");const double weight=area/distance2;out.incidence_special.push_back(special);out.incidence_neighbor.push_back(neighbor);out.incidence_weighted_displacement.insert(out.incidence_weighted_displacement.end(),{weight*displacement.x,weight*displacement.y,weight*displacement.z});for(int row=0;row<3;++row)for(int column=0;column<3;++column)normal_matrix[special*9+row*3+column]+=weight*displacement[row]*displacement[column];};auto add_connection=[&](int a,int b,double area){if(a<0||b<0||a==b||!system.active[a]||!system.active[b])throw std::runtime_error("invalid cell gradient connection");add_incidence(a,b,area);add_incidence(b,a,area);};const AmrHierarchy& hierarchy=*system.hierarchy;const int bs=system.brick_size,cells=bs*bs*bs;for(int level_index=0;level_index<static_cast<int>(hierarchy.levels().size());++level_index){const AmrLevel& level=hierarchy.levels()[level_index];const int offset=system.level_offset[level_index];const double area=level.h*level.h;for(int brick=0;brick<static_cast<int>(level.bricks.size());++brick){const BrickMetadata& record=level.bricks[brick];if(!record.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=offset+brick*cells+i+bs*(j+bs*k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if((system.cut_face_mask[a]&(1u<<axis))||(out.compact_plus[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=record.same_level_neighbor[2*axis+1];if(other_brick<0||!level.bricks[other_brick].active())continue;c[axis]=0;}const int b=offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(system.active[b]&&(out.special_mask[a]||out.special_mask[b]))add_connection(a,b,area);}}}}for(const auto& connection:compact_regular)add_connection(connection.lower_dof,connection.upper_dof,connection.open_area);auto add_special=[&](const CoarseFinePressureConnection& connection){add_connection(connection.coarse_dof,connection.fine_dof,connection.open_area);};for(const auto& connection:system.coarse_fine)add_special(connection);for(const auto& connection:system.embedded)add_special(connection);out.inverse.resize(out.special_dof.size()*9);for(int special=0;special<static_cast<int>(out.special_dof.size());++special)cell_symmetric_pseudoinverse_3x3(normal_matrix.data()+special*9,out.inverse.data()+special*9);return out;
	}

	int composite_cell_neighbor_cpu(const CompositeAmrPressureSystem& system,
		const std::vector<unsigned char>& compact_plus,int level_index,int brick,int i,int j,
		int k,int axis,int sign)
	{
		const AmrLevel& level=system.hierarchy->levels()[level_index];const int bs=system.brick_size,cells=bs*bs*bs;int c[3]={i,j,k},other_brick=brick;c[axis]+=sign;if(c[axis]<0||c[axis]>=bs){other_brick=level.bricks[brick].same_level_neighbor[2*axis+(sign>0)];if(other_brick<0||!level.bricks[other_brick].active())return -1;c[axis]=sign>0?0:bs-1;}const int node=system.level_offset[level_index]+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(!system.active[node])return -1;const int current=system.level_offset[level_index]+brick*cells+i+bs*(j+bs*k),lower=sign>0?current:node;if((system.cut_face_mask[lower]&(1u<<axis))||(compact_plus[lower]&(1u<<axis)))return -1;return node;
	}

	void conservative_composite_cell_momentum_cpu(const CompositeAmrPressureSystem& system,
		const AmrHostFields& fields,const CompositeAmrFluxes& fluxes,double dt,
		CompositeCellMomentumState& state,bool external_aero,double freestream_speed)
	{
		if(!system.hierarchy||fields.levels().size()!=system.hierarchy->levels().size()||fluxes.coarse_fine_velocity.size()!=system.coarse_fine.size()||fluxes.embedded_velocity.size()!=system.embedded.size()||!(dt>0)||!(freestream_speed>=0)||state.x.size()!=static_cast<std::size_t>(system.storage_size)||state.y.size()!=state.x.size()||state.z.size()!=state.x.size())throw std::invalid_argument("invalid composite cell momentum CPU state");const AmrHierarchy& hierarchy=*system.hierarchy;const int bs=system.brick_size,cells=bs*bs*bs;std::vector<double> dx(system.storage_size),dy(system.storage_size),dz(system.storage_size);const std::vector<CompositeEbMomentumRegularConnection> compact_regular=build_composite_eb_momentum_regular_connections(system);std::vector<unsigned char> compact_plus(system.storage_size,0);for(const auto& connection:compact_regular)if(connection.lower_dof>=0&&connection.lower_dof<system.storage_size)compact_plus[connection.lower_dof]|=static_cast<unsigned char>(1u<<connection.face.component);
		auto transfer=[&](int a,int b,double swept){if(a<0||b<0||a==b||!system.active[a]||!system.active[b])throw std::runtime_error("invalid composite cell momentum connection");const int donor=swept>=0?a:b;const double mx=swept*state.x[donor],my=swept*state.y[donor],mz=swept*state.z[donor];dx[a]-=mx;dx[b]+=mx;dy[a]-=my;dy[b]+=my;dz[a]-=mz;dz[b]+=mz;};
		auto mapped_face=[&](const AmrMacFaceAddress& face){const auto& level=fields.levels()[face.level];if(face.component==0)return static_cast<double>(level.u[level.layout.u_index(face.brick,face.i,face.j,face.k)]);if(face.component==1)return static_cast<double>(level.v[level.layout.v_index(face.brick,face.i,face.j,face.k)]);return static_cast<double>(level.w[level.layout.w_index(face.brick,face.i,face.j,face.k)]);};
		for(int level_index=0;level_index<static_cast<int>(hierarchy.levels().size());++level_index){const AmrLevel& metadata=hierarchy.levels()[level_index];const AmrHostLevelFields& values=fields.levels()[level_index];const int offset=system.level_offset[level_index];const double area=metadata.h*metadata.h;auto face_value=[&](int brick,int axis,int i,int j,int k){if(axis==0)return static_cast<double>(values.u[values.layout.u_index(brick,i+1,j,k)]);if(axis==1)return static_cast<double>(values.v[values.layout.v_index(brick,i,j+1,k)]);return static_cast<double>(values.w[values.layout.w_index(brick,i,j,k+1)]);};for(int brick=0;brick<static_cast<int>(metadata.bricks.size());++brick){const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=offset+brick*cells+i+bs*(j+bs*k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if((system.cut_face_mask[a]&(1u<<axis))||(compact_plus[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=record.same_level_neighbor[2*axis+1];if(other_brick<0||!metadata.bricks[other_brick].active())continue;c[axis]=0;}const int b=offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(system.active[b])transfer(a,b,dt*area*face_value(brick,axis,i,j,k));}if(external_aero&&(record.flags&BRICK_XMIN)&&i==0){const double swept=dt*area*static_cast<double>(values.u[values.layout.u_index(brick,0,j,k)]);if(swept>=0){dx[a]+=swept*freestream_speed;}else{dx[a]+=swept*state.x[a];dy[a]+=swept*state.y[a];dz[a]+=swept*state.z[a];}}if(external_aero&&(record.flags&BRICK_XMAX)&&i==bs-1){const double swept=dt*area*static_cast<double>(values.u[values.layout.u_index(brick,bs,j,k)]);if(swept>=0){dx[a]-=swept*state.x[a];dy[a]-=swept*state.y[a];dz[a]-=swept*state.z[a];}else dx[a]-=swept*freestream_speed;}}}}
		for(std::size_t edge=0;edge<system.coarse_fine.size();++edge){const auto& connection=system.coarse_fine[edge];const int a=connection.direction>0?connection.coarse_dof:connection.fine_dof,b=connection.direction>0?connection.fine_dof:connection.coarse_dof;transfer(a,b,dt*connection.open_area*fluxes.coarse_fine_velocity[edge]);}for(std::size_t edge=0;edge<system.embedded.size();++edge){const auto& connection=system.embedded[edge];const int a=connection.direction>0?connection.coarse_dof:connection.fine_dof,b=connection.direction>0?connection.fine_dof:connection.coarse_dof;transfer(a,b,dt*connection.open_area*fluxes.embedded_velocity[edge]);}for(const auto& connection:compact_regular)transfer(connection.lower_dof,connection.upper_dof,dt*connection.open_area*mapped_face(connection.face));for(int q=0;q<system.storage_size;++q)if(system.active[q]){const double inverse=1.0/system.volume[q];state.x[q]=static_cast<Real>(static_cast<double>(state.x[q])+dx[q]*inverse);state.y[q]=static_cast<Real>(static_cast<double>(state.y[q])+dy[q]*inverse);state.z[q]=static_cast<Real>(static_cast<double>(state.z[q])+dz[q]*inverse);}
	}

	void conservative_composite_cell_momentum_cpu(const CompositeAmrPressureSystem& system,
		const AmrHostFields& fields,const std::vector<double>& embedded_velocity,double dt,
		CompositeCellMomentumState& state,bool external_aero,double freestream_speed)
	{
		if(!system.coarse_fine.empty())throw std::invalid_argument("coarse/fine fluxes are required for multilevel composite cell momentum transport");CompositeAmrFluxes fluxes=make_zero_composite_fluxes(system);fluxes.embedded_velocity=embedded_velocity;conservative_composite_cell_momentum_cpu(system,fields,fluxes,dt,state,external_aero,freestream_speed);
	}

	void diffuse_composite_cell_momentum_cpu(const CompositeAmrPressureSystem& system,
		double kinematic_viscosity,double dt,CompositeCellMomentumState& state)
	{
		if(!system.hierarchy||kinematic_viscosity<0||!(dt>0)||state.x.size()!=static_cast<std::size_t>(system.storage_size)||state.y.size()!=state.x.size()||state.z.size()!=state.x.size())throw std::invalid_argument("invalid composite cell momentum diffusion state");const AmrHierarchy& hierarchy=*system.hierarchy;const int bs=system.brick_size,cells=bs*bs*bs;std::vector<double> dx(system.storage_size),dy(system.storage_size),dz(system.storage_size);const std::vector<CompositeEbMomentumRegularConnection> compact_regular=build_composite_eb_momentum_regular_connections(system);std::vector<unsigned char> compact_plus(system.storage_size,0);for(const auto& connection:compact_regular)compact_plus[connection.lower_dof]|=static_cast<unsigned char>(1u<<connection.face.component);auto transfer=[&](int a,int b,double conductance){if(a<0||b<0||a==b||!system.active[a]||!system.active[b]||!(conductance>0)||!std::isfinite(conductance))throw std::runtime_error("invalid composite cell momentum diffusion connection");const double scale=dt*kinematic_viscosity*conductance,mx=scale*(state.x[b]-state.x[a]),my=scale*(state.y[b]-state.y[a]),mz=scale*(state.z[b]-state.z[a]);dx[a]+=mx;dx[b]-=mx;dy[a]+=my;dy[b]-=my;dz[a]+=mz;dz[b]-=mz;};for(int level_index=0;level_index<static_cast<int>(hierarchy.levels().size());++level_index){const AmrLevel& level=hierarchy.levels()[level_index];const int offset=system.level_offset[level_index];for(int brick=0;brick<static_cast<int>(level.bricks.size());++brick){const BrickMetadata& record=level.bricks[brick];if(!record.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=offset+brick*cells+i+bs*(j+bs*k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if((system.cut_face_mask[a]&(1u<<axis))||(compact_plus[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=record.same_level_neighbor[2*axis+1];if(other_brick<0||!level.bricks[other_brick].active())continue;c[axis]=0;}const int b=offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(system.active[b])transfer(a,b,level.h);}}}}auto special=[&](const CoarseFinePressureConnection& connection){return connection.open_area*pressure_gradient_factor(connection);};for(const auto& connection:system.coarse_fine){const int a=connection.direction>0?connection.coarse_dof:connection.fine_dof,b=connection.direction>0?connection.fine_dof:connection.coarse_dof;transfer(a,b,special(connection));}for(const auto& connection:system.embedded){const int a=connection.direction>0?connection.coarse_dof:connection.fine_dof,b=connection.direction>0?connection.fine_dof:connection.coarse_dof;transfer(a,b,special(connection));}for(const auto& connection:compact_regular){const Vec3d delta=system.centroid[connection.upper_dof]-system.centroid[connection.lower_dof];const double distance2=length2(delta),normal_distance=std::abs(delta[connection.face.component]);transfer(connection.lower_dof,connection.upper_dof,connection.open_area*normal_distance/distance2);}for(int q=0;q<system.storage_size;++q)if(system.active[q]){const double inverse=1.0/system.volume[q];state.x[q]=static_cast<Real>(static_cast<double>(state.x[q])+dx[q]*inverse);state.y[q]=static_cast<Real>(static_cast<double>(state.y[q])+dy[q]*inverse);state.z[q]=static_cast<Real>(static_cast<double>(state.z[q])+dz[q]*inverse);}
	}

	std::vector<double> composite_cell_smagorinsky_viscosity_cpu(
		const CompositeAmrPressureSystem& system,const CompositeCellMomentumState& state,
		double molecular_viscosity,double smagorinsky_cs)
	{
		if(!system.hierarchy||molecular_viscosity<0||smagorinsky_cs<0||state.x.size()!=static_cast<std::size_t>(system.storage_size)||state.y.size()!=state.x.size()||state.z.size()!=state.x.size())throw std::invalid_argument("invalid composite cell Smagorinsky state");const std::vector<CompositeEbMomentumRegularConnection> compact_regular=build_composite_eb_momentum_regular_connections(system);const CompositeCellGradientTopology topology=build_composite_cell_gradient_topology(system,compact_regular);std::vector<double> viscosity(system.storage_size,molecular_viscosity);auto component=[&](int q,int c){return static_cast<double>(c==0?state.x[q]:(c==1?state.y[q]:state.z[q]));};auto strain=[&](const double* gradient){const double sxx=gradient[0],syy=gradient[4],szz=gradient[8],sxy=0.5*(gradient[1]+gradient[3]),sxz=0.5*(gradient[2]+gradient[6]),syz=0.5*(gradient[5]+gradient[7]);return std::sqrt(std::max(0.0,2*(sxx*sxx+syy*syy+szz*szz+2*(sxy*sxy+sxz*sxz+syz*syz))));};const int bs=system.brick_size,cells=bs*bs*bs;for(int level_index=0;level_index<static_cast<int>(system.hierarchy->levels().size());++level_index){const AmrLevel& level=system.hierarchy->levels()[level_index];for(int brick=0;brick<static_cast<int>(level.bricks.size());++brick){if(!level.bricks[brick].active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=system.level_offset[level_index]+brick*cells+i+bs*(j+bs*k);if(!system.active[a]||topology.special_mask[a]||smagorinsky_cs==0)continue;double gradient[9];for(int c=0;c<3;++c)for(int axis=0;axis<3;++axis){const int lower=composite_cell_neighbor_cpu(system,topology.compact_plus,level_index,brick,i,j,k,axis,-1),upper=composite_cell_neighbor_cpu(system,topology.compact_plus,level_index,brick,i,j,k,axis,1);if(lower>=0&&upper>=0)gradient[c*3+axis]=(component(upper,c)-component(lower,c))/(2*level.h);else if(upper>=0)gradient[c*3+axis]=(component(upper,c)-component(a,c))/level.h;else if(lower>=0)gradient[c*3+axis]=(component(a,c)-component(lower,c))/level.h;else gradient[c*3+axis]=0;}const double length=smagorinsky_cs*std::cbrt(system.volume[a]);viscosity[a]+=length*length*strain(gradient);}}}std::vector<double> rhs(topology.special_dof.size()*9,0);for(int edge=0;edge<static_cast<int>(topology.incidence_special.size());++edge){const int special=topology.incidence_special[edge],node=topology.special_dof[special],neighbor=topology.incidence_neighbor[edge];for(int c=0;c<3;++c){const double delta=component(neighbor,c)-component(node,c);for(int axis=0;axis<3;++axis)rhs[special*9+c*3+axis]+=topology.incidence_weighted_displacement[3*edge+axis]*delta;}}for(int special=0;special<static_cast<int>(topology.special_dof.size());++special){double gradient[9];for(int c=0;c<3;++c)for(int row=0;row<3;++row){gradient[c*3+row]=0;for(int column=0;column<3;++column)gradient[c*3+row]+=topology.inverse[special*9+row*3+column]*rhs[special*9+c*3+column];}const int node=topology.special_dof[special];const double length=smagorinsky_cs*std::cbrt(system.volume[node]);viscosity[node]+=length*length*strain(gradient);}return viscosity;
	}

	void diffuse_composite_cell_momentum_smagorinsky_cpu(
		const CompositeAmrPressureSystem& system,double molecular_viscosity,
		double smagorinsky_cs,double dt,CompositeCellMomentumState& state)
	{
		if(!(dt>0))throw std::invalid_argument("invalid composite cell Smagorinsky timestep");const std::vector<double> viscosity=composite_cell_smagorinsky_viscosity_cpu(system,state,molecular_viscosity,smagorinsky_cs);const AmrHierarchy& hierarchy=*system.hierarchy;const int bs=system.brick_size,cells=bs*bs*bs;std::vector<double> dx(system.storage_size),dy(system.storage_size),dz(system.storage_size);const std::vector<CompositeEbMomentumRegularConnection> compact_regular=build_composite_eb_momentum_regular_connections(system);std::vector<unsigned char> compact_plus(system.storage_size,0);for(const auto& connection:compact_regular)compact_plus[connection.lower_dof]|=static_cast<unsigned char>(1u<<connection.face.component);auto transfer=[&](int a,int b,double conductance){if(a<0||b<0||a==b||!system.active[a]||!system.active[b]||!(conductance>0))throw std::runtime_error("invalid composite cell Smagorinsky connection");const double scale=dt*conductance*0.5*(viscosity[a]+viscosity[b]),mx=scale*(state.x[b]-state.x[a]),my=scale*(state.y[b]-state.y[a]),mz=scale*(state.z[b]-state.z[a]);dx[a]+=mx;dx[b]-=mx;dy[a]+=my;dy[b]-=my;dz[a]+=mz;dz[b]-=mz;};for(int level_index=0;level_index<static_cast<int>(hierarchy.levels().size());++level_index){const AmrLevel& level=hierarchy.levels()[level_index];const int offset=system.level_offset[level_index];for(int brick=0;brick<static_cast<int>(level.bricks.size());++brick){const BrickMetadata& record=level.bricks[brick];if(!record.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=offset+brick*cells+i+bs*(j+bs*k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if((system.cut_face_mask[a]&(1u<<axis))||(compact_plus[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=record.same_level_neighbor[2*axis+1];if(other_brick<0||!level.bricks[other_brick].active())continue;c[axis]=0;}const int b=offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(system.active[b])transfer(a,b,level.h);}}}}auto special=[&](const CoarseFinePressureConnection& connection){const int a=connection.direction>0?connection.coarse_dof:connection.fine_dof,b=connection.direction>0?connection.fine_dof:connection.coarse_dof;transfer(a,b,connection.open_area*pressure_gradient_factor(connection));};for(const auto& connection:system.coarse_fine)special(connection);for(const auto& connection:system.embedded)special(connection);for(const auto& connection:compact_regular){const Vec3d delta=system.centroid[connection.upper_dof]-system.centroid[connection.lower_dof];transfer(connection.lower_dof,connection.upper_dof,connection.open_area*std::abs(delta[connection.face.component])/length2(delta));}for(int q=0;q<system.storage_size;++q)if(system.active[q]){const double inverse=1.0/system.volume[q];state.x[q]=static_cast<Real>(static_cast<double>(state.x[q])+dx[q]*inverse);state.y[q]=static_cast<Real>(static_cast<double>(state.y[q])+dy[q]*inverse);state.z[q]=static_cast<Real>(static_cast<double>(state.z[q])+dz[q]*inverse);}
	}

	void apply_composite_cell_pressure_impulse_cpu(const CompositeAmrPressureSystem& system,
		const std::vector<double>& pressure,double dt,double density,
		CompositeCellMomentumState& state,bool include_physical_boundaries)
	{
		if(!system.hierarchy||pressure.size()!=static_cast<std::size_t>(system.storage_size)||
			!(dt>0)||!(density>0)||!std::isfinite(dt)||!std::isfinite(density)||
			state.x.size()!=pressure.size()||state.y.size()!=pressure.size()||state.z.size()!=pressure.size())
			throw std::invalid_argument("invalid composite cell pressure impulse state");
		const AmrHierarchy& hierarchy=*system.hierarchy;const int bs=system.brick_size,cells=bs*bs*bs;
		const double scale=dt/density;std::vector<double> dx(system.storage_size),dy(system.storage_size),dz(system.storage_size);
		for(int q=0;q<system.storage_size;++q)if(system.active[q]&&!std::isfinite(pressure[q]))throw std::invalid_argument("composite cell pressure impulse contains non-finite pressure");
		const std::vector<CompositeEbMomentumRegularConnection> compact_regular=build_composite_eb_momentum_regular_connections(system);
		std::vector<unsigned char> compact_plus(system.storage_size,0);for(const auto& connection:compact_regular)compact_plus[connection.lower_dof]|=static_cast<unsigned char>(1u<<connection.face.component);
		auto axis_delta=[&](int dof,int axis,double value){if(axis==0)dx[dof]+=value;else if(axis==1)dy[dof]+=value;else dz[dof]+=value;};
		auto internal=[&](int a,int b,int axis,double area,double upper_weight){if(a<0||b<0||a==b||!system.active[a]||!system.active[b]||!(area>0)||upper_weight<0||upper_weight>1)throw std::runtime_error("invalid composite pressure impulse connection");const double face_pressure=(1-upper_weight)*pressure[a]+upper_weight*pressure[b],impulse=scale*area*face_pressure;axis_delta(a,axis,-impulse);axis_delta(b,axis,impulse);};
		auto interpolation_weight=[&](int a,int b,int axis,double face_coordinate){const double ca=system.centroid[a][axis],cb=system.centroid[b][axis],distance=cb-ca;if(!(std::abs(distance)>1e-14))throw std::runtime_error("composite pressure impulse has coincident centroids");const double weight=(face_coordinate-ca)/distance;if(weight<-1e-8||weight>1+1e-8)throw std::runtime_error("composite pressure impulse face lies outside endpoint centroids");return std::clamp(weight,0.0,1.0);};
		for(int level_index=0;level_index<static_cast<int>(hierarchy.levels().size());++level_index)
		{
			const AmrLevel& level=hierarchy.levels()[level_index];const int offset=system.level_offset[level_index];const double area=level.h*level.h;
			for(int brick=0;brick<static_cast<int>(level.bricks.size());++brick)
			{
				const BrickMetadata& record=level.bricks[brick];if(!record.active())continue;
				for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i)
				{
					const int a=offset+brick*cells+i+bs*(j+bs*k);if(!system.active[a])continue;
					for(int axis=0;axis<3;++axis)
					{
						if((system.cut_face_mask[a]&(1u<<axis))||(compact_plus[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=record.same_level_neighbor[2*axis+1];if(other_brick<0||!level.bricks[other_brick].active())continue;c[axis]=0;}const int b=offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(system.active[b])internal(a,b,axis,area,0.5);
					}
					if(include_physical_boundaries)
					{
						const double boundary=scale*area*pressure[a];if(i==0&&(record.flags&BRICK_XMIN))dx[a]+=boundary;if(i==bs-1&&(record.flags&BRICK_XMAX)&&!system.pressure_outlet_xmax)dx[a]-=boundary;if(j==0&&(record.flags&BRICK_YMIN))dy[a]+=boundary;if(j==bs-1&&(record.flags&BRICK_YMAX))dy[a]-=boundary;if(k==0&&(record.flags&BRICK_ZMIN))dz[a]+=boundary;if(k==bs-1&&(record.flags&BRICK_ZMAX))dz[a]-=boundary;
					}
				}
			}
		}
		for(const auto& connection:compact_regular){const int axis=connection.face.component;const AmrLevel& level=hierarchy.levels()[connection.face.level];const BrickMetadata& record=level.bricks[connection.face.brick];const int coordinate=axis==0?connection.face.i:(axis==1?connection.face.j:connection.face.k);const double face_coordinate=record.origin[axis]+coordinate*level.h;internal(connection.lower_dof,connection.upper_dof,axis,connection.open_area,interpolation_weight(connection.lower_dof,connection.upper_dof,axis,face_coordinate));}
		auto special=[&](const CoarseFinePressureConnection& connection){const int a=connection.direction>0?connection.coarse_dof:connection.fine_dof,b=connection.direction>0?connection.fine_dof:connection.coarse_dof;internal(a,b,connection.axis,connection.open_area,interpolation_weight(a,b,connection.axis,connection.face_centroid[connection.axis]));};
		for(const auto& connection:system.coarse_fine)special(connection);for(const auto& connection:system.embedded)special(connection);
		for(const CompositeSurfacePressurePatch& patch:system.surface_patches)
		{
			if(patch.plus_dof<0||patch.minus_dof<0||patch.plus_dof>=system.storage_size||patch.minus_dof>=system.storage_size||!system.active[patch.plus_dof]||!system.active[patch.minus_dof]||!(patch.area>=0)||!std::isfinite(patch.area)||!std::isfinite(patch.normal.x)||!std::isfinite(patch.normal.y)||!std::isfinite(patch.normal.z))throw std::runtime_error("invalid composite pressure surface patch");const Vec3d coefficient=patch.normal*patch.area;dx[patch.plus_dof]+=scale*pressure[patch.plus_dof]*coefficient.x;dy[patch.plus_dof]+=scale*pressure[patch.plus_dof]*coefficient.y;dz[patch.plus_dof]+=scale*pressure[patch.plus_dof]*coefficient.z;dx[patch.minus_dof]-=scale*pressure[patch.minus_dof]*coefficient.x;dy[patch.minus_dof]-=scale*pressure[patch.minus_dof]*coefficient.y;dz[patch.minus_dof]-=scale*pressure[patch.minus_dof]*coefficient.z;
		}
		for(int q=0;q<system.storage_size;++q)if(system.active[q]){const double inverse=1.0/system.volume[q];state.x[q]=static_cast<Real>(static_cast<double>(state.x[q])+dx[q]*inverse);state.y[q]=static_cast<Real>(static_cast<double>(state.y[q])+dy[q]*inverse);state.z[q]=static_cast<Real>(static_cast<double>(state.z[q])+dz[q]*inverse);}
	}

	void reconstruct_composite_cell_fluxes_cpu(const CompositeAmrPressureSystem& system,
		const CompositeCellMomentumState& state,AmrHostFields& fields,CompositeAmrFluxes& fluxes)
	{
		if(!system.hierarchy||fields.levels().size()!=system.hierarchy->levels().size()||state.x.size()!=static_cast<std::size_t>(system.storage_size)||state.y.size()!=state.x.size()||state.z.size()!=state.x.size())throw std::invalid_argument("invalid composite cell flux reconstruction state");const AmrHierarchy& hierarchy=*system.hierarchy;const int bs=system.brick_size,cells=bs*bs*bs;fluxes=make_zero_composite_fluxes(system);const std::vector<CompositeEbMomentumRegularConnection> compact_regular=build_composite_eb_momentum_regular_connections(system);std::vector<unsigned char> compact_plus(system.storage_size,0);for(const auto& connection:compact_regular)compact_plus[connection.lower_dof]|=static_cast<unsigned char>(1u<<connection.face.component);auto component=[&](int dof,int axis){return static_cast<double>(axis==0?state.x[dof]:(axis==1?state.y[dof]:state.z[dof]));};auto set_positive=[&](AmrHostLevelFields& value,int brick,int axis,int i,int j,int k,double face){if(axis==0)value.u[value.layout.u_index(brick,i+1,j,k)]=static_cast<Real>(face);else if(axis==1)value.v[value.layout.v_index(brick,i,j+1,k)]=static_cast<Real>(face);else value.w[value.layout.w_index(brick,i,j,k+1)]=static_cast<Real>(face);};auto set_mapped=[&](const AmrMacFaceAddress& face,double value){auto& level=fields.levels()[face.level];if(face.component==0)level.u[level.layout.u_index(face.brick,face.i,face.j,face.k)]=static_cast<Real>(value);else if(face.component==1)level.v[level.layout.v_index(face.brick,face.i,face.j,face.k)]=static_cast<Real>(value);else level.w[level.layout.w_index(face.brick,face.i,face.j,face.k)]=static_cast<Real>(value);};auto upper_weight=[&](int a,int b,int axis,double face_coordinate){const double ca=system.centroid[a][axis],cb=system.centroid[b][axis],distance=cb-ca;if(!(std::abs(distance)>1e-14))throw std::runtime_error("composite flux interpolation has coincident centroids");const double weight=(face_coordinate-ca)/distance;if(weight<-1e-8||weight>1+1e-8)throw std::runtime_error("composite flux face lies outside endpoint centroids");return std::clamp(weight,0.0,1.0);};auto connection_value=[&](const CoarseFinePressureConnection& connection){const int a=connection.direction>0?connection.coarse_dof:connection.fine_dof,b=connection.direction>0?connection.fine_dof:connection.coarse_dof,axis=connection.axis;const double weight=upper_weight(a,b,axis,connection.face_centroid[axis]);return (1-weight)*component(a,axis)+weight*component(b,axis);};for(int level_index=0;level_index<static_cast<int>(hierarchy.levels().size());++level_index){const AmrLevel& level=hierarchy.levels()[level_index];auto& values=fields.levels()[level_index];const int offset=system.level_offset[level_index];for(int brick=0;brick<static_cast<int>(level.bricks.size());++brick){const BrickMetadata& record=level.bricks[brick];if(!record.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=offset+brick*cells+i+bs*(j+bs*k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if((system.cut_face_mask[a]&(1u<<axis))||(compact_plus[a]&(1u<<axis)))continue;int c[3]={i,j,k},other_brick=brick;++c[axis];if(c[axis]>=bs){other_brick=record.same_level_neighbor[2*axis+1];if(other_brick<0||!level.bricks[other_brick].active())continue;c[axis]=0;}const int b=offset+other_brick*cells+c[0]+bs*(c[1]+bs*c[2]);if(!system.active[b])continue;const double value=0.5*(component(a,axis)+component(b,axis));set_positive(values,brick,axis,i,j,k,value);if(other_brick!=brick){int lower[3]={c[0],c[1],c[2]};lower[axis]=-1;set_positive(values,other_brick,axis,lower[0],lower[1],lower[2],value);}}}}}for(const auto& connection:compact_regular){const AmrLevel& level=hierarchy.levels()[connection.face.level];const BrickMetadata& record=level.bricks[connection.face.brick];const int coordinate[3]={connection.face.i,connection.face.j,connection.face.k},axis=connection.face.component;const double face_coordinate=record.origin[axis]+coordinate[axis]*level.h,weight=upper_weight(connection.lower_dof,connection.upper_dof,axis,face_coordinate);set_mapped(connection.face,(1-weight)*component(connection.lower_dof,axis)+weight*component(connection.upper_dof,axis));}for(std::size_t edge=0;edge<system.coarse_fine.size();++edge)fluxes.coarse_fine_velocity[edge]=connection_value(system.coarse_fine[edge]);for(std::size_t edge=0;edge<system.embedded.size();++edge)fluxes.embedded_velocity[edge]=connection_value(system.embedded[edge]);
	}

	DeviceCompositeCellMomentumTransport::DeviceCompositeCellMomentumTransport(
		const CompositeAmrPressureSystem& system,DeviceAmrFields& fields):system_(&system),fields_(&fields)
	{
		if(!system.hierarchy||fields.level_count()!=static_cast<int>(system.hierarchy->levels().size())||system.level_offset.size()!=system.hierarchy->levels().size())throw std::invalid_argument("composite cell momentum transport hierarchy mismatch");storage_size_=system.storage_size;const int bs=system.brick_size,cells=bs*bs*bs;if(storage_size_<=0||bs<=0)throw std::invalid_argument("invalid composite cell momentum storage layout");for(int level=0;level<fields.level_count();++level){const int begin=system.level_offset[level],count=static_cast<int>(system.hierarchy->levels()[level].bricks.size())*cells;if(begin<0||count<=0||begin+count>storage_size_)throw std::invalid_argument("invalid composite cell momentum level storage layout");}active_host_=system.active;volume_host_=system.volume;std::vector<Real> volume(storage_size_);for(int q=0;q<storage_size_;++q){if(system.active[q]&&!(system.volume[q]>0))throw std::invalid_argument("active composite cell momentum node has non-positive volume");volume[q]=static_cast<Real>(system.volume[q]);}
		auto interpolation_weight=[&](int a,int b,int axis,double face_coordinate){const double ca=system.centroid[a][axis],cb=system.centroid[b][axis],distance=cb-ca;if(!(std::abs(distance)>1e-14))throw std::invalid_argument("composite cell flux interpolation has coincident centroids");const double weight=(face_coordinate-ca)/distance;if(weight<-1e-8||weight>1+1e-8)throw std::invalid_argument("composite cell flux face lies outside endpoint centroids");return static_cast<Real>(std::clamp(weight,0.0,1.0));};const std::vector<CompositeEbMomentumRegularConnection> regular=build_composite_eb_momentum_regular_connections(system);const CompositeCellGradientTopology gradient_topology=build_composite_cell_gradient_topology(system,regular);gradient_special_count_=static_cast<int>(gradient_topology.special_dof.size());gradient_incidence_count_=static_cast<int>(gradient_topology.incidence_special.size());std::vector<Real> gradient_weighted_displacement(gradient_topology.incidence_weighted_displacement.size()),gradient_inverse(gradient_topology.inverse.size());for(std::size_t q=0;q<gradient_weighted_displacement.size();++q)gradient_weighted_displacement[q]=static_cast<Real>(gradient_topology.incidence_weighted_displacement[q]);for(std::size_t q=0;q<gradient_inverse.size();++q)gradient_inverse[q]=static_cast<Real>(gradient_topology.inverse[q]);regular_count_=static_cast<int>(regular.size());std::vector<int> regular_a(regular_count_),regular_b(regular_count_);std::vector<Real> regular_area(regular_count_),regular_conductance(regular_count_),regular_upper_weight(regular_count_);std::vector<std::int8_t> regular_axis(regular_count_);std::vector<AmrMacFaceAddress> regular_faces;regular_faces.reserve(regular_count_);const std::vector<unsigned char>& compact_plus=gradient_topology.compact_plus;for(int edge=0;edge<regular_count_;++edge){regular_a[edge]=regular[edge].lower_dof;regular_b[edge]=regular[edge].upper_dof;regular_axis[edge]=static_cast<std::int8_t>(regular[edge].face.component);regular_area[edge]=static_cast<Real>(regular[edge].open_area);const Vec3d delta=system.centroid[regular_b[edge]]-system.centroid[regular_a[edge]];const double distance2=length2(delta),normal_distance=std::abs(delta[regular_axis[edge]]);if(!(distance2>0&&normal_distance>0))throw std::invalid_argument("invalid composite cell momentum perimeter conductance");regular_conductance[edge]=static_cast<Real>(regular[edge].open_area*normal_distance/distance2);const AmrLevel& level=system.hierarchy->levels()[regular[edge].face.level];const double face_coordinate=level.bricks[regular[edge].face.brick].origin[regular_axis[edge]]+(regular_axis[edge]==0?regular[edge].face.i:(regular_axis[edge]==1?regular[edge].face.j:regular[edge].face.k))*level.h;regular_upper_weight[edge]=interpolation_weight(regular_a[edge],regular_b[edge],regular_axis[edge],face_coordinate);regular_faces.push_back(regular[edge].face);}regular_face_map_=std::make_unique<DeviceAmrMacFaceMap>(fields,regular_faces);
		auto connection_conductance=[](const CoarseFinePressureConnection& connection){const double value=connection.open_area*pressure_gradient_factor(connection);if(!(value>0)||!std::isfinite(value))throw std::invalid_argument("invalid composite cell momentum special conductance");return static_cast<Real>(value);};coarse_fine_count_=static_cast<int>(system.coarse_fine.size());std::vector<int> coarse_fine_a(coarse_fine_count_),coarse_fine_b(coarse_fine_count_);std::vector<Real> coarse_fine_area(coarse_fine_count_),coarse_fine_conductance(coarse_fine_count_),coarse_fine_upper_weight(coarse_fine_count_);std::vector<std::int8_t> coarse_fine_axis(coarse_fine_count_);for(int edge=0;edge<coarse_fine_count_;++edge){const auto& connection=system.coarse_fine[edge];coarse_fine_a[edge]=connection.direction>0?connection.coarse_dof:connection.fine_dof;coarse_fine_b[edge]=connection.direction>0?connection.fine_dof:connection.coarse_dof;coarse_fine_axis[edge]=connection.axis;coarse_fine_area[edge]=static_cast<Real>(connection.open_area);coarse_fine_conductance[edge]=connection_conductance(connection);coarse_fine_upper_weight[edge]=interpolation_weight(coarse_fine_a[edge],coarse_fine_b[edge],connection.axis,connection.face_centroid[connection.axis]);}
		embedded_count_=static_cast<int>(system.embedded.size());std::vector<int> embedded_a(embedded_count_),embedded_b(embedded_count_);std::vector<Real> embedded_area(embedded_count_),embedded_conductance(embedded_count_),embedded_upper_weight(embedded_count_);std::vector<std::int8_t> embedded_axis(embedded_count_);for(int edge=0;edge<embedded_count_;++edge){const auto& connection=system.embedded[edge];embedded_a[edge]=connection.direction>0?connection.coarse_dof:connection.fine_dof;embedded_b[edge]=connection.direction>0?connection.fine_dof:connection.coarse_dof;embedded_axis[edge]=connection.axis;embedded_area[edge]=static_cast<Real>(connection.open_area);embedded_conductance[edge]=connection_conductance(connection);embedded_upper_weight[edge]=interpolation_weight(embedded_a[edge],embedded_b[edge],connection.axis,connection.face_centroid[connection.axis]);}
		std::vector<int> surface_dof;std::vector<Real> surface_coefficient;surface_dof.reserve(2*system.surface_patches.size());surface_coefficient.reserve(6*system.surface_patches.size());for(const CompositeSurfacePressurePatch& patch:system.surface_patches){if(patch.plus_dof<0||patch.minus_dof<0||patch.plus_dof>=storage_size_||patch.minus_dof>=storage_size_||!system.active[patch.plus_dof]||!system.active[patch.minus_dof]||!(patch.area>=0)||!std::isfinite(patch.area)||!std::isfinite(patch.normal.x)||!std::isfinite(patch.normal.y)||!std::isfinite(patch.normal.z))throw std::invalid_argument("invalid composite cell momentum surface patch");const Vec3d coefficient=patch.normal*patch.area;surface_dof.push_back(patch.plus_dof);surface_coefficient.insert(surface_coefficient.end(),{static_cast<Real>(coefficient.x),static_cast<Real>(coefficient.y),static_cast<Real>(coefficient.z)});surface_dof.push_back(patch.minus_dof);surface_coefficient.insert(surface_coefficient.end(),{static_cast<Real>(-coefficient.x),static_cast<Real>(-coefficient.y),static_cast<Real>(-coefficient.z)});}surface_count_=static_cast<int>(surface_dof.size());
		x_=allocate<Real>(storage_size_,"allocate composite cell momentum x");y_=allocate<Real>(storage_size_,"allocate composite cell momentum y");z_=allocate<Real>(storage_size_,"allocate composite cell momentum z");delta_x_=allocate<Real>(storage_size_,"allocate composite cell momentum dx");delta_y_=allocate<Real>(storage_size_,"allocate composite cell momentum dy");delta_z_=allocate<Real>(storage_size_,"allocate composite cell momentum dz");viscosity_=allocate<Real>(storage_size_,"allocate composite cell momentum viscosity");volume_=upload(volume,"upload composite cell momentum volumes");active_=upload(system.active,"upload composite cell momentum active mask");cut_face_mask_=upload(system.cut_face_mask,"upload composite cell momentum cut-face mask");compact_plus_mask_=upload(compact_plus,"upload composite cell momentum compact ownership");gradient_special_mask_=upload(gradient_topology.special_mask,"upload composite cell momentum gradient special mask");gradient_special_dof_=upload(gradient_topology.special_dof,"upload composite cell momentum gradient nodes");gradient_incidence_special_=upload(gradient_topology.incidence_special,"upload composite cell momentum gradient incidence nodes");gradient_incidence_neighbor_=upload(gradient_topology.incidence_neighbor,"upload composite cell momentum gradient incidence neighbours");gradient_incidence_weighted_displacement_=upload(gradient_weighted_displacement,"upload composite cell momentum gradient displacements");gradient_inverse_=upload(gradient_inverse,"upload composite cell momentum gradient inverses");gradient_rhs_=allocate<Real>(static_cast<std::size_t>(gradient_special_count_)*9,"allocate composite cell momentum gradient rhs");regular_velocity_=allocate<Real>(regular_count_,"allocate composite cell momentum perimeter flux");regular_a_=upload(regular_a,"upload composite cell momentum perimeter lower nodes");regular_b_=upload(regular_b,"upload composite cell momentum perimeter upper nodes");regular_axis_=upload(regular_axis,"upload composite cell momentum perimeter axes");regular_area_=upload(regular_area,"upload composite cell momentum perimeter areas");regular_conductance_=upload(regular_conductance,"upload composite cell momentum perimeter conductances");regular_upper_weight_=upload(regular_upper_weight,"upload composite cell momentum perimeter interpolation");coarse_fine_a_=upload(coarse_fine_a,"upload composite cell momentum coarse/fine lower nodes");coarse_fine_b_=upload(coarse_fine_b,"upload composite cell momentum coarse/fine upper nodes");coarse_fine_axis_=upload(coarse_fine_axis,"upload composite cell momentum coarse/fine axes");coarse_fine_area_=upload(coarse_fine_area,"upload composite cell momentum coarse/fine areas");coarse_fine_conductance_=upload(coarse_fine_conductance,"upload composite cell momentum coarse/fine conductances");coarse_fine_upper_weight_=upload(coarse_fine_upper_weight,"upload composite cell momentum coarse/fine interpolation");embedded_a_=upload(embedded_a,"upload composite cell momentum aperture lower nodes");embedded_b_=upload(embedded_b,"upload composite cell momentum aperture upper nodes");embedded_axis_=upload(embedded_axis,"upload composite cell momentum aperture axes");embedded_area_=upload(embedded_area,"upload composite cell momentum aperture areas");embedded_conductance_=upload(embedded_conductance,"upload composite cell momentum aperture conductances");embedded_upper_weight_=upload(embedded_upper_weight,"upload composite cell momentum aperture interpolation");surface_dof_=upload(surface_dof,"upload composite cell momentum surface nodes");surface_coefficient_=upload(surface_coefficient,"upload composite cell momentum surface coefficients");bytes_=static_cast<std::size_t>(8*storage_size_+4*regular_count_+3*coarse_fine_count_+3*embedded_count_+3*surface_count_+18*gradient_special_count_+3*gradient_incidence_count_)*sizeof(Real)+static_cast<std::size_t>(4*storage_size_+regular_count_+coarse_fine_count_+embedded_count_)*sizeof(unsigned char)+static_cast<std::size_t>(2*regular_count_+2*coarse_fine_count_+2*embedded_count_+surface_count_+gradient_special_count_+2*gradient_incidence_count_)*sizeof(int);
	}

	DeviceCompositeCellMomentumTransport::~DeviceCompositeCellMomentumTransport()
	{
		for(void* pointer:{(void*)x_,(void*)y_,(void*)z_,(void*)delta_x_,(void*)delta_y_,(void*)delta_z_,(void*)viscosity_,(void*)volume_,(void*)regular_velocity_,(void*)active_,(void*)cut_face_mask_,(void*)compact_plus_mask_,(void*)gradient_special_mask_,(void*)gradient_special_dof_,(void*)gradient_incidence_special_,(void*)gradient_incidence_neighbor_,(void*)gradient_incidence_weighted_displacement_,(void*)gradient_inverse_,(void*)gradient_rhs_,(void*)coarse_fine_a_,(void*)coarse_fine_b_,(void*)coarse_fine_axis_,(void*)coarse_fine_area_,(void*)coarse_fine_conductance_,(void*)coarse_fine_upper_weight_,(void*)embedded_a_,(void*)embedded_b_,(void*)embedded_axis_,(void*)embedded_area_,(void*)embedded_conductance_,(void*)embedded_upper_weight_,(void*)regular_a_,(void*)regular_b_,(void*)regular_axis_,(void*)regular_area_,(void*)regular_conductance_,(void*)regular_upper_weight_,(void*)surface_dof_,(void*)surface_coefficient_})if(pointer)cudaFree(pointer);
	}

	void DeviceCompositeCellMomentumTransport::upload_state(const CompositeCellMomentumState& state)
	{
		if(state.x.size()!=static_cast<std::size_t>(storage_size_)||state.y.size()!=state.x.size()||state.z.size()!=state.x.size())throw std::invalid_argument("composite cell momentum upload size mismatch");check(cudaMemcpy(x_,state.x.data(),state.x.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload composite cell momentum x");check(cudaMemcpy(y_,state.y.data(),state.y.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload composite cell momentum y");check(cudaMemcpy(z_,state.z.data(),state.z.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload composite cell momentum z");
	}

	void DeviceCompositeCellMomentumTransport::download_state(CompositeCellMomentumState& state)const
	{
		state.x.resize(storage_size_);state.y.resize(storage_size_);state.z.resize(storage_size_);check(cudaMemcpy(state.x.data(),x_,state.x.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download composite cell momentum x");check(cudaMemcpy(state.y.data(),y_,state.y.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download composite cell momentum y");check(cudaMemcpy(state.z.data(),z_,state.z.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download composite cell momentum z");
	}

	void DeviceCompositeCellMomentumTransport::step(const Real* coarse_fine_velocity,
		const Real* embedded_velocity,Real dt,bool external_aero,Real freestream_speed)
	{
		if((coarse_fine_count_&&!coarse_fine_velocity)||(embedded_count_&&!embedded_velocity)||!(dt>Real(0))||!(freestream_speed>=Real(0)))throw std::invalid_argument("invalid composite cell momentum timestep");check(cudaMemset(delta_x_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell momentum dx");check(cudaMemset(delta_y_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell momentum dy");check(cudaMemset(delta_z_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell momentum dz");const int bs=system_->brick_size,cells=bs*bs*bs;for(int level=0;level<fields_->level_count();++level){const DeviceAmrFieldLevelView view=fields_->level_view(level);const int count=view.brick_count*cells;const Real h=static_cast<Real>(system_->hierarchy->levels()[level].h);transport_regular_cell_momentum_kernel<<<(count+255)/256,256>>>(view,h,system_->level_offset[level],active_,cut_face_mask_,compact_plus_mask_,x_,y_,z_,dt,external_aero,freestream_speed,delta_x_,delta_y_,delta_z_);}if(regular_count_){regular_face_map_->gather(regular_velocity_);transport_compact_cell_momentum_kernel<<<(regular_count_+255)/256,256>>>(regular_a_,regular_b_,regular_area_,regular_velocity_,x_,y_,z_,dt,delta_x_,delta_y_,delta_z_,regular_count_);}if(coarse_fine_count_)transport_compact_cell_momentum_kernel<<<(coarse_fine_count_+255)/256,256>>>(coarse_fine_a_,coarse_fine_b_,coarse_fine_area_,coarse_fine_velocity,x_,y_,z_,dt,delta_x_,delta_y_,delta_z_,coarse_fine_count_);if(embedded_count_)transport_compact_cell_momentum_kernel<<<(embedded_count_+255)/256,256>>>(embedded_a_,embedded_b_,embedded_area_,embedded_velocity,x_,y_,z_,dt,delta_x_,delta_y_,delta_z_,embedded_count_);apply_cell_momentum_kernel<<<(storage_size_+255)/256,256>>>(active_,volume_,delta_x_,delta_y_,delta_z_,x_,y_,z_,storage_size_);check(cudaDeviceSynchronize(),"composite cell momentum transport");
	}

	void DeviceCompositeCellMomentumTransport::step(const Real* embedded_velocity,Real dt,
		bool external_aero,Real freestream_speed)
	{
		if(coarse_fine_count_)throw std::invalid_argument("coarse/fine velocity is required for multilevel composite cell momentum transport");step(nullptr,embedded_velocity,dt,external_aero,freestream_speed);
	}

	void DeviceCompositeCellMomentumTransport::diffuse(Real kinematic_viscosity,Real dt)
	{
		if(kinematic_viscosity<Real(0)||!(dt>Real(0)))throw std::invalid_argument("invalid composite cell momentum diffusion timestep");if(kinematic_viscosity==Real(0))return;check(cudaMemset(delta_x_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell momentum diffusion dx");check(cudaMemset(delta_y_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell momentum diffusion dy");check(cudaMemset(delta_z_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell momentum diffusion dz");const int bs=system_->brick_size,cells=bs*bs*bs;for(int level=0;level<fields_->level_count();++level){const DeviceAmrFieldLevelView view=fields_->level_view(level);const int count=view.brick_count*cells;const Real h=static_cast<Real>(system_->hierarchy->levels()[level].h);diffuse_regular_cell_momentum_kernel<<<(count+255)/256,256>>>(view,h,system_->level_offset[level],active_,cut_face_mask_,compact_plus_mask_,x_,y_,z_,kinematic_viscosity*dt,delta_x_,delta_y_,delta_z_);}if(regular_count_)diffuse_compact_cell_momentum_kernel<<<(regular_count_+255)/256,256>>>(regular_a_,regular_b_,regular_conductance_,x_,y_,z_,kinematic_viscosity*dt,delta_x_,delta_y_,delta_z_,regular_count_);if(coarse_fine_count_)diffuse_compact_cell_momentum_kernel<<<(coarse_fine_count_+255)/256,256>>>(coarse_fine_a_,coarse_fine_b_,coarse_fine_conductance_,x_,y_,z_,kinematic_viscosity*dt,delta_x_,delta_y_,delta_z_,coarse_fine_count_);if(embedded_count_)diffuse_compact_cell_momentum_kernel<<<(embedded_count_+255)/256,256>>>(embedded_a_,embedded_b_,embedded_conductance_,x_,y_,z_,kinematic_viscosity*dt,delta_x_,delta_y_,delta_z_,embedded_count_);apply_cell_momentum_kernel<<<(storage_size_+255)/256,256>>>(active_,volume_,delta_x_,delta_y_,delta_z_,x_,y_,z_,storage_size_);check(cudaDeviceSynchronize(),"composite cell momentum diffusion");
	}

	void DeviceCompositeCellMomentumTransport::diffuse_smagorinsky(Real molecular_viscosity,
		Real smagorinsky_cs,Real dt)
	{
		if(molecular_viscosity<Real(0)||smagorinsky_cs<Real(0)||!(dt>Real(0)))throw std::invalid_argument("invalid composite cell Smagorinsky timestep");if(smagorinsky_cs==Real(0)){diffuse(molecular_viscosity,dt);return;}if(gradient_special_count_)check(cudaMemset(gradient_rhs_,0,static_cast<std::size_t>(gradient_special_count_)*9*sizeof(Real)),"clear composite cell momentum gradient rhs");const int bs=system_->brick_size,cells=bs*bs*bs;for(int level=0;level<fields_->level_count();++level){const DeviceAmrFieldLevelView view=fields_->level_view(level);const int count=view.brick_count*cells;const Real h=static_cast<Real>(system_->hierarchy->levels()[level].h);compute_regular_cell_viscosity_kernel<<<(count+255)/256,256>>>(view,h,system_->level_offset[level],active_,cut_face_mask_,compact_plus_mask_,gradient_special_mask_,volume_,x_,y_,z_,molecular_viscosity,smagorinsky_cs,viscosity_);}if(gradient_incidence_count_)accumulate_special_cell_gradient_kernel<<<(gradient_incidence_count_+255)/256,256>>>(gradient_incidence_special_,gradient_special_dof_,gradient_incidence_neighbor_,gradient_incidence_weighted_displacement_,x_,y_,z_,gradient_rhs_,gradient_incidence_count_);if(gradient_special_count_)solve_special_cell_viscosity_kernel<<<(gradient_special_count_+255)/256,256>>>(gradient_special_dof_,gradient_inverse_,gradient_rhs_,volume_,molecular_viscosity,smagorinsky_cs,viscosity_,gradient_special_count_);check(cudaMemset(delta_x_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell Smagorinsky dx");check(cudaMemset(delta_y_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell Smagorinsky dy");check(cudaMemset(delta_z_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell Smagorinsky dz");for(int level=0;level<fields_->level_count();++level){const DeviceAmrFieldLevelView view=fields_->level_view(level);const int count=view.brick_count*cells;const Real h=static_cast<Real>(system_->hierarchy->levels()[level].h);diffuse_regular_cell_momentum_variable_kernel<<<(count+255)/256,256>>>(view,h,system_->level_offset[level],active_,cut_face_mask_,compact_plus_mask_,viscosity_,x_,y_,z_,dt,delta_x_,delta_y_,delta_z_);}if(regular_count_)diffuse_compact_cell_momentum_variable_kernel<<<(regular_count_+255)/256,256>>>(regular_a_,regular_b_,regular_conductance_,viscosity_,x_,y_,z_,dt,delta_x_,delta_y_,delta_z_,regular_count_);if(coarse_fine_count_)diffuse_compact_cell_momentum_variable_kernel<<<(coarse_fine_count_+255)/256,256>>>(coarse_fine_a_,coarse_fine_b_,coarse_fine_conductance_,viscosity_,x_,y_,z_,dt,delta_x_,delta_y_,delta_z_,coarse_fine_count_);if(embedded_count_)diffuse_compact_cell_momentum_variable_kernel<<<(embedded_count_+255)/256,256>>>(embedded_a_,embedded_b_,embedded_conductance_,viscosity_,x_,y_,z_,dt,delta_x_,delta_y_,delta_z_,embedded_count_);apply_cell_momentum_kernel<<<(storage_size_+255)/256,256>>>(active_,volume_,delta_x_,delta_y_,delta_z_,x_,y_,z_,storage_size_);check(cudaDeviceSynchronize(),"composite cell Smagorinsky diffusion");
	}

	void DeviceCompositeCellMomentumTransport::apply_pressure_impulse(const Real* pressure,
		Real dt,Real density,bool include_physical_boundaries)
	{
		if(!pressure||!(dt>Real(0))||!(density>Real(0)))throw std::invalid_argument("invalid composite cell pressure impulse");check(cudaMemset(delta_x_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell pressure impulse dx");check(cudaMemset(delta_y_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell pressure impulse dy");check(cudaMemset(delta_z_,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite cell pressure impulse dz");const Real scale=dt/density;const int bs=system_->brick_size,cells=bs*bs*bs;for(int level=0;level<fields_->level_count();++level){const DeviceAmrFieldLevelView view=fields_->level_view(level);const int count=view.brick_count*cells;const Real h=static_cast<Real>(system_->hierarchy->levels()[level].h);pressure_regular_cell_impulse_kernel<<<(count+255)/256,256>>>(view,h,system_->level_offset[level],active_,cut_face_mask_,compact_plus_mask_,pressure,scale,system_->pressure_outlet_xmax,include_physical_boundaries,delta_x_,delta_y_,delta_z_);}if(regular_count_)pressure_compact_cell_impulse_kernel<<<(regular_count_+255)/256,256>>>(regular_a_,regular_b_,regular_axis_,regular_area_,regular_upper_weight_,pressure,scale,delta_x_,delta_y_,delta_z_,regular_count_);if(coarse_fine_count_)pressure_compact_cell_impulse_kernel<<<(coarse_fine_count_+255)/256,256>>>(coarse_fine_a_,coarse_fine_b_,coarse_fine_axis_,coarse_fine_area_,coarse_fine_upper_weight_,pressure,scale,delta_x_,delta_y_,delta_z_,coarse_fine_count_);if(embedded_count_)pressure_compact_cell_impulse_kernel<<<(embedded_count_+255)/256,256>>>(embedded_a_,embedded_b_,embedded_axis_,embedded_area_,embedded_upper_weight_,pressure,scale,delta_x_,delta_y_,delta_z_,embedded_count_);if(surface_count_)pressure_surface_cell_impulse_kernel<<<(surface_count_+255)/256,256>>>(surface_dof_,surface_coefficient_,pressure,scale,delta_x_,delta_y_,delta_z_,surface_count_);apply_cell_momentum_kernel<<<(storage_size_+255)/256,256>>>(active_,volume_,delta_x_,delta_y_,delta_z_,x_,y_,z_,storage_size_);check(cudaDeviceSynchronize(),"composite cell pressure impulse");
	}

	void DeviceCompositeCellMomentumTransport::reconstruct_fluxes(Real* coarse_fine_velocity,
		Real* embedded_velocity)
	{
		if((coarse_fine_count_&&!coarse_fine_velocity)||(embedded_count_&&!embedded_velocity))throw std::invalid_argument("invalid composite cell momentum flux reconstruction");const int bs=system_->brick_size,cells=bs*bs*bs;for(int level=0;level<fields_->level_count();++level){const DeviceAmrFieldLevelView view=fields_->level_view(level);const int count=view.brick_count*cells;reconstruct_regular_cell_flux_kernel<<<(count+255)/256,256>>>(view,system_->level_offset[level],active_,cut_face_mask_,compact_plus_mask_,x_,y_,z_);}if(regular_count_){reconstruct_compact_cell_flux_kernel<<<(regular_count_+255)/256,256>>>(regular_a_,regular_b_,regular_axis_,regular_upper_weight_,x_,y_,z_,regular_velocity_,regular_count_);regular_face_map_->scatter(regular_velocity_);}if(coarse_fine_count_)reconstruct_compact_cell_flux_kernel<<<(coarse_fine_count_+255)/256,256>>>(coarse_fine_a_,coarse_fine_b_,coarse_fine_axis_,coarse_fine_upper_weight_,x_,y_,z_,coarse_fine_velocity,coarse_fine_count_);if(embedded_count_)reconstruct_compact_cell_flux_kernel<<<(embedded_count_+255)/256,256>>>(embedded_a_,embedded_b_,embedded_axis_,embedded_upper_weight_,x_,y_,z_,embedded_velocity,embedded_count_);check(cudaDeviceSynchronize(),"reconstruct composite cell momentum fluxes");
	}

	std::array<double,3> DeviceCompositeCellMomentumTransport::momentum()const
	{
		CompositeCellMomentumState state;download_state(state);std::array<double,3> result{};for(int q=0;q<storage_size_;++q)if(active_host_[q]){result[0]+=volume_host_[q]*static_cast<double>(state.x[q]);result[1]+=volume_host_[q]*static_cast<double>(state.y[q]);result[2]+=volume_host_[q]*static_cast<double>(state.z[q]);}return result;
	}

	std::size_t DeviceCompositeCellMomentumTransport::bytes()const{return bytes_+(regular_face_map_?regular_face_map_->bytes():0);}

	std::vector<NormalMomentumInterfaceTile> build_normal_momentum_interface_tiles(
		const CompositeAmrPressureSystem& system)
	{
		if(!system.hierarchy||system.brick_size<=0)throw std::invalid_argument("normal momentum tiles require a composite AMR hierarchy");const int bs=system.brick_size,cells=bs*bs*bs;
		auto cell_address=[&](int dof)
		{
			AmrMacFaceAddress out;for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level){const int begin=system.level_offset[level],count=static_cast<int>(system.hierarchy->levels()[level].bricks.size())*cells;if(dof<begin||dof>=begin+count)continue;const int work=dof-begin,local=work%cells;out.level=level;out.brick=work/cells;out.i=local%bs;out.j=(local/bs)%bs;out.k=local/(bs*bs);return out;}throw std::invalid_argument("2:1 momentum tile pressure DOF is not a regular AMR cell");
		};
		std::vector<NormalMomentumInterfaceTile> tiles;tiles.reserve(system.coarse_fine.size());
		for(const CoarseFinePressureConnection& connection:system.coarse_fine)
		{
			AmrMacFaceAddress coarse=cell_address(connection.coarse_dof),fine=cell_address(connection.fine_dof);if(fine.level!=coarse.level+1||connection.axis<0||connection.axis>2||(connection.direction!=1&&connection.direction!=-1))throw std::invalid_argument("invalid 2:1 normal momentum topology");
			coarse.component=fine.component=connection.axis;NormalMomentumInterfaceTile tile;tile.coarse_interface=tile.coarse_interior=coarse;tile.fine_interface=tile.fine_interior=fine;tile.direction=connection.direction;int* coarse_interface=connection.axis==0?&tile.coarse_interface.i:(connection.axis==1?&tile.coarse_interface.j:&tile.coarse_interface.k);int* coarse_interior=connection.axis==0?&tile.coarse_interior.i:(connection.axis==1?&tile.coarse_interior.j:&tile.coarse_interior.k);int* fine_interface=connection.axis==0?&tile.fine_interface.i:(connection.axis==1?&tile.fine_interface.j:&tile.fine_interface.k);int* fine_interior=connection.axis==0?&tile.fine_interior.i:(connection.axis==1?&tile.fine_interior.j:&tile.fine_interior.k);*coarse_interface+=connection.direction>0?1:0;*coarse_interior+=connection.direction>0?0:1;*fine_interface+=connection.direction>0?0:1;*fine_interior+=connection.direction>0?1:0;const double hc=system.hierarchy->levels()[coarse.level].h,hf=system.hierarchy->levels()[fine.level].h;if(std::abs(hc-2*hf)>1e-10*hc)throw std::invalid_argument("normal momentum interface is not 2:1");tile.area=connection.open_area;tile.dual_volume=connection.open_area*0.5*(hc+hf);if(!(tile.area>0)||!(tile.dual_volume>0))throw std::invalid_argument("normal momentum interface has non-positive measure");tiles.push_back(tile);
		}
		return tiles;
	}

	std::vector<TangentialMomentumInterfaceConnection> build_tangential_momentum_interface_connections(
		const CompositeAmrPressureSystem& system)
	{
		if(!system.hierarchy||system.brick_size<=0)throw std::invalid_argument("tangential momentum connections require a composite AMR hierarchy");const AmrHierarchy& hierarchy=*system.hierarchy;const int bs=system.brick_size,cells=bs*bs*bs;
		auto cell_address=[&](int dof){AmrMacFaceAddress out;for(int level=0;level<static_cast<int>(hierarchy.levels().size());++level){const int begin=system.level_offset[level],count=static_cast<int>(hierarchy.levels()[level].bricks.size())*cells;if(dof<begin||dof>=begin+count)continue;const int work=dof-begin,local=work%cells;out.level=level;out.brick=work/cells;out.i=local%bs;out.j=(local/bs)%bs;out.k=local/(bs*bs);return out;}throw std::invalid_argument("2:1 tangential momentum pressure DOF is not a regular AMR cell");};
		auto centre=[&](const AmrMacFaceAddress& address){const AmrLevel& level=hierarchy.levels()[address.level];const BrickMetadata& brick=level.bricks[address.brick];return brick.origin+Vec3d{(address.i+(address.component==0?0.0:0.5))*level.h,(address.j+(address.component==1?0.0:0.5))*level.h,(address.k+(address.component==2?0.0:0.5))*level.h};};
		struct Key{AmrMacFaceAddress coarse,fine;bool operator==(const Key& other)const{return coarse.level==other.coarse.level&&coarse.brick==other.coarse.brick&&coarse.component==other.coarse.component&&coarse.i==other.coarse.i&&coarse.j==other.coarse.j&&coarse.k==other.coarse.k&&fine.level==other.fine.level&&fine.brick==other.fine.brick&&fine.component==other.fine.component&&fine.i==other.fine.i&&fine.j==other.fine.j&&fine.k==other.fine.k;}};
		struct KeyHash{std::size_t operator()(const Key& key)const{std::size_t h=1469598103934665603ull;auto add=[&](int value){h^=static_cast<std::uint32_t>(value);h*=1099511628211ull;};for(int value:{key.coarse.level,key.coarse.brick,key.coarse.component,key.coarse.i,key.coarse.j,key.coarse.k,key.fine.level,key.fine.brick,key.fine.component,key.fine.i,key.fine.j,key.fine.k})add(value);return h;}};
		std::vector<TangentialMomentumInterfaceConnection> output;std::unordered_map<Key,std::size_t,KeyHash> index;
		for(int pressure_index=0;pressure_index<static_cast<int>(system.coarse_fine.size());++pressure_index)
		{
			const CoarseFinePressureConnection& pressure=system.coarse_fine[pressure_index];
			const AmrMacFaceAddress coarse_cell=cell_address(pressure.coarse_dof),fine_cell=cell_address(pressure.fine_dof);if(fine_cell.level!=coarse_cell.level+1)throw std::invalid_argument("invalid tangential 2:1 level pair");const double hc=hierarchy.levels()[coarse_cell.level].h,hf=hierarchy.levels()[fine_cell.level].h;if(std::abs(hc-2*hf)>1e-10*hc)throw std::invalid_argument("tangential momentum interface is not 2:1");
			for(int component=0;component<3;++component)if(component!=pressure.axis)for(int coarse_side=0;coarse_side<2;++coarse_side)for(int fine_side=0;fine_side<2;++fine_side)
			{
				AmrMacFaceAddress coarse=coarse_cell,fine=fine_cell;coarse.component=fine.component=component;int* coarse_coordinate=component==0?&coarse.i:(component==1?&coarse.j:&coarse.k);int* fine_coordinate=component==0?&fine.i:(component==1?&fine.j:&fine.k);*coarse_coordinate+=coarse_side;*fine_coordinate+=fine_side;coarse=canonical_amr_mac_face_address(hierarchy,coarse);fine=canonical_amr_mac_face_address(hierarchy,fine);const Vec3d coarse_centre=centre(coarse),fine_centre=centre(fine);double overlap=1;
				for(int axis=0;axis<3;++axis)if(axis!=pressure.axis){const double tile_lo=pressure.face_centroid[axis]-0.5*hf,tile_hi=pressure.face_centroid[axis]+0.5*hf,coarse_lo=coarse_centre[axis]-0.5*hc,coarse_hi=coarse_centre[axis]+0.5*hc,fine_lo=fine_centre[axis]-0.5*hf,fine_hi=fine_centre[axis]+0.5*hf;overlap*=std::max(0.0,std::min({tile_hi,coarse_hi,fine_hi})-std::max({tile_lo,coarse_lo,fine_lo}));}if(overlap<=1e-14*hf*hf)continue;
				const Key key{coarse,fine};const auto found=index.find(key);if(found==index.end()){index.emplace(key,output.size());output.push_back({coarse,fine,overlap,0.5*(hc+hf),static_cast<std::int8_t>(pressure.axis),static_cast<std::int8_t>(component),pressure.direction,{{pressure_index,overlap}}});}else{TangentialMomentumInterfaceConnection& existing=output[found->second];if(existing.interface_axis!=pressure.axis||existing.component!=component||existing.direction!=pressure.direction)throw std::runtime_error("inconsistent duplicate tangential momentum overlap");existing.open_area+=overlap;existing.flux_sources.push_back({pressure_index,overlap});}
			}
		}
		return output;
	}

	std::vector<AmrMacFaceDirectionalLink> build_momentum_interface_replaced_links(
		const CompositeAmrPressureSystem& system)
	{
		const std::vector<NormalMomentumInterfaceTile> normal=build_normal_momentum_interface_tiles(system);const std::vector<TangentialMomentumInterfaceConnection> tangential=build_tangential_momentum_interface_connections(system);std::vector<AmrMacFaceDirectionalLink> output;
		auto add=[&](const AmrMacFaceAddress& face,int axis,int direction){for(const auto& old:output)if(old.face.level==face.level&&old.face.brick==face.brick&&old.face.component==face.component&&old.face.i==face.i&&old.face.j==face.j&&old.face.k==face.k&&old.transport_axis==axis&&old.direction==direction)return;output.push_back({face,static_cast<std::int8_t>(axis),static_cast<std::int8_t>(direction)});};
		for(const auto& tile:normal){const int axis=tile.fine_interface.component,direction=tile.direction;add(tile.coarse_interior,axis,direction);add(tile.coarse_interface,axis,-direction);add(tile.fine_interface,axis,direction);add(tile.fine_interior,axis,-direction);}for(const auto& connection:tangential){add(connection.coarse,connection.interface_axis,connection.direction);add(connection.fine,connection.interface_axis,-connection.direction);}return output;
	}

	CompositeAmrMomentumInterfaceTopology build_composite_amr_momentum_interface_topology(
		const CompositeAmrPressureSystem& system)
	{
		if(!system.hierarchy)throw std::invalid_argument("composite momentum topology requires a hierarchy");const AmrHierarchy& hierarchy=*system.hierarchy;CompositeAmrMomentumInterfaceTopology output;
		struct AddressKey{int level,brick,component,i,j,k;bool operator==(const AddressKey& other)const{return level==other.level&&brick==other.brick&&component==other.component&&i==other.i&&j==other.j&&k==other.k;}};struct AddressHash{std::size_t operator()(const AddressKey& key)const{std::size_t h=1469598103934665603ull;for(int value:{key.level,key.brick,key.component,key.i,key.j,key.k}){h^=static_cast<std::uint32_t>(value);h*=1099511628211ull;}return h;}};std::unordered_map<AddressKey,int,AddressHash> node_index;std::vector<unsigned char> normal_owner;
		auto add_node=[&](AmrMacFaceAddress address,double volume,bool normal)
		{
			address=canonical_amr_mac_face_address(hierarchy,address);if(!(volume>0)||!std::isfinite(volume))throw std::invalid_argument("composite momentum node has non-positive dual volume");const AddressKey key{address.level,address.brick,address.component,address.i,address.j,address.k};const auto found=node_index.find(key);if(found==node_index.end()){const int node=static_cast<int>(output.nodes.size());node_index.emplace(key,node);output.nodes.push_back(address);output.dual_volume.push_back(volume);output.interface_axis_mask.push_back(0);normal_owner.push_back(normal?1:0);return node;}const int node=found->second;if(normal&&!normal_owner[node]){output.dual_volume[node]=volume;normal_owner[node]=1;}else if(normal==static_cast<bool>(normal_owner[node])&&std::abs(output.dual_volume[node]-volume)>1e-10*std::max(output.dual_volume[node],volume))throw std::runtime_error("composite momentum node has conflicting dual-volume ownership");return node;
		};
		struct EdgeKey{int lower,upper,axis;bool operator==(const EdgeKey& other)const{return lower==other.lower&&upper==other.upper&&axis==other.axis;}};struct EdgeHash{std::size_t operator()(const EdgeKey& key)const{return (static_cast<std::size_t>(key.lower)*1315423911u^static_cast<std::size_t>(key.upper))*31u+static_cast<unsigned>(key.axis);}};std::unordered_map<EdgeKey,int,EdgeHash> edge_index;
		auto add_connection=[&](int first,int second,double area,int axis,int direction)
		{
			if(first<0||second<0||first==second||!(area>0)||axis<0||axis>2||(direction!=1&&direction!=-1))throw std::invalid_argument("invalid composite momentum interface connection");const int lower=direction>0?first:second,upper=direction>0?second:first;const int component=output.nodes[lower].component;if(component!=output.nodes[upper].component)throw std::runtime_error("composite momentum connection mixes velocity components");const EdgeKey key{lower,upper,axis};const auto found=edge_index.find(key);int edge;if(found==edge_index.end()){edge=static_cast<int>(output.connections.size());edge_index.emplace(key,edge);output.connections.push_back({lower,upper,area,static_cast<std::int8_t>(axis),static_cast<std::int8_t>(component)});}else{edge=found->second;output.connections[edge].open_area+=area;}output.interface_axis_mask[lower]|=static_cast<std::uint8_t>(1u<<axis);output.interface_axis_mask[upper]|=static_cast<std::uint8_t>(1u<<axis);return edge;
		};
		const std::vector<NormalMomentumInterfaceTile> normal=build_normal_momentum_interface_tiles(system);for(const NormalMomentumInterfaceTile& tile:normal)
		{
			const AmrMacFaceAddress coarse_interior=canonical_amr_mac_face_address(hierarchy,tile.coarse_interior),fine_interface=canonical_amr_mac_face_address(hierarchy,tile.fine_interface),fine_interior=canonical_amr_mac_face_address(hierarchy,tile.fine_interior),coarse_alias=canonical_amr_mac_face_address(hierarchy,tile.coarse_interface);const int coarse=add_node(coarse_interior,regular_mac_dual_volume(hierarchy,coarse_interior),false),interface_node=add_node(fine_interface,tile.dual_volume,true),fine=add_node(fine_interior,regular_mac_dual_volume(hierarchy,fine_interior),false),axis=tile.fine_interface.component;add_connection(coarse,interface_node,tile.area,axis,tile.direction);add_connection(interface_node,fine,tile.area,axis,tile.direction);output.coarse_aliases.push_back({coarse_alias,interface_node,tile.area});
		}
		const std::vector<TangentialMomentumInterfaceConnection> tangential=build_tangential_momentum_interface_connections(system);for(const TangentialMomentumInterfaceConnection& connection:tangential)
		{
			const AmrMacFaceAddress coarse=canonical_amr_mac_face_address(hierarchy,connection.coarse),fine=canonical_amr_mac_face_address(hierarchy,connection.fine);const int coarse_node=add_node(coarse,regular_mac_dual_volume(hierarchy,coarse),false),fine_node=add_node(fine,regular_mac_dual_volume(hierarchy,fine),false),edge=add_connection(coarse_node,fine_node,connection.open_area,connection.interface_axis,connection.direction);for(const auto& source:connection.flux_sources)output.mass_flux_sources.push_back({edge,source.coarse_fine_connection,source.overlap_area});
		}
		for(const CompositeAmrMomentumCoarseAlias& alias:output.coarse_aliases){const AmrMacFaceAddress address=canonical_amr_mac_face_address(hierarchy,alias.coarse_face);if(node_index.find({address.level,address.brick,address.component,address.i,address.j,address.k})!=node_index.end())throw std::runtime_error("coarse normal momentum alias is also independently owned at a refinement corner");}
		return output;
	}

	std::vector<double> composite_amr_momentum_connection_velocities_cpu(
		const CompositeAmrMomentumInterfaceTopology& topology,const std::vector<double>& node_velocity,
		const std::vector<double>& coarse_fine_velocity)
	{
		if(node_velocity.size()!=topology.nodes.size())throw std::invalid_argument("composite momentum node velocity size mismatch");std::vector<double> velocity(topology.connections.size(),0),source_area(topology.connections.size(),0);for(int edge=0;edge<static_cast<int>(topology.connections.size());++edge){const auto& connection=topology.connections[edge];if(connection.lower_node<0||connection.upper_node<0||connection.lower_node>=static_cast<int>(node_velocity.size())||connection.upper_node>=static_cast<int>(node_velocity.size()))throw std::invalid_argument("composite momentum connection endpoint outside node state");if(connection.component==connection.transport_axis)velocity[edge]=0.5*(node_velocity[connection.lower_node]+node_velocity[connection.upper_node]);}
		for(const auto& source:topology.mass_flux_sources){if(source.momentum_connection<0||source.momentum_connection>=static_cast<int>(topology.connections.size())||source.coarse_fine_connection<0||source.coarse_fine_connection>=static_cast<int>(coarse_fine_velocity.size())||!(source.overlap_area>0))throw std::invalid_argument("composite momentum mass-flux source outside topology");velocity[source.momentum_connection]+=source.overlap_area*coarse_fine_velocity[source.coarse_fine_connection];source_area[source.momentum_connection]+=source.overlap_area;}
		for(int edge=0;edge<static_cast<int>(topology.connections.size());++edge){const auto& connection=topology.connections[edge];if(connection.component==connection.transport_axis){if(source_area[edge]>0)throw std::runtime_error("normal momentum connection unexpectedly has pressure-tile sources");}else{if(std::abs(source_area[edge]-connection.open_area)>1e-10*connection.open_area)throw std::runtime_error("tangential momentum source areas do not tile connection");velocity[edge]/=connection.open_area;}if(!std::isfinite(velocity[edge]))throw std::runtime_error("non-finite composite momentum connection velocity");}return velocity;
	}

	DeviceAmrMacFaceMap::DeviceAmrMacFaceMap(DeviceAmrFields& fields,
		const std::vector<AmrMacFaceAddress>& addresses)
	{
		size_=static_cast<int>(addresses.size());level_count_=fields.level_count();std::vector<DeviceAmrFieldLevelView> levels(level_count_);for(int level=0;level<level_count_;++level)levels[level]=fields.level_view(level);levels_=upload(levels,"upload compact MAC field views");std::vector<int> level(size_);std::vector<std::uint64_t> index(size_);std::vector<std::int8_t> component(size_);
		struct Key{int level,brick,component,i,j,k;bool operator==(const Key& other)const{return level==other.level&&brick==other.brick&&component==other.component&&i==other.i&&j==other.j&&k==other.k;}};struct Hash{std::size_t operator()(const Key& key)const{std::size_t h=0;for(int value:{key.level,key.brick,key.component,key.i,key.j,key.k})h=(h*1315423911u)^static_cast<std::uint32_t>(value);return h;}};std::unordered_map<Key,int,Hash> unique;
		for(int q=0;q<size_;++q){const auto& address=addresses[q];if(address.level<0||address.level>=level_count_||address.brick<0||address.brick>=levels[address.level].brick_count||address.component<0||address.component>2)throw std::invalid_argument("compact MAC map address outside field hierarchy");const int bs=levels[address.level].layout.brick_size,nx=address.component==0?bs+1:bs,ny=address.component==1?bs+1:bs,nz=address.component==2?bs+1:bs;if(address.i<0||address.j<0||address.k<0||address.i>=nx||address.j>=ny||address.k>=nz)throw std::invalid_argument("compact MAC map local face index outside brick");const Key key{address.level,address.brick,address.component,address.i,address.j,address.k};if(!unique.emplace(key,q).second)throw std::invalid_argument("compact MAC map addresses must be unique");level[q]=address.level;component[q]=static_cast<std::int8_t>(address.component);const BrickFieldLayout& layout=levels[address.level].layout;index[q]=address.component==0?layout.u_index(address.brick,address.i,address.j,address.k):(address.component==1?layout.v_index(address.brick,address.i,address.j,address.k):layout.w_index(address.brick,address.i,address.j,address.k));}
		level_=upload(level,"upload compact MAC map levels");index_=upload(index,"upload compact MAC map indices");component_=upload(component,"upload compact MAC map components");bytes_=level_count_*sizeof(DeviceAmrFieldLevelView)+size_*(sizeof(int)+sizeof(std::uint64_t)+sizeof(std::int8_t));
	}

	DeviceAmrMacFaceMap::~DeviceAmrMacFaceMap(){for(void* pointer:{(void*)levels_,(void*)level_,(void*)index_,(void*)component_})if(pointer)cudaFree(pointer);}

	void DeviceAmrMacFaceMap::gather(Real* compact)const{if(size_&&!compact)throw std::invalid_argument("compact MAC gather destination is null");if(size_)gather_mac_face_map_kernel<<<(size_+255)/256,256>>>(levels_,level_,index_,component_,compact,size_);check(cudaDeviceSynchronize(),"gather compact MAC face map");}
	void DeviceAmrMacFaceMap::scatter(const Real* compact)const{if(size_&&!compact)throw std::invalid_argument("compact MAC scatter source is null");if(size_)scatter_mac_face_map_kernel<<<(size_+255)/256,256>>>(levels_,level_,index_,component_,compact,size_);check(cudaDeviceSynchronize(),"scatter compact MAC face map");}

	DevicePairwiseMomentumTransport::DevicePairwiseMomentumTransport(const std::vector<double>& dual_volume,
		const std::vector<PairwiseMomentumConnection>& connections)
	{
		node_count_=static_cast<int>(dual_volume.size());connection_count_=static_cast<int>(connections.size());if(!node_count_)throw std::invalid_argument("pairwise GPU transport requires nodes");
		std::vector<Real> volume(node_count_),area(connection_count_);std::vector<int> a(connection_count_),b(connection_count_);
		for(int node=0;node<node_count_;++node){if(!(dual_volume[node]>0)||!std::isfinite(dual_volume[node]))throw std::invalid_argument("pairwise GPU dual volume must be finite and positive");volume[node]=static_cast<Real>(dual_volume[node]);}
		for(int edge=0;edge<connection_count_;++edge){const auto& source=connections[edge];if(source.a<0||source.b<0||source.a==source.b||source.a>=node_count_||source.b>=node_count_||!(source.open_area>=0)||!std::isfinite(source.open_area))throw std::invalid_argument("invalid pairwise GPU connection");a[edge]=source.a;b[edge]=source.b;area[edge]=static_cast<Real>(source.open_area);}
		volume_=upload(volume,"upload pairwise momentum volumes");a_=upload(a,"upload pairwise momentum first nodes");b_=upload(b,"upload pairwise momentum second nodes");area_=upload(area,"upload pairwise momentum areas");delta_x_=allocate<Real>(node_count_,"allocate pairwise momentum x scratch");delta_y_=allocate<Real>(node_count_,"allocate pairwise momentum y scratch");delta_z_=allocate<Real>(node_count_,"allocate pairwise momentum z scratch");bytes_=node_count_*4*sizeof(Real)+connection_count_*(2*sizeof(int)+sizeof(Real));
	}

	DevicePairwiseMomentumTransport::~DevicePairwiseMomentumTransport(){for(void* pointer:{(void*)a_,(void*)b_,(void*)volume_,(void*)area_,(void*)delta_x_,(void*)delta_y_,(void*)delta_z_})if(pointer)cudaFree(pointer);}

	void DevicePairwiseMomentumTransport::step(Real* velocity_x,Real* velocity_y,Real* velocity_z,
		const Real* connection_normal_velocity,Real dt)
	{
		if(!velocity_x||!velocity_y||!velocity_z||(!connection_normal_velocity&&connection_count_))throw std::invalid_argument("pairwise GPU transport has null state");if(!(dt>Real(0)))throw std::invalid_argument("pairwise GPU transport timestep must be positive");
		check(cudaMemset(delta_x_,0,node_count_*sizeof(Real)),"clear pairwise momentum x scratch");check(cudaMemset(delta_y_,0,node_count_*sizeof(Real)),"clear pairwise momentum y scratch");check(cudaMemset(delta_z_,0,node_count_*sizeof(Real)),"clear pairwise momentum z scratch");
		if(connection_count_)accumulate_pairwise_momentum_kernel<<<(connection_count_+255)/256,256>>>(a_,b_,area_,connection_normal_velocity,velocity_x,velocity_y,velocity_z,dt,delta_x_,delta_y_,delta_z_,connection_count_);
		apply_pairwise_momentum_kernel<<<(node_count_+255)/256,256>>>(volume_,delta_x_,delta_y_,delta_z_,velocity_x,velocity_y,velocity_z,node_count_);check(cudaDeviceSynchronize(),"pairwise conservative momentum transport");
	}

	void DevicePairwiseMomentumTransport::step_scalar(Real* velocity,const Real* connection_normal_velocity,Real dt)
	{
		if(!velocity||(!connection_normal_velocity&&connection_count_))throw std::invalid_argument("pairwise GPU scalar transport has null state");if(!(dt>Real(0)))throw std::invalid_argument("pairwise GPU scalar transport timestep must be positive");check(cudaMemset(delta_x_,0,node_count_*sizeof(Real)),"clear pairwise scalar scratch");if(connection_count_)accumulate_pairwise_scalar_kernel<<<(connection_count_+255)/256,256>>>(a_,b_,area_,connection_normal_velocity,velocity,dt,delta_x_,connection_count_);apply_pairwise_scalar_kernel<<<(node_count_+255)/256,256>>>(volume_,delta_x_,velocity,node_count_);check(cudaDeviceSynchronize(),"pairwise conservative scalar transport");
	}

	DeviceCompositeEbMomentumTransport::DeviceCompositeEbMomentumTransport(
		const CompositeAmrPressureSystem& system,DeviceAmrFields& fields)
	{
		if(!system.hierarchy||system.brick_size<=0)throw std::invalid_argument("compact/perimeter momentum transport requires a composite hierarchy");const AmrHierarchy& hierarchy=*system.hierarchy;if(fields.level_count()!=static_cast<int>(hierarchy.levels().size()))throw std::invalid_argument("compact/perimeter momentum field hierarchy mismatch");embedded_count_=static_cast<int>(system.embedded.size());if(!embedded_count_)throw std::invalid_argument("compact/perimeter momentum transport requires embedded apertures");
		std::unordered_map<int,int> node_index;std::vector<int> node_dof;auto add_node=[&](int dof){if(dof<0||dof>=system.storage_size||!system.active[dof])throw std::invalid_argument("compact/perimeter momentum node is inactive");const auto found=node_index.find(dof);if(found!=node_index.end())return found->second;const int id=static_cast<int>(node_dof.size());node_index.emplace(dof,id);node_dof.push_back(dof);return id;};
		std::vector<int> embedded_a(embedded_count_),embedded_b(embedded_count_);std::vector<std::int8_t> embedded_axis(embedded_count_);std::vector<Real> embedded_area(embedded_count_),embedded_mass(embedded_count_);embedded_mass_host_.resize(embedded_count_);embedded_axis_host_.resize(embedded_count_);for(int edge=0;edge<embedded_count_;++edge){const CoarseFinePressureConnection& connection=system.embedded[edge];const int first=add_node(connection.coarse_dof),second=add_node(connection.fine_dof);embedded_a[edge]=connection.direction>0?first:second;embedded_b[edge]=connection.direction>0?second:first;embedded_axis[edge]=connection.axis;embedded_area[edge]=static_cast<Real>(connection.open_area);const double mass=connection.open_area/pressure_gradient_factor(connection);if(!(mass>0)||!std::isfinite(mass))throw std::runtime_error("compact aperture momentum mass is invalid");embedded_mass[edge]=static_cast<Real>(mass);embedded_mass_host_[edge]=mass;embedded_axis_host_[edge]=connection.axis;}
		const std::vector<CompositeEbMomentumRegularConnection> regular=build_composite_eb_momentum_regular_connections(system);regular_count_=static_cast<int>(regular.size());std::vector<int> regular_a(regular_count_),regular_b(regular_count_);std::vector<Real> regular_area(regular_count_);for(int edge=0;edge<regular_count_;++edge){regular_a[edge]=add_node(regular[edge].lower_dof);regular_b[edge]=add_node(regular[edge].upper_dof);regular_area[edge]=static_cast<Real>(regular[edge].open_area);}node_count_=static_cast<int>(node_dof.size());
		struct CellAddress{int level=-1,brick=-1,i=-1,j=-1,k=-1;};const int bs=system.brick_size,cells=bs*bs*bs;auto cell_address=[&](int dof){CellAddress out;for(int level=0;level<static_cast<int>(hierarchy.levels().size());++level){const int begin=system.level_offset[level],count=static_cast<int>(hierarchy.levels()[level].bricks.size())*cells;if(dof<begin||dof>=begin+count)continue;const int work=dof-begin,local=work%cells;out.level=level;out.brick=work/cells;out.i=local%bs;out.j=(local/bs)%bs;out.k=local/(bs*bs);break;}return out;};
		struct FaceKey{int level,brick,component,i,j,k;bool operator==(const FaceKey& other)const{return level==other.level&&brick==other.brick&&component==other.component&&i==other.i&&j==other.j&&k==other.k;}};struct FaceHash{std::size_t operator()(const FaceKey& key)const{std::size_t h=0;for(int value:{key.level,key.brick,key.component,key.i,key.j,key.k})h=(h*1315423911u)^static_cast<std::uint32_t>(value);return h;}};std::unordered_map<FaceKey,int,FaceHash> carrier_index;std::vector<AmrMacFaceAddress> carrier_address;std::vector<int> incidence_node,incidence_carrier;std::vector<std::int8_t> incidence_component;std::vector<Real> incidence_half_mass;auto add_carrier=[&](AmrMacFaceAddress face,double mass){face=canonical_amr_mac_face_address(hierarchy,face);const FaceKey key{face.level,face.brick,face.component,face.i,face.j,face.k};const auto found=carrier_index.find(key);if(found!=carrier_index.end()){const int id=found->second;if(std::abs(carrier_mass_host_[id]-mass)>1e-10*std::max(carrier_mass_host_[id],mass))throw std::runtime_error("compact/perimeter carrier has inconsistent physical mass");return id;}const int id=static_cast<int>(carrier_address.size());carrier_index.emplace(key,id);carrier_address.push_back(face);carrier_mass_host_.push_back(mass);carrier_component_host_.push_back(static_cast<std::int8_t>(face.component));return id;};
		for(int node=0;node<node_count_;++node){const int dof=node_dof[node];const CellAddress cell=cell_address(dof);if(cell.level<0)continue;const AmrLevel& level=hierarchy.levels()[cell.level];const BrickMetadata& metadata=level.bricks[cell.brick];if(!metadata.active())continue;for(int axis=0;axis<3;++axis)for(int sign=-1;sign<=1;sign+=2){int coordinate[3]={cell.i,cell.j,cell.k},other_brick=cell.brick;coordinate[axis]+=sign;if(coordinate[axis]<0||coordinate[axis]>=bs){other_brick=metadata.same_level_neighbor[2*axis+(sign>0)];if(other_brick<0||!level.bricks[other_brick].active())continue;coordinate[axis]=sign>0?0:bs-1;}const int other=system.dof(cell.level,other_brick,coordinate[0],coordinate[1],coordinate[2]);if(other<0||other>=system.storage_size||!system.active[other])continue;const int lower=sign>0?dof:other;if(system.cut_face_mask[lower]&(1u<<axis))continue;AmrMacFaceAddress face{cell.level,cell.brick,axis,cell.i,cell.j,cell.k};if(axis==0)face.i+=sign>0;else if(axis==1)face.j+=sign>0;else face.k+=sign>0;const double mass=composite_mac_carrier_volume(system,dof,other);if(!(mass>0)||!std::isfinite(mass))throw std::runtime_error("compact/perimeter regular carrier mass is invalid");const int carrier=add_carrier(face,mass);incidence_node.push_back(node);incidence_carrier.push_back(carrier);incidence_component.push_back(static_cast<std::int8_t>(axis));incidence_half_mass.push_back(static_cast<Real>(0.5*mass));}}
		carrier_count_=static_cast<int>(carrier_address.size());incidence_count_=static_cast<int>(incidence_node.size());std::vector<int> regular_carrier(regular_count_);for(int edge=0;edge<regular_count_;++edge){const AmrMacFaceAddress face=canonical_amr_mac_face_address(hierarchy,regular[edge].face);const FaceKey key{face.level,face.brick,face.component,face.i,face.j,face.k};const auto found=carrier_index.find(key);if(found==carrier_index.end())throw std::runtime_error("compact/perimeter normal carrier is absent from node incidence");regular_carrier[edge]=found->second;}
		carrier_map_=std::make_unique<DeviceAmrMacFaceMap>(fields,carrier_address);std::vector<AmrMacFaceAddress> alias_address;std::vector<int> alias_source;for(int carrier=0;carrier<carrier_count_;++carrier){const AmrMacFaceAddress& face=carrier_address[carrier];const int normal_coordinate=face.component==0?face.i:(face.component==1?face.j:face.k);if(normal_coordinate!=bs)continue;const AmrLevel& level=hierarchy.levels()[face.level];const int neighbour=level.bricks[face.brick].same_level_neighbor[2*face.component+1];if(neighbour<0||!level.bricks[neighbour].active())continue;AmrMacFaceAddress alias=face;alias.brick=neighbour;if(alias.component==0)alias.i=0;else if(alias.component==1)alias.j=0;else alias.k=0;alias_address.push_back(alias);alias_source.push_back(carrier);}alias_count_=static_cast<int>(alias_address.size());if(alias_count_){alias_map_=std::make_unique<DeviceAmrMacFaceMap>(fields,alias_address);alias_source_=upload(alias_source,"upload compact/perimeter same-level alias sources");alias_value_=allocate<Real>(alias_count_,"allocate compact/perimeter same-level alias values");}
		carrier_state_=allocate<Real>(carrier_count_,"allocate compact/perimeter carrier state");carrier_delta_=allocate<Real>(carrier_count_,"allocate compact/perimeter carrier delta");incidence_node_=upload(incidence_node,"upload compact/perimeter incidence nodes");incidence_carrier_=upload(incidence_carrier,"upload compact/perimeter incidence carriers");incidence_component_=upload(incidence_component,"upload compact/perimeter incidence components");incidence_half_mass_=upload(incidence_half_mass,"upload compact/perimeter incidence masses");embedded_a_=upload(embedded_a,"upload compact momentum lower aperture nodes");embedded_b_=upload(embedded_b,"upload compact momentum upper aperture nodes");embedded_axis_=upload(embedded_axis,"upload compact momentum aperture axes");embedded_area_=upload(embedded_area,"upload compact momentum aperture areas");embedded_mass_=upload(embedded_mass,"upload compact momentum aperture masses");regular_a_=upload(regular_a,"upload compact momentum lower perimeter nodes");regular_b_=upload(regular_b,"upload compact momentum upper perimeter nodes");regular_carrier_=upload(regular_carrier,"upload compact momentum perimeter carriers");regular_area_=upload(regular_area,"upload compact momentum perimeter areas");node_sum_=allocate<Real>(static_cast<std::size_t>(node_count_)*3,"allocate compact/perimeter node momentum");node_weight_=allocate<Real>(static_cast<std::size_t>(node_count_)*3,"allocate compact/perimeter node mass");node_delta_=allocate<Real>(static_cast<std::size_t>(node_count_)*3,"allocate compact/perimeter node increments");bytes_=static_cast<std::size_t>(2*carrier_count_+9*node_count_+alias_count_)*sizeof(Real)+static_cast<std::size_t>(alias_count_)*sizeof(int)+static_cast<std::size_t>(incidence_count_)*(2*sizeof(int)+sizeof(std::int8_t)+sizeof(Real))+static_cast<std::size_t>(embedded_count_)*(2*sizeof(int)+sizeof(std::int8_t)+2*sizeof(Real))+static_cast<std::size_t>(regular_count_)*(3*sizeof(int)+sizeof(Real));
	}

	DeviceCompositeEbMomentumTransport::~DeviceCompositeEbMomentumTransport()
	{
		for(void* pointer:{(void*)carrier_state_,(void*)carrier_delta_,(void*)alias_source_,(void*)alias_value_,(void*)incidence_node_,(void*)incidence_carrier_,(void*)incidence_component_,(void*)incidence_half_mass_,(void*)embedded_a_,(void*)embedded_b_,(void*)embedded_axis_,(void*)embedded_area_,(void*)embedded_mass_,(void*)regular_a_,(void*)regular_b_,(void*)regular_carrier_,(void*)regular_area_,(void*)node_sum_,(void*)node_weight_,(void*)node_delta_})if(pointer)cudaFree(pointer);
	}

	void DeviceCompositeEbMomentumTransport::step(Real* embedded_velocity,Real dt)
	{
		if(!embedded_velocity)throw std::invalid_argument("compact/perimeter embedded velocity is null");if(!(dt>Real(0)))throw std::invalid_argument("compact/perimeter momentum timestep must be positive");carrier_map_->gather(carrier_state_);const std::size_t components=static_cast<std::size_t>(node_count_)*3;check(cudaMemset(node_sum_,0,components*sizeof(Real)),"clear compact/perimeter node momentum");check(cudaMemset(node_weight_,0,components*sizeof(Real)),"clear compact/perimeter node mass");check(cudaMemset(node_delta_,0,components*sizeof(Real)),"clear compact/perimeter node increments");check(cudaMemset(carrier_delta_,0,static_cast<std::size_t>(carrier_count_)*sizeof(Real)),"clear compact/perimeter carrier increments");if(incidence_count_)accumulate_eb_momentum_carrier_state_kernel<<<(incidence_count_+255)/256,256>>>(incidence_node_,incidence_carrier_,incidence_component_,incidence_half_mass_,carrier_state_,node_sum_,node_weight_,incidence_count_);accumulate_eb_momentum_aperture_state_kernel<<<(embedded_count_+255)/256,256>>>(embedded_a_,embedded_b_,embedded_axis_,embedded_mass_,embedded_velocity,node_sum_,node_weight_,embedded_count_);transport_eb_momentum_aperture_kernel<<<(embedded_count_+255)/256,256>>>(embedded_a_,embedded_b_,embedded_area_,embedded_velocity,node_sum_,node_weight_,dt,node_delta_,embedded_count_);if(regular_count_)transport_eb_momentum_regular_kernel<<<(regular_count_+255)/256,256>>>(regular_a_,regular_b_,regular_carrier_,regular_area_,carrier_state_,node_sum_,node_weight_,dt,node_delta_,regular_count_);if(incidence_count_)scatter_eb_momentum_carrier_delta_kernel<<<(incidence_count_+255)/256,256>>>(incidence_node_,incidence_carrier_,incidence_component_,node_weight_,node_delta_,carrier_delta_,incidence_count_);if(carrier_count_)apply_eb_momentum_carrier_delta_kernel<<<(carrier_count_+255)/256,256>>>(carrier_state_,carrier_delta_,carrier_count_);apply_eb_momentum_aperture_delta_kernel<<<(embedded_count_+255)/256,256>>>(embedded_a_,embedded_b_,embedded_axis_,node_weight_,node_delta_,embedded_velocity,embedded_count_);carrier_map_->scatter(carrier_state_);if(alias_count_){gather_eb_momentum_alias_kernel<<<(alias_count_+255)/256,256>>>(alias_source_,carrier_state_,alias_value_,alias_count_);alias_map_->scatter(alias_value_);}
	}

	std::array<double,3> DeviceCompositeEbMomentumTransport::momentum(const Real* embedded_velocity)const
	{
		if(!embedded_velocity)throw std::invalid_argument("compact/perimeter momentum diagnostic velocity is null");carrier_map_->gather(carrier_state_);std::vector<Real> carriers(carrier_count_),apertures(embedded_count_);if(carrier_count_)check(cudaMemcpy(carriers.data(),carrier_state_,static_cast<std::size_t>(carrier_count_)*sizeof(Real),cudaMemcpyDeviceToHost),"download compact/perimeter carriers");if(embedded_count_)check(cudaMemcpy(apertures.data(),embedded_velocity,static_cast<std::size_t>(embedded_count_)*sizeof(Real),cudaMemcpyDeviceToHost),"download compact/perimeter apertures");std::array<double,3> result{};for(int q=0;q<carrier_count_;++q)result[carrier_component_host_[q]]+=carrier_mass_host_[q]*static_cast<double>(carriers[q]);for(int q=0;q<embedded_count_;++q)result[embedded_axis_host_[q]]+=embedded_mass_host_[q]*static_cast<double>(apertures[q]);return result;
	}

	std::size_t DeviceCompositeEbMomentumTransport::bytes()const{return bytes_+(carrier_map_?carrier_map_->bytes():0)+(alias_map_?alias_map_->bytes():0);}

	DeviceCompositeAmrMomentumInterfaceTransport::DeviceCompositeAmrMomentumInterfaceTransport(
		const CompositeAmrPressureSystem& system,DeviceAmrFields& fields)
	{
		const CompositeAmrMomentumInterfaceTopology topology=build_composite_amr_momentum_interface_topology(system);node_count_=static_cast<int>(topology.nodes.size());connection_count_=static_cast<int>(topology.connections.size());alias_count_=static_cast<int>(topology.coarse_aliases.size());mass_flux_source_count_=static_cast<int>(topology.mass_flux_sources.size());if(!node_count_||!connection_count_)throw std::invalid_argument("composite AMR momentum transport requires a non-empty 2:1 interface");std::vector<PairwiseMomentumConnection> connections;connections.reserve(connection_count_);std::vector<int> connection_lower(connection_count_),connection_upper(connection_count_);std::vector<std::int8_t> connection_axis(connection_count_),connection_component(connection_count_);std::vector<Real> connection_area(connection_count_);for(int edge_index=0;edge_index<connection_count_;++edge_index){const auto& edge=topology.connections[edge_index];connections.push_back({edge.lower_node,edge.upper_node,edge.open_area,0});connection_lower[edge_index]=edge.lower_node;connection_upper[edge_index]=edge.upper_node;connection_axis[edge_index]=edge.transport_axis;connection_component[edge_index]=edge.component;connection_area[edge_index]=static_cast<Real>(edge.open_area);}std::vector<int> mass_flux_connection(mass_flux_source_count_),mass_flux_source(mass_flux_source_count_);std::vector<Real> mass_flux_area(mass_flux_source_count_);for(int q=0;q<mass_flux_source_count_;++q){const auto& source=topology.mass_flux_sources[q];mass_flux_connection[q]=source.momentum_connection;mass_flux_source[q]=source.coarse_fine_connection;mass_flux_area[q]=static_cast<Real>(source.overlap_area);}node_map_=std::make_unique<DeviceAmrMacFaceMap>(fields,topology.nodes);transport_=std::make_unique<DevicePairwiseMomentumTransport>(topology.dual_volume,connections);node_state_=allocate<Real>(node_count_,"allocate composite momentum node state");connection_normal_velocity_=allocate<Real>(connection_count_,"allocate composite momentum connection velocities");connection_lower_=upload(connection_lower,"upload composite momentum lower nodes");connection_upper_=upload(connection_upper,"upload composite momentum upper nodes");connection_axis_=upload(connection_axis,"upload composite momentum transport axes");connection_component_=upload(connection_component,"upload composite momentum components");connection_area_=upload(connection_area,"upload composite momentum connection areas");mass_flux_connection_=upload(mass_flux_connection,"upload composite momentum flux destinations");mass_flux_source_=upload(mass_flux_source,"upload composite momentum pressure-flux sources");mass_flux_area_=upload(mass_flux_area,"upload composite momentum pressure-flux areas");
		struct Key{int level,brick,component,i,j,k;bool operator==(const Key& other)const{return level==other.level&&brick==other.brick&&component==other.component&&i==other.i&&j==other.j&&k==other.k;}};struct Hash{std::size_t operator()(const Key& key)const{std::size_t h=0;for(int value:{key.level,key.brick,key.component,key.i,key.j,key.k})h=(h*1315423911u)^static_cast<std::uint32_t>(value);return h;}};std::unordered_map<Key,int,Hash> group_index;std::vector<AmrMacFaceAddress> group_address;std::vector<double> group_area;std::vector<int> alias_node(alias_count_),alias_group(alias_count_);std::vector<Real> alias_area(alias_count_);
		for(int q=0;q<alias_count_;++q){const CompositeAmrMomentumCoarseAlias& alias=topology.coarse_aliases[q];const Key key{alias.coarse_face.level,alias.coarse_face.brick,alias.coarse_face.component,alias.coarse_face.i,alias.coarse_face.j,alias.coarse_face.k};const auto found=group_index.find(key);int group;if(found==group_index.end()){group=static_cast<int>(group_address.size());group_index.emplace(key,group);group_address.push_back(alias.coarse_face);group_area.push_back(0);}else group=found->second;alias_node[q]=alias.fine_owned_node;alias_group[q]=group;alias_area[q]=static_cast<Real>(alias.area);group_area[group]+=alias.area;}alias_group_count_=static_cast<int>(group_address.size());if(!alias_group_count_)throw std::runtime_error("composite AMR momentum topology has no coarse alias groups");std::vector<Real> group_area_real(alias_group_count_);for(int q=0;q<alias_group_count_;++q){if(!(group_area[q]>0))throw std::runtime_error("composite momentum coarse alias has zero area");group_area_real[q]=static_cast<Real>(group_area[q]);}alias_map_=std::make_unique<DeviceAmrMacFaceMap>(fields,group_address);alias_node_=upload(alias_node,"upload momentum alias nodes");alias_group_=upload(alias_group,"upload momentum alias groups");alias_area_=upload(alias_area,"upload momentum alias areas");alias_group_area_=upload(group_area_real,"upload momentum alias group areas");alias_sum_=allocate<Real>(alias_group_count_,"allocate momentum alias sum");alias_value_=allocate<Real>(alias_group_count_,"allocate momentum alias values");bytes_=static_cast<std::size_t>(node_count_+connection_count_)*sizeof(Real)+static_cast<std::size_t>(connection_count_)*(2*sizeof(int)+2*sizeof(std::int8_t)+sizeof(Real))+static_cast<std::size_t>(mass_flux_source_count_)*(2*sizeof(int)+sizeof(Real))+static_cast<std::size_t>(alias_count_)*(2*sizeof(int)+sizeof(Real))+static_cast<std::size_t>(alias_group_count_)*3*sizeof(Real);
	}

	DeviceCompositeAmrMomentumInterfaceTransport::~DeviceCompositeAmrMomentumInterfaceTransport()
	{
		for(void* pointer:{(void*)node_state_,(void*)connection_normal_velocity_,(void*)alias_sum_,(void*)alias_value_,(void*)alias_node_,(void*)alias_group_,(void*)alias_area_,(void*)alias_group_area_,(void*)connection_lower_,(void*)connection_upper_,(void*)connection_axis_,(void*)connection_component_,(void*)connection_area_,(void*)mass_flux_connection_,(void*)mass_flux_source_,(void*)mass_flux_area_})if(pointer)cudaFree(pointer);
	}

	void DeviceCompositeAmrMomentumInterfaceTransport::step(const Real* connection_normal_velocity,Real dt)
	{
		if(!connection_normal_velocity)throw std::invalid_argument("composite momentum connection velocity is null");if(!(dt>Real(0)))throw std::invalid_argument("composite momentum timestep must be positive");node_map_->gather(node_state_);transport_->step_scalar(node_state_,connection_normal_velocity,dt);node_map_->scatter(node_state_);check(cudaMemset(alias_sum_,0,static_cast<std::size_t>(alias_group_count_)*sizeof(Real)),"clear momentum coarse alias sums");accumulate_momentum_alias_kernel<<<(alias_count_+255)/256,256>>>(alias_node_,alias_group_,alias_area_,node_state_,alias_sum_,alias_count_);normalize_momentum_alias_kernel<<<(alias_group_count_+255)/256,256>>>(alias_sum_,alias_group_area_,alias_value_,alias_group_count_);alias_map_->scatter(alias_value_);
	}

	void DeviceCompositeAmrMomentumInterfaceTransport::step_from_composite_flux(const Real* coarse_fine_velocity,Real dt)
	{
		if(!coarse_fine_velocity)throw std::invalid_argument("composite pressure-interface velocity is null");if(!(dt>Real(0)))throw std::invalid_argument("composite momentum timestep must be positive");node_map_->gather(node_state_);initialize_momentum_connection_velocity_kernel<<<(connection_count_+255)/256,256>>>(connection_lower_,connection_upper_,connection_axis_,connection_component_,node_state_,connection_normal_velocity_,connection_count_);if(mass_flux_source_count_)accumulate_momentum_mass_flux_kernel<<<(mass_flux_source_count_+255)/256,256>>>(mass_flux_connection_,mass_flux_source_,mass_flux_area_,coarse_fine_velocity,connection_normal_velocity_,mass_flux_source_count_);normalize_tangential_momentum_velocity_kernel<<<(connection_count_+255)/256,256>>>(connection_axis_,connection_component_,connection_area_,connection_normal_velocity_,connection_count_);transport_->step_scalar(node_state_,connection_normal_velocity_,dt);node_map_->scatter(node_state_);check(cudaMemset(alias_sum_,0,static_cast<std::size_t>(alias_group_count_)*sizeof(Real)),"clear momentum coarse alias sums");accumulate_momentum_alias_kernel<<<(alias_count_+255)/256,256>>>(alias_node_,alias_group_,alias_area_,node_state_,alias_sum_,alias_count_);normalize_momentum_alias_kernel<<<(alias_group_count_+255)/256,256>>>(alias_sum_,alias_group_area_,alias_value_,alias_group_count_);alias_map_->scatter(alias_value_);
	}

	std::size_t DeviceCompositeAmrMomentumInterfaceTransport::bytes()const
	{
		return bytes_+(node_map_?node_map_->bytes():0)+(alias_map_?alias_map_->bytes():0)+(transport_?transport_->bytes():0);
	}

	DeviceAmrAdvection::DeviceAmrAdvection(const AmrHierarchy& hierarchy,const TriangleBvh& fabric,double protection_cells,const CompositeAmrPressureSystem* composite_system):locator_(hierarchy)
	{
		if(!(protection_cells>=1.0))throw std::invalid_argument("AMR fabric advection protection must cover at least one cell");levels_.resize(hierarchy.levels().size());std::vector<DeviceAmrFieldLevelView> empty_views(levels_.size());const std::vector<AmrMacFaceDirectionalLink> replaced=composite_system?build_momentum_interface_replaced_links(*composite_system):std::vector<AmrMacFaceDirectionalLink>{};replaced_interface_links_=replaced.size();
		for(std::size_t level=0;level<levels_.size();++level){const AmrLevel& metadata=hierarchy.levels()[level];Level& allocation=levels_[level];allocation.layout=BrickFieldLayout::make(hierarchy.brick_size(),hierarchy.ghost_cells());allocation.brick_count=static_cast<int>(metadata.bricks.size());if(!allocation.brick_count)continue;const std::size_t un=static_cast<std::size_t>(allocation.brick_count)*allocation.layout.u_stride,vn=static_cast<std::size_t>(allocation.brick_count)*allocation.layout.v_stride,wn=static_cast<std::size_t>(allocation.brick_count)*allocation.layout.w_stride;const int bs=hierarchy.brick_size();const double radius=protection_cells*metadata.h;std::vector<unsigned char> lu(un,0x3f),lv(vn,0x3f),lw(wn,0x3f);
			for(int brick=0;brick<allocation.brick_count;++brick){const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;Aabb3d expanded{{record.origin.x-radius,record.origin.y-radius,record.origin.z-radius},{record.origin.x+bs*metadata.h+radius,record.origin.y+bs*metadata.h+radius,record.origin.z+bs*metadata.h+radius}};if(fabric.query_aabb(expanded).empty()){active_faces_+=static_cast<std::size_t>(3)*bs*bs*(bs+1);continue;}for(int component=0;component<3;++component){const int ni=component==0?bs+1:bs,nj=component==1?bs+1:bs,nk=component==2?bs+1:bs;for(int k=0;k<nk;++k)for(int j=0;j<nj;++j)for(int i=0;i<ni;++i){const Vec3d point=face_centre(record,metadata.h,component,i,j,k);const double queried=fabric.distance(point,radius);const std::size_t index=component==0?allocation.layout.u_index(brick,i,j,k):(component==1?allocation.layout.v_index(brick,i,j,k):allocation.layout.w_index(brick,i,j,k));if(std::isfinite(queried)){unsigned char links=0x40;if(queried>1e-7*metadata.h)for(int axis=0;axis<3;++axis)for(int direction=-1;direction<=1;direction+=2){Vec3d other=point;other[axis]+=direction*metadata.h;if(!fabric.intersect_segment(point,other,1e-7,1.0-1e-7).hit)links|=static_cast<unsigned char>(1u<<(2*axis+(direction>0)));}if(component==0)lu[index]=links;else if(component==1)lv[index]=links;else lw[index]=links;++protected_faces_;}++active_faces_;}}}for(const auto& link:replaced)if(link.face.level==static_cast<int>(level)){const std::size_t index=link.face.component==0?allocation.layout.u_index(link.face.brick,link.face.i,link.face.j,link.face.k):(link.face.component==1?allocation.layout.v_index(link.face.brick,link.face.i,link.face.j,link.face.k):allocation.layout.w_index(link.face.brick,link.face.i,link.face.j,link.face.k));unsigned char& mask=link.face.component==0?lu[index]:(link.face.component==1?lv[index]:lw[index]);mask&=static_cast<unsigned char>(~(1u<<(2*link.transport_axis+(link.direction>0))));}
			allocation.u=allocate<Real>(un,"allocate AMR advection u scratch");allocation.v=allocate<Real>(vn,"allocate AMR advection v scratch");allocation.w=allocate<Real>(wn,"allocate AMR advection w scratch");allocation.forward_u=allocate<Real>(un,"allocate AMR forward u scratch");allocation.forward_v=allocate<Real>(vn,"allocate AMR forward v scratch");allocation.forward_w=allocate<Real>(wn,"allocate AMR forward w scratch");allocation.links_u=allocate<unsigned char>(un,"allocate AMR advection u links");allocation.links_v=allocate<unsigned char>(vn,"allocate AMR advection v links");allocation.links_w=allocate<unsigned char>(wn,"allocate AMR advection w links");check(cudaMemset(allocation.u,0,un*sizeof(Real)),"clear AMR advection u scratch");check(cudaMemset(allocation.v,0,vn*sizeof(Real)),"clear AMR advection v scratch");check(cudaMemset(allocation.w,0,wn*sizeof(Real)),"clear AMR advection w scratch");check(cudaMemset(allocation.forward_u,0,un*sizeof(Real)),"clear AMR forward u scratch");check(cudaMemset(allocation.forward_v,0,vn*sizeof(Real)),"clear AMR forward v scratch");check(cudaMemset(allocation.forward_w,0,wn*sizeof(Real)),"clear AMR forward w scratch");check(cudaMemcpy(allocation.links_u,lu.data(),un,cudaMemcpyHostToDevice),"upload AMR u same-side links");check(cudaMemcpy(allocation.links_v,lv.data(),vn,cudaMemcpyHostToDevice),"upload AMR v same-side links");check(cudaMemcpy(allocation.links_w,lw.data(),wn,cudaMemcpyHostToDevice),"upload AMR w same-side links");bytes_+=(un+vn+wn)*(2*sizeof(Real)+sizeof(unsigned char));}
		device_views_=allocate<DeviceAmrFieldLevelView>(levels_.size(),"allocate AMR advection field views");device_forward_views_=allocate<DeviceAmrFieldLevelView>(levels_.size(),"allocate AMR forward field views");bytes_+=2*levels_.size()*sizeof(DeviceAmrFieldLevelView);
	}

	DeviceAmrAdvection::~DeviceAmrAdvection(){for(Level& level:levels_)for(void* pointer:{(void*)level.u,(void*)level.v,(void*)level.w,(void*)level.forward_u,(void*)level.forward_v,(void*)level.forward_w,(void*)level.links_u,(void*)level.links_v,(void*)level.links_w})if(pointer)cudaFree(pointer);if(device_views_)cudaFree(device_views_);if(device_forward_views_)cudaFree(device_forward_views_);}

	void DeviceAmrAdvection::advect(DeviceAmrFields& fields,Real dt)
	{
		if(!(dt>Real(0)))throw std::invalid_argument("AMR advection timestep must be positive");if(fields.level_count()!=static_cast<int>(levels_.size()))throw std::invalid_argument("AMR advection field hierarchy mismatch");std::vector<DeviceAmrFieldLevelView> views(levels_.size()),forward_views(levels_.size());for(int level=0;level<fields.level_count();++level){views[level]=fields.level_view(level);forward_views[level]=views[level];forward_views[level].u=levels_[level].forward_u;forward_views[level].v=levels_[level].forward_v;forward_views[level].w=levels_[level].forward_w;}check(cudaMemcpy(device_views_,views.data(),views.size()*sizeof(DeviceAmrFieldLevelView),cudaMemcpyHostToDevice),"upload AMR advection field views");check(cudaMemcpy(device_forward_views_,forward_views.data(),forward_views.size()*sizeof(DeviceAmrFieldLevelView),cudaMemcpyHostToDevice),"upload AMR forward field views");const GpuAmrHierarchyView hierarchy=locator_.view();
		// Build the forward state on every level before the correction sweep. Both sweeps
		// therefore sample one consistent hierarchy rather than a mixture of time levels.
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const int work=output.brick_count*3*output.layout.brick_size*output.layout.brick_size*(output.layout.brick_size+1);advect_forward_level_kernel<<<(work+255)/256,256>>>(hierarchy,device_views_,level,dt,output.links_u,output.links_v,output.links_w,output.forward_u,output.forward_v,output.forward_w);}
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const int work=output.brick_count*3*output.layout.brick_size*output.layout.brick_size*(output.layout.brick_size+1);advect_correct_level_kernel<<<(work+255)/256,256>>>(hierarchy,device_views_,device_forward_views_,level,dt,output.links_u,output.links_v,output.links_w,output.u,output.v,output.w);}
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const DeviceAmrFieldLevelView view=views[level];check(cudaMemcpy(view.u,output.u,static_cast<std::size_t>(output.brick_count)*output.layout.u_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR advected u");check(cudaMemcpy(view.v,output.v,static_cast<std::size_t>(output.brick_count)*output.layout.v_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR advected v");check(cudaMemcpy(view.w,output.w,static_cast<std::size_t>(output.brick_count)*output.layout.w_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR advected w");}check(cudaDeviceSynchronize(),"AMR advection");
	}

	void DeviceAmrAdvection::advect_conservative_uniform(DeviceAmrFields& fields,Real dt)
	{
		if(!(dt>Real(0)))throw std::invalid_argument("conservative uniform advection timestep must be positive");if(fields.level_count()!=1||levels_.size()!=1)throw std::invalid_argument("conservative uniform advection currently requires exactly one AMR level");const DeviceAmrFieldLevelView view=fields.level_view(0);check(cudaMemcpy(device_views_,&view,sizeof(view),cudaMemcpyHostToDevice),"upload conservative uniform field view");Level& output=levels_[0];const int work=output.brick_count*3*output.layout.brick_size*output.layout.brick_size*(output.layout.brick_size+1);if(work)advect_conservative_uniform_kernel<<<(work+255)/256,256>>>(locator_.view(),device_views_,dt,output.links_u,output.links_v,output.links_w,output.u,output.v,output.w,false,Real(0));check(cudaMemcpy(view.u,output.u,static_cast<std::size_t>(output.brick_count)*output.layout.u_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit conservative uniform u");check(cudaMemcpy(view.v,output.v,static_cast<std::size_t>(output.brick_count)*output.layout.v_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit conservative uniform v");check(cudaMemcpy(view.w,output.w,static_cast<std::size_t>(output.brick_count)*output.layout.w_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit conservative uniform w");check(cudaDeviceSynchronize(),"conservative uniform MAC advection");
	}

	void DeviceAmrAdvection::advect_conservative_uniform_external(DeviceAmrFields& fields,Real dt,Real freestream_speed)
	{
		if(!(dt>Real(0))||!(freestream_speed>=Real(0)))throw std::invalid_argument("invalid conservative external-aero timestep or freestream");if(fields.level_count()!=1||levels_.size()!=1)throw std::invalid_argument("conservative external-aero advection currently requires one AMR level");const DeviceAmrFieldLevelView view=fields.level_view(0);check(cudaMemcpy(device_views_,&view,sizeof(view),cudaMemcpyHostToDevice),"upload conservative external field view");Level& output=levels_[0];const int work=output.brick_count*3*output.layout.brick_size*output.layout.brick_size*(output.layout.brick_size+1);if(work)advect_conservative_uniform_kernel<<<(work+255)/256,256>>>(locator_.view(),device_views_,dt,output.links_u,output.links_v,output.links_w,output.u,output.v,output.w,true,freestream_speed);check(cudaMemcpy(view.u,output.u,static_cast<std::size_t>(output.brick_count)*output.layout.u_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit conservative external u");check(cudaMemcpy(view.v,output.v,static_cast<std::size_t>(output.brick_count)*output.layout.v_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit conservative external v");check(cudaMemcpy(view.w,output.w,static_cast<std::size_t>(output.brick_count)*output.layout.w_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit conservative external w");check(cudaDeviceSynchronize(),"conservative uniform external-aero advection");
	}

	void DeviceAmrAdvection::diffuse_smagorinsky(DeviceAmrFields& fields,Real molecular_nu,Real cs,Real dt)
	{
		if(molecular_nu<Real(0)||cs<Real(0)||!(dt>Real(0)))throw std::invalid_argument("invalid AMR LES/diffusion parameters");if(fields.level_count()!=static_cast<int>(levels_.size()))throw std::invalid_argument("AMR LES field hierarchy mismatch");std::vector<DeviceAmrFieldLevelView> views(levels_.size());for(int level=0;level<fields.level_count();++level)views[level]=fields.level_view(level);check(cudaMemcpy(device_views_,views.data(),views.size()*sizeof(DeviceAmrFieldLevelView),cudaMemcpyHostToDevice),"upload AMR LES field views");const GpuAmrHierarchyView hierarchy=locator_.view();
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const int cells=output.brick_count*output.layout.brick_size*output.layout.brick_size*output.layout.brick_size;smagorinsky_kernel<<<(cells+255)/256,256>>>(hierarchy,device_views_,level,cs);}
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const int work=output.brick_count*3*output.layout.brick_size*output.layout.brick_size*(output.layout.brick_size+1);diffuse_level_kernel<<<(work+255)/256,256>>>(hierarchy,device_views_,level,molecular_nu,dt,output.links_u,output.links_v,output.links_w,output.u,output.v,output.w);}
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const DeviceAmrFieldLevelView view=views[level];check(cudaMemcpy(view.u,output.u,static_cast<std::size_t>(output.brick_count)*output.layout.u_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR diffused u");check(cudaMemcpy(view.v,output.v,static_cast<std::size_t>(output.brick_count)*output.layout.v_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR diffused v");check(cudaMemcpy(view.w,output.w,static_cast<std::size_t>(output.brick_count)*output.layout.w_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR diffused w");}check(cudaDeviceSynchronize(),"AMR Smagorinsky diffusion");
	}
}
