#include "core/fluid/amr_advection.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace paracfd::core
{
	namespace
	{
		void check(cudaError_t error,const char* operation){if(error!=cudaSuccess)throw std::runtime_error(std::string(operation)+": "+cudaGetErrorString(error));}
		template<class T>T* allocate(std::size_t count,const char* operation){if(!count)return nullptr;T* pointer=nullptr;check(cudaMalloc(&pointer,count*sizeof(T)),operation);return pointer;}

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
		__device__ Real sample_local_clamped(const DeviceAmrFieldLevelView& view,const GpuBrickRecord& record,int component,int brick,GpuAmrPoint point)
		{
			Real x=(Real(point.x)-Real(record.origin_x))/Real(record.h)-(component==0?Real(0):Real(0.5));
			Real y=(Real(point.y)-Real(record.origin_y))/Real(record.h)-(component==1?Real(0):Real(0.5));
			Real z=(Real(point.z)-Real(record.origin_z))/Real(record.h)-(component==2?Real(0):Real(0.5));
			const int i0=static_cast<int>(floor(x)),j0=static_cast<int>(floor(y)),k0=static_cast<int>(floor(z));x-=i0;y-=j0;z-=k0;Real result=0;
			for(int dz=0;dz<2;++dz)for(int dy=0;dy<2;++dy)for(int dx=0;dx<2;++dx){const Real wx=dx?x:Real(1)-x,wy=dy?y:Real(1)-y,wz=dz?z:Real(1)-z;result+=wx*wy*wz*value_at_clamped(view,component,brick,i0+dx,j0+dy,k0+dz);}return result;
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
			GpuAmrPoint probe=node;const float epsilon=1e-5f*source_record.h;probe.x=fminf(hierarchy.domain_hi.x-epsilon,fmaxf(hierarchy.domain_lo.x+epsilon,probe.x));probe.y=fminf(hierarchy.domain_hi.y-epsilon,fmaxf(hierarchy.domain_lo.y+epsilon,probe.y));probe.z=fminf(hierarchy.domain_hi.z-epsilon,fmaxf(hierarchy.domain_lo.z+epsilon,probe.z));const GpuBrickLocation location=locate_finest(hierarchy,probe);if(location.brick<0)return value_at_clamped(source_view,component,source_brick,i,j,k);const GpuBrickRecord record=hierarchy.levels[location.level].bricks[location.brick];const DeviceAmrFieldLevelView view=views[location.level];if(location.level!=source_level)return sample_local_clamped(view,record,component,location.brick,node);const int ni=nearest_lattice_index(node.x,record.origin_x,record.h,component==0?0.0f:0.5f),nj=nearest_lattice_index(node.y,record.origin_y,record.h,component==1?0.0f:0.5f),nk=nearest_lattice_index(node.z,record.origin_z,record.h,component==2?0.0f:0.5f);return value_at_clamped(view,component,location.brick,ni,nj,nk);
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
		__global__ void advect_forward_level_kernel(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,Real dt,const unsigned char* protect_u,const unsigned char* protect_v,const unsigned char* protect_w,Real* u_out,Real* v_out,Real* w_out)
		{
			const DeviceAmrFieldLevelView view=views[level];const int bs=view.layout.brick_size,per_component=bs*bs*(bs+1),q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*3*per_component;if(q>=total)return;int component,brick,i,j,k;decode_face_work(view,q,component,brick,i,j,k);if(view.flags[brick]&BRICK_COVERED)return;const GpuBrickRecord record=hierarchy.levels[level].bricks[brick];const GpuAmrPoint point=face_node_point(record,component,i,j,k);const std::size_t index=component==0?view.layout.u_index(brick,i,j,k):(component==1?view.layout.v_index(brick,i,j,k):view.layout.w_index(brick,i,j,k));const Real current=component==0?view.u[index]:(component==1?view.v[index]:view.w[index]);const unsigned char protected_face=component==0?protect_u[index]:(component==1?protect_v[index]:protect_w[index]);const Real result=protected_face?current:sample_finest(hierarchy,views,component,rk2_trace(hierarchy,views,level,record,brick,point,dt),current);if(component==0)u_out[index]=result;else if(component==1)v_out[index]=result;else w_out[index]=result;
		}
		__global__ void advect_correct_level_kernel(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,const DeviceAmrFieldLevelView* forward_views,int level,Real dt,const unsigned char* protect_u,const unsigned char* protect_v,const unsigned char* protect_w,Real* u_out,Real* v_out,Real* w_out)
		{
			const DeviceAmrFieldLevelView view=views[level],forward=forward_views[level];const int bs=view.layout.brick_size,per_component=bs*bs*(bs+1),q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*3*per_component;if(q>=total)return;int component,brick,i,j,k;decode_face_work(view,q,component,brick,i,j,k);if(view.flags[brick]&BRICK_COVERED)return;const GpuBrickRecord record=hierarchy.levels[level].bricks[brick];const GpuAmrPoint point=face_node_point(record,component,i,j,k);const std::size_t index=component==0?view.layout.u_index(brick,i,j,k):(component==1?view.layout.v_index(brick,i,j,k):view.layout.w_index(brick,i,j,k));const Real current=component==0?view.u[index]:(component==1?view.v[index]:view.w[index]),phf=component==0?forward.u[index]:(component==1?forward.v[index]:forward.w[index]);const unsigned char protected_face=component==0?protect_u[index]:(component==1?protect_v[index]:protect_w[index]);Real result=current;if(!protected_face){Real mn,mx;const GpuAmrPoint back=rk2_trace(hierarchy,views,level,record,brick,point,dt);sample_finest_bounded(hierarchy,views,component,back,current,&mn,&mx);const GpuAmrPoint forward_point=rk2_trace(hierarchy,views,level,record,brick,point,-dt);const Real phn=sample_finest(hierarchy,forward_views,component,forward_point,phf);result=phf+Real(0.5)*(current-phn);result=max(mn,min(mx,result));}if(component==0)u_out[index]=result;else if(component==1)v_out[index]=result;else w_out[index]=result;
		}
		__device__ Real cell_value_clamped(const DeviceAmrFieldLevelView& view,int brick,int i,int j,int k)
		{
			const int bs=view.layout.brick_size;i=max(0,min(bs-1,i));j=max(0,min(bs-1,j));k=max(0,min(bs-1,k));return view.nut[view.layout.cell_index(brick,i,j,k)];
		}
		__device__ Real sample_cell_local_clamped(const DeviceAmrFieldLevelView& view,const GpuBrickRecord& record,int brick,GpuAmrPoint point)
		{
			Real x=(Real(point.x)-Real(record.origin_x))/Real(record.h)-Real(0.5),y=(Real(point.y)-Real(record.origin_y))/Real(record.h)-Real(0.5),z=(Real(point.z)-Real(record.origin_z))/Real(record.h)-Real(0.5);const int i0=static_cast<int>(floor(x)),j0=static_cast<int>(floor(y)),k0=static_cast<int>(floor(z));x-=i0;y-=j0;z-=k0;Real result=0;for(int dz=0;dz<2;++dz)for(int dy=0;dy<2;++dy)for(int dx=0;dx<2;++dx){const Real wx=dx?x:Real(1)-x,wy=dy?y:Real(1)-y,wz=dz?z:Real(1)-z;result+=wx*wy*wz*cell_value_clamped(view,brick,i0+dx,j0+dy,k0+dz);}return result;
		}
		__device__ __noinline__ Real cell_value_at_node(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int source_level,const GpuBrickRecord& source_record,int source_brick,int i,int j,int k)
		{
			const DeviceAmrFieldLevelView source_view=views[source_level];const int bs=source_view.layout.brick_size;if(i>=0&&j>=0&&k>=0&&i<bs&&j<bs&&k<bs)return cell_value_clamped(source_view,source_brick,i,j,k);int face=-1,outside=0;if(i<0){face=0;++outside;}else if(i>=bs){face=1;++outside;}if(j<0){face=2;++outside;}else if(j>=bs){face=3;++outside;}if(k<0){face=4;++outside;}else if(k>=bs){face=5;++outside;}if(outside==1){const int neighbor=source_view.neighbors[source_brick*6+face];if(neighbor>=0&&!(source_view.flags[neighbor]&BRICK_COVERED)){if(face==0)i+=bs;else if(face==1)i-=bs;else if(face==2)j+=bs;else if(face==3)j-=bs;else if(face==4)k+=bs;else k-=bs;return cell_value_clamped(source_view,neighbor,i,j,k);}}const GpuAmrPoint node{source_record.origin_x+(i+0.5f)*source_record.h,source_record.origin_y+(j+0.5f)*source_record.h,source_record.origin_z+(k+0.5f)*source_record.h};if(node.x<hierarchy.domain_lo.x||node.y<hierarchy.domain_lo.y||node.z<hierarchy.domain_lo.z||node.x>=hierarchy.domain_hi.x||node.y>=hierarchy.domain_hi.y||node.z>=hierarchy.domain_hi.z)return cell_value_clamped(source_view,source_brick,i,j,k);const GpuBrickLocation location=locate_finest(hierarchy,node);if(location.brick<0)return cell_value_clamped(source_view,source_brick,i,j,k);const GpuBrickRecord record=hierarchy.levels[location.level].bricks[location.brick];const DeviceAmrFieldLevelView view=views[location.level];if(location.level!=source_level)return sample_cell_local_clamped(view,record,location.brick,node);const int ni=nearest_lattice_index(node.x,record.origin_x,record.h,0.5f),nj=nearest_lattice_index(node.y,record.origin_y,record.h,0.5f),nk=nearest_lattice_index(node.z,record.origin_z,record.h,0.5f);return cell_value_clamped(view,location.brick,ni,nj,nk);
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
		__global__ void diffuse_level_kernel(GpuAmrHierarchyView hierarchy,const DeviceAmrFieldLevelView* views,int level,Real molecular_nu,Real dt,const unsigned char* protect_u,const unsigned char* protect_v,const unsigned char* protect_w,Real* u_out,Real* v_out,Real* w_out)
		{
			const DeviceAmrFieldLevelView view=views[level];const int bs=view.layout.brick_size,per_component=bs*bs*(bs+1),q=blockIdx.x*blockDim.x+threadIdx.x,total=view.brick_count*3*per_component;if(q>=total)return;int r=q,local=r%per_component;r/=per_component;const int component=r%3,brick=r/3;if(view.flags[brick]&BRICK_COVERED)return;int i=0,j=0,k=0;if(component==0){i=local%(bs+1);local/=bs+1;j=local%bs;k=local/bs;}else if(component==1){i=local%bs;local/=bs;j=local%(bs+1);k=local/(bs+1);}else{i=local%bs;local/=bs;j=local%bs;k=local/bs;}const GpuBrickRecord record=hierarchy.levels[level].bricks[brick];const std::size_t index=component==0?view.layout.u_index(brick,i,j,k):(component==1?view.layout.v_index(brick,i,j,k):view.layout.w_index(brick,i,j,k));const Real current=component==0?view.u[index]:(component==1?view.v[index]:view.w[index]);const unsigned char protected_face=component==0?protect_u[index]:(component==1?protect_v[index]:protect_w[index]);Real result=current;
			if(!protected_face){const Real neighbours=stencil_face_value(hierarchy,views,level,record,component,brick,i-1,j,k)+stencil_face_value(hierarchy,views,level,record,component,brick,i+1,j,k)+stencil_face_value(hierarchy,views,level,record,component,brick,i,j-1,k)+stencil_face_value(hierarchy,views,level,record,component,brick,i,j+1,k)+stencil_face_value(hierarchy,views,level,record,component,brick,i,j,k-1)+stencil_face_value(hierarchy,views,level,record,component,brick,i,j,k+1);int ci=i,cj=j,ck=k;if(component==0)--ci;else if(component==1)--cj;else --ck;const int ci2=ci+(component==0),cj2=cj+(component==1),ck2=ck+(component==2);const Real nut=Real(0.5)*(stencil_cell_value(hierarchy,views,level,record,brick,ci,cj,ck)+stencil_cell_value(hierarchy,views,level,record,brick,ci2,cj2,ck2));const Real h=hierarchy.levels[level].h;result=current+dt*(molecular_nu+nut)*(neighbours-Real(6)*current)/(h*h);}
			if(component==0)u_out[index]=result;else if(component==1)v_out[index]=result;else w_out[index]=result;
		}

		Vec3d face_centre(const BrickMetadata& brick,double h,int component,int i,int j,int k)
		{
			return brick.origin+Vec3d{(i+(component==0?0.0:0.5))*h,(j+(component==1?0.0:0.5))*h,(k+(component==2?0.0:0.5))*h};
		}
	}

	DeviceAmrAdvection::DeviceAmrAdvection(const AmrHierarchy& hierarchy,const TriangleBvh& fabric,double protection_cells):locator_(hierarchy)
	{
		if(!(protection_cells>=1.0))throw std::invalid_argument("AMR fabric advection protection must cover at least one cell");levels_.resize(hierarchy.levels().size());std::vector<DeviceAmrFieldLevelView> empty_views(levels_.size());
		for(std::size_t level=0;level<levels_.size();++level){const AmrLevel& metadata=hierarchy.levels()[level];Level& allocation=levels_[level];allocation.layout=BrickFieldLayout::make(hierarchy.brick_size(),hierarchy.ghost_cells());allocation.brick_count=static_cast<int>(metadata.bricks.size());if(!allocation.brick_count)continue;const std::size_t un=static_cast<std::size_t>(allocation.brick_count)*allocation.layout.u_stride,vn=static_cast<std::size_t>(allocation.brick_count)*allocation.layout.v_stride,wn=static_cast<std::size_t>(allocation.brick_count)*allocation.layout.w_stride;std::vector<unsigned char> pu(un,0),pv(vn,0),pw(wn,0);const int bs=hierarchy.brick_size();const double radius=protection_cells*metadata.h;
			for(int brick=0;brick<allocation.brick_count;++brick){const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;Aabb3d expanded{{record.origin.x-radius,record.origin.y-radius,record.origin.z-radius},{record.origin.x+bs*metadata.h+radius,record.origin.y+bs*metadata.h+radius,record.origin.z+bs*metadata.h+radius}};if(fabric.query_aabb(expanded).empty()){active_faces_+=static_cast<std::size_t>(3)*bs*bs*(bs+1);continue;}for(int component=0;component<3;++component){const int ni=component==0?bs+1:bs,nj=component==1?bs+1:bs,nk=component==2?bs+1:bs;for(int k=0;k<nk;++k)for(int j=0;j<nj;++j)for(int i=0;i<ni;++i){const Vec3d point=face_centre(record,metadata.h,component,i,j,k);const bool protect=fabric.distance(point,radius)<=radius;const std::size_t index=component==0?allocation.layout.u_index(brick,i,j,k):(component==1?allocation.layout.v_index(brick,i,j,k):allocation.layout.w_index(brick,i,j,k));if(component==0)pu[index]=protect;else if(component==1)pv[index]=protect;else pw[index]=protect;++active_faces_;protected_faces_+=protect;}}}
			allocation.u=allocate<Real>(un,"allocate AMR advection u scratch");allocation.v=allocate<Real>(vn,"allocate AMR advection v scratch");allocation.w=allocate<Real>(wn,"allocate AMR advection w scratch");allocation.forward_u=allocate<Real>(un,"allocate AMR forward u scratch");allocation.forward_v=allocate<Real>(vn,"allocate AMR forward v scratch");allocation.forward_w=allocate<Real>(wn,"allocate AMR forward w scratch");allocation.protect_u=allocate<unsigned char>(un,"allocate AMR advection u protection");allocation.protect_v=allocate<unsigned char>(vn,"allocate AMR advection v protection");allocation.protect_w=allocate<unsigned char>(wn,"allocate AMR advection w protection");check(cudaMemset(allocation.u,0,un*sizeof(Real)),"clear AMR advection u scratch");check(cudaMemset(allocation.v,0,vn*sizeof(Real)),"clear AMR advection v scratch");check(cudaMemset(allocation.w,0,wn*sizeof(Real)),"clear AMR advection w scratch");check(cudaMemset(allocation.forward_u,0,un*sizeof(Real)),"clear AMR forward u scratch");check(cudaMemset(allocation.forward_v,0,vn*sizeof(Real)),"clear AMR forward v scratch");check(cudaMemset(allocation.forward_w,0,wn*sizeof(Real)),"clear AMR forward w scratch");check(cudaMemcpy(allocation.protect_u,pu.data(),un,cudaMemcpyHostToDevice),"upload AMR u fabric protection");check(cudaMemcpy(allocation.protect_v,pv.data(),vn,cudaMemcpyHostToDevice),"upload AMR v fabric protection");check(cudaMemcpy(allocation.protect_w,pw.data(),wn,cudaMemcpyHostToDevice),"upload AMR w fabric protection");bytes_+=(un+vn+wn)*(2*sizeof(Real)+sizeof(unsigned char));}
		device_views_=allocate<DeviceAmrFieldLevelView>(levels_.size(),"allocate AMR advection field views");device_forward_views_=allocate<DeviceAmrFieldLevelView>(levels_.size(),"allocate AMR forward field views");bytes_+=2*levels_.size()*sizeof(DeviceAmrFieldLevelView);
	}

	DeviceAmrAdvection::~DeviceAmrAdvection(){for(Level& level:levels_)for(void* pointer:{(void*)level.u,(void*)level.v,(void*)level.w,(void*)level.forward_u,(void*)level.forward_v,(void*)level.forward_w,(void*)level.protect_u,(void*)level.protect_v,(void*)level.protect_w})if(pointer)cudaFree(pointer);if(device_views_)cudaFree(device_views_);if(device_forward_views_)cudaFree(device_forward_views_);}

	void DeviceAmrAdvection::advect(DeviceAmrFields& fields,Real dt)
	{
		if(!(dt>Real(0)))throw std::invalid_argument("AMR advection timestep must be positive");if(fields.level_count()!=static_cast<int>(levels_.size()))throw std::invalid_argument("AMR advection field hierarchy mismatch");std::vector<DeviceAmrFieldLevelView> views(levels_.size()),forward_views(levels_.size());for(int level=0;level<fields.level_count();++level){views[level]=fields.level_view(level);forward_views[level]=views[level];forward_views[level].u=levels_[level].forward_u;forward_views[level].v=levels_[level].forward_v;forward_views[level].w=levels_[level].forward_w;}check(cudaMemcpy(device_views_,views.data(),views.size()*sizeof(DeviceAmrFieldLevelView),cudaMemcpyHostToDevice),"upload AMR advection field views");check(cudaMemcpy(device_forward_views_,forward_views.data(),forward_views.size()*sizeof(DeviceAmrFieldLevelView),cudaMemcpyHostToDevice),"upload AMR forward field views");const GpuAmrHierarchyView hierarchy=locator_.view();
		// Build the forward state on every level before the correction sweep. Both sweeps
		// therefore sample one consistent hierarchy rather than a mixture of time levels.
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const int work=output.brick_count*3*output.layout.brick_size*output.layout.brick_size*(output.layout.brick_size+1);advect_forward_level_kernel<<<(work+255)/256,256>>>(hierarchy,device_views_,level,dt,output.protect_u,output.protect_v,output.protect_w,output.forward_u,output.forward_v,output.forward_w);}
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const int work=output.brick_count*3*output.layout.brick_size*output.layout.brick_size*(output.layout.brick_size+1);advect_correct_level_kernel<<<(work+255)/256,256>>>(hierarchy,device_views_,device_forward_views_,level,dt,output.protect_u,output.protect_v,output.protect_w,output.u,output.v,output.w);}
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const DeviceAmrFieldLevelView view=views[level];check(cudaMemcpy(view.u,output.u,static_cast<std::size_t>(output.brick_count)*output.layout.u_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR advected u");check(cudaMemcpy(view.v,output.v,static_cast<std::size_t>(output.brick_count)*output.layout.v_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR advected v");check(cudaMemcpy(view.w,output.w,static_cast<std::size_t>(output.brick_count)*output.layout.w_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR advected w");}check(cudaDeviceSynchronize(),"AMR advection");
	}

	void DeviceAmrAdvection::diffuse_smagorinsky(DeviceAmrFields& fields,Real molecular_nu,Real cs,Real dt)
	{
		if(molecular_nu<Real(0)||cs<Real(0)||!(dt>Real(0)))throw std::invalid_argument("invalid AMR LES/diffusion parameters");if(fields.level_count()!=static_cast<int>(levels_.size()))throw std::invalid_argument("AMR LES field hierarchy mismatch");std::vector<DeviceAmrFieldLevelView> views(levels_.size());for(int level=0;level<fields.level_count();++level)views[level]=fields.level_view(level);check(cudaMemcpy(device_views_,views.data(),views.size()*sizeof(DeviceAmrFieldLevelView),cudaMemcpyHostToDevice),"upload AMR LES field views");const GpuAmrHierarchyView hierarchy=locator_.view();
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const int cells=output.brick_count*output.layout.brick_size*output.layout.brick_size*output.layout.brick_size;smagorinsky_kernel<<<(cells+255)/256,256>>>(hierarchy,device_views_,level,cs);}
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const int work=output.brick_count*3*output.layout.brick_size*output.layout.brick_size*(output.layout.brick_size+1);diffuse_level_kernel<<<(work+255)/256,256>>>(hierarchy,device_views_,level,molecular_nu,dt,output.protect_u,output.protect_v,output.protect_w,output.u,output.v,output.w);}
		for(int level=0;level<fields.level_count();++level){Level& output=levels_[level];if(!output.brick_count)continue;const DeviceAmrFieldLevelView view=views[level];check(cudaMemcpy(view.u,output.u,static_cast<std::size_t>(output.brick_count)*output.layout.u_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR diffused u");check(cudaMemcpy(view.v,output.v,static_cast<std::size_t>(output.brick_count)*output.layout.v_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR diffused v");check(cudaMemcpy(view.w,output.w,static_cast<std::size_t>(output.brick_count)*output.layout.w_stride*sizeof(Real),cudaMemcpyDeviceToDevice),"commit AMR diffused w");}check(cudaDeviceSynchronize(),"AMR Smagorinsky diffusion");
	}
}
