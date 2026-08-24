#include "core/fluid/amr_fields.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <stdexcept>

namespace paracfd::core
{
	BrickFieldLayout BrickFieldLayout::make(int bs, int g)
	{
		BrickFieldLayout l; l.brick_size=bs; l.ghost=g; l.cell_n=bs+2*g; l.u_nx=bs+1+2*g; l.v_ny=bs+1+2*g; l.w_nz=bs+1+2*g;
		l.cell_stride=static_cast<std::size_t>(l.cell_n)*l.cell_n*l.cell_n;
		l.u_stride=static_cast<std::size_t>(l.u_nx)*l.cell_n*l.cell_n; l.v_stride=static_cast<std::size_t>(l.cell_n)*l.v_ny*l.cell_n; l.w_stride=static_cast<std::size_t>(l.cell_n)*l.cell_n*l.w_nz; return l;
	}
	AmrHostFields::AmrHostFields(const AmrHierarchy& h)
	{
		levels_.resize(h.levels().size());
		for (std::size_t l=0;l<levels_.size();++l)
		{
			auto& f=levels_[l]; f.layout=BrickFieldLayout::make(h.brick_size(),h.ghost_cells()); const std::size_t n=h.levels()[l].bricks.size();
			f.u.assign(n*f.layout.u_stride,Real{}); f.v.assign(n*f.layout.v_stride,Real{}); f.w.assign(n*f.layout.w_stride,Real{}); f.p.assign(n*f.layout.cell_stride,Real{}); f.nut.assign(n*f.layout.cell_stride,Real{}); f.temp.assign(n*f.layout.cell_stride,Real{});
		}
	}
	std::size_t AmrHostFields::bytes() const { std::size_t n=0; for(const auto& f:levels_) n+=(f.u.size()+f.v.size()+f.w.size()+f.p.size()+f.nut.size()+f.temp.size())*sizeof(Real); return n; }

	void AmrHostFields::exchange_same_level_pressure_halos(const AmrHierarchy& h)
	{
		for(std::size_t l=0;l<levels_.size();++l)
		{
			auto& f=levels_[l]; const int bs=f.layout.brick_size;
			for(int b=0;b<static_cast<int>(h.levels()[l].bricks.size());++b)
			{
				const auto& meta=h.levels()[l].bricks[b]; if(!meta.active()) continue;
				for(int face=0;face<6;++face)
				{
					const int nb=meta.same_level_neighbor[face]; if(nb<0 || !h.levels()[l].bricks[nb].active()) continue;
					for(int a=0;a<bs;++a) for(int c=0;c<bs;++c)
					{
						int di=0,dj=a,dk=c,si=0,sj=a,sk=c;
						if(face==0){di=-1;si=bs-1;} else if(face==1){di=bs;si=0;}
						else if(face==2){di=a;dj=-1;dk=c;si=a;sj=bs-1;sk=c;} else if(face==3){di=a;dj=bs;dk=c;si=a;sj=0;sk=c;}
						else if(face==4){di=a;dj=c;dk=-1;si=a;sj=c;sk=bs-1;} else {di=a;dj=c;dk=bs;si=a;sj=c;sk=0;}
						f.p[f.layout.cell_index(b,di,dj,dk)]=f.p[f.layout.cell_index(nb,si,sj,sk)];
					}
				}
			}
		}
	}

	namespace
	{
		void check(cudaError_t e,const char* what){if(e!=cudaSuccess) throw std::runtime_error(std::string(what)+": "+cudaGetErrorString(e));}
		template<class T> T* alloc(std::size_t n){T* p=nullptr;check(cudaMalloc(&p,n*sizeof(T)),"cudaMalloc AMR pool");check(cudaMemset(p,0,n*sizeof(T)),"cudaMemset AMR pool");return p;}
		__global__ void fill_kernel(Real* p,std::size_t n,Real v){std::size_t i=static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;if(i<n)p[i]=v;}
		__device__ void atomic_max_positive(Real* destination,Real value)
		{
			if constexpr(sizeof(Real)==sizeof(float))atomicMax(reinterpret_cast<unsigned int*>(destination),__float_as_uint(static_cast<float>(value)));
			else atomicMax(reinterpret_cast<unsigned long long*>(destination),static_cast<unsigned long long>(__double_as_longlong(static_cast<double>(value))));
		}
		__global__ void max_abs_kernel(const Real* values,std::size_t count,Real* maximum)
		{
			const std::size_t q=static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;if(q<count)atomic_max_positive(maximum,abs(values[q]));
		}
		__global__ void max_abs_active_velocity_kernel(const Real* u,const Real* v,const Real* w,BrickFieldLayout layout,const std::uint32_t* flags,int bricks,Real* maximum)
		{
			const int bs=layout.brick_size,per_component=bs*bs*(bs+1),q=blockIdx.x*blockDim.x+threadIdx.x,total=bricks*3*per_component;if(q>=total)return;int local=q%per_component,r=q/per_component,component=r%3,brick=r/3;if(flags[brick]&BRICK_COVERED)return;int i=0,j=0,k=0;if(component==0){i=local%(bs+1);local/=bs+1;j=local%bs;k=local/bs;}else if(component==1){i=local%bs;local/=bs;j=local%(bs+1);k=local/(bs+1);}else{i=local%bs;local/=bs;j=local%bs;k=local/bs;}const std::size_t index=component==0?layout.u_index(brick,i,j,k):(component==1?layout.v_index(brick,i,j,k):layout.w_index(brick,i,j,k));atomic_max_positive(maximum,abs(component==0?u[index]:(component==1?v[index]:w[index])));
		}
		__global__ void initialize_velocity_kernel(Real* u,Real* v,Real* w,BrickFieldLayout layout,int bricks,Real speed)
		{
			std::size_t q=static_cast<std::size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
			const std::size_t un=static_cast<std::size_t>(bricks)*layout.u_stride,vn=static_cast<std::size_t>(bricks)*layout.v_stride,wn=static_cast<std::size_t>(bricks)*layout.w_stride,total=un+vn+wn;
			if(q<un)u[q]=speed;else if(q<un+vn)v[q-un]=Real(0);else if(q<total)w[q-un-vn]=Real(0);
		}
		__global__ void halo_kernel(Real* p,BrickFieldLayout layout,const int* neighbors,int bricks)
		{
			const int bs=layout.brick_size; const int per_face=bs*bs; const int q=blockIdx.x*blockDim.x+threadIdx.x; const int total=bricks*6*per_face; if(q>=total)return;
			int r=q; const int ac=r%per_face; r/=per_face; const int face=r%6; const int b=r/6; const int nb=neighbors[b*6+face]; if(nb<0)return; const int a=ac%bs,c=ac/bs;
			int di=0,dj=a,dk=c,si=0,sj=a,sk=c;
			if(face==0){di=-1;si=bs-1;} else if(face==1){di=bs;si=0;}
			else if(face==2){di=a;dj=-1;dk=c;si=a;sj=bs-1;sk=c;} else if(face==3){di=a;dj=bs;dk=c;si=a;sj=0;sk=c;}
			else if(face==4){di=a;dj=c;dk=-1;si=a;sj=c;sk=bs-1;} else {di=a;dj=c;dk=bs;si=a;sj=c;sk=0;}
			p[layout.cell_index(b,di,dj,dk)]=p[layout.cell_index(nb,si,sj,sk)];
		}

		__global__ void external_boundary_kernel(Real* u,Real* v,Real* w,Real* p,BrickFieldLayout layout,
			const std::uint32_t* flags,int bricks,Real speed)
		{
			const int span=layout.brick_size+1,per_face=span*span,q=blockIdx.x*blockDim.x+threadIdx.x,total=bricks*6*per_face;if(q>=total)return;
			int r=q,a=r%span;r/=span;int c=r%span;r/=span;int face=r%6,brick=r/6;const int bs=layout.brick_size;const std::uint32_t flag=flags[brick];if(flag&BRICK_COVERED)return;
			const std::uint32_t required=1u<<face;if(!(flag&required))return;
			if(face==0) // X-min prescribed freestream
			{
				if(a<bs&&c<bs){p[layout.cell_index(brick,-1,a,c)]=p[layout.cell_index(brick,0,a,c)];u[layout.u_index(brick,0,a,c)]=speed;u[layout.u_index(brick,-1,a,c)]=speed;}
				if(c<bs){v[layout.v_index(brick,-1,a,c)]=-v[layout.v_index(brick,0,a,c)];}
				if(a<bs){w[layout.w_index(brick,-1,a,c)]=-w[layout.w_index(brick,0,a,c)];}
			}
			else if(face==1) // X-max convective/zero-gradient velocity, p=0 at boundary
			{
				if(a<bs&&c<bs){p[layout.cell_index(brick,bs,a,c)]=-p[layout.cell_index(brick,bs-1,a,c)];u[layout.u_index(brick,bs+1,a,c)]=u[layout.u_index(brick,bs,a,c)];}
				if(c<bs)v[layout.v_index(brick,bs,a,c)]=v[layout.v_index(brick,bs-1,a,c)];
				if(a<bs)w[layout.w_index(brick,bs,a,c)]=w[layout.w_index(brick,bs-1,a,c)];
			}
			else if(face==2||face==3) // Y far field: free slip
			{
				const bool upper=face==3;const int ghost=upper?bs:-1,inner=upper?bs-1:0,boundary=upper?bs:0;
				if(a<bs&&c<bs){p[layout.cell_index(brick,a,ghost,c)]=p[layout.cell_index(brick,a,inner,c)];v[layout.v_index(brick,a,boundary,c)]=Real(0);v[layout.v_index(brick,a,upper?bs+1:-1,c)]=Real(0);}
				if(c<bs)u[layout.u_index(brick,a,ghost,c)]=u[layout.u_index(brick,a,inner,c)];
				if(a<bs)w[layout.w_index(brick,a,ghost,c)]=w[layout.w_index(brick,a,inner,c)];
			}
			else // Z far field: free slip, identical at Z-min and Z-max (no ground)
			{
				const bool upper=face==5;const int ghost=upper?bs:-1,inner=upper?bs-1:0,boundary=upper?bs:0;
				if(a<bs&&c<bs){p[layout.cell_index(brick,a,c,ghost)]=p[layout.cell_index(brick,a,c,inner)];w[layout.w_index(brick,a,c,boundary)]=Real(0);w[layout.w_index(brick,a,c,upper?bs+1:-1)]=Real(0);}
				if(c<bs)u[layout.u_index(brick,a,c,ghost)]=u[layout.u_index(brick,a,c,inner)];
				if(a<bs)v[layout.v_index(brick,a,c,ghost)]=v[layout.v_index(brick,a,c,inner)];
			}
		}

		__host__ __device__ std::uint64_t gpu_coord_hash(int x, int y, int z)
		{
			std::uint64_t key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(x)) << 42) ^
				(static_cast<std::uint64_t>(static_cast<std::uint32_t>(y)) << 21) ^
				static_cast<std::uint32_t>(z);
			std::uint64_t q = key + 0x9e3779b97f4a7c15ull;
			q = (q ^ (q >> 30)) * 0xbf58476d1ce4e5b9ull;
			q = (q ^ (q >> 27)) * 0x94d049bb133111ebull;
			return q ^ (q >> 31);
		}

		__device__ GpuBrickLocation locate_finest(GpuAmrHierarchyView hierarchy, GpuAmrPoint point)
		{
			GpuBrickLocation out;
			if (point.x < hierarchy.domain_lo.x || point.y < hierarchy.domain_lo.y || point.z < hierarchy.domain_lo.z ||
				point.x >= hierarchy.domain_hi.x || point.y >= hierarchy.domain_hi.y || point.z >= hierarchy.domain_hi.z) return out;
			for (int level = hierarchy.level_count - 1; level >= 0; --level)
			{
				const GpuAmrLevelView lev = hierarchy.levels[level];
				const float brick_width = hierarchy.brick_size * lev.h;
				const int bx = static_cast<int>(floorf((point.x - hierarchy.domain_lo.x) / brick_width));
				const int by = static_cast<int>(floorf((point.y - hierarchy.domain_lo.y) / brick_width));
				const int bz = static_cast<int>(floorf((point.z - hierarchy.domain_lo.z) / brick_width));
				std::uint32_t slot = static_cast<std::uint32_t>(gpu_coord_hash(bx, by, bz)) & lev.lookup_mask;
				int brick = -1;
				for (int probe = 0; probe < lev.lookup_size; ++probe)
				{
					const GpuBrickLookupEntry entry = lev.lookup[slot];
					if (entry.brick < 0) break;
					if (entry.x == bx && entry.y == by && entry.z == bz) { brick = entry.brick; break; }
					slot = (slot + 1) & lev.lookup_mask;
				}
				if (brick < 0 || brick >= lev.brick_count || (lev.bricks[brick].flags & BRICK_COVERED)) continue;
				const GpuBrickRecord record = lev.bricks[brick];
				out.level = level; out.brick = brick;
				out.i = max(0, min(hierarchy.brick_size - 1, static_cast<int>(floorf((point.x - record.origin_x) / record.h))));
				out.j = max(0, min(hierarchy.brick_size - 1, static_cast<int>(floorf((point.y - record.origin_y) / record.h))));
				out.k = max(0, min(hierarchy.brick_size - 1, static_cast<int>(floorf((point.z - record.origin_z) / record.h))));
				return out;
			}
			return out;
		}

		__global__ void locate_points_kernel(GpuAmrHierarchyView hierarchy, const GpuAmrPoint* points,
			GpuBrickLocation* locations, int count)
		{
			const int point = blockIdx.x * blockDim.x + threadIdx.x;
			if (point < count) locations[point] = locate_finest(hierarchy, points[point]);
		}
	}

	DeviceAmrFields::DeviceAmrFields(const AmrHierarchy& h)
	{
		levels_.resize(h.levels().size());
		for(std::size_t l=0;l<levels_.size();++l)
		{
			auto& d=levels_[l]; d.layout=BrickFieldLayout::make(h.brick_size(),h.ghost_cells()); d.brick_count=static_cast<int>(h.levels()[l].bricks.size()); if(d.brick_count==0)continue;
			d.u=alloc<Real>(d.brick_count*d.layout.u_stride);d.v=alloc<Real>(d.brick_count*d.layout.v_stride);d.w=alloc<Real>(d.brick_count*d.layout.w_stride);d.p=alloc<Real>(d.brick_count*d.layout.cell_stride);d.nut=alloc<Real>(d.brick_count*d.layout.cell_stride);
			std::vector<int> nb(static_cast<std::size_t>(d.brick_count)*6,-1); for(int b=0;b<d.brick_count;++b)for(int f=0;f<6;++f){int n=h.levels()[l].bricks[b].same_level_neighbor[f];if(n>=0&&h.levels()[l].bricks[n].active())nb[b*6+f]=n;}
			d.neighbors=alloc<int>(nb.size());check(cudaMemcpy(d.neighbors,nb.data(),nb.size()*sizeof(int),cudaMemcpyHostToDevice),"upload AMR neighbors");std::vector<std::uint32_t> flags(d.brick_count);for(int b=0;b<d.brick_count;++b)flags[b]=h.levels()[l].bricks[b].flags;d.flags=alloc<std::uint32_t>(flags.size());check(cudaMemcpy(d.flags,flags.data(),flags.size()*sizeof(std::uint32_t),cudaMemcpyHostToDevice),"upload AMR brick flags"); bytes_+=(d.brick_count*(d.layout.u_stride+d.layout.v_stride+d.layout.w_stride+2*d.layout.cell_stride))*sizeof(Real)+nb.size()*sizeof(int)+flags.size()*sizeof(std::uint32_t);
		}
		max_abs_scratch_=alloc<Real>(1);bytes_+=sizeof(Real);
	}
	DeviceAmrFields::~DeviceAmrFields(){for(auto& d:levels_)for(void* p:{(void*)d.u,(void*)d.v,(void*)d.w,(void*)d.p,(void*)d.nut,(void*)d.temp,(void*)d.neighbors,(void*)d.flags})if(p)cudaFree(p);if(max_abs_scratch_)cudaFree(max_abs_scratch_);}
	void DeviceAmrFields::upload(const AmrHostFields& h){for(std::size_t l=0;l<levels_.size();++l){auto& d=levels_[l];if(!d.brick_count)continue;const auto& s=h.levels()[l];check(cudaMemcpy(d.u,s.u.data(),s.u.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload u");check(cudaMemcpy(d.v,s.v.data(),s.v.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload v");check(cudaMemcpy(d.w,s.w.data(),s.w.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload w");check(cudaMemcpy(d.p,s.p.data(),s.p.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload p");check(cudaMemcpy(d.nut,s.nut.data(),s.nut.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload nut");}}
	void DeviceAmrFields::download(AmrHostFields& h)const{for(std::size_t l=0;l<levels_.size();++l){const auto& d=levels_[l];if(!d.brick_count)continue;auto& s=h.levels()[l];check(cudaMemcpy(s.u.data(),d.u,s.u.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download u");check(cudaMemcpy(s.v.data(),d.v,s.v.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download v");check(cudaMemcpy(s.w.data(),d.w,s.w.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download w");check(cudaMemcpy(s.p.data(),d.p,s.p.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download p");check(cudaMemcpy(s.nut.data(),d.nut,s.nut.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download nut");}}
	void DeviceAmrFields::download_pressure(AmrHostFields& h)const{for(std::size_t l=0;l<levels_.size();++l){const auto& d=levels_[l];if(d.brick_count)check(cudaMemcpy(h.levels()[l].p.data(),d.p,h.levels()[l].p.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download p");}}
	void DeviceAmrFields::fill(Real v){for(auto& d:levels_)if(d.brick_count){auto run=[&](Real* p,std::size_t n){int block=256;fill_kernel<<<static_cast<unsigned>((n+block-1)/block),block>>>(p,n,v);};run(d.u,d.brick_count*d.layout.u_stride);run(d.v,d.brick_count*d.layout.v_stride);run(d.w,d.brick_count*d.layout.w_stride);run(d.p,d.brick_count*d.layout.cell_stride);run(d.nut,d.brick_count*d.layout.cell_stride);}check(cudaDeviceSynchronize(),"fill AMR fields");}
	void DeviceAmrFields::initialize_freestream(Real speed){for(auto& d:levels_)if(d.brick_count){const std::size_t n=static_cast<std::size_t>(d.brick_count)*(d.layout.u_stride+d.layout.v_stride+d.layout.w_stride);initialize_velocity_kernel<<<static_cast<unsigned>((n+255)/256),256>>>(d.u,d.v,d.w,d.layout,d.brick_count,speed);fill_kernel<<<static_cast<unsigned>((static_cast<std::size_t>(d.brick_count)*d.layout.cell_stride+255)/256),256>>>(d.p,static_cast<std::size_t>(d.brick_count)*d.layout.cell_stride,Real(0));}check(cudaDeviceSynchronize(),"initialize AMR freestream");}
	void DeviceAmrFields::apply_external_aero_boundaries(Real speed){for(auto& d:levels_)if(d.brick_count){const int n=d.brick_count*6*(d.layout.brick_size+1)*(d.layout.brick_size+1);external_boundary_kernel<<<(n+255)/256,256>>>(d.u,d.v,d.w,d.p,d.layout,d.flags,d.brick_count,speed);}check(cudaDeviceSynchronize(),"external aerodynamic boundaries");}
	void DeviceAmrFields::exchange_same_level_pressure_halos(){for(auto& d:levels_)if(d.brick_count){int n=d.brick_count*6*d.layout.brick_size*d.layout.brick_size;halo_kernel<<<(n+255)/256,256>>>(d.p,d.layout,d.neighbors,d.brick_count);}check(cudaDeviceSynchronize(),"AMR halo exchange");}
	double max_abs_device_values(const Real* values,std::size_t count,Real* scratch)
	{
		if(!scratch)throw std::invalid_argument("max-absolute reduction has no scratch scalar");check(cudaMemset(scratch,0,sizeof(Real)),"clear max-absolute scalar");if(values&&count)max_abs_kernel<<<static_cast<unsigned>((count+255)/256),256>>>(values,count,scratch);Real host=0;check(cudaMemcpy(&host,scratch,sizeof(Real),cudaMemcpyDeviceToHost),"download max-absolute scalar");return static_cast<double>(host);
	}
	double DeviceAmrFields::max_abs_velocity() const
	{
		check(cudaMemset(max_abs_scratch_,0,sizeof(Real)),"clear pooled velocity maximum");for(const auto& d:levels_)if(d.brick_count){const int work=d.brick_count*3*d.layout.brick_size*d.layout.brick_size*(d.layout.brick_size+1);max_abs_active_velocity_kernel<<<(work+255)/256,256>>>(d.u,d.v,d.w,d.layout,d.flags,d.brick_count,max_abs_scratch_);}Real host=0;check(cudaMemcpy(&host,max_abs_scratch_,sizeof(Real),cudaMemcpyDeviceToHost),"download active velocity maximum");return static_cast<double>(host);
	}

	DeviceAmrLocator::DeviceAmrLocator(const AmrHierarchy& hierarchy)
	{
		allocations_.resize(hierarchy.levels().size());
		std::vector<GpuAmrLevelView> views(hierarchy.levels().size());
		for (std::size_t level = 0; level < hierarchy.levels().size(); ++level)
		{
			const AmrLevel& source = hierarchy.levels()[level];
			std::vector<GpuBrickLookupEntry> lookup(source.lookup.size());
			for (std::size_t i = 0; i < lookup.size(); ++i)
				lookup[i] = {source.lookup[i].coord.x, source.lookup[i].coord.y, source.lookup[i].coord.z, source.lookup[i].brick_id};
			std::vector<GpuBrickRecord> bricks(source.bricks.size());
			for (std::size_t i = 0; i < bricks.size(); ++i)
			{
				const BrickMetadata& brick = source.bricks[i];
				bricks[i] = {static_cast<float>(brick.origin.x), static_cast<float>(brick.origin.y),
					static_cast<float>(brick.origin.z), static_cast<float>(brick.h), brick.flags};
			}
			Allocation& allocation = allocations_[level];
			allocation.lookup = alloc<GpuBrickLookupEntry>(lookup.size());
			allocation.bricks = bricks.empty() ? nullptr : alloc<GpuBrickRecord>(bricks.size());
			check(cudaMemcpy(allocation.lookup, lookup.data(), lookup.size() * sizeof(GpuBrickLookupEntry), cudaMemcpyHostToDevice), "upload AMR lookup");
			if (!bricks.empty()) check(cudaMemcpy(allocation.bricks, bricks.data(), bricks.size() * sizeof(GpuBrickRecord), cudaMemcpyHostToDevice), "upload AMR brick records");
			views[level] = {allocation.lookup, allocation.bricks, source.lookup_mask, static_cast<int>(source.lookup.size()), static_cast<int>(source.bricks.size()), static_cast<float>(source.h)};
			bytes_ += lookup.size() * sizeof(GpuBrickLookupEntry) + bricks.size() * sizeof(GpuBrickRecord);
		}
		device_levels_ = alloc<GpuAmrLevelView>(views.size());
		check(cudaMemcpy(device_levels_, views.data(), views.size() * sizeof(GpuAmrLevelView), cudaMemcpyHostToDevice), "upload AMR level views");
		const Aabb3d& domain = hierarchy.domain();
		view_ = {device_levels_, static_cast<int>(views.size()), hierarchy.brick_size(),
			{static_cast<float>(domain.lo.x), static_cast<float>(domain.lo.y), static_cast<float>(domain.lo.z)},
			{static_cast<float>(domain.hi.x), static_cast<float>(domain.hi.y), static_cast<float>(domain.hi.z)}};
		bytes_ += views.size() * sizeof(GpuAmrLevelView);
	}

	DeviceAmrLocator::~DeviceAmrLocator()
	{
		for (Allocation& allocation : allocations_)
		{
			if (allocation.lookup) cudaFree(allocation.lookup);
			if (allocation.bricks) cudaFree(allocation.bricks);
		}
		if (device_levels_) cudaFree(device_levels_);
	}

	void DeviceAmrLocator::locate_points(const GpuAmrPoint* points, GpuBrickLocation* locations, int count) const
	{
		if (count <= 0) return;
		locate_points_kernel<<<(count + 255) / 256, 256>>>(view_, points, locations, count);
		check(cudaDeviceSynchronize(), "GPU AMR point location");
	}
}
