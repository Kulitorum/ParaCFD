#include "core/fluid/amr_pressure.h"

#include <cuda_runtime.h>

#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform_reduce.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace paracfd::core
{
	struct DeviceCompositeAmrLevelView
	{
		int level_offset = 0;
		int brick_count = 0;
		float h = 0;
		const int* neighbors = nullptr;
		const std::uint32_t* flags = nullptr;
	};
	struct DeviceCompositeAmrFluxLevelView
	{
		int level_offset = 0;
		float h = 0;
		DeviceAmrFieldLevelView fields;
	};

	namespace
	{
		void check(cudaError_t error, const char* operation)
		{
			if (error != cudaSuccess) throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(error));
		}
		template<class T> T* upload(const std::vector<T>& host, const char* operation)
		{
			if(host.empty())return nullptr;T* device=nullptr;check(cudaMalloc(&device,host.size()*sizeof(T)),operation);check(cudaMemcpy(device,host.data(),host.size()*sizeof(T),cudaMemcpyHostToDevice),operation);return device;
		}
		template<class T> T* allocate_zero(std::size_t count,const char* operation)
		{
			if(!count)return nullptr;T* device=nullptr;check(cudaMalloc(&device,count*sizeof(T)),operation);check(cudaMemset(device,0,count*sizeof(T)),operation);return device;
		}
		__host__ __device__ int local_index(int bs, int i, int j, int k) { return (k * bs + j) * bs + i; }
		__global__ void structured_apply_kernel(const DeviceCompositeAmrLevelView* levels, int level, int bs,
			bool outlet, const unsigned char* active, const std::uint8_t* cut_face_mask, const Real* pressure, Real* output)
		{
			const DeviceCompositeAmrLevelView source = levels[level]; const int cells_per_brick = bs * bs * bs;
			const int work = blockIdx.x * blockDim.x + threadIdx.x; if (work >= source.brick_count * cells_per_brick) return;
			const int brick = work / cells_per_brick, local = work % cells_per_brick; if (source.flags[brick] & BRICK_COVERED) return;
			const int i = local % bs, j = (local / bs) % bs, k = local / (bs * bs); const int a = source.level_offset + work;if(!active[a])return; Real value = Real(0);
			const int coordinate[3] = {i,j,k};
			for (int axis = 0; axis < 3; ++axis) for (int sign = -1; sign <= 1; sign += 2)
			{
				int q[3] = {i,j,k}; q[axis] += sign; int b = -1;
				if (q[axis] >= 0 && q[axis] < bs) b = source.level_offset + brick * cells_per_brick + local_index(bs, q[0], q[1], q[2]);
				else
				{
					const int neighbour = source.neighbors[brick * 6 + 2 * axis + (sign > 0)];
					if (neighbour >= 0 && !(source.flags[neighbour] & BRICK_COVERED)) { q[axis] = sign > 0 ? 0 : bs - 1; b = source.level_offset + neighbour * cells_per_brick + local_index(bs, q[0], q[1], q[2]); }
				}
				const int lower=sign>0?a:b;if(b>=0&&active[b]&&!(cut_face_mask[lower]&(1u<<axis))) value += Real(source.h) * (pressure[a] - pressure[b]);
			}
			if (outlet && (source.flags[brick] & BRICK_XMAX) && coordinate[0] == bs - 1) value += Real(2) * Real(source.h) * pressure[a]; output[a] = value;
		}

		__global__ void coarse_fine_apply_kernel(const int* coarse, const int* fine, const Real* coefficient,
			int count, const unsigned char* active, const Real* pressure, Real* output)
		{
			const int edge = blockIdx.x * blockDim.x + threadIdx.x; if (edge >= count) return; const int a = coarse[edge], b = fine[edge];if(!active[a]||!active[b])return; const Real flux = coefficient[edge] * (pressure[a] - pressure[b]); atomicAdd(output + a, flux); atomicAdd(output + b, -flux);
		}
		__global__ void gauge_apply_kernel(const int* dof,const Real* coefficient,int count,const Real* pressure,Real* output){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<count)output[dof[q]]+=coefficient[q]*pressure[dof[q]];}
		__global__ void subtract_kernel(Real* residual,const Real* value,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)residual[q]=active[q]?residual[q]-value[q]:Real(0);}
		__global__ void update_pressure_residual_kernel(Real* pressure,Real* residual,const Real* direction,const Real* Ad,Real alpha,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n&&active[q]){pressure[q]+=alpha*direction[q];residual[q]-=alpha*Ad[q];}}
		__global__ void precondition_kernel(const Real* residual,Real* z,const Real* diagonal,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)z[q]=active[q]&&diagonal[q]>Real(0)?residual[q]/diagonal[q]:Real(0);}
		__global__ void update_direction_kernel(const Real* z,Real* direction,Real beta,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)direction[q]=active[q]?z[q]+beta*direction[q]:Real(0);}
		__global__ void copy_active_kernel(const Real* source,Real* destination,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)destination[q]=active[q]?source[q]:Real(0);}
		__global__ void restrict_aggregate_kernel(const Real* residual,const unsigned char* active,const int* aggregate,int n,int base_offset,int base_cells,Real* coarse_rhs)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=n||!active[q])return;const int coarse=aggregate[q]-base_offset;if(coarse>=0&&coarse<base_cells)atomicAdd(coarse_rhs+coarse,residual[q]);
		}
		__global__ void coarse_jacobi_kernel(const Real* rhs,const Real* x,Real* next,const Real* positive,const Real* diagonal,const int* neighbors,int bricks,int bs,Real omega)
		{
			const int cells_per_brick=bs*bs*bs,work=blockIdx.x*blockDim.x+threadIdx.x;if(work>=bricks*cells_per_brick)return;if(!(diagonal[work]>Real(0))){next[work]=Real(0);return;}const int brick=work/cells_per_brick,local=work%cells_per_brick,i=local%bs,j=(local/bs)%bs,k=local/(bs*bs);Real Ax=diagonal[work]*x[work];for(int axis=0;axis<3;++axis){int c[3]={i,j,k},other=-1;if(++c[axis]<bs)other=brick*cells_per_brick+local_index(bs,c[0],c[1],c[2]);else{const int neighbour=neighbors[brick*6+2*axis+1];if(neighbour>=0){c[axis]=0;other=neighbour*cells_per_brick+local_index(bs,c[0],c[1],c[2]);}}if(other>=0)Ax-=positive[work*3+axis]*x[other];c[0]=i;c[1]=j;c[2]=k;int lower=brick;if(--c[axis]<0){lower=neighbors[brick*6+2*axis];if(lower>=0)c[axis]=bs-1;}if(lower>=0){const int lower_cell=lower*cells_per_brick+local_index(bs,c[0],c[1],c[2]);Ax-=positive[lower_cell*3+axis]*x[lower_cell];}}next[work]=x[work]+omega*(rhs[work]-Ax)/diagonal[work];
		}
		__global__ void coarse_extra_jacobi_kernel(const int* a,const int* b,const Real* coefficient,int count,const Real* x,const Real* diagonal,Real omega,Real* next)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int qa=a[edge],qb=b[edge];const Real c=omega*coefficient[edge];atomicAdd(next+qa,c*x[qb]/diagonal[qa]);atomicAdd(next+qb,c*x[qa]/diagonal[qb]);
		}
		__global__ void prolong_aggregate_kernel(const int* aggregate,const unsigned char* active,int n,int base_offset,int base_cells,const Real* coarse_x,Real* output)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=n||!active[q])return;const int coarse=aggregate[q]-base_offset;if(coarse>=0&&coarse<base_cells)output[q]+=coarse_x[coarse];
		}
		__device__ Real positive_face_value(const DeviceAmrFieldLevelView& fields,int brick,int axis,int i,int j,int k)
		{
			if(axis==0)return fields.u[fields.layout.u_index(brick,i+1,j,k)];if(axis==1)return fields.v[fields.layout.v_index(brick,i,j+1,k)];return fields.w[fields.layout.w_index(brick,i,j,k+1)];
		}
		__device__ void set_positive_face_value(const DeviceAmrFieldLevelView& fields,int brick,int axis,int i,int j,int k,Real value)
		{
			if(axis==0)fields.u[fields.layout.u_index(brick,i+1,j,k)]=value;else if(axis==1)fields.v[fields.layout.v_index(brick,i,j+1,k)]=value;else fields.w[fields.layout.w_index(brick,i,j,k+1)]=value;
		}
		__global__ void structured_integrated_flux_kernel(const DeviceCompositeAmrFluxLevelView* levels,int level,int bs,
			const unsigned char* active,const std::uint8_t* cut_face_mask,Real* integrated)
		{
			const DeviceCompositeAmrFluxLevelView source=levels[level];const DeviceAmrFieldLevelView fields=source.fields;const int cells_per_brick=bs*bs*bs,work=blockIdx.x*blockDim.x+threadIdx.x;if(work>=fields.brick_count*cells_per_brick)return;const int brick=work/cells_per_brick,local=work%cells_per_brick;if(fields.flags[brick]&BRICK_COVERED)return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs),a=source.level_offset+work;if(!active[a])return;const Real area=Real(source.h)*Real(source.h);Real net=Real(0);const int coordinate[3]={i,j,k};
			for(int axis=0;axis<3;++axis)
			{
				int q[3]={i,j,k},upper=-1;if(++q[axis]<bs)upper=source.level_offset+brick*cells_per_brick+local_index(bs,q[0],q[1],q[2]);else{const int neighbour=fields.neighbors[brick*6+2*axis+1];if(neighbour>=0&&!(fields.flags[neighbour]&BRICK_COVERED)){q[axis]=0;upper=source.level_offset+neighbour*cells_per_brick+local_index(bs,q[0],q[1],q[2]);}}
				if(upper>=0&&active[upper]&&!(cut_face_mask[a]&(1u<<axis)))net+=positive_face_value(fields,brick,axis,i,j,k)*area;
				q[0]=i;q[1]=j;q[2]=k;int lower=-1,lower_brick=brick;if(--q[axis]>=0)lower=source.level_offset+brick*cells_per_brick+local_index(bs,q[0],q[1],q[2]);else{lower_brick=fields.neighbors[brick*6+2*axis];if(lower_brick>=0&&!(fields.flags[lower_brick]&BRICK_COVERED)){q[axis]=bs-1;lower=source.level_offset+lower_brick*cells_per_brick+local_index(bs,q[0],q[1],q[2]);}}
				if(lower>=0&&active[lower]&&!(cut_face_mask[lower]&(1u<<axis)))net-=positive_face_value(fields,lower_brick,axis,q[0],q[1],q[2])*area;
			}
			if((fields.flags[brick]&BRICK_XMIN)&&coordinate[0]==0)net-=fields.u[fields.layout.u_index(brick,0,j,k)]*area;if((fields.flags[brick]&BRICK_XMAX)&&coordinate[0]==bs-1)net+=fields.u[fields.layout.u_index(brick,bs,j,k)]*area;integrated[a]=net;
		}
		__global__ void special_integrated_flux_kernel(const int* first,const int* second,const std::int8_t* direction,
			const Real* area,const Real* velocity,int count,const unsigned char* active,Real* integrated)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int lower=direction[edge]>0?first[edge]:second[edge],upper=direction[edge]>0?second[edge]:first[edge];if(!active[lower]||!active[upper])return;const Real flux=area[edge]*velocity[edge];atomicAdd(integrated+lower,flux);atomicAdd(integrated+upper,-flux);
		}
		__global__ void normalize_divergence_rhs_kernel(const Real* integrated,const Real* volume,const unsigned char* active,
			Real rho_over_dt,Real* divergence,Real* rhs,int n)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=n)return;if(active[q]&&volume[q]>Real(0)){divergence[q]=integrated[q]/volume[q];if(rhs)rhs[q]=-rho_over_dt*integrated[q];}else{divergence[q]=Real(0);if(rhs)rhs[q]=Real(0);}
		}
		__global__ void structured_flux_correction_kernel(const DeviceCompositeAmrFluxLevelView* levels,int level,int bs,bool outlet,
			const unsigned char* active,const std::uint8_t* cut_face_mask,const Real* pressure,Real scale)
		{
			const DeviceCompositeAmrFluxLevelView source=levels[level];const DeviceAmrFieldLevelView fields=source.fields;const int cells_per_brick=bs*bs*bs,work=blockIdx.x*blockDim.x+threadIdx.x;if(work>=fields.brick_count*cells_per_brick)return;const int brick=work/cells_per_brick,local=work%cells_per_brick;if(fields.flags[brick]&BRICK_COVERED)return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs),a=source.level_offset+work;if(!active[a])return;
			for(int axis=0;axis<3;++axis)
			{
				if(cut_face_mask[a]&(1u<<axis))continue;int q[3]={i,j,k},upper=-1,neighbour=-1;if(++q[axis]<bs)upper=source.level_offset+brick*cells_per_brick+local_index(bs,q[0],q[1],q[2]);else{neighbour=fields.neighbors[brick*6+2*axis+1];if(neighbour>=0&&!(fields.flags[neighbour]&BRICK_COVERED)){q[axis]=0;upper=source.level_offset+neighbour*cells_per_brick+local_index(bs,q[0],q[1],q[2]);}}
				if(upper>=0&&active[upper]){const Real corrected=positive_face_value(fields,brick,axis,i,j,k)-scale*(pressure[upper]-pressure[a])/Real(source.h);set_positive_face_value(fields,brick,axis,i,j,k,corrected);if(neighbour>=0){int opposite[3]={i,j,k};opposite[axis]=-1;set_positive_face_value(fields,neighbour,axis,opposite[0],opposite[1],opposite[2],corrected);}}
			}
			if(outlet&&(fields.flags[brick]&BRICK_XMAX)&&i==bs-1)fields.u[fields.layout.u_index(brick,bs,j,k)]+=scale*pressure[a]/(Real(0.5)*Real(source.h));
		}
		__global__ void special_flux_correction_kernel(const int* first,const int* second,const std::int8_t* direction,
			const Real* distance,Real* velocity,int count,const unsigned char* active,const Real* pressure,Real scale)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int lower=direction[edge]>0?first[edge]:second[edge],upper=direction[edge]>0?second[edge]:first[edge];if(active[lower]&&active[upper])velocity[edge]-=scale*(pressure[upper]-pressure[lower])/distance[edge];
		}
		__global__ void initialize_special_freestream_kernel(const std::int8_t* axis,Real* velocity,int count,Real speed){const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count)velocity[edge]=axis[edge]==0?speed:Real(0);}
		__device__ Real indexed_face_value(const DeviceAmrFieldLevelView& fields,std::int8_t axis,std::uint64_t index)
		{
			return axis==0?fields.u[index]:(axis==1?fields.v[index]:fields.w[index]);
		}
		__device__ void set_indexed_face_value(const DeviceAmrFieldLevelView& fields,std::int8_t axis,std::uint64_t index,Real value)
		{
			if(axis==0)fields.u[index]=value;else if(axis==1)fields.v[index]=value;else fields.w[index]=value;
		}
		__global__ void gather_coarse_fine_velocity_kernel(const DeviceCompositeAmrFluxLevelView* levels,const int* fine_level,const std::uint64_t* fine_index,const std::int8_t* axis,Real* velocity,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count)velocity[edge]=indexed_face_value(levels[fine_level[edge]].fields,axis[edge],fine_index[edge]);
		}
		__global__ void scatter_coarse_fine_tiles_kernel(const DeviceCompositeAmrFluxLevelView* levels,const int* fine_level,const std::uint64_t* fine_index,const std::int8_t* axis,const int* group,const Real* area,const Real* velocity,Real* group_sum,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;set_indexed_face_value(levels[fine_level[edge]].fields,axis[edge],fine_index[edge],velocity[edge]);atomicAdd(group_sum+group[edge],area[edge]*velocity[edge]);
		}
		__global__ void scatter_coarse_fine_groups_kernel(const DeviceCompositeAmrFluxLevelView* levels,const int* group_level,const std::uint64_t* group_index,const std::int8_t* group_axis,const Real* group_area,const Real* group_sum,int count)
		{
			const int group=blockIdx.x*blockDim.x+threadIdx.x;if(group<count&&group_area[group]>Real(0))set_indexed_face_value(levels[group_level[group]].fields,group_axis[group],group_index[group],group_sum[group]/group_area[group]);
		}
		__global__ void accumulate_eb_carrier_kernel(const DeviceCompositeAmrFluxLevelView* levels,const int* node,const int* level,const std::uint64_t* index,const std::int8_t* axis,const Real* area,Real* sum,Real* weight,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;const int slot=node[q]*3+axis[q];atomicAdd(sum+slot,area[q]*indexed_face_value(levels[level[q]].fields,axis[q],index[q]));atomicAdd(weight+slot,area[q]);
		}
		__global__ void accumulate_eb_aperture_kernel(const int* node_a,const int* node_b,const std::int8_t* axis,const Real* area,const Real* velocity,Real* sum,Real* weight,int special_offset,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int global=special_offset+edge,component=axis[global],a=node_a[edge]*3+component,b=node_b[edge]*3+component;const Real weighted=area[global]*velocity[global];atomicAdd(sum+a,weighted);atomicAdd(sum+b,weighted);atomicAdd(weight+a,area[global]);atomicAdd(weight+b,area[global]);
		}
		__global__ void transport_eb_aperture_kernel(const int* node_a,const int* node_b,const int* negative_edge,const int* positive_edge,const std::int8_t* direction,const std::int8_t* axis,const Real* transport_length,const Real* velocity,const Real* sum,const Real* weight,Real dt,Real molecular_nu,Real smagorinsky_cs,Real* output,int special_offset,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int global=special_offset+edge,component=axis[global],first=node_a[edge],second=node_b[edge],lower=direction[global]>0?first:second,upper=direction[global]>0?second:first;const Real current=velocity[global],lower_weight=weight[lower*3+component],upper_weight=weight[upper*3+component],lower_value=lower_weight>Real(0)?sum[lower*3+component]/lower_weight:current,upper_value=upper_weight>Real(0)?sum[upper*3+component]/upper_weight:current,d=max(transport_length[edge],Real(1e-12));const Real courant=min(Real(1),abs(current)*dt/d),normal_jump=abs(upper_value-lower_value),eddy_nu=smagorinsky_cs*smagorinsky_cs*d*normal_jump,diffusion=min(Real(0.25),(molecular_nu+eddy_nu)*dt/(d*d));const bool positive=current>=Real(0);const int upstream_edge=positive?negative_edge[edge]:positive_edge[edge],downstream_edge=positive?positive_edge[edge]:negative_edge[edge],upstream2_edge=upstream_edge>=0?(positive?negative_edge[upstream_edge]:positive_edge[upstream_edge]):-1;const bool complete=upstream_edge>=0&&upstream2_edge>=0&&downstream_edge>=0;const Real upstream=complete?velocity[special_offset+upstream_edge]:current,upstream2=complete?velocity[special_offset+upstream2_edge]:current,downstream=complete?velocity[special_offset+downstream_edge]:current;output[edge]=detail::bounded_compact_transport_update(current,lower_value,upper_value,upstream,upstream2,downstream,courant,diffusion,complete);
		}
		__global__ void commit_eb_aperture_kernel(Real* velocity,const Real* source,int special_offset,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count)velocity[special_offset+edge]=source[edge];
		}
		__global__ void embedded_cfl_rate_kernel(const Real* velocity,const Real* transport_length,
			Real* rate,int special_offset,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge<count)rate[edge]=abs(velocity[special_offset+edge])/max(transport_length[edge],Real(1e-12));
		}
		__global__ void scatter_regular_pressure_kernel(const DeviceCompositeAmrFluxLevelView* levels,int level,int bs,const unsigned char* active,const Real* pressure)
		{
			const DeviceCompositeAmrFluxLevelView source=levels[level];const DeviceAmrFieldLevelView fields=source.fields;const int cells_per_brick=bs*bs*bs,work=blockIdx.x*blockDim.x+threadIdx.x;if(work>=fields.brick_count*cells_per_brick)return;const int brick=work/cells_per_brick,local=work%cells_per_brick,a=source.level_offset+work;if(!active[a]||(fields.flags[brick]&BRICK_COVERED))return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs);fields.p[fields.layout.cell_index(brick,i,j,k)]=pressure[a];
		}
		struct DotActive{const Real* a;const Real* b;const unsigned char* active;__host__ __device__ double operator()(int q)const{return active[q]?static_cast<double>(a[q])*static_cast<double>(b[q]):0.0;}};
		double dot_active(const Real* a,const Real* b,const unsigned char* active,int n){auto begin=thrust::counting_iterator<int>(0);return thrust::transform_reduce(thrust::device,begin,begin+n,DotActive{a,b,active},0.0,thrust::plus<double>());}
	}

	DeviceCompositeAmrPressureOperator::DeviceCompositeAmrPressureOperator(const CompositeAmrPressureSystem& system)
	{
		if (!system.hierarchy) throw std::invalid_argument("composite AMR pressure system has no hierarchy");
		storage_size_ = system.storage_size; brick_size_ = system.brick_size; outlet_ = system.pressure_outlet_xmax; level_count_ = static_cast<int>(system.hierarchy->levels().size()); allocations_.resize(level_count_); brick_counts_.resize(level_count_); std::vector<DeviceCompositeAmrLevelView> views(level_count_);
		for (int level = 0; level < level_count_; ++level)
		{
			const AmrLevel& source = system.hierarchy->levels()[level]; Allocation& allocation = allocations_[level]; brick_counts_[level] = static_cast<int>(source.bricks.size()); std::vector<int> neighbors(source.bricks.size() * 6); std::vector<std::uint32_t> flags(source.bricks.size());
			for (int brick = 0; brick < static_cast<int>(source.bricks.size()); ++brick) { flags[brick] = source.bricks[brick].flags; for (int face=0;face<6;++face) neighbors[brick*6+face] = source.bricks[brick].same_level_neighbor[face]; }
			if (!neighbors.empty()) { check(cudaMalloc(&allocation.neighbors, neighbors.size()*sizeof(int)), "cudaMalloc composite AMR neighbors"); check(cudaMemcpy(allocation.neighbors, neighbors.data(), neighbors.size()*sizeof(int), cudaMemcpyHostToDevice), "upload composite AMR neighbors"); check(cudaMalloc(&allocation.flags, flags.size()*sizeof(std::uint32_t)), "cudaMalloc composite AMR flags"); check(cudaMemcpy(allocation.flags, flags.data(), flags.size()*sizeof(std::uint32_t), cudaMemcpyHostToDevice), "upload composite AMR flags"); }
			views[level] = {system.level_offset[level], static_cast<int>(source.bricks.size()), source.h, allocation.neighbors, allocation.flags}; bytes_ += neighbors.size()*sizeof(int)+flags.size()*sizeof(std::uint32_t);
		}
		if (!views.empty()) { check(cudaMalloc(&levels_, views.size()*sizeof(DeviceCompositeAmrLevelView)), "cudaMalloc composite AMR views"); check(cudaMemcpy(levels_, views.data(), views.size()*sizeof(DeviceCompositeAmrLevelView), cudaMemcpyHostToDevice), "upload composite AMR views"); bytes_ += views.size()*sizeof(DeviceCompositeAmrLevelView); }
		active_=upload(system.active,"upload composite AMR operator active mask");cut_face_mask_=upload(system.cut_face_mask,"upload composite AMR cut-face mask");bytes_+=system.active.size()+system.cut_face_mask.size()*sizeof(std::uint8_t);std::vector<CoarseFinePressureConnection> special=system.coarse_fine;special.insert(special.end(),system.embedded.begin(),system.embedded.end());connection_count_ = static_cast<int>(special.size()); if (connection_count_)
		{
			std::vector<int> coarse(connection_count_), fine(connection_count_); std::vector<Real> coefficient(connection_count_); for(int edge=0;edge<connection_count_;++edge){coarse[edge]=special[edge].coarse_dof;fine[edge]=special[edge].fine_dof;coefficient[edge]=static_cast<Real>(special[edge].open_area/special[edge].centre_distance);}
			check(cudaMalloc(&coarse_dof_,coarse.size()*sizeof(int)),"cudaMalloc coarse/fine coarse DOFs");check(cudaMalloc(&fine_dof_,fine.size()*sizeof(int)),"cudaMalloc coarse/fine fine DOFs");check(cudaMalloc(&coefficient_,coefficient.size()*sizeof(Real)),"cudaMalloc coarse/fine coefficients");check(cudaMemcpy(coarse_dof_,coarse.data(),coarse.size()*sizeof(int),cudaMemcpyHostToDevice),"upload coarse DOFs");check(cudaMemcpy(fine_dof_,fine.data(),fine.size()*sizeof(int),cudaMemcpyHostToDevice),"upload fine DOFs");check(cudaMemcpy(coefficient_,coefficient.data(),coefficient.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload coarse/fine coefficients");bytes_+=coarse.size()*sizeof(int)+fine.size()*sizeof(int)+coefficient.size()*sizeof(Real);
		}
		gauge_count_=static_cast<int>(system.gauges.size());if(gauge_count_){std::vector<int> dof(gauge_count_);std::vector<Real> coefficient(gauge_count_);for(int q=0;q<gauge_count_;++q){dof[q]=system.gauges[q].dof;coefficient[q]=static_cast<Real>(system.gauges[q].coefficient);}gauge_dof_=upload(dof,"upload composite pressure gauge DOFs");gauge_coefficient_=upload(coefficient,"upload composite pressure gauge coefficients");bytes_+=dof.size()*sizeof(int)+coefficient.size()*sizeof(Real);}
	}

	DeviceCompositeAmrPressureOperator::~DeviceCompositeAmrPressureOperator()
	{
		for (Allocation& allocation : allocations_) { if(allocation.neighbors)cudaFree(allocation.neighbors);if(allocation.flags)cudaFree(allocation.flags); } if(levels_)cudaFree(levels_);if(coarse_dof_)cudaFree(coarse_dof_);if(fine_dof_)cudaFree(fine_dof_);if(coefficient_)cudaFree(coefficient_);if(gauge_dof_)cudaFree(gauge_dof_);if(gauge_coefficient_)cudaFree(gauge_coefficient_);if(active_)cudaFree(active_);if(cut_face_mask_)cudaFree(cut_face_mask_);
	}

	void DeviceCompositeAmrPressureOperator::apply(const Real* pressure, Real* output) const
	{
		check(cudaMemset(output,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite AMR output");const int cells_per_brick=brick_size_*brick_size_*brick_size_;
		for(int level=0;level<level_count_;++level){const int work=brick_counts_[level]*cells_per_brick;if(work)structured_apply_kernel<<<(work+255)/256,256>>>(levels_,level,brick_size_,outlet_,active_,cut_face_mask_,pressure,output);}if(connection_count_)coarse_fine_apply_kernel<<<(connection_count_+255)/256,256>>>(coarse_dof_,fine_dof_,coefficient_,connection_count_,active_,pressure,output);if(gauge_count_)gauge_apply_kernel<<<(gauge_count_+255)/256,256>>>(gauge_dof_,gauge_coefficient_,gauge_count_,pressure,output);check(cudaDeviceSynchronize(),"apply composite AMR pressure operator");
	}

	DeviceCompositeAmrPressureSolver::DeviceCompositeAmrPressureSolver(const CompositeAmrPressureSystem& system):op_(system),n_(system.storage_size)
	{
		if(!system.hierarchy||system.hierarchy->levels().empty()||system.preconditioner_aggregate.size()!=static_cast<std::size_t>(n_))throw std::invalid_argument("composite AMR solver geometric aggregation metadata");std::vector<double> diagonal_double;system.diagonal_cpu(diagonal_double);std::vector<Real> diagonal(diagonal_double.size());for(std::size_t q=0;q<diagonal.size();++q)diagonal[q]=static_cast<Real>(diagonal_double[q]);
		r_=allocate_zero<Real>(n_,"allocate AMR residual");z_=allocate_zero<Real>(n_,"allocate AMR preconditioned residual");direction_=allocate_zero<Real>(n_,"allocate AMR direction");Ad_=allocate_zero<Real>(n_,"allocate AMR operator temporary");diagonal_=upload(diagonal,"upload AMR diagonal");active_=upload(system.active,"upload AMR active mask");aggregate_=upload(system.preconditioner_aggregate,"upload AMR coarse aggregate map");
		const AmrLevel& base=system.hierarchy->levels().front();brick_size_=system.brick_size;base_offset_=system.level_offset.front();base_bricks_=static_cast<int>(base.bricks.size());base_cells_=base_bricks_*brick_size_*brick_size_*brick_size_;const int cells_per_brick=brick_size_*brick_size_*brick_size_;std::vector<int> neighbors(static_cast<std::size_t>(base_bricks_)*6,-1);for(int brick=0;brick<base_bricks_;++brick)for(int face=0;face<6;++face)neighbors[brick*6+face]=base.bricks[brick].same_level_neighbor[face];std::vector<double> coarse_positive_double(static_cast<std::size_t>(base_cells_)*3,0),coarse_diagonal_double(base_cells_,0);
		std::vector<int> base_extra_a,base_extra_b;std::vector<double> base_extra_coefficient_double;auto global_coord=[&](int coarse){const int brick=coarse/cells_per_brick,local=coarse%cells_per_brick;return std::array<int,3>{base.bricks[brick].coord.x*brick_size_+local%brick_size_,base.bricks[brick].coord.y*brick_size_+(local/brick_size_)%brick_size_,base.bricks[brick].coord.z*brick_size_+local/(brick_size_*brick_size_)};};auto add_coarse_edge=[&](int a,int b,double coefficient){if(a<0||b<0||a==b||!system.active[a]||!system.active[b])return;const int ca=system.preconditioner_aggregate[a]-base_offset_,cb=system.preconditioner_aggregate[b]-base_offset_;if(ca<0||cb<0||ca>=base_cells_||cb>=base_cells_||ca==cb)return;const auto ga=global_coord(ca),gb=global_coord(cb);int axis=-1,manhattan=0;for(int d=0;d<3;++d){const int delta=gb[d]-ga[d];manhattan+=std::abs(delta);if(delta)axis=d;}if(manhattan==1&&axis>=0){const int lower=ga[axis]<gb[axis]?ca:cb;coarse_positive_double[static_cast<std::size_t>(lower)*3+axis]+=coefficient;}else{base_extra_a.push_back(ca);base_extra_b.push_back(cb);base_extra_coefficient_double.push_back(coefficient);}coarse_diagonal_double[ca]+=coefficient;coarse_diagonal_double[cb]+=coefficient;};
		for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level){const AmrLevel& source=system.hierarchy->levels()[level];const double coefficient=source.h;for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick){const BrickMetadata& meta=source.bricks[brick];if(!meta.active())continue;for(int k=0;k<brick_size_;++k)for(int j=0;j<brick_size_;++j)for(int i=0;i<brick_size_;++i){const int a=system.dof(level,brick,i,j,k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if(system.cut_face_mask[a]&(1u<<axis))continue;int c[3]={i,j,k},b=-1;if(++c[axis]<brick_size_)b=system.dof(level,brick,c[0],c[1],c[2]);else{const int neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&source.bricks[neighbour].active()){c[axis]=0;b=system.dof(level,neighbour,c[0],c[1],c[2]);}}add_coarse_edge(a,b,coefficient);}if(system.pressure_outlet_xmax&&(meta.flags&BRICK_XMAX)&&i==brick_size_-1){const int coarse=system.preconditioner_aggregate[a]-base_offset_;if(coarse>=0&&coarse<base_cells_)coarse_diagonal_double[coarse]+=2.0*source.h;}}}}
		for(const CoarseFinePressureConnection& edge:system.coarse_fine)add_coarse_edge(edge.coarse_dof,edge.fine_dof,edge.open_area/edge.centre_distance);for(const CoarseFinePressureConnection& edge:system.embedded)add_coarse_edge(edge.coarse_dof,edge.fine_dof,edge.open_area/edge.centre_distance);for(const CompositePressureGauge& gauge:system.gauges){const int coarse=system.preconditioner_aggregate[gauge.dof]-base_offset_;if(coarse>=0&&coarse<base_cells_)coarse_diagonal_double[coarse]+=gauge.coefficient;}
		std::vector<Real> base_positive(coarse_positive_double.size()),base_diagonal(base_cells_),base_extra_coefficient(base_extra_coefficient_double.size());for(std::size_t q=0;q<base_positive.size();++q)base_positive[q]=static_cast<Real>(coarse_positive_double[q]);for(int q=0;q<base_cells_;++q)base_diagonal[q]=static_cast<Real>(coarse_diagonal_double[q]);for(std::size_t q=0;q<base_extra_coefficient.size();++q)base_extra_coefficient[q]=static_cast<Real>(base_extra_coefficient_double[q]);base_extra_count_=static_cast<int>(base_extra_a.size());base_neighbors_=upload(neighbors,"upload AMR base neighbors");base_positive_=upload(base_positive,"upload AMR Galerkin base faces");base_diagonal_=upload(base_diagonal,"upload AMR base diagonal");base_extra_a_=upload(base_extra_a,"upload AMR nonlocal aggregate edge A");base_extra_b_=upload(base_extra_b,"upload AMR nonlocal aggregate edge B");base_extra_coefficient_=upload(base_extra_coefficient,"upload AMR nonlocal aggregate coefficient");base_rhs_=allocate_zero<Real>(base_cells_,"allocate AMR base RHS");base_x_=allocate_zero<Real>(base_cells_,"allocate AMR base correction");base_tmp_=allocate_zero<Real>(base_cells_,"allocate AMR base temporary");bytes_=op_.bytes()+static_cast<std::size_t>(5)*n_*sizeof(Real)+system.active.size()+system.preconditioner_aggregate.size()*sizeof(int)+neighbors.size()*sizeof(int)+base_positive.size()*sizeof(Real)+static_cast<std::size_t>(4)*base_cells_*sizeof(Real)+base_extra_a.size()*static_cast<std::size_t>(2*sizeof(int)+sizeof(Real));
	}
	DeviceCompositeAmrPressureSolver::~DeviceCompositeAmrPressureSolver(){for(void* pointer:{(void*)r_,(void*)z_,(void*)direction_,(void*)Ad_,(void*)diagonal_,(void*)active_,(void*)aggregate_,(void*)base_neighbors_,(void*)base_extra_a_,(void*)base_extra_b_,(void*)base_positive_,(void*)base_diagonal_,(void*)base_extra_coefficient_,(void*)base_rhs_,(void*)base_x_,(void*)base_tmp_})if(pointer)cudaFree(pointer);}
	void DeviceCompositeAmrPressureSolver::apply_preconditioner(const Real* residual,Real* output)
	{
		precondition_kernel<<<(n_+255)/256,256>>>(residual,output,diagonal_,active_,n_);check(cudaMemset(base_rhs_,0,static_cast<std::size_t>(base_cells_)*sizeof(Real)),"clear AMR base RHS");check(cudaMemset(base_x_,0,static_cast<std::size_t>(base_cells_)*sizeof(Real)),"clear AMR base correction");restrict_aggregate_kernel<<<(n_+255)/256,256>>>(residual,active_,aggregate_,n_,base_offset_,base_cells_,base_rhs_);Real* current=base_x_;Real* next=base_tmp_;constexpr int sweeps=16;for(int sweep=0;sweep<sweeps;++sweep){coarse_jacobi_kernel<<<(base_cells_+255)/256,256>>>(base_rhs_,current,next,base_positive_,base_diagonal_,base_neighbors_,base_bricks_,brick_size_,Real(0.7));if(base_extra_count_)coarse_extra_jacobi_kernel<<<(base_extra_count_+255)/256,256>>>(base_extra_a_,base_extra_b_,base_extra_coefficient_,base_extra_count_,current,base_diagonal_,Real(0.7),next);Real* swap=current;current=next;next=swap;}prolong_aggregate_kernel<<<(n_+255)/256,256>>>(aggregate_,active_,n_,base_offset_,base_cells_,current,output);
	}
	AmrGpuSolveResult DeviceCompositeAmrPressureSolver::solve(Real* pressure,const Real* rhs,double tolerance,int max_iterations,bool warm_start)
	{
		AmrGpuSolveResult result;if(!warm_start)check(cudaMemset(pressure,0,static_cast<std::size_t>(n_)*sizeof(Real)),"clear composite AMR pressure");if(warm_start){op_.apply(pressure,Ad_);check(cudaMemcpy(r_,rhs,static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),"copy composite AMR rhs");subtract_kernel<<<(n_+255)/256,256>>>(r_,Ad_,active_,n_);}else check(cudaMemcpy(r_,rhs,static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),"initial composite AMR residual");
		const double rhs2=dot_active(rhs,rhs,active_,n_);if(!(rhs2>0)){result.converged=true;return result;}double residual2=dot_active(r_,r_,active_,n_);result.relative_residual=std::sqrt(residual2/rhs2);if(result.relative_residual<=tolerance){result.converged=true;return result;}apply_preconditioner(r_,z_);copy_active_kernel<<<(n_+255)/256,256>>>(z_,direction_,active_,n_);double rz=dot_active(r_,z_,active_,n_);
		for(int iteration=0;iteration<max_iterations;++iteration){op_.apply(direction_,Ad_);const double dAd=dot_active(direction_,Ad_,active_,n_);if(!(dAd>0))break;const Real alpha=static_cast<Real>(rz/dAd);update_pressure_residual_kernel<<<(n_+255)/256,256>>>(pressure,r_,direction_,Ad_,alpha,active_,n_);residual2=dot_active(r_,r_,active_,n_);result.iterations=iteration+1;result.relative_residual=std::sqrt(residual2/rhs2);if(result.relative_residual<=tolerance){result.converged=true;break;}apply_preconditioner(r_,z_);const double next_rz=dot_active(r_,z_,active_,n_);const Real beta=static_cast<Real>(next_rz/rz);update_direction_kernel<<<(n_+255)/256,256>>>(z_,direction_,beta,active_,n_);rz=next_rz;}check(cudaDeviceSynchronize(),"solve composite AMR pressure");return result;
	}

	DeviceCompositeAmrProjection::DeviceCompositeAmrProjection(const CompositeAmrPressureSystem& system,DeviceAmrFields& fields)
		:solver_(system),fields_(&fields)
	{
		if(!system.hierarchy)throw std::invalid_argument("composite AMR projection has no hierarchy");level_count_=static_cast<int>(system.hierarchy->levels().size());if(fields.level_count()!=level_count_)throw std::invalid_argument("composite AMR projection field hierarchy mismatch");storage_size_=system.storage_size;brick_size_=system.brick_size;outlet_=system.pressure_outlet_xmax;coarse_fine_count_=static_cast<int>(system.coarse_fine.size());embedded_count_=static_cast<int>(system.embedded.size());special_count_=coarse_fine_count_+embedded_count_;
		std::vector<DeviceCompositeAmrFluxLevelView> views(level_count_);brick_counts_.resize(level_count_);for(int level=0;level<level_count_;++level){const DeviceAmrFieldLevelView field_view=fields.level_view(level);views[level]={system.level_offset[level],system.hierarchy->levels()[level].h,field_view};brick_counts_[level]=field_view.brick_count;}levels_=upload(views,"upload composite AMR flux level views");
		active_=upload(system.active,"upload composite AMR projection active mask");cut_face_mask_=upload(system.cut_face_mask,"upload composite AMR projection cut-face mask");std::vector<Real> volume(system.volume.size());for(std::size_t q=0;q<volume.size();++q)volume[q]=static_cast<Real>(system.volume[q]);volume_=upload(volume,"upload composite AMR projection volumes");
		std::vector<CoarseFinePressureConnection> special=system.coarse_fine;special.insert(special.end(),system.embedded.begin(),system.embedded.end());std::vector<int> first(special_count_),second(special_count_);std::vector<std::int8_t> direction(special_count_),axis(special_count_);std::vector<Real> area(special_count_),distance(special_count_);for(int edge=0;edge<special_count_;++edge){first[edge]=special[edge].coarse_dof;second[edge]=special[edge].fine_dof;direction[edge]=special[edge].direction;axis[edge]=special[edge].axis;area[edge]=static_cast<Real>(special[edge].open_area);distance[edge]=static_cast<Real>(special[edge].centre_distance);}first_dof_=upload(first,"upload composite flux first DOFs");second_dof_=upload(second,"upload composite flux second DOFs");direction_=upload(direction,"upload composite flux direction");axis_=upload(axis,"upload composite flux axes");open_area_=upload(area,"upload composite flux areas");centre_distance_=upload(distance,"upload composite flux distances");special_velocity_=allocate_zero<Real>(special_count_,"allocate composite special velocity");max_abs_scratch_=allocate_zero<Real>(1,"allocate special-velocity maximum");
		// A coarse/fine connection is collocated with one unique fine MAC face. Keep
		// that compact flux state synchronized with transport, then scatter the four
		// corrected fine tiles back to their fine faces and their area-mean coarse face.
		if(coarse_fine_count_)
		{
			std::vector<int> fine_level(coarse_fine_count_),group(coarse_fine_count_),group_level;std::vector<std::uint64_t> fine_index(coarse_fine_count_),group_index;std::vector<std::int8_t> group_axis;std::vector<Real> group_area;std::unordered_map<std::uint64_t,int> groups;
			auto face_address=[&](int dof,int face_axis,bool positive_face)->std::pair<int,std::uint64_t>
			{
				const int cells_per_brick=brick_size_*brick_size_*brick_size_;for(int level=0;level<level_count_;++level){const int begin=system.level_offset[level],end=begin+brick_counts_[level]*cells_per_brick;if(dof<begin||dof>=end)continue;const int work=dof-begin,brick=work/cells_per_brick,local=work%cells_per_brick;int i=local%brick_size_,j=(local/brick_size_)%brick_size_,k=local/(brick_size_*brick_size_);if(face_axis==0)i+=positive_face;else if(face_axis==1)j+=positive_face;else k+=positive_face;const BrickFieldLayout layout=views[level].fields.layout;const std::uint64_t index=face_axis==0?layout.u_index(brick,i,j,k):(face_axis==1?layout.v_index(brick,i,j,k):layout.w_index(brick,i,j,k));return {level,index};}throw std::runtime_error("coarse/fine pressure DOF has no regular MAC field address");
			};
			for(int edge=0;edge<coarse_fine_count_;++edge)
			{
				const CoarseFinePressureConnection& connection=system.coarse_fine[edge];const auto fine=face_address(connection.fine_dof,connection.axis,connection.direction<0);fine_level[edge]=fine.first;fine_index[edge]=fine.second;const std::uint64_t key=(static_cast<std::uint64_t>(static_cast<std::uint32_t>(connection.coarse_dof))<<3)|(static_cast<std::uint64_t>(connection.axis)<<1)|(connection.direction>0);auto found=groups.find(key);int id;if(found==groups.end()){id=static_cast<int>(group_level.size());groups.emplace(key,id);const auto coarse=face_address(connection.coarse_dof,connection.axis,connection.direction>0);group_level.push_back(coarse.first);group_index.push_back(coarse.second);group_axis.push_back(connection.axis);group_area.push_back(Real(0));}else id=found->second;group[edge]=id;group_area[id]+=static_cast<Real>(connection.open_area);
			}
			coarse_fine_group_count_=static_cast<int>(group_level.size());cf_fine_level_=upload(fine_level,"upload coarse/fine fine levels");cf_fine_index_=upload(fine_index,"upload coarse/fine fine face indices");cf_group_=upload(group,"upload coarse/fine group IDs");cf_group_level_=upload(group_level,"upload coarse/fine coarse levels");cf_group_index_=upload(group_index,"upload coarse/fine coarse face indices");cf_group_axis_=upload(group_axis,"upload coarse/fine group axes");cf_group_area_=upload(group_area,"upload coarse/fine group areas");cf_group_sum_=allocate_zero<Real>(coarse_fine_group_count_,"allocate coarse/fine group flux sums");bytes_+=fine_level.size()*sizeof(int)+fine_index.size()*sizeof(std::uint64_t)+group.size()*sizeof(int)+group_level.size()*sizeof(int)+group_index.size()*sizeof(std::uint64_t)+group_axis.size()*sizeof(std::int8_t)+2*group_area.size()*sizeof(Real);
		}
		if(embedded_count_)
		{
			std::unordered_map<int,int> node_map;std::vector<int> node_dof,node_a(embedded_count_),node_b(embedded_count_);std::vector<Real> transport_length(embedded_count_);auto node=[&](int dof){auto found=node_map.find(dof);if(found!=node_map.end())return found->second;const int id=static_cast<int>(node_dof.size());node_map.emplace(dof,id);node_dof.push_back(dof);return id;};for(int edge=0;edge<embedded_count_;++edge){const auto& connection=system.embedded[edge];node_a[edge]=node(connection.coarse_dof);node_b[edge]=node(connection.fine_dof);const double a_scale=std::cbrt(std::max(0.0,system.volume[connection.coarse_dof])),b_scale=std::cbrt(std::max(0.0,system.volume[connection.fine_dof]));transport_length[edge]=static_cast<Real>(std::max(connection.centre_distance,0.5*(a_scale+b_scale)));}embedded_node_count_=static_cast<int>(node_dof.size());
			// Build a compact axis-aligned aperture chain at each graph node. A candidate must
			// extend through the shared fluid node in the required direction and use the same
			// velocity component; nearest face centroids resolve junction branches deterministically.
			std::vector<std::vector<int>> incident(embedded_node_count_);for(int edge=0;edge<embedded_count_;++edge){incident[node_a[edge]].push_back(edge);incident[node_b[edge]].push_back(edge);}std::vector<int> negative_edge(embedded_count_,-1),positive_edge(embedded_count_,-1);auto lower_node=[&](int edge){return system.embedded[edge].direction>0?node_a[edge]:node_b[edge];};auto upper_node=[&](int edge){return system.embedded[edge].direction>0?node_b[edge]:node_a[edge];};auto coordinate=[](Vec3d point,int component){return component==0?point.x:(component==1?point.y:point.z);};auto nearest_extension=[&](int edge,int shared,bool negative){int best=-1;double best_distance=std::numeric_limits<double>::infinity();const auto& source=system.embedded[edge];for(int candidate:incident[shared]){if(candidate==edge||system.embedded[candidate].axis!=source.axis)continue;if(negative?upper_node(candidate)!=shared:lower_node(candidate)!=shared)continue;const double delta=coordinate(system.embedded[candidate].face_centroid,source.axis)-coordinate(source.face_centroid,source.axis);if(negative?delta>=-1e-12:delta<=1e-12)continue;const double distance=length2(system.embedded[candidate].face_centroid-source.face_centroid);if(distance<best_distance){best_distance=distance;best=candidate;}}return best;};for(int edge=0;edge<embedded_count_;++edge){negative_edge[edge]=nearest_extension(edge,lower_node(edge),true);positive_edge[edge]=nearest_extension(edge,upper_node(edge),false);}for(int edge=0;edge<embedded_count_;++edge)if(negative_edge[edge]>=0&&positive_edge[edge]>=0&&negative_edge[negative_edge[edge]]>=0&&positive_edge[positive_edge[edge]]>=0)++embedded_high_order_stencil_count_;
			std::vector<int> carrier_node,carrier_level;std::vector<std::uint64_t> carrier_index;std::vector<std::int8_t> carrier_axis;std::vector<Real> carrier_area;const int cells_per_brick=brick_size_*brick_size_*brick_size_;
			for(int node_id=0;node_id<embedded_node_count_;++node_id)
			{
				const int dof=node_dof[node_id];int level=-1,brick=-1,local=-1;for(int candidate=0;candidate<level_count_;++candidate){const int begin=system.level_offset[candidate],end=begin+brick_counts_[candidate]*cells_per_brick;if(dof>=begin&&dof<end){level=candidate;const int work=dof-begin;brick=work/cells_per_brick;local=work%cells_per_brick;break;}}if(level<0)continue;const AmrLevel& metadata=system.hierarchy->levels()[level];const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;const int i=local%brick_size_,j=(local/brick_size_)%brick_size_,k=local/(brick_size_*brick_size_);const BrickFieldLayout layout=views[level].fields.layout;
				for(int face_axis=0;face_axis<3;++face_axis)for(int sign=-1;sign<=1;sign+=2)
				{
					int c[3]={i,j,k},other=-1,other_brick=brick;c[face_axis]+=sign;if(c[face_axis]>=0&&c[face_axis]<brick_size_)other=system.dof(level,brick,c[0],c[1],c[2]);else{other_brick=record.same_level_neighbor[2*face_axis+(sign>0)];if(other_brick>=0&&metadata.bricks[other_brick].active()){c[face_axis]=sign>0?0:brick_size_-1;other=system.dof(level,other_brick,c[0],c[1],c[2]);}}if(other<0||!system.active[other])continue;const int lower=sign>0?dof:other;if(system.cut_face_mask[lower]&(1u<<face_axis))continue;int fi=i,fj=j,fk=k;if(face_axis==0)fi+=sign>0;else if(face_axis==1)fj+=sign>0;else fk+=sign>0;const std::uint64_t index=face_axis==0?layout.u_index(brick,fi,fj,fk):(face_axis==1?layout.v_index(brick,fi,fj,fk):layout.w_index(brick,fi,fj,fk));carrier_node.push_back(node_id);carrier_level.push_back(level);carrier_index.push_back(index);carrier_axis.push_back(static_cast<std::int8_t>(face_axis));carrier_area.push_back(static_cast<Real>(metadata.h*metadata.h));
				}
			}
			embedded_carrier_count_=static_cast<int>(carrier_node.size());eb_node_a_=upload(node_a,"upload EB transport first nodes");eb_node_b_=upload(node_b,"upload EB transport second nodes");eb_negative_edge_=upload(negative_edge,"upload EB negative graph neighbours");eb_positive_edge_=upload(positive_edge,"upload EB positive graph neighbours");eb_carrier_node_=upload(carrier_node,"upload EB transport carrier nodes");eb_carrier_level_=upload(carrier_level,"upload EB transport carrier levels");eb_carrier_index_=upload(carrier_index,"upload EB transport carrier indices");eb_carrier_axis_=upload(carrier_axis,"upload EB transport carrier axes");eb_carrier_area_=upload(carrier_area,"upload EB transport carrier areas");eb_node_axis_sum_=allocate_zero<Real>(static_cast<std::size_t>(embedded_node_count_)*3,"allocate EB node component sums");eb_node_axis_weight_=allocate_zero<Real>(static_cast<std::size_t>(embedded_node_count_)*3,"allocate EB node component weights");eb_transport_length_=upload(transport_length,"upload EB aperture transport lengths");eb_transport_scratch_=allocate_zero<Real>(embedded_count_,"allocate EB aperture transport scratch");bytes_+=4*node_a.size()*sizeof(int)+2*carrier_node.size()*sizeof(int)+carrier_index.size()*sizeof(std::uint64_t)+carrier_axis.size()*sizeof(std::int8_t)+carrier_area.size()*sizeof(Real)+static_cast<std::size_t>(embedded_node_count_)*6*sizeof(Real)+static_cast<std::size_t>(embedded_count_)*2*sizeof(Real);
		}
		integrated_=allocate_zero<Real>(storage_size_,"allocate composite integrated flux");divergence_=allocate_zero<Real>(storage_size_,"allocate composite divergence");rhs_=allocate_zero<Real>(storage_size_,"allocate composite RHS");pressure_=allocate_zero<Real>(storage_size_,"allocate composite pressure");
		bytes_+=solver_.bytes()+views.size()*sizeof(DeviceCompositeAmrFluxLevelView)+system.active.size()+system.cut_face_mask.size()*sizeof(std::uint8_t)+volume.size()*sizeof(Real)+static_cast<std::size_t>(special_count_)*(2*sizeof(int)+2*sizeof(std::int8_t)+3*sizeof(Real))+static_cast<std::size_t>(4)*storage_size_*sizeof(Real)+sizeof(Real);
	}

	DeviceCompositeAmrProjection::~DeviceCompositeAmrProjection()
	{
		for(void* pointer:{(void*)levels_,(void*)active_,(void*)cut_face_mask_,(void*)volume_,(void*)integrated_,(void*)divergence_,(void*)rhs_,(void*)pressure_,(void*)first_dof_,(void*)second_dof_,(void*)direction_,(void*)axis_,(void*)open_area_,(void*)centre_distance_,(void*)special_velocity_,(void*)max_abs_scratch_,(void*)cf_fine_level_,(void*)cf_fine_index_,(void*)cf_group_,(void*)cf_group_level_,(void*)cf_group_index_,(void*)cf_group_axis_,(void*)cf_group_area_,(void*)cf_group_sum_,(void*)eb_node_a_,(void*)eb_node_b_,(void*)eb_negative_edge_,(void*)eb_positive_edge_,(void*)eb_carrier_node_,(void*)eb_carrier_level_,(void*)eb_carrier_index_,(void*)eb_carrier_axis_,(void*)eb_carrier_area_,(void*)eb_node_axis_sum_,(void*)eb_node_axis_weight_,(void*)eb_transport_length_,(void*)eb_transport_scratch_})if(pointer)cudaFree(pointer);
	}

	void DeviceCompositeAmrProjection::clear_special_fluxes(){if(special_count_)check(cudaMemset(special_velocity_,0,static_cast<std::size_t>(special_count_)*sizeof(Real)),"clear composite special fluxes");}
	void DeviceCompositeAmrProjection::initialize_special_freestream(Real speed){if(special_count_)initialize_special_freestream_kernel<<<(special_count_+255)/256,256>>>(axis_,special_velocity_,special_count_,speed);check(cudaDeviceSynchronize(),"initialize composite special freestream");}
	void DeviceCompositeAmrProjection::sync_coarse_fine_from_fields()
	{
		if(coarse_fine_count_)gather_coarse_fine_velocity_kernel<<<(coarse_fine_count_+255)/256,256>>>(levels_,cf_fine_level_,cf_fine_index_,axis_,special_velocity_,coarse_fine_count_);check(cudaDeviceSynchronize(),"gather transported coarse/fine velocities");
	}
	void DeviceCompositeAmrProjection::transport_embedded_apertures(Real dt,Real molecular_nu,Real smagorinsky_cs)
	{
		if(!(dt>Real(0))||molecular_nu<Real(0)||smagorinsky_cs<Real(0))throw std::invalid_argument("embedded aperture transport timestep/viscosity");if(!embedded_count_)return;const std::size_t node_components=static_cast<std::size_t>(embedded_node_count_)*3;check(cudaMemset(eb_node_axis_sum_,0,node_components*sizeof(Real)),"clear EB transport component sums");check(cudaMemset(eb_node_axis_weight_,0,node_components*sizeof(Real)),"clear EB transport component weights");if(embedded_carrier_count_)accumulate_eb_carrier_kernel<<<(embedded_carrier_count_+255)/256,256>>>(levels_,eb_carrier_node_,eb_carrier_level_,eb_carrier_index_,eb_carrier_axis_,eb_carrier_area_,eb_node_axis_sum_,eb_node_axis_weight_,embedded_carrier_count_);accumulate_eb_aperture_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,axis_,open_area_,special_velocity_,eb_node_axis_sum_,eb_node_axis_weight_,coarse_fine_count_,embedded_count_);transport_eb_aperture_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,eb_negative_edge_,eb_positive_edge_,direction_,axis_,eb_transport_length_,special_velocity_,eb_node_axis_sum_,eb_node_axis_weight_,dt,molecular_nu,smagorinsky_cs,eb_transport_scratch_,coarse_fine_count_,embedded_count_);commit_eb_aperture_kernel<<<(embedded_count_+255)/256,256>>>(special_velocity_,eb_transport_scratch_,coarse_fine_count_,embedded_count_);check(cudaDeviceSynchronize(),"transport compact EB aperture velocities");
	}
	void DeviceCompositeAmrProjection::upload_special_fluxes(const CompositeAmrFluxes& host)
	{
		if(host.coarse_fine_velocity.size()!=static_cast<std::size_t>(coarse_fine_count_)||host.embedded_velocity.size()!=static_cast<std::size_t>(special_count_-coarse_fine_count_))throw std::invalid_argument("composite special flux upload size");std::vector<Real> values(special_count_);for(int q=0;q<coarse_fine_count_;++q)values[q]=static_cast<Real>(host.coarse_fine_velocity[q]);for(int q=coarse_fine_count_;q<special_count_;++q)values[q]=static_cast<Real>(host.embedded_velocity[q-coarse_fine_count_]);if(special_count_)check(cudaMemcpy(special_velocity_,values.data(),values.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload composite special fluxes");
	}
	void DeviceCompositeAmrProjection::download_special_fluxes(CompositeAmrFluxes& host)const
	{
		std::vector<Real> values(special_count_);if(special_count_)check(cudaMemcpy(values.data(),special_velocity_,values.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download composite special fluxes");host.coarse_fine_velocity.resize(coarse_fine_count_);host.embedded_velocity.resize(special_count_-coarse_fine_count_);for(int q=0;q<coarse_fine_count_;++q)host.coarse_fine_velocity[q]=values[q];for(int q=coarse_fine_count_;q<special_count_;++q)host.embedded_velocity[q-coarse_fine_count_]=values[q];
	}
	double DeviceCompositeAmrProjection::max_abs_special_velocity()const{return max_abs_device_values(special_velocity_,special_count_,max_abs_scratch_);}
	double DeviceCompositeAmrProjection::max_embedded_cfl_rate()const
	{
		if(!embedded_count_)return 0;embedded_cfl_rate_kernel<<<(embedded_count_+255)/256,256>>>(special_velocity_,eb_transport_length_,eb_transport_scratch_,coarse_fine_count_,embedded_count_);return max_abs_device_values(eb_transport_scratch_,embedded_count_,max_abs_scratch_);
	}

	namespace
	{
		void launch_integrated_flux(const DeviceCompositeAmrFluxLevelView* levels,const std::vector<int>& brick_counts,int level_count,int brick_size,const unsigned char* active,const std::uint8_t* cut_face_mask,int storage_size,const int* first,const int* second,const std::int8_t* direction,const Real* area,const Real* velocity,int special_count,Real* integrated)
		{
			check(cudaMemset(integrated,0,static_cast<std::size_t>(storage_size)*sizeof(Real)),"clear composite integrated flux");for(int level=0;level<level_count;++level){const int work=brick_counts[level]*brick_size*brick_size*brick_size;if(work)structured_integrated_flux_kernel<<<(work+255)/256,256>>>(levels,level,brick_size,active,cut_face_mask,integrated);}if(special_count)special_integrated_flux_kernel<<<(special_count+255)/256,256>>>(first,second,direction,area,velocity,special_count,active,integrated);
		}
	}

	void DeviceCompositeAmrProjection::compute_divergence()
	{
		launch_integrated_flux(levels_,brick_counts_,level_count_,brick_size_,active_,cut_face_mask_,storage_size_,first_dof_,second_dof_,direction_,open_area_,special_velocity_,special_count_,integrated_);normalize_divergence_rhs_kernel<<<(storage_size_+255)/256,256>>>(integrated_,volume_,active_,Real(0),divergence_,nullptr,storage_size_);check(cudaDeviceSynchronize(),"compute composite AMR divergence");
	}
	void DeviceCompositeAmrProjection::build_projection_rhs(Real rho,Real dt)
	{
		if(!(rho>Real(0))||!(dt>Real(0)))throw std::invalid_argument("composite projection rho/dt");launch_integrated_flux(levels_,brick_counts_,level_count_,brick_size_,active_,cut_face_mask_,storage_size_,first_dof_,second_dof_,direction_,open_area_,special_velocity_,special_count_,integrated_);normalize_divergence_rhs_kernel<<<(storage_size_+255)/256,256>>>(integrated_,volume_,active_,rho/dt,divergence_,rhs_,storage_size_);check(cudaDeviceSynchronize(),"build composite AMR projection RHS");
	}
	void DeviceCompositeAmrProjection::correct_fluxes(Real rho,Real dt)
	{
		if(!(rho>Real(0))||!(dt>Real(0)))throw std::invalid_argument("composite projection correction rho/dt");const Real scale=dt/rho;for(int level=0;level<level_count_;++level){const int work=brick_counts_[level]*brick_size_*brick_size_*brick_size_;if(work)structured_flux_correction_kernel<<<(work+255)/256,256>>>(levels_,level,brick_size_,outlet_,active_,cut_face_mask_,pressure_,scale);}if(special_count_)special_flux_correction_kernel<<<(special_count_+255)/256,256>>>(first_dof_,second_dof_,direction_,centre_distance_,special_velocity_,special_count_,active_,pressure_,scale);if(coarse_fine_count_){check(cudaMemset(cf_group_sum_,0,static_cast<std::size_t>(coarse_fine_group_count_)*sizeof(Real)),"clear coarse/fine group flux sums");scatter_coarse_fine_tiles_kernel<<<(coarse_fine_count_+255)/256,256>>>(levels_,cf_fine_level_,cf_fine_index_,axis_,cf_group_,open_area_,special_velocity_,cf_group_sum_,coarse_fine_count_);scatter_coarse_fine_groups_kernel<<<(coarse_fine_group_count_+255)/256,256>>>(levels_,cf_group_level_,cf_group_index_,cf_group_axis_,cf_group_area_,cf_group_sum_,coarse_fine_group_count_);}for(int level=0;level<level_count_;++level){const int work=brick_counts_[level]*brick_size_*brick_size_*brick_size_;if(work)scatter_regular_pressure_kernel<<<(work+255)/256,256>>>(levels_,level,brick_size_,active_,pressure_);}check(cudaDeviceSynchronize(),"correct composite AMR fluxes");
	}
	AmrGpuSolveResult DeviceCompositeAmrProjection::project(Real rho,Real dt,double tolerance,int max_iterations,bool warm_start)
	{
		build_projection_rhs(rho,dt);AmrGpuSolveResult result=solver_.solve(pressure_,rhs_,tolerance,max_iterations,warm_start);if(result.converged){correct_fluxes(rho,dt);compute_divergence();}return result;
	}
	void DeviceCompositeAmrProjection::download_divergence(std::vector<Real>& host)const{host.resize(storage_size_);check(cudaMemcpy(host.data(),divergence_,host.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download composite divergence");}
	void DeviceCompositeAmrProjection::download_pressure(std::vector<Real>& host)const{host.resize(storage_size_);check(cudaMemcpy(host.data(),pressure_,host.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download composite pressure");}
}
