#include "core/fluid/amr_pressure.h"

#include <cuda_runtime.h>

#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform_reduce.h>

#include <cmath>
#include <stdexcept>
#include <string>
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
		__device__ int local_index(int bs, int i, int j, int k) { return (k * bs + j) * bs + i; }
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
		__global__ void initial_precondition_kernel(const Real* residual, Real* z, Real* direction, const Real* diagonal, const unsigned char* active, int n)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n){const Real v=active[q]&&diagonal[q]>Real(0)?residual[q]/diagonal[q]:Real(0);z[q]=v;direction[q]=v;}
		}
		__global__ void subtract_kernel(Real* residual,const Real* value,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)residual[q]=active[q]?residual[q]-value[q]:Real(0);}
		__global__ void update_pressure_residual_kernel(Real* pressure,Real* residual,const Real* direction,const Real* Ad,Real alpha,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n&&active[q]){pressure[q]+=alpha*direction[q];residual[q]-=alpha*Ad[q];}}
		__global__ void precondition_kernel(const Real* residual,Real* z,const Real* diagonal,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)z[q]=active[q]&&diagonal[q]>Real(0)?residual[q]/diagonal[q]:Real(0);}
		__global__ void update_direction_kernel(const Real* z,Real* direction,Real beta,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)direction[q]=active[q]?z[q]+beta*direction[q]:Real(0);}
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
	}

	DeviceCompositeAmrPressureOperator::~DeviceCompositeAmrPressureOperator()
	{
		for (Allocation& allocation : allocations_) { if(allocation.neighbors)cudaFree(allocation.neighbors);if(allocation.flags)cudaFree(allocation.flags); } if(levels_)cudaFree(levels_);if(coarse_dof_)cudaFree(coarse_dof_);if(fine_dof_)cudaFree(fine_dof_);if(coefficient_)cudaFree(coefficient_);if(active_)cudaFree(active_);if(cut_face_mask_)cudaFree(cut_face_mask_);
	}

	void DeviceCompositeAmrPressureOperator::apply(const Real* pressure, Real* output) const
	{
		check(cudaMemset(output,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite AMR output");const int cells_per_brick=brick_size_*brick_size_*brick_size_;
		for(int level=0;level<level_count_;++level){const int work=brick_counts_[level]*cells_per_brick;if(work)structured_apply_kernel<<<(work+255)/256,256>>>(levels_,level,brick_size_,outlet_,active_,cut_face_mask_,pressure,output);}if(connection_count_)coarse_fine_apply_kernel<<<(connection_count_+255)/256,256>>>(coarse_dof_,fine_dof_,coefficient_,connection_count_,active_,pressure,output);check(cudaDeviceSynchronize(),"apply composite AMR pressure operator");
	}

	DeviceCompositeAmrPressureSolver::DeviceCompositeAmrPressureSolver(const CompositeAmrPressureSystem& system):op_(system),n_(system.storage_size)
	{
		std::vector<double> diagonal_double;system.diagonal_cpu(diagonal_double);std::vector<Real> diagonal(diagonal_double.size());for(std::size_t q=0;q<diagonal.size();++q)diagonal[q]=static_cast<Real>(diagonal_double[q]);std::vector<Real> zero(n_);
		r_=upload(zero,"allocate AMR residual");z_=upload(zero,"allocate AMR preconditioned residual");direction_=upload(zero,"allocate AMR direction");Ad_=upload(zero,"allocate AMR operator temporary");diagonal_=upload(diagonal,"upload AMR diagonal");active_=upload(system.active,"upload AMR active mask");bytes_=op_.bytes()+static_cast<std::size_t>(5)*n_*sizeof(Real)+system.active.size();
	}
	DeviceCompositeAmrPressureSolver::~DeviceCompositeAmrPressureSolver(){for(void* pointer:{(void*)r_,(void*)z_,(void*)direction_,(void*)Ad_,(void*)diagonal_,(void*)active_})if(pointer)cudaFree(pointer);}
	AmrGpuSolveResult DeviceCompositeAmrPressureSolver::solve(Real* pressure,const Real* rhs,double tolerance,int max_iterations,bool warm_start)
	{
		AmrGpuSolveResult result;if(!warm_start)check(cudaMemset(pressure,0,static_cast<std::size_t>(n_)*sizeof(Real)),"clear composite AMR pressure");if(warm_start){op_.apply(pressure,Ad_);check(cudaMemcpy(r_,rhs,static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),"copy composite AMR rhs");subtract_kernel<<<(n_+255)/256,256>>>(r_,Ad_,active_,n_);}else check(cudaMemcpy(r_,rhs,static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),"initial composite AMR residual");
		const double rhs2=dot_active(rhs,rhs,active_,n_);if(!(rhs2>0)){result.converged=true;return result;}double residual2=dot_active(r_,r_,active_,n_);result.relative_residual=std::sqrt(residual2/rhs2);if(result.relative_residual<=tolerance){result.converged=true;return result;}initial_precondition_kernel<<<(n_+255)/256,256>>>(r_,z_,direction_,diagonal_,active_,n_);double rz=dot_active(r_,z_,active_,n_);
		for(int iteration=0;iteration<max_iterations;++iteration){op_.apply(direction_,Ad_);const double dAd=dot_active(direction_,Ad_,active_,n_);if(!(dAd>0))break;const Real alpha=static_cast<Real>(rz/dAd);update_pressure_residual_kernel<<<(n_+255)/256,256>>>(pressure,r_,direction_,Ad_,alpha,active_,n_);residual2=dot_active(r_,r_,active_,n_);result.iterations=iteration+1;result.relative_residual=std::sqrt(residual2/rhs2);if(result.relative_residual<=tolerance){result.converged=true;break;}precondition_kernel<<<(n_+255)/256,256>>>(r_,z_,diagonal_,active_,n_);const double next_rz=dot_active(r_,z_,active_,n_);const Real beta=static_cast<Real>(next_rz/rz);update_direction_kernel<<<(n_+255)/256,256>>>(z_,direction_,beta,active_,n_);rz=next_rz;}check(cudaDeviceSynchronize(),"solve composite AMR pressure");return result;
	}
}
