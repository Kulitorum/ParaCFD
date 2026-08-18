#include "core/fluid/amr_pressure.h"
#include "core/fluid/amr_advection.h"

#include <cuda_runtime.h>

#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform_reduce.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
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
		int symmetric_pseudoinverse_3x3(const double* input,double* inverse)
		{
			double a[9],vectors[9]={1,0,0,0,1,0,0,0,1};for(int q=0;q<9;++q)a[q]=input[q];
			for(int sweep=0;sweep<16;++sweep)
			{
				bool changed=false;for(int pair=0;pair<3;++pair)
				{
					const int p=pair==0?0:(pair==1?0:1),q=pair==0?1:(pair==1?2:2);const double apq=a[p*3+q],scale=std::max({std::abs(a[p*3+p]),std::abs(a[q*3+q]),1e-300});if(std::abs(apq)<=1e-14*scale)continue;changed=true;const double phi=0.5*std::atan2(2*apq,a[q*3+q]-a[p*3+p]),c=std::cos(phi),s=std::sin(phi),app=a[p*3+p],aqq=a[q*3+q];for(int k=0;k<3;++k)if(k!=p&&k!=q){const double akp=a[k*3+p],akq=a[k*3+q];a[k*3+p]=a[p*3+k]=c*akp-s*akq;a[k*3+q]=a[q*3+k]=s*akp+c*akq;}a[p*3+p]=c*c*app-2*s*c*apq+s*s*aqq;a[q*3+q]=s*s*app+2*s*c*apq+c*c*aqq;a[p*3+q]=a[q*3+p]=0;for(int k=0;k<3;++k){const double vkp=vectors[k*3+p],vkq=vectors[k*3+q];vectors[k*3+p]=c*vkp-s*vkq;vectors[k*3+q]=s*vkp+c*vkq;}
				}
				if(!changed)break;
			}
			// The geometry solve runs in double, but the inverse and its RHS are production
			// FP32 data. Retaining a mode with condition number above 1e4 would amplify
			// float roundoff into a fictitious derivative at nearly coplanar seams. Treat
			// that direction as geometrically unsupported. FP64 validation builds retain
			// the much smaller double-precision cutoff.
			for(int q=0;q<9;++q)inverse[q]=0;const double maximum=std::max({std::abs(a[0]),std::abs(a[4]),std::abs(a[8])}),relative_cutoff=sizeof(Real)==4?1e-4:1e-10,threshold=maximum*relative_cutoff;int rank=0;for(int mode=0;mode<3;++mode){const double eigenvalue=a[mode*3+mode];if(!(eigenvalue>threshold))continue;++rank;for(int row=0;row<3;++row)for(int column=0;column<3;++column)inverse[row*3+column]+=vectors[row*3+mode]*vectors[column*3+mode]/eigenvalue;}return rank;
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

		__global__ void reconstruct_pressure_gradient_kernel(const int* node_dof,const int* offset,
			const int* neighbour,const Real* weight,const Real* pressure,Real* gradient,int node_count)
		{
			const int node=blockIdx.x*blockDim.x+threadIdx.x;if(node>=node_count)return;const Real centre=pressure[node_dof[node]];Real gx=0,gy=0,gz=0;for(int q=offset[node];q<offset[node+1];++q){const Real delta=pressure[neighbour[q]]-centre;gx+=weight[3*q]*delta;gy+=weight[3*q+1]*delta;gz+=weight[3*q+2]*delta;}gradient[3*node]=gx;gradient[3*node+1]=gy;gradient[3*node+2]=gz;
		}
		// Irregular pressure edges touch only a compact subset of the grid.  A thread per
		// edge would require unordered atomic additions at both endpoints, making the
		// FP32 operator depend on warp scheduling.  The CPU-built incidence table gives
		// one thread exclusive ownership of each touched DOF and a stable edge order.
		__global__ void irregular_orthogonal_gather_kernel(
			const int* incident_dof,const int* incident_offset,const int* incident_edge,
			const std::int8_t* incident_sign,int incident_dof_count,int connection_count,
			const int* connection_lower,const int* connection_upper,const Real* connection_coefficient,
			const int* regular_lower,const int* regular_upper,const Real* regular_delta,
			const Real* regular_area,const unsigned char* active,const Real* pressure,Real* output)
		{
			const int node=blockIdx.x*blockDim.x+threadIdx.x;if(node>=incident_dof_count)return;
			double sum=0.0;
			for(int q=incident_offset[node];q<incident_offset[node+1];++q)
			{
				const int encoded=incident_edge[q];int a,b;Real coefficient;
				if(encoded<connection_count)
				{
					a=connection_lower[encoded];b=connection_upper[encoded];coefficient=connection_coefficient[encoded];
				}
				else
				{
					const int edge=encoded-connection_count;a=regular_lower[edge];b=regular_upper[edge];
					coefficient=regular_area[edge]*regular_delta[edge];
				}
				if(active[a]&&active[b])
				{
					const Real flux=coefficient*(pressure[a]-pressure[b]);
					sum+=static_cast<double>(incident_sign[q])*static_cast<double>(flux);
				}
			}
			output[incident_dof[node]]+=static_cast<Real>(sum);
		}
		__global__ void irregular_nonorthogonal_gather_kernel(
			const int* incident_dof,const int* incident_offset,const int* incident_edge,
			const std::int8_t* incident_sign,int incident_dof_count,int connection_count,
			const int* connection_lower,const int* connection_upper,const int* connection_lower_node,
			const int* connection_upper_node,const Real* connection_correction,
			const Real* connection_upper_weight,const Real* connection_area,
			const int* regular_lower,const int* regular_upper,const int* regular_lower_node,
			const int* regular_upper_node,const Real* regular_correction,
			const Real* regular_upper_weight,const Real* regular_area,
			const Real* gradient,const unsigned char* active,Real* output)
		{
			const int node=blockIdx.x*blockDim.x+threadIdx.x;if(node>=incident_dof_count)return;
			double sum=0.0;
			for(int q=incident_offset[node];q<incident_offset[node+1];++q)
			{
				const int encoded=incident_edge[q];int a,b,lower_node,upper_node;const Real* correction;
				Real upper_weight,area;
				if(encoded<connection_count)
				{
					a=connection_lower[encoded];b=connection_upper[encoded];
					lower_node=connection_lower_node[encoded];upper_node=connection_upper_node[encoded];
					correction=connection_correction+3*encoded;
					upper_weight=connection_upper_weight[encoded];area=connection_area[encoded];
				}
				else
				{
					const int edge=encoded-connection_count;a=regular_lower[edge];b=regular_upper[edge];
					lower_node=regular_lower_node[edge];upper_node=regular_upper_node[edge];
					correction=regular_correction+3*edge;
					upper_weight=regular_upper_weight[edge];area=regular_area[edge];
				}
				if(!active[a]||!active[b]||lower_node<0||upper_node<0)continue;
				const Real lower_weight=Real(1)-upper_weight;Real derivative=Real(0);
				for(int component=0;component<3;++component)
				{
					const Real face_gradient=lower_weight*gradient[3*lower_node+component]+
						upper_weight*gradient[3*upper_node+component];
					derivative+=correction[component]*face_gradient;
				}
				const Real flux=-area*derivative;
				sum+=static_cast<double>(incident_sign[q])*static_cast<double>(flux);
			}
			output[incident_dof[node]]+=static_cast<Real>(sum);
		}
		__global__ void gauge_apply_kernel(const int* dof,const Real* coefficient,int count,const Real* pressure,Real* output){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<count)output[dof[q]]+=coefficient[q]*pressure[dof[q]];}
		__global__ void subtract_kernel(Real* residual,const Real* value,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)residual[q]=active[q]?residual[q]-value[q]:Real(0);}
		__global__ void update_pressure_residual_kernel(Real* pressure,Real* residual,const Real* direction,const Real* Ad,Real alpha,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n&&active[q]){pressure[q]+=alpha*direction[q];residual[q]-=alpha*Ad[q];}}
		__global__ void precondition_kernel(const Real* residual,Real* z,const Real* diagonal,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)z[q]=active[q]&&diagonal[q]>Real(0)?residual[q]/diagonal[q]:Real(0);}
		__global__ void update_direction_kernel(const Real* z,Real* direction,Real beta,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)direction[q]=active[q]?z[q]+beta*direction[q]:Real(0);}
		__global__ void copy_active_kernel(const Real* source,Real* destination,const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)destination[q]=active[q]?source[q]:Real(0);}
		__global__ void axpy_active_kernel(Real* value,const Real* increment,Real scale,
			const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n&&active[q])value[q]+=scale*increment[q];}
		__global__ void scaled_copy_active_kernel(Real* destination,const Real* source,Real scale,
			const unsigned char* active,int n){const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q<n)destination[q]=active[q]?scale*source[q]:Real(0);}
		__global__ void restrict_aggregate_gather_kernel(const int* offset,const int* dof,
			const Real* residual,int base_cells,Real* coarse_rhs)
		{
			const int coarse=blockIdx.x*blockDim.x+threadIdx.x;if(coarse>=base_cells)return;
			double sum=0.0;for(int q=offset[coarse];q<offset[coarse+1];++q)sum+=static_cast<double>(residual[dof[q]]);
			coarse_rhs[coarse]=static_cast<Real>(sum);
		}
		__global__ void coarse_jacobi_kernel(const Real* rhs,const Real* x,Real* next,const Real* positive,const Real* diagonal,const int* neighbors,int bricks,int bs,Real omega)
		{
			const int cells_per_brick=bs*bs*bs,work=blockIdx.x*blockDim.x+threadIdx.x;if(work>=bricks*cells_per_brick)return;if(!(diagonal[work]>Real(0))){next[work]=Real(0);return;}const int brick=work/cells_per_brick,local=work%cells_per_brick,i=local%bs,j=(local/bs)%bs,k=local/(bs*bs);Real Ax=diagonal[work]*x[work];for(int axis=0;axis<3;++axis){int c[3]={i,j,k},other=-1;if(++c[axis]<bs)other=brick*cells_per_brick+local_index(bs,c[0],c[1],c[2]);else{const int neighbour=neighbors[brick*6+2*axis+1];if(neighbour>=0){c[axis]=0;other=neighbour*cells_per_brick+local_index(bs,c[0],c[1],c[2]);}}if(other>=0)Ax-=positive[work*3+axis]*x[other];c[0]=i;c[1]=j;c[2]=k;int lower=brick;if(--c[axis]<0){lower=neighbors[brick*6+2*axis];if(lower>=0)c[axis]=bs-1;}if(lower>=0){const int lower_cell=lower*cells_per_brick+local_index(bs,c[0],c[1],c[2]);Ax-=positive[lower_cell*3+axis]*x[lower_cell];}}next[work]=x[work]+omega*(rhs[work]-Ax)/diagonal[work];
		}
		__global__ void coarse_extra_jacobi_gather_kernel(const int* incident_dof,
			const int* incident_offset,const int* incident_edge,int incident_dof_count,
			const int* a,const int* b,const Real* coefficient,const Real* x,
			const Real* diagonal,Real omega,Real* next)
		{
			const int node=blockIdx.x*blockDim.x+threadIdx.x;if(node>=incident_dof_count)return;
			const int cell=incident_dof[node];if(!(diagonal[cell]>Real(0)))return;double sum=0.0;
			for(int q=incident_offset[node];q<incident_offset[node+1];++q)
			{
				const int edge=incident_edge[q],other=a[edge]==cell?b[edge]:a[edge];
				const Real contribution=omega*coefficient[edge]*x[other]/diagonal[cell];
				sum+=static_cast<double>(contribution);
			}
			next[cell]+=static_cast<Real>(sum);
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
		__global__ void special_integrated_flux_gather_kernel(const int* incident_dof,
			const int* incident_offset,const int* incident_edge,const std::int8_t* incident_sign,
			int incident_dof_count,const Real* area,const Real* velocity,Real* integrated)
		{
			const int node=blockIdx.x*blockDim.x+threadIdx.x;if(node>=incident_dof_count)return;
			double sum=0.0;for(int q=incident_offset[node];q<incident_offset[node+1];++q)
			{
				const int edge=incident_edge[q];const Real flux=area[edge]*velocity[edge];
				sum+=static_cast<double>(incident_sign[q])*static_cast<double>(flux);
			}
			integrated[incident_dof[node]]+=static_cast<Real>(sum);
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
			const Real* gradient_factor,Real* velocity,int count,const unsigned char* active,const Real* pressure,Real scale)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int lower=direction[edge]>0?first[edge]:second[edge],upper=direction[edge]>0?second[edge]:first[edge];if(active[lower]&&active[upper])velocity[edge]-=scale*(pressure[upper]-pressure[lower])*gradient_factor[edge];
		}
		__global__ void special_nonorthogonal_flux_correction_kernel(const int* first,const int* second,
			const std::int8_t* direction,const int* lower_node,const int* upper_node,const Real* correction,
			const Real* upper_weight,const Real* gradient,Real* velocity,int count,const unsigned char* active,Real scale)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count||lower_node[edge]<0||upper_node[edge]<0)return;const int lower=direction[edge]>0?first[edge]:second[edge],upper=direction[edge]>0?second[edge]:first[edge];if(!active[lower]||!active[upper])return;const Real w=upper_weight[edge],lw=Real(1)-w;Real derivative=0;for(int component=0;component<3;++component){const Real face_gradient=lw*gradient[3*lower_node[edge]+component]+w*gradient[3*upper_node[edge]+component];derivative+=correction[3*edge+component]*face_gradient;}velocity[edge]-=scale*derivative;
		}
		__global__ void regular_nonorthogonal_flux_correction_kernel(
			const DeviceCompositeAmrFluxLevelView* levels,const int* level,
			const std::uint64_t* index,const std::uint64_t* mirror_index,const std::int8_t* axis,
			const int* lower,const int* upper,const int* lower_node,const int* upper_node,
			const Real* two_point_delta,const Real* correction,const Real* upper_weight,
			const Real* gradient,const Real* pressure,Real scale,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;
			Real derivative=two_point_delta[edge]*(pressure[upper[edge]]-pressure[lower[edge]]);
			if(lower_node[edge]>=0&&upper_node[edge]>=0)
			{
				const Real w=upper_weight[edge],lw=Real(1)-w;
				for(int component=0;component<3;++component)
				{
					const Real face_gradient=lw*gradient[3*lower_node[edge]+component]+
						w*gradient[3*upper_node[edge]+component];
					derivative+=correction[3*edge+component]*face_gradient;
				}
			}
			const DeviceAmrFieldLevelView view=levels[level[edge]].fields;
			Real* values=axis[edge]==0?view.u:(axis[edge]==1?view.v:view.w);
			values[index[edge]]-=scale*derivative;
			if(mirror_index[edge]!=~std::uint64_t(0))values[mirror_index[edge]]-=scale*derivative;
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
		__global__ void accumulate_eb_gradient_rhs_kernel(const int* node_a,const int* node_b,const Real* displacement,const Real* ls_weight,const Real* sum,const Real* weight,Real* gradient_rhs,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int a=node_a[edge],b=node_b[edge];for(int component=0;component<3;++component){const Real wa=weight[a*3+component],wb=weight[b*3+component],sample_weight=ls_weight[edge*3+component];if(!(wa>Real(0)&&wb>Real(0)&&sample_weight>Real(0)))continue;const Real delta=sum[b*3+component]/wb-sum[a*3+component]/wa;for(int derivative=0;derivative<3;++derivative){const Real contribution=sample_weight*displacement[edge*9+component*3+derivative]*delta;atomicAdd(gradient_rhs+a*9+component*3+derivative,contribution);atomicAdd(gradient_rhs+b*9+component*3+derivative,contribution);}}
		}
		__global__ void solve_eb_gradient_kernel(Real* gradient,const Real* inverse,int node_count)
		{
			const int node=blockIdx.x*blockDim.x+threadIdx.x;if(node>=node_count)return;for(int component=0;component<3;++component){Real rhs[3];for(int q=0;q<3;++q)rhs[q]=gradient[node*9+component*3+q];for(int row=0;row<3;++row){Real value=Real(0);for(int column=0;column<3;++column)value+=inverse[node*27+component*9+row*3+column]*rhs[column];gradient[node*9+component*3+row]=value;}}
		}
		__global__ void transport_eb_aperture_kernel(const int* node_a,const int* node_b,const int* negative_edge,const int* positive_edge,const std::int8_t* direction,const std::int8_t* axis,const Real* transport_length,const Real* velocity,const Real* sum,const Real* weight,const Real* gradient_sum,Real dt,Real molecular_nu,Real smagorinsky_cs,Real* output,Real* diffusion_rate,int special_offset,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int global=special_offset+edge,component=axis[global],first=node_a[edge],second=node_b[edge],lower=direction[global]>0?first:second,upper=direction[global]>0?second:first;const Real current=velocity[global],lower_weight=weight[lower*3+component],upper_weight=weight[upper*3+component],lower_value=lower_weight>Real(0)?sum[lower*3+component]/lower_weight:current,upper_value=upper_weight>Real(0)?sum[upper*3+component]/upper_weight:current,d=max(transport_length[edge],Real(1e-12));Real gradient[9],advector[3];for(int velocity_component=0;velocity_component<3;++velocity_component){const Real a_weight=weight[lower*3+velocity_component],b_weight=weight[upper*3+velocity_component];Real value=Real(0),samples=Real(0);if(a_weight>Real(0)){value+=sum[lower*3+velocity_component]/a_weight;samples+=Real(1);}if(b_weight>Real(0)){value+=sum[upper*3+velocity_component]/b_weight;samples+=Real(1);}advector[velocity_component]=samples>Real(0)?value/samples:Real(0);}advector[component]=current;for(int slot=0;slot<9;++slot)gradient[slot]=Real(0.5)*(gradient_sum[lower*9+slot]+gradient_sum[upper*9+slot]);const Real normal_rate=abs(upper_value-lower_value)/d,strain=max(sqrt(max(Real(0),detail::compact_strain_magnitude_squared(gradient))),normal_rate),length=smagorinsky_cs*d,eddy_nu=length*length*strain;diffusion_rate[edge]=(molecular_nu+eddy_nu)/(d*d);const Real courant=min(Real(1),abs(current)*dt/d);const bool positive=current>=Real(0);const int upstream_edge=positive?negative_edge[edge]:positive_edge[edge],downstream_edge=positive?positive_edge[edge]:negative_edge[edge],upstream2_edge=upstream_edge>=0?(positive?negative_edge[upstream_edge]:positive_edge[upstream_edge]):-1;const bool complete=upstream_edge>=0&&upstream2_edge>=0&&downstream_edge>=0;const Real upstream=complete?velocity[special_offset+upstream_edge]:current,upstream2=complete?velocity[special_offset+upstream2_edge]:current,downstream=complete?velocity[special_offset+downstream_edge]:current;const Real normal_update=detail::bounded_compact_transport_update(current,lower_value,upper_value,upstream,upstream2,downstream,courant,Real(0),complete);output[edge]=detail::bounded_compact_tangential_update(normal_update,current,lower_value,upper_value,advector,gradient,component,dt);
		}
		__global__ void accumulate_eb_dual_momentum_kernel(const int* node_a,const int* node_b,const std::int8_t* axis,const Real* area,const Real* gradient_factor,const Real* velocity,Real* sum,Real* weight,int special_offset,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int global=special_offset+edge,component=axis[global],a=node_a[edge]*3+component,b=node_b[edge]*3+component;const Real mass=Real(0.5)*area[global]/max(gradient_factor[global],Real(1e-12)),weighted=mass*velocity[edge];atomicAdd(sum+a,weighted);atomicAdd(sum+b,weighted);atomicAdd(weight+a,mass);atomicAdd(weight+b,mass);
		}
		__global__ void accumulate_conservative_carrier_kernel(const DeviceCompositeAmrFluxLevelView* levels,const int* node,const int* level,const std::uint64_t* index,const std::int8_t* axis,const Real* mass,Real* sum,Real* weight,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;const int slot=node[q]*3+axis[q];const Real m=Real(0.5)*mass[q],u=indexed_face_value(levels[level[q]].fields,axis[q],index[q]);atomicAdd(sum+slot,m*u);atomicAdd(weight+slot,m);
		}
		__global__ void accumulate_unmapped_carrier_momentum_kernel(const DeviceCompositeAmrFluxLevelView* levels,const int* node,const int* level,const std::uint64_t* index,const std::int8_t* axis,const Real* unmapped_mass,Real* sum,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count||!(unmapped_mass[q]>Real(0)))return;const int slot=node[q]*3+axis[q],component=axis[q];atomicAdd(sum+slot,unmapped_mass[q]*indexed_face_value(levels[level[q]].fields,static_cast<std::int8_t>(component),index[q]));
		}
		__global__ void advect_eb_node_momentum_kernel(const int* node_a,const int* node_b,
			const std::int8_t* direction,const Real* area,const Real* normal_velocity,
			const Real* sum,const Real* weight,Real dt,Real* momentum_delta,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int lower=direction[edge]>0?node_a[edge]:node_b[edge],upper=direction[edge]>0?node_b[edge]:node_a[edge];const Real swept=dt*area[edge]*normal_velocity[edge],donor=swept>=Real(0)?lower:upper;for(int component=0;component<3;++component){const int source=donor*3+component,lower_slot=lower*3+component,upper_slot=upper*3+component;if(!(weight[source]>Real(0)&&weight[lower_slot]>Real(0)&&weight[upper_slot]>Real(0)))continue;const Real transfer=swept*sum[source]/weight[source];atomicAdd(momentum_delta+lower_slot,-transfer);atomicAdd(momentum_delta+upper_slot,transfer);}
		}
		__global__ void reconstruct_compatible_eb_kernel(const int* node_a,const int* node_b,const std::int8_t* axis,const Real* sum,const Real* weight,const Real* advected,Real* output,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int component=axis[edge],sa=node_a[edge]*3+component,sb=node_b[edge]*3+component;const Real wa=weight[sa],wb=weight[sb],current=advected[edge],ua=wa>Real(0)?sum[sa]/wa:current,ub=wb>Real(0)?sum[sb]/wb:current;output[edge]=Real(0.5)*(ua+ub);
		}
		__global__ void reconstruct_compatible_carrier_kernel(const DeviceCompositeAmrFluxLevelView* levels,const int* node,const int* level,const std::uint64_t* index,const std::int8_t* axis,const Real* sum,const Real* weight,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;const int slot=node[q]*3+axis[q];if(!(weight[slot]>Real(0)))return;Real* destination=axis[q]==0?levels[level[q]].fields.u:(axis[q]==1?levels[level[q]].fields.v:levels[level[q]].fields.w);const Real current=destination[index[q]],mean=sum[slot]/weight[slot];atomicAdd(destination+index[q],Real(0.5)*(mean-current));
		}
		__global__ void count_node_diffusion_neighbours_kernel(const int* node_a,const int* node_b,const Real* weight,Real* neighbour_count,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int a=node_a[edge],b=node_b[edge];for(int component=0;component<3;++component){const int sa=a*3+component,sb=b*3+component;if(weight[sa]>Real(0)&&weight[sb]>Real(0)){atomicAdd(neighbour_count+sa,Real(1));atomicAdd(neighbour_count+sb,Real(1));}}
		}
		__global__ void exchange_node_momentum_kernel(const int* node_a,const int* node_b,const Real* diffusion_rate,const Real* sum,const Real* weight,const Real* neighbour_count,Real dt,Real* momentum_delta,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int a=node_a[edge],b=node_b[edge];for(int component=0;component<3;++component){const int sa=a*3+component,sb=b*3+component;const Real ma=weight[sa],mb=weight[sb];if(!(ma>Real(0)&&mb>Real(0)))continue;const Real degree=max(Real(1),max(neighbour_count[sa],neighbour_count[sb])),theta=min(Real(0.25),dt*diffusion_rate[edge])/degree,ua=sum[sa]/ma,ub=sum[sb]/mb,transfer=detail::compact_diffusive_momentum_transfer(ma,mb,ua,ub,theta);atomicAdd(momentum_delta+sa,transfer);atomicAdd(momentum_delta+sb,-transfer);}
		}
		__global__ void apply_node_momentum_to_eb_kernel(const int* node_a,const int* node_b,const std::int8_t* axis,const Real* weight,const Real* momentum_delta,const Real* advected,Real* output,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int component=axis[edge],sa=node_a[edge]*3+component,sb=node_b[edge]*3+component;const Real wa=weight[sa],wb=weight[sb],du_a=wa>Real(0)?momentum_delta[sa]/wa:Real(0),du_b=wb>Real(0)?momentum_delta[sb]/wb:Real(0);output[edge]=advected[edge]+Real(0.5)*(du_a+du_b);
		}
		__global__ void apply_node_momentum_to_carrier_kernel(const DeviceCompositeAmrFluxLevelView* levels,const int* node,const int* level,const std::uint64_t* index,const std::int8_t* axis,const Real* weight,const Real* momentum_delta,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;const int slot=node[q]*3+axis[q];if(!(weight[slot]>Real(0)))return;Real* destination=axis[q]==0?levels[level[q]].fields.u:(axis[q]==1?levels[level[q]].fields.v:levels[level[q]].fields.w);atomicAdd(destination+index[q],Real(0.5)*momentum_delta[slot]/weight[slot]);
		}
		__global__ void accumulate_wall_aperture_velocity_kernel(const int* node_a,const int* node_b,const std::int8_t* axis,const Real* area,const Real* gradient_factor,const Real* velocity,Real* sum,Real* weight,int special_offset,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int global=special_offset+edge,component=axis[global],a=node_a[edge],b=node_b[edge];const Real mass=Real(0.5)*area[global]/max(gradient_factor[global],Real(1e-12)),weighted=mass*velocity[global];if(a>=0){atomicAdd(sum+a*3+component,weighted);atomicAdd(weight+a*3+component,mass);}if(b>=0){atomicAdd(sum+b*3+component,weighted);atomicAdd(weight+b*3+component,mass);}
		}
		__global__ void accumulate_wall_carrier_velocity_kernel(const DeviceCompositeAmrFluxLevelView* levels,const int* node,const int* level,const std::uint64_t* index,const std::int8_t* axis,const Real* mass,Real* sum,Real* weight,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;const int slot=node[q]*3+axis[q];const Real m=Real(0.5)*mass[q],u=indexed_face_value(levels[level[q]].fields,axis[q],index[q]);atomicAdd(sum+slot,m*u);atomicAdd(weight+slot,m);
		}
		__global__ void compute_no_slip_wall_delta_kernel(const Real* sum,const Real* weight,const Real* wall_rate,Real dt,Real molecular_nu,Real* velocity_delta,int node_count)
		{
			const int node=blockIdx.x*blockDim.x+threadIdx.x;if(node>=node_count)return;const Real theta=max(Real(0),dt*molecular_nu*wall_rate[node]),fraction=theta/(Real(1)+theta);for(int component=0;component<3;++component){const int slot=node*3+component;velocity_delta[slot]=weight[slot]>Real(0)?-fraction*sum[slot]/weight[slot]:Real(0);}
		}
		__global__ void accumulate_smooth_wall_node_matrix_kernel(const int* node,
			const Real* normal,const Real* area,const Real* distance,const Real* sum,
			const Real* weight,Real molecular_nu,Real* matrix,Real* coefficient,int count)
		{
			const int patch=blockIdx.x*blockDim.x+threadIdx.x;
			if(patch>=count)return;
			const int target=node[patch];
			Real velocity[3];
			for(int component=0;component<3;++component)
			{
				const int slot=target*3+component;
				velocity[component]=weight[slot]>Real(0)?sum[slot]/weight[slot]:Real(0);
			}
			const Real* n=normal+3*patch;
			const Real un=velocity[0]*n[0]+velocity[1]*n[1]+velocity[2]*n[2];
			const Real tangent[3]={velocity[0]-un*n[0],velocity[1]-un*n[1],
				velocity[2]-un*n[2]};
			const Real speed=sqrt(max(Real(0),tangent[0]*tangent[0]+
				tangent[1]*tangent[1]+tangent[2]*tangent[2]));
			Real k=Real(0);
			if(speed>Real(0))
			{
				const Real friction=detail::smooth_wall_friction_velocity(
					speed,distance[patch],molecular_nu);
				k=area[patch]*friction*friction/speed;
			}
			coefficient[patch]=k;
			if(!(k>Real(0)))return;
			Real* a=matrix+6*target;
			atomicAdd(a,k*(Real(1)-n[0]*n[0]));
			atomicAdd(a+1,-k*n[0]*n[1]);
			atomicAdd(a+2,-k*n[0]*n[2]);
			atomicAdd(a+3,k*(Real(1)-n[1]*n[1]));
			atomicAdd(a+4,-k*n[1]*n[2]);
			atomicAdd(a+5,k*(Real(1)-n[2]*n[2]));
		}
		__global__ void solve_implicit_smooth_wall_node_kernel(const Real* sum,
			const Real* weight,const Real* matrix,Real dt,Real* velocity_delta,int count)
		{
			const int node=blockIdx.x*blockDim.x+threadIdx.x;
			if(node>=count)return;
			const Real* k=matrix+6*node;
			Real root_mass[3],old_velocity[3],rhs[3];
			for(int component=0;component<3;++component)
			{
				const int slot=node*3+component;
				root_mass[component]=sqrt(max(Real(0),weight[slot]));
				old_velocity[component]=weight[slot]>Real(0)?sum[slot]/weight[slot]:Real(0);
				rhs[component]=root_mass[component]*old_velocity[component];
			}
			const Real k00=k[0],k01=k[1],k02=k[2],k11=k[3],k12=k[4],k22=k[5];
			auto scaled=[&](int a,int b,Real value)
			{
				return root_mass[a]>Real(0)&&root_mass[b]>Real(0)?
					dt*value/(root_mass[a]*root_mass[b]):Real(0);
			};
			const Real a00=Real(1)+scaled(0,0,k00),a01=scaled(0,1,k01),
				a02=scaled(0,2,k02),a11=Real(1)+scaled(1,1,k11),
				a12=scaled(1,2,k12),a22=Real(1)+scaled(2,2,k22);
			// I + dt M^-1/2 K M^-1/2 is symmetric positive definite because
			// every patch contributes k(I-nn^T), k >= 0.  Cholesky therefore
			// gives an unconditionally dissipative wall impulse without a wall CFL.
			const Real l00=sqrt(a00),l10=a01/l00,l20=a02/l00;
			const Real l11=sqrt(a11-l10*l10),l21=(a12-l20*l10)/l11;
			const Real l22=sqrt(a22-l20*l20-l21*l21);
			const Real y0=rhs[0]/l00,y1=(rhs[1]-l10*y0)/l11,
				y2=(rhs[2]-l20*y0-l21*y1)/l22;
			const Real solved2=y2/l22,solved1=(y1-l21*solved2)/l11,
				solved0=(y0-l10*solved1-l20*solved2)/l00;
			const Real solved[3]={solved0,solved1,solved2};
			for(int component=0;component<3;++component)
			{
				const int slot=node*3+component;
				velocity_delta[slot]=root_mass[component]>Real(0)?
					solved[component]/root_mass[component]-old_velocity[component]:Real(0);
			}
		}
		__global__ void finalize_smooth_wall_patch_force_kernel(const int* node,
			const Real* normal,const Real* sum,const Real* weight,const Real* velocity_delta,
			const Real* coefficient,Real* force_per_density,int count)
		{
			const int patch=blockIdx.x*blockDim.x+threadIdx.x;
			if(patch>=count)return;
			const int target=node[patch];
			Real velocity[3];
			for(int component=0;component<3;++component)
			{
				const int slot=target*3+component;
				velocity[component]=(weight[slot]>Real(0)?sum[slot]/weight[slot]:Real(0))+
					velocity_delta[slot];
			}
			const Real* n=normal+3*patch;
			const Real un=velocity[0]*n[0]+velocity[1]*n[1]+velocity[2]*n[2],
				k=coefficient[patch];
			force_per_density[3*patch]=k*(velocity[0]-un*n[0]);
			force_per_density[3*patch+1]=k*(velocity[1]-un*n[1]);
			force_per_density[3*patch+2]=k*(velocity[2]-un*n[2]);
		}
		__global__ void accumulate_wall_aperture_momentum_kernel(const int* node_a,const int* node_b,
			const std::int8_t* axis,const Real* area,const Real* gradient_factor,const Real* velocity,
			double* momentum,int special_offset,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count||(node_a[edge]<0&&node_b[edge]<0))return;const int global=special_offset+edge,component=axis[global];const Real mass=area[global]/max(gradient_factor[global],Real(1e-12));atomicAdd(momentum+component,static_cast<double>(mass)*velocity[global]);
		}
		__global__ void accumulate_wall_carrier_momentum_kernel(const DeviceCompositeAmrFluxLevelView* levels,
			const int* level,const std::uint64_t* index,const std::int8_t* axis,const Real* mass,
			double* momentum,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;atomicAdd(momentum+axis[q],static_cast<double>(mass[q])*indexed_face_value(levels[level[q]].fields,axis[q],index[q]));
		}
		__global__ void apply_no_slip_wall_to_aperture_kernel(const int* node_a,const int* node_b,const std::int8_t* axis,const Real* velocity_delta,Real* velocity,int special_offset,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;const int component=axis[special_offset+edge],a=node_a[edge],b=node_b[edge];Real delta=Real(0);if(a>=0)delta+=Real(0.5)*velocity_delta[a*3+component];if(b>=0)delta+=Real(0.5)*velocity_delta[b*3+component];velocity[special_offset+edge]+=delta;
		}
		__global__ void apply_no_slip_wall_to_carrier_kernel(const DeviceCompositeAmrFluxLevelView* levels,const int* node,const int* level,const std::uint64_t* index,const std::int8_t* axis,const Real* velocity_delta,int count)
		{
			const int q=blockIdx.x*blockDim.x+threadIdx.x;if(q>=count)return;Real* destination=axis[q]==0?levels[level[q]].fields.u:(axis[q]==1?levels[level[q]].fields.v:levels[level[q]].fields.w);atomicAdd(destination+index[q],Real(0.5)*velocity_delta[node[q]*3+axis[q]]);
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
		__global__ void accumulate_embedded_node_flux_kernel(const int* node_a,const int* node_b,
			const std::int8_t* direction,const Real* area,const Real* velocity,Real* balance,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;
			const int lower=direction[edge]>0?node_a[edge]:node_b[edge],upper=direction[edge]>0?node_b[edge]:node_a[edge];
			const Real flux=area[edge]*velocity[edge];atomicAdd(balance+lower,flux);atomicAdd(balance+upper,-flux);
		}
		__global__ void accumulate_regular_perimeter_flux_kernel(const DeviceCompositeAmrFluxLevelView* levels,
			const int* lower_node,const int* upper_node,const int* level,const std::uint64_t* index,
			const std::int8_t* axis,const Real* area,Real* balance,int count)
		{
			const int edge=blockIdx.x*blockDim.x+threadIdx.x;if(edge>=count)return;
			const Real flux=area[edge]*indexed_face_value(levels[level[edge]].fields,axis[edge],index[edge]);
			if(lower_node[edge]>=0)atomicAdd(balance+lower_node[edge],flux);
			if(upper_node[edge]>=0)atomicAdd(balance+upper_node[edge],-flux);
		}
		__global__ void scatter_regular_pressure_kernel(const DeviceCompositeAmrFluxLevelView* levels,int level,int bs,const unsigned char* active,const Real* pressure)
		{
			const DeviceCompositeAmrFluxLevelView source=levels[level];const DeviceAmrFieldLevelView fields=source.fields;const int cells_per_brick=bs*bs*bs,work=blockIdx.x*blockDim.x+threadIdx.x;if(work>=fields.brick_count*cells_per_brick)return;const int brick=work/cells_per_brick,local=work%cells_per_brick,a=source.level_offset+work;if(!active[a]||(fields.flags[brick]&BRICK_COVERED))return;const int i=local%bs,j=(local/bs)%bs,k=local/(bs*bs);fields.p[fields.layout.cell_index(brick,i,j,k)]=pressure[a];
		}
		__global__ void irregular_schwarz_lu_kernel(const int* block_offset,const int* factor_offset,
			const int* dof,const int* pivot,const Real* lu,const Real* residual,Real* scratch,
			Real* output,int block_count)
		{
			const int block=blockIdx.x*blockDim.x+threadIdx.x;if(block>=block_count)return;
			const int begin=block_offset[block],end=block_offset[block+1],count=end-begin;
			const int matrix=factor_offset[block];Real* rhs=scratch+begin;
			for(int row=0;row<count;++row)rhs[row]=residual[dof[begin+row]];
			for(int column=0;column<count;++column)
			{
				const int other=pivot[begin+column];if(other!=column){const Real swap=rhs[column];rhs[column]=rhs[other];rhs[other]=swap;}
			}
			for(int row=0;row<count;++row)
			{
				Real value=rhs[row];for(int column=0;column<row;++column)value-=lu[matrix+row*count+column]*rhs[column];rhs[row]=value;
			}
			for(int row=count-1;row>=0;--row)
			{
				Real value=rhs[row];for(int column=row+1;column<count;++column)value-=lu[matrix+row*count+column]*rhs[column];rhs[row]=value/lu[matrix+row*count+row];
			}
			for(int row=0;row<count;++row)output[dof[begin+row]]=rhs[row];
		}

		struct HostIrregularSchwarz
		{
			std::vector<int> block_offset{0},factor_offset{0},dof,pivot;
			std::vector<Real> lu;
			std::size_t rejected_blocks=0;
		};

		HostIrregularSchwarz build_irregular_schwarz(const CompositeAmrPressureSystem& system)
		{
			HostIrregularSchwarz result;if(!system.hierarchy||system.hierarchy->levels().empty())return result;
			const int base_offset=system.level_offset.front(),bs=system.brick_size;
			const int base_cells=static_cast<int>(system.hierarchy->levels().front().bricks.size())*bs*bs*bs;
			std::vector<unsigned char> selected(base_cells,0);
			auto select=[&](int pressure_dof)
			{
				if(pressure_dof<0||pressure_dof>=system.storage_size||!system.active[pressure_dof])return;
				const int aggregate=system.preconditioner_aggregate[pressure_dof]-base_offset;
				if(aggregate>=0&&aggregate<base_cells)selected[aggregate]=1;
			};
			auto select_gradient_node=[&](int node)
			{
				if(node<0||node>=static_cast<int>(system.pressure_gradient_dof.size()))return;
				select(system.pressure_gradient_dof[node]);
				for(int q=system.pressure_gradient_offset[node];q<system.pressure_gradient_offset[node+1];++q)
					select(system.pressure_gradient_neighbour[q]);
			};
			auto select_connection=[&](const CoarseFinePressureConnection& edge,bool always)
			{
				if(!always&&length2(edge.nonorthogonal_correction)<=1e-24)return;
				select(edge.coarse_dof);select(edge.fine_dof);
				select_gradient_node(edge.lower_gradient_node);select_gradient_node(edge.upper_gradient_node);
			};
			for(const CoarseFinePressureConnection& edge:system.coarse_fine)select_connection(edge,false);
			for(const CoarseFinePressureConnection& edge:system.embedded)select_connection(edge,true);
			for(const RegularPressureCorrection& edge:system.regular_pressure_corrections)
			{
				select(edge.lower_dof);select(edge.upper_dof);
				select_gradient_node(edge.lower_gradient_node);select_gradient_node(edge.upper_gradient_node);
			}

			std::vector<std::vector<int>> candidate_dof(base_cells);
			for(int q=0;q<system.storage_size;++q)if(system.active[q])
			{
				const int aggregate=system.preconditioner_aggregate[q]-base_offset;
				if(aggregate>=0&&aggregate<base_cells&&selected[aggregate])candidate_dof[aggregate].push_back(q);
			}
			std::vector<int> block_of(system.storage_size,-1),local_of(system.storage_size,-1),candidate_aggregate;
			for(int aggregate=0;aggregate<base_cells;++aggregate)if(!candidate_dof[aggregate].empty())
			{
				const int block=static_cast<int>(candidate_aggregate.size());candidate_aggregate.push_back(aggregate);
				for(int local=0;local<static_cast<int>(candidate_dof[aggregate].size());++local)
				{
					const int q=candidate_dof[aggregate][local];block_of[q]=block;local_of[q]=local;
				}
			}
			std::vector<std::vector<double>> matrix(candidate_aggregate.size());
			for(int block=0;block<static_cast<int>(candidate_aggregate.size());++block)
			{
				const std::size_t count=candidate_dof[candidate_aggregate[block]].size();matrix[block].assign(count*count,0);
			}
			auto add=[&](int row,int column,double value)
			{
				if(value==0||row<0||column<0||row>=system.storage_size||column>=system.storage_size)return;
				const int block=block_of[row];if(block<0||block_of[column]!=block)return;
				const int count=static_cast<int>(candidate_dof[candidate_aggregate[block]].size());
				matrix[block][static_cast<std::size_t>(local_of[row])*count+local_of[column]]+=value;
			};
			auto add_edge=[&](int a,int b,double coefficient)
			{
				if(a<0||b<0||a==b||a>=system.storage_size||b>=system.storage_size||!system.active[a]||!system.active[b]||coefficient==0)return;
				add(a,a,coefficient);add(a,b,-coefficient);add(b,b,coefficient);add(b,a,-coefficient);
			};
			for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level)
			{
				const AmrLevel& source=system.hierarchy->levels()[level];const double coefficient=source.h;
				for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick)
				{
					const BrickMetadata& meta=source.bricks[brick];if(!meta.active())continue;
					for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i)
					{
						const int a=system.dof(level,brick,i,j,k);if(!system.active[a])continue;
						for(int axis=0;axis<3;++axis)
						{
							if(system.cut_face_mask[a]&(1u<<axis))continue;int coordinate[3]={i,j,k},b=-1;
							if(++coordinate[axis]<bs)b=system.dof(level,brick,coordinate[0],coordinate[1],coordinate[2]);
							else {const int neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&source.bricks[neighbour].active()){coordinate[axis]=0;b=system.dof(level,neighbour,coordinate[0],coordinate[1],coordinate[2]);}}
							add_edge(a,b,coefficient);
						}
						if(system.pressure_outlet_xmax&&(meta.flags&BRICK_XMAX)&&i==bs-1)add(a,a,2.0*source.h);
					}
				}
			}
			for(const CoarseFinePressureConnection& edge:system.coarse_fine)add_edge(edge.coarse_dof,edge.fine_dof,edge.open_area*pressure_gradient_factor(edge));
			for(const CoarseFinePressureConnection& edge:system.embedded)add_edge(edge.coarse_dof,edge.fine_dof,edge.open_area*pressure_gradient_factor(edge));
			for(const RegularPressureCorrection& edge:system.regular_pressure_corrections)add_edge(edge.lower_dof,edge.upper_dof,edge.open_area*edge.two_point_delta);
			for(const CompositePressureGauge& gauge:system.gauges)if(gauge.dof>=0&&system.active[gauge.dof])add(gauge.dof,gauge.dof,gauge.coefficient);

			auto add_gradient_flux=[&](int lower,int upper,int lower_node,int upper_node,
				double upper_weight,double area,Vec3d correction)
			{
				if(lower_node<0||upper_node<0||length2(correction)<=1e-24)return;
				auto add_node=[&](int node,double interpolation)
				{
					const int centre=system.pressure_gradient_dof[node];
					for(int q=system.pressure_gradient_offset[node];q<system.pressure_gradient_offset[node+1];++q)
					{
						const int neighbour=system.pressure_gradient_neighbour[q];
						const double coefficient=-area*interpolation*dot(correction,system.pressure_gradient_weight[q]);
						add(lower,neighbour,coefficient);add(lower,centre,-coefficient);
						add(upper,neighbour,-coefficient);add(upper,centre,coefficient);
					}
				};
				add_node(lower_node,1-upper_weight);add_node(upper_node,upper_weight);
			};
			auto add_connection_correction=[&](const CoarseFinePressureConnection& edge)
			{
				const int lower=edge.direction>0?edge.coarse_dof:edge.fine_dof;
				const int upper=edge.direction>0?edge.fine_dof:edge.coarse_dof;
				add_gradient_flux(lower,upper,edge.lower_gradient_node,edge.upper_gradient_node,
					edge.upper_gradient_weight,edge.open_area,edge.nonorthogonal_correction);
			};
			for(const CoarseFinePressureConnection& edge:system.coarse_fine)add_connection_correction(edge);
			for(const CoarseFinePressureConnection& edge:system.embedded)add_connection_correction(edge);
			for(const RegularPressureCorrection& edge:system.regular_pressure_corrections)
				add_gradient_flux(edge.lower_dof,edge.upper_dof,edge.lower_gradient_node,edge.upper_gradient_node,
					edge.upper_gradient_weight,edge.open_area,edge.nonorthogonal_correction);

			for(int block=0;block<static_cast<int>(candidate_aggregate.size());++block)
			{
				std::vector<double> factor=std::move(matrix[block]);const int count=static_cast<int>(candidate_dof[candidate_aggregate[block]].size());
				double scale=0;for(double value:factor)scale=std::max(scale,std::abs(value));
				const double pivot_floor=64*std::numeric_limits<double>::epsilon()*
					std::max(scale,std::numeric_limits<double>::min())*count;
				std::vector<int> pivot(count);bool valid=scale>0;
				for(int column=0;column<count&&valid;++column)
				{
					int selected_row=column;double magnitude=std::abs(factor[static_cast<std::size_t>(column)*count+column]);
					for(int row=column+1;row<count;++row){const double candidate=std::abs(factor[static_cast<std::size_t>(row)*count+column]);if(candidate>magnitude){magnitude=candidate;selected_row=row;}}
					if(!(magnitude>pivot_floor)||!std::isfinite(magnitude)){valid=false;break;}pivot[column]=selected_row;
					if(selected_row!=column)for(int q=0;q<count;++q)std::swap(factor[static_cast<std::size_t>(column)*count+q],factor[static_cast<std::size_t>(selected_row)*count+q]);
					const double diagonal=factor[static_cast<std::size_t>(column)*count+column];
					for(int row=column+1;row<count;++row){double& multiplier=factor[static_cast<std::size_t>(row)*count+column];multiplier/=diagonal;for(int q=column+1;q<count;++q)factor[static_cast<std::size_t>(row)*count+q]-=multiplier*factor[static_cast<std::size_t>(column)*count+q];}
				}
				if(valid)for(double value:factor)if(!std::isfinite(value)||!std::isfinite(static_cast<Real>(value))){valid=false;break;}
				if(!valid){++result.rejected_blocks;continue;}
				const std::vector<int>& dofs=candidate_dof[candidate_aggregate[block]];
				result.dof.insert(result.dof.end(),dofs.begin(),dofs.end());result.pivot.insert(result.pivot.end(),pivot.begin(),pivot.end());
				for(double value:factor)result.lu.push_back(static_cast<Real>(value));
				result.block_offset.push_back(static_cast<int>(result.dof.size()));result.factor_offset.push_back(static_cast<int>(result.lu.size()));
			}
			return result;
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
			std::vector<int> coarse(connection_count_), fine(connection_count_); std::vector<Real> coefficient(connection_count_); for(int edge=0;edge<connection_count_;++edge){coarse[edge]=special[edge].coarse_dof;fine[edge]=special[edge].fine_dof;coefficient[edge]=static_cast<Real>(special[edge].open_area*pressure_gradient_factor(special[edge]));}
			check(cudaMalloc(&coarse_dof_,coarse.size()*sizeof(int)),"cudaMalloc coarse/fine coarse DOFs");check(cudaMalloc(&fine_dof_,fine.size()*sizeof(int)),"cudaMalloc coarse/fine fine DOFs");check(cudaMalloc(&coefficient_,coefficient.size()*sizeof(Real)),"cudaMalloc coarse/fine coefficients");check(cudaMemcpy(coarse_dof_,coarse.data(),coarse.size()*sizeof(int),cudaMemcpyHostToDevice),"upload coarse DOFs");check(cudaMemcpy(fine_dof_,fine.data(),fine.size()*sizeof(int),cudaMemcpyHostToDevice),"upload fine DOFs");check(cudaMemcpy(coefficient_,coefficient.data(),coefficient.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload coarse/fine coefficients");bytes_+=coarse.size()*sizeof(int)+fine.size()*sizeof(int)+coefficient.size()*sizeof(Real);
			std::vector<int> lower(connection_count_),upper(connection_count_),lower_node(connection_count_),upper_node(connection_count_);std::vector<Real> correction(static_cast<std::size_t>(connection_count_)*3),upper_weight(connection_count_),area(connection_count_);for(int edge=0;edge<connection_count_;++edge){lower[edge]=special[edge].direction>0?special[edge].coarse_dof:special[edge].fine_dof;upper[edge]=special[edge].direction>0?special[edge].fine_dof:special[edge].coarse_dof;lower_node[edge]=special[edge].lower_gradient_node;upper_node[edge]=special[edge].upper_gradient_node;for(int component=0;component<3;++component)correction[3*edge+component]=static_cast<Real>(special[edge].nonorthogonal_correction[component]);upper_weight[edge]=static_cast<Real>(special[edge].upper_gradient_weight);area[edge]=static_cast<Real>(special[edge].open_area);}connection_lower_node_=upload(lower_node,"upload pressure correction lower nodes");connection_upper_node_=upload(upper_node,"upload pressure correction upper nodes");connection_correction_=upload(correction,"upload nonorthogonal pressure corrections");connection_upper_weight_=upload(upper_weight,"upload pressure-gradient interpolation weights");connection_area_=upload(area,"upload pressure correction areas");
			// The correction kernel needs geometrically ordered endpoints; reuse the
			// base endpoint arrays only for its symmetric two-point contribution.
			cudaFree(coarse_dof_);cudaFree(fine_dof_);coarse_dof_=upload(lower,"upload ordered pressure lower DOFs");fine_dof_=upload(upper,"upload ordered pressure upper DOFs");bytes_+=lower_node.size()*2*sizeof(int)+correction.size()*sizeof(Real)+upper_weight.size()*sizeof(Real)+area.size()*sizeof(Real);
		}
		gradient_node_count_=static_cast<int>(system.pressure_gradient_dof.size());if(gradient_node_count_){std::vector<Real> weight(static_cast<std::size_t>(system.pressure_gradient_weight.size())*3);for(std::size_t q=0;q<system.pressure_gradient_weight.size();++q)for(int component=0;component<3;++component)weight[3*q+component]=static_cast<Real>(system.pressure_gradient_weight[q][component]);gradient_node_dof_=upload(system.pressure_gradient_dof,"upload compact pressure-gradient DOFs");gradient_offset_=upload(system.pressure_gradient_offset,"upload compact pressure-gradient offsets");gradient_neighbour_=upload(system.pressure_gradient_neighbour,"upload compact pressure-gradient neighbours");gradient_weight_=upload(weight,"upload compact pressure-gradient weights");gradient_=allocate_zero<Real>(static_cast<std::size_t>(gradient_node_count_)*3,"allocate compact pressure gradients");bytes_+=system.pressure_gradient_dof.size()*sizeof(int)+system.pressure_gradient_offset.size()*sizeof(int)+system.pressure_gradient_neighbour.size()*sizeof(int)+weight.size()*sizeof(Real)+static_cast<std::size_t>(gradient_node_count_)*3*sizeof(Real);}
		regular_correction_count_=static_cast<int>(system.regular_pressure_corrections.size());if(regular_correction_count_){std::vector<int> lower(regular_correction_count_),upper(regular_correction_count_),lower_node(regular_correction_count_),upper_node(regular_correction_count_);std::vector<Real> delta(regular_correction_count_),correction(static_cast<std::size_t>(regular_correction_count_)*3),weight(regular_correction_count_),area(regular_correction_count_);for(int edge=0;edge<regular_correction_count_;++edge){const auto& source=system.regular_pressure_corrections[edge];lower[edge]=source.lower_dof;upper[edge]=source.upper_dof;lower_node[edge]=source.lower_gradient_node;upper_node[edge]=source.upper_gradient_node;delta[edge]=static_cast<Real>(source.two_point_delta);weight[edge]=static_cast<Real>(source.upper_gradient_weight);area[edge]=static_cast<Real>(source.open_area);for(int component=0;component<3;++component)correction[3*edge+component]=static_cast<Real>(source.nonorthogonal_correction[component]);}regular_correction_lower_=upload(lower,"upload regular pressure-correction lower DOFs");regular_correction_upper_=upload(upper,"upload regular pressure-correction upper DOFs");regular_correction_lower_node_=upload(lower_node,"upload regular pressure-correction lower gradient nodes");regular_correction_upper_node_=upload(upper_node,"upload regular pressure-correction upper gradient nodes");regular_correction_two_point_delta_=upload(delta,"upload regular pressure two-point deltas");regular_correction_vector_=upload(correction,"upload regular nonorthogonal pressure corrections");regular_correction_upper_weight_=upload(weight,"upload regular pressure-gradient interpolation");regular_correction_area_=upload(area,"upload regular pressure-correction areas");bytes_+=static_cast<std::size_t>(4*regular_correction_count_)*sizeof(int)+static_cast<std::size_t>(6*regular_correction_count_)*sizeof(Real);}
		struct IrregularIncident { int dof=0,edge=0;std::int8_t sign=0; };
		std::vector<IrregularIncident> incidents;incidents.reserve(static_cast<std::size_t>(2)*(connection_count_+regular_correction_count_));
		auto add_incident_edge=[&](int lower,int upper,int encoded)
		{
			if(lower<0||upper<0||lower>=storage_size_||upper>=storage_size_||lower==upper||
				!system.active[lower]||!system.active[upper])return;
			incidents.push_back({lower,encoded,1});incidents.push_back({upper,encoded,-1});
		};
		for(int edge=0;edge<connection_count_;++edge)
		{
			const CoarseFinePressureConnection& source=special[edge];
			const int lower=source.direction>0?source.coarse_dof:source.fine_dof;
			const int upper=source.direction>0?source.fine_dof:source.coarse_dof;
			add_incident_edge(lower,upper,edge);
		}
		for(int edge=0;edge<regular_correction_count_;++edge)
		{
			const RegularPressureCorrection& source=system.regular_pressure_corrections[edge];
			add_incident_edge(source.lower_dof,source.upper_dof,connection_count_+edge);
		}
		std::sort(incidents.begin(),incidents.end(),[](const IrregularIncident& a,const IrregularIncident& b)
		{
			if(a.dof!=b.dof)return a.dof<b.dof;if(a.edge!=b.edge)return a.edge<b.edge;return a.sign<b.sign;
		});
		if(!incidents.empty())
		{
			std::vector<int> dof,offset,edge;std::vector<std::int8_t> sign;dof.reserve(incidents.size());
			offset.reserve(incidents.size()+1);edge.reserve(incidents.size());sign.reserve(incidents.size());
			for(std::size_t q=0;q<incidents.size();++q)
			{
				if(q==0||incidents[q].dof!=incidents[q-1].dof){dof.push_back(incidents[q].dof);offset.push_back(static_cast<int>(q));}
				edge.push_back(incidents[q].edge);sign.push_back(incidents[q].sign);
			}
			offset.push_back(static_cast<int>(incidents.size()));irregular_incident_dof_count_=static_cast<int>(dof.size());
			irregular_incident_dof_=upload(dof,"upload deterministic irregular pressure DOFs");
			irregular_incident_offset_=upload(offset,"upload deterministic irregular pressure offsets");
			irregular_incident_edge_=upload(edge,"upload deterministic irregular pressure edges");
			irregular_incident_sign_=upload(sign,"upload deterministic irregular pressure signs");
			bytes_+=dof.size()*sizeof(int)+offset.size()*sizeof(int)+edge.size()*sizeof(int)+sign.size()*sizeof(std::int8_t);
		}
		gauge_count_=static_cast<int>(system.gauges.size());if(gauge_count_){std::vector<int> dof(gauge_count_);std::vector<Real> coefficient(gauge_count_);for(int q=0;q<gauge_count_;++q){dof[q]=system.gauges[q].dof;coefficient[q]=static_cast<Real>(system.gauges[q].coefficient);}gauge_dof_=upload(dof,"upload composite pressure gauge DOFs");gauge_coefficient_=upload(coefficient,"upload composite pressure gauge coefficients");bytes_+=dof.size()*sizeof(int)+coefficient.size()*sizeof(Real);}
	}

	DeviceCompositeAmrPressureOperator::~DeviceCompositeAmrPressureOperator()
	{
		for (Allocation& allocation : allocations_) { if(allocation.neighbors)cudaFree(allocation.neighbors);if(allocation.flags)cudaFree(allocation.flags); }for(void* pointer:{(void*)levels_,(void*)coarse_dof_,(void*)fine_dof_,(void*)coefficient_,(void*)gradient_node_dof_,(void*)gradient_offset_,(void*)gradient_neighbour_,(void*)gradient_weight_,(void*)gradient_,(void*)connection_lower_node_,(void*)connection_upper_node_,(void*)connection_correction_,(void*)connection_upper_weight_,(void*)connection_area_,(void*)regular_correction_lower_,(void*)regular_correction_upper_,(void*)regular_correction_lower_node_,(void*)regular_correction_upper_node_,(void*)regular_correction_two_point_delta_,(void*)regular_correction_vector_,(void*)regular_correction_upper_weight_,(void*)regular_correction_area_,(void*)irregular_incident_dof_,(void*)irregular_incident_offset_,(void*)irregular_incident_edge_,(void*)irregular_incident_sign_,(void*)gauge_dof_,(void*)gauge_coefficient_,(void*)active_,(void*)cut_face_mask_})if(pointer)cudaFree(pointer);
	}

	void DeviceCompositeAmrPressureOperator::apply_orthogonal(const Real* pressure, Real* output) const
	{
		check(cudaMemset(output,0,static_cast<std::size_t>(storage_size_)*sizeof(Real)),"clear composite AMR output");const int cells_per_brick=brick_size_*brick_size_*brick_size_;
		for(int level=0;level<level_count_;++level){const int work=brick_counts_[level]*cells_per_brick;if(work)structured_apply_kernel<<<(work+255)/256,256>>>(levels_,level,brick_size_,outlet_,active_,cut_face_mask_,pressure,output);}if(irregular_incident_dof_count_)irregular_orthogonal_gather_kernel<<<(irregular_incident_dof_count_+255)/256,256>>>(irregular_incident_dof_,irregular_incident_offset_,irregular_incident_edge_,irregular_incident_sign_,irregular_incident_dof_count_,connection_count_,coarse_dof_,fine_dof_,coefficient_,regular_correction_lower_,regular_correction_upper_,regular_correction_two_point_delta_,regular_correction_area_,active_,pressure,output);if(gauge_count_)gauge_apply_kernel<<<(gauge_count_+255)/256,256>>>(gauge_dof_,gauge_coefficient_,gauge_count_,pressure,output);check(cudaDeviceSynchronize(),"apply orthogonal composite AMR pressure operator");
	}

	void DeviceCompositeAmrPressureOperator::apply(const Real* pressure, Real* output) const
	{
		apply_orthogonal(pressure,output);if(gradient_node_count_)reconstruct_pressure_gradient_kernel<<<(gradient_node_count_+255)/256,256>>>(gradient_node_dof_,gradient_offset_,gradient_neighbour_,gradient_weight_,pressure,gradient_,gradient_node_count_);if(irregular_incident_dof_count_&&gradient_node_count_)irregular_nonorthogonal_gather_kernel<<<(irregular_incident_dof_count_+255)/256,256>>>(irregular_incident_dof_,irregular_incident_offset_,irregular_incident_edge_,irregular_incident_sign_,irregular_incident_dof_count_,connection_count_,coarse_dof_,fine_dof_,connection_lower_node_,connection_upper_node_,connection_correction_,connection_upper_weight_,connection_area_,regular_correction_lower_,regular_correction_upper_,regular_correction_lower_node_,regular_correction_upper_node_,regular_correction_vector_,regular_correction_upper_weight_,regular_correction_area_,gradient_,active_,output);check(cudaDeviceSynchronize(),"apply nonorthogonal composite AMR pressure correction");
	}

	DeviceCompositeAmrPressureSolver::DeviceCompositeAmrPressureSolver(const CompositeAmrPressureSystem& system):op_(system),n_(system.storage_size)
	{
		if(std::getenv("PARACFD_PRESSURE_TRACE"))trace_system_=&system;
		if(!system.hierarchy||system.hierarchy->levels().empty()||system.preconditioner_aggregate.size()!=static_cast<std::size_t>(n_))throw std::invalid_argument("composite AMR solver geometric aggregation metadata");std::vector<double> diagonal_double;system.diagonal_cpu(diagonal_double);std::vector<Real> diagonal(diagonal_double.size());for(std::size_t q=0;q<diagonal.size();++q)diagonal[q]=static_cast<Real>(diagonal_double[q]);
		r_=allocate_zero<Real>(n_,"allocate AMR residual");z_=allocate_zero<Real>(n_,"allocate AMR preconditioned direction");direction_=allocate_zero<Real>(n_,"allocate AMR Krylov direction");Ad_=allocate_zero<Real>(n_,"allocate AMR operator direction");t_=allocate_zero<Real>(n_,"allocate AMR true-residual operator value");correction_=allocate_zero<Real>(n_,"allocate AMR GMRES cycle correction");defect_rhs_=allocate_zero<Real>(n_,"allocate AMR nonorthogonal defect RHS");gmres_v_=allocate_zero<Real>(static_cast<std::size_t>(gmres_restart_+1)*n_,"allocate AMR GMRES basis");gmres_z_=allocate_zero<Real>(static_cast<std::size_t>(gmres_restart_)*n_,"allocate AMR flexible GMRES basis");gmres_recycle_u_=allocate_zero<Real>(static_cast<std::size_t>(gmres_recycle_capacity_)*n_,"allocate AMR GCRO solution recycle basis");gmres_recycle_c_=allocate_zero<Real>(static_cast<std::size_t>(gmres_recycle_capacity_)*n_,"allocate AMR GCRO image recycle basis");diagonal_=upload(diagonal,"upload AMR base diagonal");active_=upload(system.active,"upload AMR active mask");aggregate_=upload(system.preconditioner_aggregate,"upload AMR coarse aggregate map");
		const AmrLevel& base=system.hierarchy->levels().front();brick_size_=system.brick_size;base_offset_=system.level_offset.front();base_bricks_=static_cast<int>(base.bricks.size());base_cells_=base_bricks_*brick_size_*brick_size_*brick_size_;const int cells_per_brick=brick_size_*brick_size_*brick_size_;std::vector<int> neighbors(static_cast<std::size_t>(base_bricks_)*6,-1);for(int brick=0;brick<base_bricks_;++brick)for(int face=0;face<6;++face)neighbors[brick*6+face]=base.bricks[brick].same_level_neighbor[face];std::vector<double> coarse_positive_double(static_cast<std::size_t>(base_cells_)*3,0),coarse_diagonal_double(base_cells_,0);
		std::vector<int> aggregate_restrict_offset(static_cast<std::size_t>(base_cells_)+1,0);
		for(int q=0;q<n_;++q)if(system.active[q])
		{
			const int coarse=system.preconditioner_aggregate[q]-base_offset_;
			if(coarse>=0&&coarse<base_cells_)++aggregate_restrict_offset[coarse+1];
		}
		for(int coarse=0;coarse<base_cells_;++coarse)aggregate_restrict_offset[coarse+1]+=aggregate_restrict_offset[coarse];
		std::vector<int> aggregate_restrict_dof(aggregate_restrict_offset.back()),aggregate_cursor=aggregate_restrict_offset;
		for(int q=0;q<n_;++q)if(system.active[q])
		{
			const int coarse=system.preconditioner_aggregate[q]-base_offset_;
			if(coarse>=0&&coarse<base_cells_)aggregate_restrict_dof[aggregate_cursor[coarse]++]=q;
		}
		const HostIrregularSchwarz schwarz=build_irregular_schwarz(system);
		schwarz_block_count_=static_cast<int>(schwarz.block_offset.size())-1;
		schwarz_dof_count_=static_cast<int>(schwarz.dof.size());
		for(int block=0;block<schwarz_block_count_;++block)schwarz_max_block_size_=std::max(
			schwarz_max_block_size_,schwarz.block_offset[block+1]-schwarz.block_offset[block]);
		schwarz_block_offset_=upload(schwarz.block_offset,"upload irregular Schwarz block offsets");
		schwarz_factor_offset_=upload(schwarz.factor_offset,"upload irregular Schwarz factor offsets");
		schwarz_dof_=upload(schwarz.dof,"upload irregular Schwarz DOFs");
		schwarz_pivot_=upload(schwarz.pivot,"upload irregular Schwarz pivots");
		schwarz_lu_=upload(schwarz.lu,"upload irregular Schwarz LU factors");
		schwarz_rhs_=allocate_zero<Real>(schwarz.dof.size(),"allocate irregular Schwarz RHS");
		if(trace_system_)
		{
			std::array<std::size_t,6> histogram{};
			for(int block=0;block<schwarz_block_count_;++block)
			{
				const int count=schwarz.block_offset[block+1]-schwarz.block_offset[block];
				const int bucket=count<=32?0:(count<=64?1:(count<=96?2:(count<=128?3:(count<=192?4:5))));++histogram[bucket];
			}
			std::fprintf(stderr,"[pressure-schwarz] factored aggregates=%d rejected=%zu dofs=%zu dense-entries=%zu (%.3f MiB FP32) max-size=%d histogram <=32/64/96/128/192/>192=[%zu/%zu/%zu/%zu/%zu/%zu]\n",
				schwarz_block_count_,schwarz.rejected_blocks,schwarz.dof.size(),schwarz.lu.size(),
				schwarz.lu.size()*sizeof(Real)/(1024.0*1024.0),schwarz_max_block_size_,
				histogram[0],histogram[1],histogram[2],histogram[3],histogram[4],histogram[5]);
		}
		aggregate_restrict_offset_=upload(aggregate_restrict_offset,"upload deterministic aggregate offsets");
		aggregate_restrict_dof_=upload(aggregate_restrict_dof,"upload deterministic aggregate DOFs");
		std::vector<int> base_extra_a,base_extra_b;std::vector<double> base_extra_coefficient_double;auto global_coord=[&](int coarse){const int brick=coarse/cells_per_brick,local=coarse%cells_per_brick;return std::array<int,3>{base.bricks[brick].coord.x*brick_size_+local%brick_size_,base.bricks[brick].coord.y*brick_size_+(local/brick_size_)%brick_size_,base.bricks[brick].coord.z*brick_size_+local/(brick_size_*brick_size_)};};auto add_coarse_edge=[&](int a,int b,double coefficient){if(a<0||b<0||a==b||!system.active[a]||!system.active[b])return;const int ca=system.preconditioner_aggregate[a]-base_offset_,cb=system.preconditioner_aggregate[b]-base_offset_;if(ca<0||cb<0||ca>=base_cells_||cb>=base_cells_||ca==cb)return;const auto ga=global_coord(ca),gb=global_coord(cb);int axis=-1,manhattan=0;for(int d=0;d<3;++d){const int delta=gb[d]-ga[d];manhattan+=std::abs(delta);if(delta)axis=d;}if(manhattan==1&&axis>=0){const int lower=ga[axis]<gb[axis]?ca:cb;coarse_positive_double[static_cast<std::size_t>(lower)*3+axis]+=coefficient;}else{base_extra_a.push_back(ca);base_extra_b.push_back(cb);base_extra_coefficient_double.push_back(coefficient);}coarse_diagonal_double[ca]+=coefficient;coarse_diagonal_double[cb]+=coefficient;};
		for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level){const AmrLevel& source=system.hierarchy->levels()[level];const double coefficient=source.h;for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick){const BrickMetadata& meta=source.bricks[brick];if(!meta.active())continue;for(int k=0;k<brick_size_;++k)for(int j=0;j<brick_size_;++j)for(int i=0;i<brick_size_;++i){const int a=system.dof(level,brick,i,j,k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if(system.cut_face_mask[a]&(1u<<axis))continue;int c[3]={i,j,k},b=-1;if(++c[axis]<brick_size_)b=system.dof(level,brick,c[0],c[1],c[2]);else{const int neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&source.bricks[neighbour].active()){c[axis]=0;b=system.dof(level,neighbour,c[0],c[1],c[2]);}}add_coarse_edge(a,b,coefficient);}if(system.pressure_outlet_xmax&&(meta.flags&BRICK_XMAX)&&i==brick_size_-1){const int coarse=system.preconditioner_aggregate[a]-base_offset_;if(coarse>=0&&coarse<base_cells_)coarse_diagonal_double[coarse]+=2.0*source.h;}}}}
		for(const CoarseFinePressureConnection& edge:system.coarse_fine)add_coarse_edge(edge.coarse_dof,edge.fine_dof,edge.open_area*pressure_gradient_factor(edge));for(const CoarseFinePressureConnection& edge:system.embedded)add_coarse_edge(edge.coarse_dof,edge.fine_dof,edge.open_area*pressure_gradient_factor(edge));for(const RegularPressureCorrection& edge:system.regular_pressure_corrections)add_coarse_edge(edge.lower_dof,edge.upper_dof,edge.open_area*edge.two_point_delta);for(const CompositePressureGauge& gauge:system.gauges){const int coarse=system.preconditioner_aggregate[gauge.dof]-base_offset_;if(coarse>=0&&coarse<base_cells_)coarse_diagonal_double[coarse]+=gauge.coefficient;}
		struct BaseExtraIncident { int dof=0,edge=0; };
		std::vector<BaseExtraIncident> base_extra_incidents;base_extra_incidents.reserve(2*base_extra_a.size());
		for(int edge=0;edge<static_cast<int>(base_extra_a.size());++edge)
		{
			base_extra_incidents.push_back({base_extra_a[edge],edge});
			base_extra_incidents.push_back({base_extra_b[edge],edge});
		}
		std::sort(base_extra_incidents.begin(),base_extra_incidents.end(),[](const BaseExtraIncident& a,const BaseExtraIncident& b)
		{
			if(a.dof!=b.dof)return a.dof<b.dof;return a.edge<b.edge;
		});
		std::vector<int> base_extra_incident_dof,base_extra_incident_offset,base_extra_incident_edge;
		for(std::size_t q=0;q<base_extra_incidents.size();++q)
		{
			if(q==0||base_extra_incidents[q].dof!=base_extra_incidents[q-1].dof)
			{
				base_extra_incident_dof.push_back(base_extra_incidents[q].dof);
				base_extra_incident_offset.push_back(static_cast<int>(q));
			}
			base_extra_incident_edge.push_back(base_extra_incidents[q].edge);
		}
		if(!base_extra_incidents.empty())base_extra_incident_offset.push_back(static_cast<int>(base_extra_incidents.size()));
		base_extra_incident_dof_count_=static_cast<int>(base_extra_incident_dof.size());
		base_extra_incident_dof_=upload(base_extra_incident_dof,"upload deterministic base-extra DOFs");
		base_extra_incident_offset_=upload(base_extra_incident_offset,"upload deterministic base-extra offsets");
		base_extra_incident_edge_=upload(base_extra_incident_edge,"upload deterministic base-extra edges");
		std::vector<Real> base_positive(coarse_positive_double.size()),base_diagonal(base_cells_),base_extra_coefficient(base_extra_coefficient_double.size());for(std::size_t q=0;q<base_positive.size();++q)base_positive[q]=static_cast<Real>(coarse_positive_double[q]);for(int q=0;q<base_cells_;++q)base_diagonal[q]=static_cast<Real>(coarse_diagonal_double[q]);for(std::size_t q=0;q<base_extra_coefficient.size();++q)base_extra_coefficient[q]=static_cast<Real>(base_extra_coefficient_double[q]);base_extra_count_=static_cast<int>(base_extra_a.size());base_neighbors_=upload(neighbors,"upload AMR base neighbors");base_positive_=upload(base_positive,"upload AMR Galerkin base faces");base_diagonal_=upload(base_diagonal,"upload AMR base diagonal");base_extra_a_=upload(base_extra_a,"upload AMR nonlocal aggregate edge A");base_extra_b_=upload(base_extra_b,"upload AMR nonlocal aggregate edge B");base_extra_coefficient_=upload(base_extra_coefficient,"upload AMR nonlocal aggregate coefficient");base_rhs_=allocate_zero<Real>(base_cells_,"allocate AMR base RHS");base_x_=allocate_zero<Real>(base_cells_,"allocate AMR base correction");base_tmp_=allocate_zero<Real>(base_cells_,"allocate AMR base temporary");const std::size_t full_vectors=static_cast<std::size_t>(2*gmres_restart_+2*gmres_recycle_capacity_+9);bytes_=op_.bytes()+full_vectors*n_*sizeof(Real)+system.active.size()+system.preconditioner_aggregate.size()*sizeof(int)+neighbors.size()*sizeof(int)+base_positive.size()*sizeof(Real)+static_cast<std::size_t>(4)*base_cells_*sizeof(Real)+base_extra_a.size()*static_cast<std::size_t>(2*sizeof(int)+sizeof(Real))+aggregate_restrict_offset.size()*sizeof(int)+aggregate_restrict_dof.size()*sizeof(int)+base_extra_incident_dof.size()*sizeof(int)+base_extra_incident_offset.size()*sizeof(int)+base_extra_incident_edge.size()*sizeof(int)+schwarz.block_offset.size()*sizeof(int)+schwarz.factor_offset.size()*sizeof(int)+schwarz.dof.size()*static_cast<std::size_t>(2*sizeof(int)+sizeof(Real))+schwarz.lu.size()*sizeof(Real);
	}
	DeviceCompositeAmrPressureSolver::~DeviceCompositeAmrPressureSolver(){for(void* pointer:{(void*)r_,(void*)z_,(void*)direction_,(void*)Ad_,(void*)t_,(void*)correction_,(void*)defect_rhs_,(void*)gmres_v_,(void*)gmres_z_,(void*)gmres_recycle_u_,(void*)gmres_recycle_c_,(void*)diagonal_,(void*)active_,(void*)aggregate_,(void*)aggregate_restrict_offset_,(void*)aggregate_restrict_dof_,(void*)base_neighbors_,(void*)base_extra_a_,(void*)base_extra_b_,(void*)base_extra_incident_dof_,(void*)base_extra_incident_offset_,(void*)base_extra_incident_edge_,(void*)base_positive_,(void*)base_diagonal_,(void*)base_extra_coefficient_,(void*)base_rhs_,(void*)base_x_,(void*)base_tmp_,(void*)schwarz_block_offset_,(void*)schwarz_factor_offset_,(void*)schwarz_dof_,(void*)schwarz_pivot_,(void*)schwarz_lu_,(void*)schwarz_rhs_})if(pointer)cudaFree(pointer);}
	void DeviceCompositeAmrPressureSolver::apply_preconditioner(const Real* residual,Real* output)
	{
		precondition_kernel<<<(n_+255)/256,256>>>(residual,output,diagonal_,active_,n_);check(cudaMemset(base_x_,0,static_cast<std::size_t>(base_cells_)*sizeof(Real)),"clear AMR base correction");restrict_aggregate_gather_kernel<<<(base_cells_+255)/256,256>>>(aggregate_restrict_offset_,aggregate_restrict_dof_,residual,base_cells_,base_rhs_);Real* current=base_x_;Real* next=base_tmp_;constexpr int sweeps=16;for(int sweep=0;sweep<sweeps;++sweep){coarse_jacobi_kernel<<<(base_cells_+255)/256,256>>>(base_rhs_,current,next,base_positive_,base_diagonal_,base_neighbors_,base_bricks_,brick_size_,Real(0.7));if(base_extra_incident_dof_count_)coarse_extra_jacobi_gather_kernel<<<(base_extra_incident_dof_count_+255)/256,256>>>(base_extra_incident_dof_,base_extra_incident_offset_,base_extra_incident_edge_,base_extra_incident_dof_count_,base_extra_a_,base_extra_b_,base_extra_coefficient_,current,base_diagonal_,Real(0.7),next);Real* swap=current;current=next;next=swap;}prolong_aggregate_kernel<<<(n_+255)/256,256>>>(aggregate_,active_,n_,base_offset_,base_cells_,current,output);
	}
	void DeviceCompositeAmrPressureSolver::apply_irregular_schwarz(const Real* residual,Real* output)
	{
		check(cudaMemset(output,0,static_cast<std::size_t>(n_)*sizeof(Real)),"clear irregular Schwarz correction");
		if(schwarz_block_count_)irregular_schwarz_lu_kernel<<<(schwarz_block_count_+127)/128,128>>>(
			schwarz_block_offset_,schwarz_factor_offset_,schwarz_dof_,schwarz_pivot_,schwarz_lu_,
			residual,schwarz_rhs_,output,schwarz_block_count_);
	}
	void DeviceCompositeAmrPressureSolver::trace_failure_residual(const Real* residual) const
	{
		if(!trace_system_)return;
		const CompositeAmrPressureSystem& system=*trace_system_;
		std::vector<Real> host(n_);check(cudaMemcpy(host.data(),residual,
			static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToHost),
			"download failed pressure residual diagnostic");
		std::vector<unsigned char> eb(n_,0),gauge(n_,0);
		for(const CoarseFinePressureConnection& edge:system.embedded)
		{
			if(edge.coarse_dof>=0&&edge.coarse_dof<n_)eb[edge.coarse_dof]=1;
			if(edge.fine_dof>=0&&edge.fine_dof<n_)eb[edge.fine_dof]=1;
		}
		for(const CompositeSurfacePressurePatch& patch:system.surface_patches)
		{
			if(patch.plus_dof>=0&&patch.plus_dof<n_)eb[patch.plus_dof]=1;
			if(patch.minus_dof>=0&&patch.minus_dof<n_)eb[patch.minus_dof]=1;
		}
		for(const CompositePressureGauge& value:system.gauges)
			if(value.dof>=0&&value.dof<n_)gauge[value.dof]=1;

		// Reconstruct exactly the same open pressure graph used when components and
		// gauges are finalized. This runs only after an explicitly traced failure.
		std::vector<int> parent(n_,-1);for(int q=0;q<n_;++q)if(system.active[q])parent[q]=q;
		auto root=[&](int q)
		{
			while(parent[q]!=q){parent[q]=parent[parent[q]];q=parent[q];}return q;
		};
		auto join=[&](int a,int b)
		{
			if(a<0||b<0||a>=n_||b>=n_||!system.active[a]||!system.active[b])return;
			a=root(a);b=root(b);if(a!=b)parent[std::max(a,b)]=std::min(a,b);
		};
		const int bs=system.brick_size;
		for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level)
		{
			const AmrLevel& source=system.hierarchy->levels()[level];
			for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick)
			{
				const BrickMetadata& meta=source.bricks[brick];if(!meta.active())continue;
				for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i)
				{
					const int a=system.dof(level,brick,i,j,k);if(!system.active[a])continue;
					for(int axis=0;axis<3;++axis)
					{
						if(system.cut_face_mask[a]&(1u<<axis))continue;
						int c[3]={i,j,k},b=-1;if(++c[axis]<bs)b=system.dof(level,brick,c[0],c[1],c[2]);
						else {const int neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&source.bricks[neighbour].active()){c[axis]=0;b=system.dof(level,neighbour,c[0],c[1],c[2]);}}
						join(a,b);
					}
				}
			}
		}
		for(const CoarseFinePressureConnection& edge:system.coarse_fine)join(edge.coarse_dof,edge.fine_dof);
		for(const CoarseFinePressureConnection& edge:system.embedded)join(edge.coarse_dof,edge.fine_dof);

		struct ComponentStats
		{
			double residual2=0,maximum=0,volume=0;
			std::size_t count=0,eb_count=0,gauge_count=0;
			bool freestream=false;
		};
		std::vector<int> root_component(n_,-1),component_of(n_,-1),top;
		std::vector<ComponentStats> components;
		double total2=0,regular2=0,eb2=0,gauge2=0,external2=0,sealed2=0;
		std::size_t active_count=0;
		for(int q=0;q<n_;++q)if(system.active[q])
		{
			const int r=root(q);if(root_component[r]<0){root_component[r]=static_cast<int>(components.size());components.push_back({});}
			const int component=root_component[r];component_of[q]=component;ComponentStats& stats=components[component];
			const double value=static_cast<double>(host[q]),square=value*value,absolute=std::abs(value);
			total2+=square;if(eb[q])eb2+=square;else regular2+=square;if(gauge[q])gauge2+=square;
			const bool external=q<static_cast<int>(system.freestream_connected.size())&&system.freestream_connected[q];
			if(external)external2+=square;else sealed2+=square;
			stats.residual2+=square;stats.maximum=std::max(stats.maximum,absolute);stats.volume+=system.volume[q];
			++stats.count;stats.eb_count+=eb[q]!=0;stats.gauge_count+=gauge[q]!=0;stats.freestream=stats.freestream||external;++active_count;
			if(top.size()<12)top.push_back(q);else
			{
				auto smallest=std::min_element(top.begin(),top.end(),[&](int a,int b){return std::abs(static_cast<double>(host[a]))<std::abs(static_cast<double>(host[b]));});
				if(absolute>std::abs(static_cast<double>(host[*smallest]))) *smallest=q;
			}
		}
		std::sort(top.begin(),top.end(),[&](int a,int b){return std::abs(static_cast<double>(host[a]))>std::abs(static_cast<double>(host[b]));});
		std::vector<int> component_order(components.size());for(int q=0;q<static_cast<int>(components.size());++q)component_order[q]=q;
		std::sort(component_order.begin(),component_order.end(),[&](int a,int b){return components[a].residual2>components[b].residual2;});
		auto share=[&](double value){return total2>0?100*value/total2:0;};
		std::fprintf(stderr,"[pressure-residual] active=%zu components=%zu L2=%.9e max=%.9e regular-L2=%.9e/share2=%.3f%% EB-incident-L2=%.9e/share2=%.3f%% gauge-L2=%.9e/share2=%.3f%% external/share2=%.3f%% sealed/share2=%.3f%%\n",
			active_count,components.size(),std::sqrt(total2),top.empty()?0:std::abs(static_cast<double>(host[top.front()])),
			std::sqrt(regular2),share(regular2),std::sqrt(eb2),share(eb2),std::sqrt(gauge2),share(gauge2),share(external2),share(sealed2));
		for(int order=0;order<std::min<int>(16,component_order.size());++order)
		{
			const int component=component_order[order];const ComponentStats& stats=components[component];
			std::fprintf(stderr,"[pressure-residual] component=%d count=%zu volume=%.9e L2=%.9e/share2=%.3f%% max=%.9e EB=%zu gauges=%zu freestream=%d\n",
				component,stats.count,stats.volume,std::sqrt(stats.residual2),share(stats.residual2),stats.maximum,
				stats.eb_count,stats.gauge_count,stats.freestream?1:0);
		}
		for(int rank=0;rank<static_cast<int>(top.size());++rank)
		{
			const int q=top[rank];const Vec3d centre=system.centroid[q];
			std::fprintf(stderr,"[pressure-residual] top=%d dof=%d r=%+.9e |r|/V=%.9e V=%.9e centroid=[%.9g %.9g %.9g] kind=%s gauge=%d freestream=%d component=%d\n",
				rank+1,q,static_cast<double>(host[q]),system.volume[q]>0?std::abs(static_cast<double>(host[q]))/system.volume[q]:std::numeric_limits<double>::infinity(),
				system.volume[q],centre.x,centre.y,centre.z,eb[q]?"EB-incident":"regular",gauge[q]?1:0,
				(q<static_cast<int>(system.freestream_connected.size())&&system.freestream_connected[q])?1:0,component_of[q]);
		}
	}
	AmrGpuSolveResult DeviceCompositeAmrPressureSolver::solve_orthogonal(Real* pressure,const Real* rhs,double tolerance,int max_iterations,bool warm_start)
	{
		AmrGpuSolveResult result;
		if(!warm_start)check(cudaMemset(pressure,0,static_cast<std::size_t>(n_)*sizeof(Real)),"clear composite AMR pressure");
		const double rhs2=dot_active(rhs,rhs,active_,n_);
		if(rhs2==0)
		{
			check(cudaMemset(pressure,0,static_cast<std::size_t>(n_)*sizeof(Real)),"clear zero-RHS composite AMR pressure");
			result.relative_residual=0;result.converged=true;return result;
		}
		if(!std::isfinite(rhs2)||rhs2<0){result.relative_residual=std::numeric_limits<double>::infinity();return result;}
		auto true_residual=[&]()
		{
			op_.apply_orthogonal(pressure,t_);
			check(cudaMemcpy(r_,rhs,static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),"reset true composite residual");
			subtract_kernel<<<(n_+255)/256,256>>>(r_,t_,active_,n_);
			return dot_active(r_,r_,active_,n_);
		};
		double residual2;
		if(warm_start)residual2=true_residual();
		else
		{
			check(cudaMemcpy(r_,rhs,static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),"initial orthogonal AMR residual");
			residual2=rhs2;
		}
		result.relative_residual=std::sqrt(residual2/rhs2);if(!std::isfinite(result.relative_residual))return result;if(result.relative_residual<=tolerance){result.converged=true;return result;}

		// The orthogonal finite-volume operator is SPD.  Keep its inner inverse on
		// PCG: the full affine-corrected operator is handled by the outer flexible
		// GMRES, so exposing PCG to that nonsymmetric correction is neither required
		// nor mathematically valid.  Periodic true-residual restarts bound FP32
		// recurrence drift without interpreting a small positive p.A.p as breakdown.
		auto restart=[&]()->double
		{
			apply_preconditioner(r_,z_);
			copy_active_kernel<<<(n_+255)/256,256>>>(z_,direction_,active_,n_);
			return dot_active(r_,z_,active_,n_);
		};
		double rz=restart();if(!(rz>0)||!std::isfinite(rz))return result;
		constexpr int reliable_update_interval=32;
		for(int iteration=0;iteration<max_iterations;++iteration)
		{
			op_.apply_orthogonal(direction_,Ad_);const double pAp=dot_active(direction_,Ad_,active_,n_);if(!(pAp>0)||!std::isfinite(pAp))break;const double alpha=rz/pAp;if(!std::isfinite(alpha)||!std::isfinite(static_cast<Real>(alpha)))break;update_pressure_residual_kernel<<<(n_+255)/256,256>>>(pressure,r_,direction_,Ad_,static_cast<Real>(alpha),active_,n_);
			residual2=dot_active(r_,r_,active_,n_);
			result.iterations=iteration+1;result.relative_residual=std::sqrt(residual2/rhs2);if(!std::isfinite(result.relative_residual))break;const bool reliable=(iteration+1)%reliable_update_interval==0;
			if(result.relative_residual<=tolerance||reliable)
			{
				residual2=true_residual();result.relative_residual=std::sqrt(residual2/rhs2);
				if(result.relative_residual<=tolerance){result.converged=true;break;}
				rz=restart();if(!(rz>0)||!std::isfinite(rz))break;continue;
			}
			apply_preconditioner(r_,z_);const double next_rz=dot_active(r_,z_,active_,n_);if(!(next_rz>0)||!std::isfinite(next_rz))break;const double beta=next_rz/rz;if(!std::isfinite(beta)||!std::isfinite(static_cast<Real>(beta)))break;update_direction_kernel<<<(n_+255)/256,256>>>(z_,direction_,static_cast<Real>(beta),active_,n_);rz=next_rz;
		}
		if(!result.converged)
		{
			residual2=true_residual();result.relative_residual=std::sqrt(residual2/rhs2);
			result.converged=std::isfinite(result.relative_residual)&&result.relative_residual<=tolerance;
		}
		check(cudaDeviceSynchronize(),"solve orthogonal composite AMR pressure with PCG");return result;
	}

	AmrGpuSolveResult DeviceCompositeAmrPressureSolver::solve(Real* pressure,const Real* rhs,
		double tolerance,int max_iterations,bool warm_start)
	{
		// Flexible right-preconditioned GCRO/LGMRES. U contains realized solution-space
		// corrections and C=A*U contains their orthonormal images. Keeping both is what
		// makes recycling valid even though the inner A0 solve is a variable/nonlinear
		// right preconditioner.
		AmrGpuSolveResult result;result.irregular_schwarz_block_count=schwarz_block_count_;
		result.irregular_schwarz_max_block_size=schwarz_max_block_size_;
		const bool trace=std::getenv("PARACFD_PRESSURE_TRACE")!=nullptr;
		const bool recycle_enabled=std::getenv("PARACFD_PRESSURE_DISABLE_RECYCLE")==nullptr;
		int restart_index=0;
		if(!recycle_enabled)gmres_recycle_count_=0;
		if(!warm_start)
		{
			gmres_recycle_count_=0;
			if(max_iterations>0)check(cudaMemset(pressure,0,static_cast<std::size_t>(n_)*sizeof(Real)),
				"clear GCRO composite pressure");
		}
		if(max_iterations<=0)return result;
		const double rhs2=dot_active(rhs,rhs,active_,n_);
		if(rhs2==0)
		{
			check(cudaMemset(pressure,0,static_cast<std::size_t>(n_)*sizeof(Real)),
				"clear zero-RHS GCRO pressure");
			result.converged=true;
			return result;
		}
		if(!std::isfinite(rhs2)||rhs2<0)
		{
			result.relative_residual=std::numeric_limits<double>::infinity();
			return result;
		}
		auto true_residual=[&]()
		{
			op_.apply(pressure,Ad_);
			check(cudaMemcpy(defect_rhs_,rhs,static_cast<std::size_t>(n_)*sizeof(Real),
				cudaMemcpyDeviceToDevice),"copy GCRO true residual");
			subtract_kernel<<<(n_+255)/256,256>>>(defect_rhs_,Ad_,active_,n_);
			return dot_active(defect_rhs_,defect_rhs_,active_,n_);
		};
		auto finite_real=[](double value)
		{
			return std::isfinite(value)&&std::isfinite(static_cast<Real>(value));
		};
		auto recycle_u=[&](int column)
		{
			return gmres_recycle_u_+static_cast<std::size_t>(column)*n_;
		};
		auto recycle_c=[&](int column)
		{
			return gmres_recycle_c_+static_cast<std::size_t>(column)*n_;
		};

		double residual2=true_residual();
		result.relative_residual=std::sqrt(residual2/rhs2);
		if(result.relative_residual<=tolerance){result.converged=true;return result;}
		while(result.iterations<max_iterations)
		{
			++restart_index;
			// Minimize over the retained space first. The paired U/C updates preserve
			// r=b-Ax in exact arithmetic.  In FP32, however, adding a small U update to
			// an already large pressure can round differently from subtracting the paired
			// C update from the residual.  Always measure that consistency against the
			// true operator before trusting the recycled residual.  A pair whose gap is
			// not comfortably below the requested solve accuracy is discarded; the true
			// residual remains the starting vector for the ordinary flexible cycle.
			const double residual_before_recycle2=residual2;
			const int recycle_attempted=gmres_recycle_count_;
			bool recycle_valid=true;
			if(recycle_attempted)
			{
				copy_active_kernel<<<(n_+255)/256,256>>>(pressure,correction_,active_,n_);
				copy_active_kernel<<<(n_+255)/256,256>>>(defect_rhs_,t_,active_,n_);
			}
			for(int pass=0;pass<2&&recycle_valid;++pass)
				for(int column=0;column<gmres_recycle_count_;++column)
				{
					const double value=dot_active(recycle_c(column),defect_rhs_,active_,n_);
					if(!finite_real(value)){recycle_valid=false;break;}
					axpy_active_kernel<<<(n_+255)/256,256>>>(correction_,recycle_u(column),
						static_cast<Real>(value),active_,n_);
					axpy_active_kernel<<<(n_+255)/256,256>>>(defect_rhs_,recycle_c(column),
						static_cast<Real>(-value),active_,n_);
				}
			const double predicted_residual2=dot_active(defect_rhs_,defect_rhs_,active_,n_);
			double recycle_gap2=0,true_recycled2=residual_before_recycle2;
			if(recycle_attempted&&recycle_valid)
			{
				copy_active_kernel<<<(n_+255)/256,256>>>(defect_rhs_,r_,active_,n_);
				op_.apply(correction_,Ad_);
				check(cudaMemcpy(defect_rhs_,rhs,static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),"copy transactional recycled true residual");
				subtract_kernel<<<(n_+255)/256,256>>>(defect_rhs_,Ad_,active_,n_);
				true_recycled2=dot_active(defect_rhs_,defect_rhs_,active_,n_);
				axpy_active_kernel<<<(n_+255)/256,256>>>(r_,defect_rhs_,Real(-1),active_,n_);
				recycle_gap2=dot_active(r_,r_,active_,n_);
			}
			const double requested_gap=Real(0.25)*tolerance*std::sqrt(rhs2);
			const bool recycle_accurate=!recycle_attempted||
				(std::isfinite(recycle_gap2)&&std::sqrt(recycle_gap2)<=requested_gap&&
				 true_recycled2<residual_before_recycle2);
			if(recycle_attempted&&recycle_valid&&recycle_accurate)
			{
				copy_active_kernel<<<(n_+255)/256,256>>>(correction_,pressure,active_,n_);
				residual2=true_recycled2;
			}
			else if(recycle_attempted)
			{
				gmres_recycle_count_=0;copy_active_kernel<<<(n_+255)/256,256>>>(t_,defect_rhs_,active_,n_);
				residual2=residual_before_recycle2;
			}
			if(recycle_attempted&&(!recycle_valid||!recycle_accurate))++result.recycle_rejections;
			result.relative_residual=std::sqrt(residual2/rhs2);
			if(trace&&(recycle_attempted||!recycle_valid||!recycle_accurate))
			{
				std::fprintf(stderr,"[pressure-gcro] restart=%d recycled residual-pair gap=%.9e relative-to-rhs=%.9e predicted/true=%.9e/%.9e retained=%d\n",
					restart_index,true_recycled2>0?std::sqrt(recycle_gap2/true_recycled2):0,
					std::sqrt(recycle_gap2/rhs2),std::sqrt(predicted_residual2/rhs2),
					std::sqrt(true_recycled2/rhs2),gmres_recycle_count_);
			}
			if(trace)std::fprintf(stderr,"[pressure-gcro] restart=%d begin-it=%d recycle=%d residual=%.9e\n",
				restart_index,result.iterations,gmres_recycle_count_,result.relative_residual);
			if(!std::isfinite(result.relative_residual))break;
			if(result.relative_residual<=tolerance)
			{
				residual2=true_residual();result.relative_residual=std::sqrt(residual2/rhs2);
				if(result.relative_residual<=tolerance){result.converged=true;break;}
			}
			const double beta=std::sqrt(residual2);
			if(!(beta>0)||!std::isfinite(beta)||!finite_real(1.0/beta))break;
			const int cycle=std::min(gmres_restart_,max_iterations-result.iterations);
			std::vector<double> h(static_cast<std::size_t>(cycle+1)*cycle,0),
				cosine(cycle),sine(cycle),g(cycle+1),y(cycle);
			g[0]=beta;
			scaled_copy_active_kernel<<<(n_+255)/256,256>>>(gmres_v_,defect_rhs_,
				static_cast<Real>(1.0/beta),active_,n_);
			int used=0;
			for(int column=0;column<cycle;++column)
			{
				Real* const v=gmres_v_+static_cast<std::size_t>(column)*n_;
				Real* const z=gmres_z_+static_cast<std::size_t>(column)*n_;
				// The orthogonal finite-volume operator is the coercive part of the full
				// non-orthogonal system.  Approximately invert it for each outer basis
				// vector; retaining z_j makes this a flexible right-preconditioner even
				// when the inner Krylov iteration takes a different path per column.
				auto record_inner=[&](const AmrGpuSolveResult& inner)
				{
					++result.preconditioner_applications;
					result.inner_iterations_total+=inner.iterations;
					result.inner_iterations_maximum=std::max(result.inner_iterations_maximum,inner.iterations);
					result.inner_relative_residual_maximum=std::max(
						result.inner_relative_residual_maximum,inner.relative_residual);
				};
				const AmrGpuSolveResult inner=solve_orthogonal(z,v,2e-2,80,false);record_inner(inner);
				op_.apply(z,Ad_);
				check(cudaMemcpy(defect_rhs_,v,static_cast<std::size_t>(n_)*sizeof(Real),
					cudaMemcpyDeviceToDevice),"copy first-stage full-operator defect");
				subtract_kernel<<<(n_+255)/256,256>>>(defect_rhs_,Ad_,active_,n_);
				const double source2=dot_active(v,v,active_,n_);
				double defect2=dot_active(defect_rhs_,defect_rhs_,active_,n_);
				const double first_eta=source2>0?std::sqrt(defect2/source2):
					std::numeric_limits<double>::infinity();
				double schwarz_omega=0,schwarz_ratio=1;
				if(schwarz_block_count_&&defect2>0&&std::isfinite(defect2))
				{
					apply_irregular_schwarz(defect_rhs_,correction_);op_.apply(correction_,t_);
					const double image2=dot_active(t_,t_,active_,n_);
					const double numerator=dot_active(defect_rhs_,t_,active_,n_);
					if(image2>0&&std::isfinite(image2)&&std::isfinite(numerator)&&finite_real(numerator/image2))
					{
						schwarz_omega=numerator/image2;
						copy_active_kernel<<<(n_+255)/256,256>>>(z,direction_,active_,n_);
						axpy_active_kernel<<<(n_+255)/256,256>>>(direction_,correction_,static_cast<Real>(schwarz_omega),active_,n_);
						op_.apply(direction_,Ad_);
						check(cudaMemcpy(defect_rhs_,v,static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),"copy Schwarz-corrected full-operator defect");
						subtract_kernel<<<(n_+255)/256,256>>>(defect_rhs_,Ad_,active_,n_);
						const double corrected2=dot_active(defect_rhs_,defect_rhs_,active_,n_);
						++result.irregular_schwarz_applications;
						schwarz_ratio=corrected2>=0&&std::isfinite(corrected2)?std::sqrt(corrected2/defect2):1;
						result.irregular_schwarz_best_defect_ratio=std::min(result.irregular_schwarz_best_defect_ratio,schwarz_ratio);
						if(corrected2<defect2){copy_active_kernel<<<(n_+255)/256,256>>>(direction_,z,active_,n_);defect2=corrected2;}
						else
						{
							op_.apply(z,Ad_);check(cudaMemcpy(defect_rhs_,v,static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),"restore pre-Schwarz full-operator defect");subtract_kernel<<<(n_+255)/256,256>>>(defect_rhs_,Ad_,active_,n_);schwarz_omega=0;schwarz_ratio=1;
						}
					}
				}
				double defect_omega=0;
				int second_iterations=0;double second_residual=0;
				// A0 is an excellent inverse in regular regions, but a large non-orthogonal
				// EB correction can make B=A0^-1 a poor full-A right preconditioner. Apply
				// one residual-correction stage only when the measured full defect warrants
				// it. Flexible GMRES permits this variable preconditioner. The scalar omega
				// is the minimum-residual line search for d-A*(omega*A0^-1*d).
				const double post_schwarz_eta=source2>0?std::sqrt(defect2/source2):std::numeric_limits<double>::infinity();
				if(std::isfinite(post_schwarz_eta)&&post_schwarz_eta>0.25&&defect2>0)
				{
					Real* const defect_correction=gmres_v_+static_cast<std::size_t>(column+1)*n_;
					const AmrGpuSolveResult second=solve_orthogonal(
						defect_correction,defect_rhs_,2e-2,80,false);record_inner(second);
					second_iterations=second.iterations;second_residual=second.relative_residual;
					op_.apply(defect_correction,t_);
					const double image2=dot_active(t_,t_,active_,n_);
					const double numerator=dot_active(defect_rhs_,t_,active_,n_);
					if(image2>0&&std::isfinite(image2)&&std::isfinite(numerator)&&
						finite_real(numerator/image2))defect_omega=numerator/image2;
					if(defect_omega!=0)
						axpy_active_kernel<<<(n_+255)/256,256>>>(z,defect_correction,
							static_cast<Real>(defect_omega),active_,n_);
					// Flexible Arnoldi must use the true image of the FP32 vector actually
					// stored in z. q1+omega*q2 is only the nominal image: the rounded z update
					// can differ by pressure-vector ULPs and would violate A*z=Ad.
					op_.apply(z,Ad_);
					check(cudaMemcpy(defect_rhs_,v,static_cast<std::size_t>(n_)*sizeof(Real),
						cudaMemcpyDeviceToDevice),"copy corrected full-operator defect");
					subtract_kernel<<<(n_+255)/256,256>>>(defect_rhs_,Ad_,active_,n_);
					defect2=dot_active(defect_rhs_,defect_rhs_,active_,n_);
				}
				if(trace&&(column==0||column+1==cycle))
				{
					std::fprintf(stderr,"[pressure-gcro] restart=%d column=%d eta1/schwarz/final=%.9e/%.9e/%.9e schwarz-omega/ratio=%.9e/%.9e defect-omega=%.9e inner1=%d/%.9e inner2=%d/%.9e\n",
						restart_index,column+1,first_eta,
						post_schwarz_eta,
						source2>0?std::sqrt(defect2/source2):std::numeric_limits<double>::infinity(),
						schwarz_omega,schwarz_ratio,defect_omega,inner.iterations,inner.relative_residual,
						second_iterations,second_residual);
				}

				// Project the true image against C and apply the identical coefficients
				// to z against U. Therefore A*z=Ad remains true and the current flexible
				// direction cannot reintroduce a recycled operator mode.
				bool column_valid=true;
				for(int pass=0;pass<2&&column_valid;++pass)
					for(int recycled=0;recycled<gmres_recycle_count_;++recycled)
					{
						const double value=dot_active(recycle_c(recycled),Ad_,active_,n_);
						if(!finite_real(value)){column_valid=false;break;}
						axpy_active_kernel<<<(n_+255)/256,256>>>(Ad_,recycle_c(recycled),
							static_cast<Real>(-value),active_,n_);
						axpy_active_kernel<<<(n_+255)/256,256>>>(z,recycle_u(recycled),
							static_cast<Real>(-value),active_,n_);
					}
				if(trace&&gmres_recycle_count_&&(column==0||column+1==cycle))
				{
					op_.apply(z,t_);
					check(cudaMemcpy(r_,t_,static_cast<std::size_t>(n_)*sizeof(Real),
						cudaMemcpyDeviceToDevice),"copy traced recycled direction image");
					axpy_active_kernel<<<(n_+255)/256,256>>>(r_,Ad_,Real(-1),active_,n_);
					const double gap2=dot_active(r_,r_,active_,n_),image2=dot_active(t_,t_,active_,n_);
					std::fprintf(stderr,"[pressure-gcro] restart=%d column=%d recycled A*z gap=%.9e\n",
						restart_index,column+1,image2>0?std::sqrt(gap2/image2):0);
				}
				if(!column_valid)break;
				const double arnoldi_scale=std::sqrt(std::max(0.0,
					dot_active(Ad_,Ad_,active_,n_)));
				if(!(arnoldi_scale>0)||!std::isfinite(arnoldi_scale))break;
				// Twice-modified Gram-Schmidt keeps an FP32 basis usable while all scalar
				// orthogonalization and the small Hessenberg solve remain FP64 on the host.
				for(int pass=0;pass<2&&column_valid;++pass)
					for(int row=0;row<=column;++row)
					{
						Real* const basis=gmres_v_+static_cast<std::size_t>(row)*n_;
						const double value=dot_active(basis,Ad_,active_,n_);
						if(!finite_real(value)){column_valid=false;break;}
						h[static_cast<std::size_t>(row)*cycle+column]+=value;
						axpy_active_kernel<<<(n_+255)/256,256>>>(Ad_,basis,
							static_cast<Real>(-value),active_,n_);
					}
				if(!column_valid)break;
				double next_norm=std::sqrt(std::max(0.0,dot_active(Ad_,Ad_,active_,n_)));
				const double basis_floor=64*static_cast<double>(std::numeric_limits<Real>::epsilon())*
					arnoldi_scale;
				if(next_norm>basis_floor&&finite_real(1.0/next_norm))
					scaled_copy_active_kernel<<<(n_+255)/256,256>>>(
						gmres_v_+static_cast<std::size_t>(column+1)*n_,Ad_,
						static_cast<Real>(1.0/next_norm),active_,n_);
				else next_norm=0;
				h[static_cast<std::size_t>(column+1)*cycle+column]=next_norm;
				for(int row=0;row<column;++row)
				{
					double& a=h[static_cast<std::size_t>(row)*cycle+column];
					double& b=h[static_cast<std::size_t>(row+1)*cycle+column];
					const double rotated=cosine[row]*a+sine[row]*b;
					b=-sine[row]*a+cosine[row]*b;a=rotated;
				}
				double& diagonal=h[static_cast<std::size_t>(column)*cycle+column];
				double& below=h[static_cast<std::size_t>(column+1)*cycle+column];
				const double norm=std::hypot(diagonal,below);
				if(!(norm>0)||!std::isfinite(norm))break;
				cosine[column]=diagonal/norm;sine[column]=below/norm;
				diagonal=norm;below=0;
				const double old_g=g[column];g[column]=cosine[column]*old_g;
				g[column+1]=-sine[column]*old_g;
				used=column+1;++result.iterations;
				if(std::abs(g[column+1])/std::sqrt(rhs2)<=tolerance||next_norm==0)break;
			}
			if(used==0)break;
			for(int row=used-1;row>=0;--row)
			{
				double value=g[row];
				for(int column=row+1;column<used;++column)
					value-=h[static_cast<std::size_t>(row)*cycle+column]*y[column];
				const double diagonal=h[static_cast<std::size_t>(row)*cycle+row];
				if(!(std::abs(diagonal)>0)||!std::isfinite(diagonal)){used=0;break;}
				y[row]=value/diagonal;
				if(!finite_real(y[row])){used=0;break;}
			}
			if(used==0)break;

			check(cudaMemset(correction_,0,static_cast<std::size_t>(n_)*sizeof(Real)),
				"clear GCRO cycle correction");
			for(int column=0;column<used;++column)
				axpy_active_kernel<<<(n_+255)/256,256>>>(correction_,
					gmres_z_+static_cast<std::size_t>(column)*n_,static_cast<Real>(y[column]),
					active_,n_);
			axpy_active_kernel<<<(n_+255)/256,256>>>(pressure,correction_,Real(1),active_,n_);
			residual2=true_residual();result.relative_residual=std::sqrt(residual2/rhs2);
			const double estimated_relative_residual=std::abs(g[used])/std::sqrt(rhs2);
			if(trace)std::fprintf(stderr,"[pressure-gcro] restart=%d used=%d end-it=%d estimated=%.9e true=%.9e recycle-before=%d\n",
				restart_index,used,result.iterations,estimated_relative_residual,
				result.relative_residual,gmres_recycle_count_);
			if(!std::isfinite(result.relative_residual))break;

			// Store a normalized pair only when its true image contains a resolved direction
			// outside the existing recycle space. Pairwise MGS updates preserve A*d=q.
			op_.apply(correction_,Ad_);
			const double candidate_before=std::sqrt(std::max(0.0,
				dot_active(Ad_,Ad_,active_,n_)));
			bool candidate_valid=candidate_before>0&&std::isfinite(candidate_before);
			for(int pass=0;pass<2&&candidate_valid;++pass)
				for(int recycled=0;recycled<gmres_recycle_count_;++recycled)
				{
					const double value=dot_active(recycle_c(recycled),Ad_,active_,n_);
					if(!finite_real(value)){candidate_valid=false;break;}
					axpy_active_kernel<<<(n_+255)/256,256>>>(Ad_,recycle_c(recycled),
						static_cast<Real>(-value),active_,n_);
					axpy_active_kernel<<<(n_+255)/256,256>>>(correction_,recycle_u(recycled),
						static_cast<Real>(-value),active_,n_);
				}
			const double candidate_after=candidate_valid?
				std::sqrt(std::max(0.0,dot_active(Ad_,Ad_,active_,n_))):0;
			const double novelty_floor=sizeof(Real)==4?
				std::max(1e-5,256*static_cast<double>(std::numeric_limits<Real>::epsilon())):1e-12;
			candidate_valid=candidate_valid&&std::isfinite(candidate_after)&&
				candidate_after>novelty_floor*candidate_before&&finite_real(1.0/candidate_after);
			if(candidate_valid&&recycle_enabled)
			{
				if(gmres_recycle_count_==gmres_recycle_capacity_)
				{
					for(int recycled=1;recycled<gmres_recycle_capacity_;++recycled)
					{
						check(cudaMemcpy(recycle_u(recycled-1),recycle_u(recycled),
							static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),
							"shift GCRO solution recycle basis");
						check(cudaMemcpy(recycle_c(recycled-1),recycle_c(recycled),
							static_cast<std::size_t>(n_)*sizeof(Real),cudaMemcpyDeviceToDevice),
							"shift GCRO image recycle basis");
					}
					--gmres_recycle_count_;
				}
				scaled_copy_active_kernel<<<(n_+255)/256,256>>>(recycle_u(gmres_recycle_count_),
					correction_,static_cast<Real>(1.0/candidate_after),active_,n_);
				scaled_copy_active_kernel<<<(n_+255)/256,256>>>(recycle_c(gmres_recycle_count_),
					Ad_,static_cast<Real>(1.0/candidate_after),active_,n_);
				++gmres_recycle_count_;
				if(trace)
				{
					op_.apply(recycle_u(gmres_recycle_count_-1),t_);
					check(cudaMemcpy(r_,t_,static_cast<std::size_t>(n_)*sizeof(Real),
						cudaMemcpyDeviceToDevice),"copy traced stored recycle image");
					axpy_active_kernel<<<(n_+255)/256,256>>>(r_,recycle_c(gmres_recycle_count_-1),Real(-1),active_,n_);
					const double gap2=dot_active(r_,r_,active_,n_),image2=dot_active(t_,t_,active_,n_);
					std::fprintf(stderr,"[pressure-gcro] restart=%d stored A*U-C gap=%.9e\n",
						restart_index,image2>0?std::sqrt(gap2/image2):0);
				}
				if(trace)std::fprintf(stderr,"[pressure-gcro] restart=%d stored recycle=%d novelty=%.9e\n",
					restart_index,gmres_recycle_count_,candidate_after/candidate_before);
			}
			if(result.relative_residual<=tolerance){result.converged=true;break;}
		}
		if(!result.converged)
		{
			residual2=true_residual();result.relative_residual=std::sqrt(residual2/rhs2);
			result.converged=std::isfinite(result.relative_residual)&&
				result.relative_residual<=tolerance;
		}
		if(trace&&!result.converged)trace_failure_residual(defect_rhs_);
		if(trace)std::fprintf(stderr,"[pressure-solver] converged=%d iterations=%d residual=%.9e recycle-rejections=%d schwarz-blocks/max=%d/%d applications=%d best-defect-ratio=%.9e\n",
			result.converged?1:0,result.iterations,result.relative_residual,result.recycle_rejections,
			result.irregular_schwarz_block_count,result.irregular_schwarz_max_block_size,
			result.irregular_schwarz_applications,result.irregular_schwarz_best_defect_ratio);
		check(cudaDeviceSynchronize(),"solve flexible GCRO/LGMRES composite pressure");
		return result;
	}

	DeviceCompositeAmrProjection::DeviceCompositeAmrProjection(const CompositeAmrPressureSystem& system,DeviceAmrFields& fields)
		:solver_(system),fields_(&fields)
	{
		if(!system.hierarchy)throw std::invalid_argument("composite AMR projection has no hierarchy");level_count_=static_cast<int>(system.hierarchy->levels().size());if(fields.level_count()!=level_count_)throw std::invalid_argument("composite AMR projection field hierarchy mismatch");storage_size_=system.storage_size;brick_size_=system.brick_size;outlet_=system.pressure_outlet_xmax;coarse_fine_count_=static_cast<int>(system.coarse_fine.size());embedded_count_=static_cast<int>(system.embedded.size());special_count_=coarse_fine_count_+embedded_count_;
		std::vector<DeviceCompositeAmrFluxLevelView> views(level_count_);brick_counts_.resize(level_count_);for(int level=0;level<level_count_;++level){const DeviceAmrFieldLevelView field_view=fields.level_view(level);views[level]={system.level_offset[level],system.hierarchy->levels()[level].h,field_view};brick_counts_[level]=field_view.brick_count;}levels_=upload(views,"upload composite AMR flux level views");
		active_=upload(system.active,"upload composite AMR projection active mask");cut_face_mask_=upload(system.cut_face_mask,"upload composite AMR projection cut-face mask");std::vector<Real> volume(system.volume.size());for(std::size_t q=0;q<volume.size();++q)volume[q]=static_cast<Real>(system.volume[q]);volume_=upload(volume,"upload composite AMR projection volumes");
		std::vector<CoarseFinePressureConnection> special=system.coarse_fine;special.insert(special.end(),system.embedded.begin(),system.embedded.end());std::vector<int> first(special_count_),second(special_count_),lower_gradient_node(special_count_),upper_gradient_node(special_count_);std::vector<std::int8_t> direction(special_count_),axis(special_count_);std::vector<Real> area(special_count_),distance(special_count_),gradient_factor(special_count_),nonorthogonal_correction(static_cast<std::size_t>(special_count_)*3),upper_gradient_weight(special_count_);for(int edge=0;edge<special_count_;++edge){first[edge]=special[edge].coarse_dof;second[edge]=special[edge].fine_dof;direction[edge]=special[edge].direction;axis[edge]=special[edge].axis;area[edge]=static_cast<Real>(special[edge].open_area);distance[edge]=static_cast<Real>(special[edge].centre_distance);gradient_factor[edge]=static_cast<Real>(pressure_gradient_factor(special[edge]));lower_gradient_node[edge]=special[edge].lower_gradient_node;upper_gradient_node[edge]=special[edge].upper_gradient_node;for(int component=0;component<3;++component)nonorthogonal_correction[3*edge+component]=static_cast<Real>(special[edge].nonorthogonal_correction[component]);upper_gradient_weight[edge]=static_cast<Real>(special[edge].upper_gradient_weight);}first_dof_=upload(first,"upload composite flux first DOFs");second_dof_=upload(second,"upload composite flux second DOFs");direction_=upload(direction,"upload composite flux direction");axis_=upload(axis,"upload composite flux axes");open_area_=upload(area,"upload composite flux areas");centre_distance_=upload(distance,"upload composite flux distances");pressure_gradient_factor_=upload(gradient_factor,"upload composite flux pressure-gradient factors");special_lower_gradient_node_=upload(lower_gradient_node,"upload special lower pressure-gradient nodes");special_upper_gradient_node_=upload(upper_gradient_node,"upload special upper pressure-gradient nodes");special_nonorthogonal_correction_=upload(nonorthogonal_correction,"upload special nonorthogonal corrections");special_upper_gradient_weight_=upload(upper_gradient_weight,"upload special pressure-gradient interpolation weights");special_velocity_=allocate_zero<Real>(special_count_,"allocate composite special velocity");max_abs_scratch_=allocate_zero<Real>(1,"allocate special-velocity maximum");
		struct SpecialIncident { int dof=0,edge=0;std::int8_t sign=0; };
		std::vector<SpecialIncident> special_incidents;special_incidents.reserve(static_cast<std::size_t>(2)*special_count_);
		for(int edge=0;edge<special_count_;++edge)
		{
			const int lower=direction[edge]>0?first[edge]:second[edge],upper=direction[edge]>0?second[edge]:first[edge];
			if(lower<0||upper<0||lower>=storage_size_||upper>=storage_size_||lower==upper||
				!system.active[lower]||!system.active[upper])continue;
			special_incidents.push_back({lower,edge,1});special_incidents.push_back({upper,edge,-1});
		}
		std::sort(special_incidents.begin(),special_incidents.end(),[](const SpecialIncident& a,const SpecialIncident& b)
		{
			if(a.dof!=b.dof)return a.dof<b.dof;if(a.edge!=b.edge)return a.edge<b.edge;return a.sign<b.sign;
		});
		std::vector<int> special_incident_dof,special_incident_offset,special_incident_edge;
		std::vector<std::int8_t> special_incident_sign;
		for(std::size_t q=0;q<special_incidents.size();++q)
		{
			if(q==0||special_incidents[q].dof!=special_incidents[q-1].dof)
			{
				special_incident_dof.push_back(special_incidents[q].dof);
				special_incident_offset.push_back(static_cast<int>(q));
			}
			special_incident_edge.push_back(special_incidents[q].edge);
			special_incident_sign.push_back(special_incidents[q].sign);
		}
		if(!special_incidents.empty())special_incident_offset.push_back(static_cast<int>(special_incidents.size()));
		special_incident_dof_count_=static_cast<int>(special_incident_dof.size());
		special_incident_dof_=upload(special_incident_dof,"upload deterministic special-flux DOFs");
		special_incident_offset_=upload(special_incident_offset,"upload deterministic special-flux offsets");
		special_incident_edge_=upload(special_incident_edge,"upload deterministic special-flux edges");
		special_incident_sign_=upload(special_incident_sign,"upload deterministic special-flux signs");
		bytes_+=special_incident_dof.size()*sizeof(int)+special_incident_offset.size()*sizeof(int)+
			special_incident_edge.size()*sizeof(int)+special_incident_sign.size()*sizeof(std::int8_t);
		pressure_gradient_node_count_=static_cast<int>(system.pressure_gradient_dof.size());if(pressure_gradient_node_count_){std::vector<Real> weight(static_cast<std::size_t>(system.pressure_gradient_weight.size())*3);for(std::size_t q=0;q<system.pressure_gradient_weight.size();++q)for(int component=0;component<3;++component)weight[3*q+component]=static_cast<Real>(system.pressure_gradient_weight[q][component]);pressure_gradient_node_dof_=upload(system.pressure_gradient_dof,"upload projection pressure-gradient DOFs");pressure_gradient_offset_=upload(system.pressure_gradient_offset,"upload projection pressure-gradient offsets");pressure_gradient_neighbour_=upload(system.pressure_gradient_neighbour,"upload projection pressure-gradient neighbours");pressure_gradient_weight_=upload(weight,"upload projection pressure-gradient weights");pressure_gradient_=allocate_zero<Real>(static_cast<std::size_t>(pressure_gradient_node_count_)*3,"allocate projection pressure gradients");bytes_+=system.pressure_gradient_dof.size()*sizeof(int)+system.pressure_gradient_offset.size()*sizeof(int)+system.pressure_gradient_neighbour.size()*sizeof(int)+weight.size()*sizeof(Real)+static_cast<std::size_t>(pressure_gradient_node_count_)*3*sizeof(Real);}
		regular_pressure_correction_count_=static_cast<int>(system.regular_pressure_corrections.size());if(regular_pressure_correction_count_){std::vector<int> lower(regular_pressure_correction_count_),upper(regular_pressure_correction_count_),lower_node(regular_pressure_correction_count_),upper_node(regular_pressure_correction_count_),level(regular_pressure_correction_count_);std::vector<std::uint64_t> index(regular_pressure_correction_count_),mirror(regular_pressure_correction_count_,~std::uint64_t(0));std::vector<std::int8_t> correction_axis(regular_pressure_correction_count_);std::vector<Real> delta(regular_pressure_correction_count_),correction(static_cast<std::size_t>(regular_pressure_correction_count_)*3),weight(regular_pressure_correction_count_);for(int edge=0;edge<regular_pressure_correction_count_;++edge){const auto& source=system.regular_pressure_corrections[edge];lower[edge]=source.lower_dof;upper[edge]=source.upper_dof;lower_node[edge]=source.lower_gradient_node;upper_node[edge]=source.upper_gradient_node;level[edge]=source.level;correction_axis[edge]=source.axis;delta[edge]=static_cast<Real>(source.two_point_delta);weight[edge]=static_cast<Real>(source.upper_gradient_weight);for(int component=0;component<3;++component)correction[3*edge+component]=static_cast<Real>(source.nonorthogonal_correction[component]);const BrickFieldLayout layout=views[source.level].fields.layout;index[edge]=source.axis==0?layout.u_index(source.brick,source.i,source.j,source.k):(source.axis==1?layout.v_index(source.brick,source.i,source.j,source.k):layout.w_index(source.brick,source.i,source.j,source.k));const int face_coordinate=source.axis==0?source.i:(source.axis==1?source.j:source.k);if(face_coordinate==brick_size_){const int neighbour=system.hierarchy->levels()[source.level].bricks[source.brick].same_level_neighbor[2*source.axis+1];if(neighbour>=0){if(source.axis==0)mirror[edge]=layout.u_index(neighbour,0,source.j,source.k);else if(source.axis==1)mirror[edge]=layout.v_index(neighbour,source.i,0,source.k);else mirror[edge]=layout.w_index(neighbour,source.i,source.j,0);}}}regular_correction_lower_=upload(lower,"upload projection regular-correction lower DOFs");regular_correction_upper_=upload(upper,"upload projection regular-correction upper DOFs");regular_correction_lower_node_=upload(lower_node,"upload projection regular-correction lower gradient nodes");regular_correction_upper_node_=upload(upper_node,"upload projection regular-correction upper gradient nodes");regular_correction_level_=upload(level,"upload projection regular-correction levels");regular_correction_index_=upload(index,"upload projection regular-correction face indices");regular_correction_mirror_index_=upload(mirror,"upload projection regular-correction mirror indices");regular_correction_axis_=upload(correction_axis,"upload projection regular-correction axes");regular_correction_two_point_delta_=upload(delta,"upload projection regular two-point deltas");regular_correction_vector_=upload(correction,"upload projection regular nonorthogonal corrections");regular_correction_upper_weight_=upload(weight,"upload projection regular pressure-gradient interpolation");bytes_+=static_cast<std::size_t>(5*regular_pressure_correction_count_)*sizeof(int)+static_cast<std::size_t>(2*regular_pressure_correction_count_)*sizeof(std::uint64_t)+static_cast<std::size_t>(regular_pressure_correction_count_)*sizeof(std::int8_t)+static_cast<std::size_t>(5*regular_pressure_correction_count_)*sizeof(Real);}
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
			std::vector<Vec3d> node_component_position_sum(static_cast<std::size_t>(embedded_node_count_)*3);std::vector<double> node_component_sample_weight(static_cast<std::size_t>(embedded_node_count_)*3,0);auto add_component_sample=[&](int node_id,int component,Vec3d position,double sample_weight){const int slot=node_id*3+component;node_component_position_sum[slot]=node_component_position_sum[slot]+position*sample_weight;node_component_sample_weight[slot]+=sample_weight;};for(int edge=0;edge<embedded_count_;++edge){const auto& connection=system.embedded[edge];add_component_sample(node_a[edge],connection.axis,connection.face_centroid,connection.open_area);add_component_sample(node_b[edge],connection.axis,connection.face_centroid,connection.open_area);}
			// Build a compact axis-aligned aperture chain at each graph node. A candidate must
			// extend through the shared fluid node in the required direction and use the same
			// velocity component; nearest face centroids resolve junction branches deterministically.
			std::vector<std::vector<int>> incident(embedded_node_count_);for(int edge=0;edge<embedded_count_;++edge){incident[node_a[edge]].push_back(edge);incident[node_b[edge]].push_back(edge);}std::vector<int> negative_edge(embedded_count_,-1),positive_edge(embedded_count_,-1);auto lower_node=[&](int edge){return system.embedded[edge].direction>0?node_a[edge]:node_b[edge];};auto upper_node=[&](int edge){return system.embedded[edge].direction>0?node_b[edge]:node_a[edge];};auto coordinate=[](Vec3d point,int component){return component==0?point.x:(component==1?point.y:point.z);};auto nearest_extension=[&](int edge,int shared,bool negative){int best=-1;double best_distance=std::numeric_limits<double>::infinity();const auto& source=system.embedded[edge];for(int candidate:incident[shared]){if(candidate==edge||system.embedded[candidate].axis!=source.axis)continue;if(negative?upper_node(candidate)!=shared:lower_node(candidate)!=shared)continue;const double delta=coordinate(system.embedded[candidate].face_centroid,source.axis)-coordinate(source.face_centroid,source.axis);if(negative?delta>=-1e-12:delta<=1e-12)continue;const double distance=length2(system.embedded[candidate].face_centroid-source.face_centroid);if(distance<best_distance){best_distance=distance;best=candidate;}}return best;};for(int edge=0;edge<embedded_count_;++edge){negative_edge[edge]=nearest_extension(edge,lower_node(edge),true);positive_edge[edge]=nearest_extension(edge,upper_node(edge),false);}for(int edge=0;edge<embedded_count_;++edge)if(negative_edge[edge]>=0&&positive_edge[edge]>=0&&negative_edge[negative_edge[edge]]>=0&&positive_edge[positive_edge[edge]]>=0)++embedded_high_order_stencil_count_;
			std::vector<int> carrier_node,carrier_level;std::vector<std::uint64_t> carrier_index;std::vector<std::int8_t> carrier_axis;std::vector<Real> carrier_area,carrier_mass;const int cells_per_brick=brick_size_*brick_size_*brick_size_;
			for(int node_id=0;node_id<embedded_node_count_;++node_id)
			{
				const int dof=node_dof[node_id];if(dof<0||dof>=storage_size_||!system.active[dof])continue;int level=-1,brick=-1,local=-1;for(int candidate=0;candidate<level_count_;++candidate){const int begin=system.level_offset[candidate],end=begin+brick_counts_[candidate]*cells_per_brick;if(dof>=begin&&dof<end){level=candidate;const int work=dof-begin;brick=work/cells_per_brick;local=work%cells_per_brick;break;}}if(level<0)continue;const AmrLevel& metadata=system.hierarchy->levels()[level];const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;const int i=local%brick_size_,j=(local/brick_size_)%brick_size_,k=local/(brick_size_*brick_size_);const BrickFieldLayout layout=views[level].fields.layout;
				for(int face_axis=0;face_axis<3;++face_axis)for(int sign=-1;sign<=1;sign+=2)
				{
					int c[3]={i,j,k},other=-1,other_brick=brick;c[face_axis]+=sign;if(c[face_axis]>=0&&c[face_axis]<brick_size_)other=system.dof(level,brick,c[0],c[1],c[2]);else{other_brick=record.same_level_neighbor[2*face_axis+(sign>0)];if(other_brick>=0&&metadata.bricks[other_brick].active()){c[face_axis]=sign>0?0:brick_size_-1;other=system.dof(level,other_brick,c[0],c[1],c[2]);}}if(other<0||!system.active[other])continue;const int lower=sign>0?dof:other;if(system.cut_face_mask[lower]&(1u<<face_axis))continue;int fi=i,fj=j,fk=k;if(face_axis==0)fi+=sign>0;else if(face_axis==1)fj+=sign>0;else fk+=sign>0;const std::uint64_t index=face_axis==0?layout.u_index(brick,fi,fj,fk):(face_axis==1?layout.v_index(brick,fi,fj,fk):layout.w_index(brick,fi,fj,fk));const Vec3d face_centroid=record.origin+Vec3d{(i+0.5)*metadata.h,(j+0.5)*metadata.h,(k+0.5)*metadata.h};Vec3d sample_position=face_centroid;sample_position[face_axis]=record.origin[face_axis]+((face_axis==0?i:(face_axis==1?j:k))+(sign>0?1:0))*metadata.h;carrier_node.push_back(node_id);carrier_level.push_back(level);carrier_index.push_back(index);carrier_axis.push_back(static_cast<std::int8_t>(face_axis));carrier_area.push_back(static_cast<Real>(metadata.h*metadata.h));carrier_mass.push_back(static_cast<Real>(composite_mac_carrier_volume(system,dof,other)));add_component_sample(node_id,face_axis,sample_position,metadata.h*metadata.h);
				}
			}
			const std::vector<CompositeEbMomentumRegularConnection> regular_connections=build_composite_eb_momentum_regular_connections(system);std::vector<int> regular_lower_node,regular_upper_node,regular_level;std::vector<std::uint64_t> regular_index;std::vector<std::int8_t> regular_axis;std::vector<Real> regular_area;regular_lower_node.reserve(regular_connections.size());regular_upper_node.reserve(regular_connections.size());regular_level.reserve(regular_connections.size());regular_index.reserve(regular_connections.size());regular_axis.reserve(regular_connections.size());regular_area.reserve(regular_connections.size());auto compact_node=[&](int dof){const auto found=node_map.find(dof);return found==node_map.end()?-1:found->second;};for(const CompositeEbMomentumRegularConnection& connection:regular_connections){const AmrMacFaceAddress& face=connection.face;if(face.level<0||face.level>=level_count_)throw std::runtime_error("compact perimeter MAC face has invalid level");const BrickFieldLayout layout=views[face.level].fields.layout;const std::uint64_t index=face.component==0?layout.u_index(face.brick,face.i,face.j,face.k):(face.component==1?layout.v_index(face.brick,face.i,face.j,face.k):layout.w_index(face.brick,face.i,face.j,face.k));const int lower=compact_node(connection.lower_dof),upper=compact_node(connection.upper_dof);if(lower<0&&upper<0)throw std::runtime_error("compact perimeter connection has no compact endpoint");regular_lower_node.push_back(lower);regular_upper_node.push_back(upper);regular_level.push_back(face.level);regular_index.push_back(index);regular_axis.push_back(static_cast<std::int8_t>(face.component));regular_area.push_back(static_cast<Real>(connection.open_area));}embedded_regular_connection_count_=static_cast<int>(regular_connections.size());eb_regular_lower_node_=upload(regular_lower_node,"upload compact perimeter lower nodes");eb_regular_upper_node_=upload(regular_upper_node,"upload compact perimeter upper nodes");eb_regular_level_=upload(regular_level,"upload compact perimeter levels");eb_regular_index_=upload(regular_index,"upload compact perimeter MAC indices");eb_regular_axis_=upload(regular_axis,"upload compact perimeter axes");eb_regular_area_=upload(regular_area,"upload compact perimeter areas");eb_node_flux_balance_=allocate_zero<Real>(embedded_node_count_,"allocate compact perimeter flux balance");bytes_+=regular_connections.size()*(3*sizeof(int)+sizeof(std::uint64_t)+sizeof(std::int8_t)+sizeof(Real))+static_cast<std::size_t>(embedded_node_count_)*sizeof(Real);
			std::unordered_map<std::uint64_t,int> carrier_occurrence;for(std::size_t q=0;q<carrier_node.size();++q){const std::uint64_t key=(static_cast<std::uint64_t>(static_cast<std::uint32_t>(carrier_level[q]))<<56)^(static_cast<std::uint64_t>(static_cast<std::uint8_t>(carrier_axis[q]))<<54)^carrier_index[q];++carrier_occurrence[key];}std::vector<Real> carrier_unmapped_mass(carrier_node.size(),Real(0));for(std::size_t q=0;q<carrier_node.size();++q){const std::uint64_t key=(static_cast<std::uint64_t>(static_cast<std::uint32_t>(carrier_level[q]))<<56)^(static_cast<std::uint64_t>(static_cast<std::uint8_t>(carrier_axis[q]))<<54)^carrier_index[q];const int occurrence=carrier_occurrence[key];if(occurrence<1||occurrence>2)throw std::runtime_error("compact EB carrier has invalid fragment multiplicity");if(occurrence==1)carrier_unmapped_mass[q]=Real(0.5)*carrier_mass[q];}
			for(int node_id=0;node_id<embedded_node_count_;++node_id)for(int component=0;component<3;++component){const int slot=node_id*3+component;if(node_component_sample_weight[slot]>0)node_component_position_sum[slot]=node_component_position_sum[slot]/node_component_sample_weight[slot];}
			std::vector<Real> edge_ls_displacement(static_cast<std::size_t>(embedded_count_)*9,Real(0)),edge_ls_weight(static_cast<std::size_t>(embedded_count_)*3,Real(0)),node_gradient_inverse(static_cast<std::size_t>(embedded_node_count_)*27,Real(0));std::vector<double> normal_matrix(static_cast<std::size_t>(embedded_node_count_)*27,0);for(int edge=0;edge<embedded_count_;++edge)for(int component=0;component<3;++component){const int a=node_a[edge],b=node_b[edge],sa=a*3+component,sb=b*3+component;if(!(node_component_sample_weight[sa]>0&&node_component_sample_weight[sb]>0))continue;const Vec3d displacement=node_component_position_sum[sb]-node_component_position_sum[sa];const double distance2=length2(displacement);if(!(distance2>1e-20))continue;const double sample_weight=system.embedded[edge].open_area/distance2;edge_ls_weight[edge*3+component]=static_cast<Real>(sample_weight);for(int derivative=0;derivative<3;++derivative)edge_ls_displacement[edge*9+component*3+derivative]=static_cast<Real>(displacement[derivative]);for(int endpoint:{a,b})for(int row=0;row<3;++row)for(int column=0;column<3;++column)normal_matrix[endpoint*27+component*9+row*3+column]+=sample_weight*displacement[row]*displacement[column];}
			embedded_gradient_rank_.assign(static_cast<std::size_t>(embedded_node_count_)*3,0);for(int node_id=0;node_id<embedded_node_count_;++node_id)for(int component=0;component<3;++component){double inverse[9];const int rank=symmetric_pseudoinverse_3x3(normal_matrix.data()+node_id*27+component*9,inverse);embedded_gradient_rank_[node_id*3+component]=static_cast<std::uint8_t>(rank);if(rank==3)++embedded_least_squares_full_rank_count_;for(int q=0;q<9;++q)node_gradient_inverse[node_id*27+component*9+q]=static_cast<Real>(inverse[q]);}
			embedded_carrier_count_=static_cast<int>(carrier_node.size());eb_node_a_=upload(node_a,"upload EB transport first nodes");eb_node_b_=upload(node_b,"upload EB transport second nodes");eb_negative_edge_=upload(negative_edge,"upload EB negative graph neighbours");eb_positive_edge_=upload(positive_edge,"upload EB positive graph neighbours");eb_carrier_node_=upload(carrier_node,"upload EB transport carrier nodes");eb_carrier_level_=upload(carrier_level,"upload EB transport carrier levels");eb_carrier_index_=upload(carrier_index,"upload EB transport carrier indices");eb_carrier_axis_=upload(carrier_axis,"upload EB transport carrier axes");eb_carrier_area_=upload(carrier_area,"upload EB transport carrier areas");eb_carrier_mass_=upload(carrier_mass,"upload EB transport carrier masses");eb_node_axis_sum_=allocate_zero<Real>(static_cast<std::size_t>(embedded_node_count_)*3,"allocate EB node component sums");eb_node_axis_weight_=allocate_zero<Real>(static_cast<std::size_t>(embedded_node_count_)*3,"allocate EB node component weights");eb_node_gradient_sum_=allocate_zero<Real>(static_cast<std::size_t>(embedded_node_count_)*9,"allocate EB node gradient sums");eb_node_gradient_inverse_=upload(node_gradient_inverse,"upload EB least-squares gradient inverses");eb_edge_ls_displacement_=upload(edge_ls_displacement,"upload EB least-squares displacements");eb_edge_ls_weight_=upload(edge_ls_weight,"upload EB least-squares weights");eb_transport_length_=upload(transport_length,"upload EB aperture transport lengths");eb_transport_scratch_=allocate_zero<Real>(embedded_count_,"allocate EB aperture transport scratch");eb_diffusion_rate_=allocate_zero<Real>(embedded_count_,"allocate EB aperture diffusion rates");eb_node_diffusion_sum_=allocate_zero<Real>(static_cast<std::size_t>(embedded_node_count_)*3,"allocate EB node diffusion sums");eb_node_neighbor_count_=allocate_zero<Real>(static_cast<std::size_t>(embedded_node_count_)*3,"allocate EB node neighbour counts");eb_diffusion_scratch_=allocate_zero<Real>(embedded_count_,"allocate conservative EB diffusion scratch");bytes_+=4*node_a.size()*sizeof(int)+2*carrier_node.size()*sizeof(int)+carrier_index.size()*sizeof(std::uint64_t)+carrier_axis.size()*sizeof(std::int8_t)+2*carrier_area.size()*sizeof(Real)+static_cast<std::size_t>(embedded_node_count_)*48*sizeof(Real)+static_cast<std::size_t>(embedded_count_)*16*sizeof(Real);
			eb_carrier_unmapped_mass_=upload(carrier_unmapped_mass,"upload EB unmapped carrier masses");bytes_+=carrier_unmapped_mass.size()*sizeof(Real);
		}
		if(!system.surface_patches.empty())
		{
			std::unordered_map<int,int> wall_node_map;std::vector<int> wall_node_dof;auto wall_node=[&](int dof){if(dof<0||dof>=storage_size_||!system.active[dof])return -1;auto found=wall_node_map.find(dof);if(found!=wall_node_map.end())return found->second;const int id=static_cast<int>(wall_node_dof.size());wall_node_map.emplace(dof,id);wall_node_dof.push_back(dof);return id;};for(const auto& patch:system.surface_patches){wall_node(patch.plus_dof);wall_node(patch.minus_dof);}fabric_wall_node_count_=static_cast<int>(wall_node_dof.size());
			std::vector<Real> wall_rate(fabric_wall_node_count_,Real(0)),wall_patch_normal,wall_patch_area,wall_patch_distance;std::vector<int> wall_patch_node;for(const auto& patch:system.surface_patches){const double normal_length=std::sqrt(length2(patch.normal));if(!(patch.area>=0&&normal_length>0))throw std::runtime_error("invalid fabric wall patch geometry");const Vec3d unit_normal=patch.normal/normal_length;for(int dof:{patch.plus_dof,patch.minus_dof}){auto found=wall_node_map.find(dof);if(found==wall_node_map.end())continue;const double volume=system.volume[dof],distance=std::abs(dot(system.centroid[dof]-patch.centroid,unit_normal)),scale=std::cbrt(std::max(0.0,volume));if(!(volume>0&&distance>std::max(1e-14,scale*1e-10)))throw std::runtime_error("fabric wall patch has zero same-side centroid distance");wall_rate[found->second]+=static_cast<Real>(patch.area/(volume*distance));wall_patch_node.push_back(found->second);wall_patch_normal.insert(wall_patch_normal.end(),{static_cast<Real>(unit_normal.x),static_cast<Real>(unit_normal.y),static_cast<Real>(unit_normal.z)});wall_patch_area.push_back(static_cast<Real>(patch.area));wall_patch_distance.push_back(static_cast<Real>(distance));wall_patch_source_triangle_id_.push_back(patch.source_triangle_id);wall_patch_source_face_id_.push_back(patch.source_face_id);wall_patch_centroid_.push_back(patch.centroid);}}
			std::vector<int> wall_edge_a(embedded_count_,-1),wall_edge_b(embedded_count_,-1);for(int edge=0;edge<embedded_count_;++edge){auto a=wall_node_map.find(system.embedded[edge].coarse_dof),b=wall_node_map.find(system.embedded[edge].fine_dof);if(a!=wall_node_map.end())wall_edge_a[edge]=a->second;if(b!=wall_node_map.end())wall_edge_b[edge]=b->second;}
			std::vector<int> carrier_node,carrier_level;std::vector<std::uint64_t> carrier_index;std::vector<std::int8_t> carrier_axis;std::vector<Real> carrier_mass;const int cells_per_brick=brick_size_*brick_size_*brick_size_;for(int node_id=0;node_id<fabric_wall_node_count_;++node_id)
			{
				const int dof=wall_node_dof[node_id];int level=-1,brick=-1,local=-1;for(int candidate=0;candidate<level_count_;++candidate){const int begin=system.level_offset[candidate],end=begin+brick_counts_[candidate]*cells_per_brick;if(dof>=begin&&dof<end){level=candidate;const int work=dof-begin;brick=work/cells_per_brick;local=work%cells_per_brick;break;}}if(level<0)continue;const AmrLevel& metadata=system.hierarchy->levels()[level];const BrickMetadata& record=metadata.bricks[brick];if(!record.active())continue;const int i=local%brick_size_,j=(local/brick_size_)%brick_size_,k=local/(brick_size_*brick_size_);const BrickFieldLayout layout=views[level].fields.layout;
				for(int face_axis=0;face_axis<3;++face_axis)for(int sign=-1;sign<=1;sign+=2){int c[3]={i,j,k},other=-1,other_brick=brick;c[face_axis]+=sign;if(c[face_axis]>=0&&c[face_axis]<brick_size_)other=system.dof(level,brick,c[0],c[1],c[2]);else{other_brick=record.same_level_neighbor[2*face_axis+(sign>0)];if(other_brick>=0&&metadata.bricks[other_brick].active()){c[face_axis]=sign>0?0:brick_size_-1;other=system.dof(level,other_brick,c[0],c[1],c[2]);}}if(other<0||!system.active[other])continue;const int lower=sign>0?dof:other;if(system.cut_face_mask[lower]&(1u<<face_axis))continue;int fi=i,fj=j,fk=k;if(face_axis==0)fi+=sign>0;else if(face_axis==1)fj+=sign>0;else fk+=sign>0;const std::uint64_t index=face_axis==0?layout.u_index(brick,fi,fj,fk):(face_axis==1?layout.v_index(brick,fi,fj,fk):layout.w_index(brick,fi,fj,fk));carrier_node.push_back(node_id);carrier_level.push_back(level);carrier_index.push_back(index);carrier_axis.push_back(static_cast<std::int8_t>(face_axis));carrier_mass.push_back(static_cast<Real>(composite_mac_carrier_volume(system,dof,other)));}
			}
			std::vector<int> unique_carrier_level;std::vector<std::uint64_t> unique_carrier_index;std::vector<std::int8_t> unique_carrier_axis;std::vector<Real> unique_carrier_mass;std::unordered_map<std::uint64_t,int> unique_carrier;for(std::size_t q=0;q<carrier_node.size();++q){const std::uint64_t key=(static_cast<std::uint64_t>(static_cast<std::uint32_t>(carrier_level[q]))<<56)^(static_cast<std::uint64_t>(static_cast<std::uint8_t>(carrier_axis[q]))<<54)^carrier_index[q];auto inserted=unique_carrier.emplace(key,static_cast<int>(unique_carrier_level.size()));if(inserted.second){unique_carrier_level.push_back(carrier_level[q]);unique_carrier_index.push_back(carrier_index[q]);unique_carrier_axis.push_back(carrier_axis[q]);unique_carrier_mass.push_back(carrier_mass[q]);}else if(std::abs(static_cast<double>(unique_carrier_mass[inserted.first->second]-carrier_mass[q]))>1e-6*std::max(1.0,std::abs(static_cast<double>(carrier_mass[q]))))throw std::runtime_error("fabric wall carrier aliases disagree on physical volume");}
			fabric_wall_carrier_count_=static_cast<int>(carrier_node.size());
			fabric_wall_patch_count_=static_cast<int>(wall_patch_node.size());
			fabric_wall_unique_carrier_count_=static_cast<int>(unique_carrier_level.size());
			wall_edge_node_a_=upload(wall_edge_a,"upload fabric wall first nodes");
			wall_edge_node_b_=upload(wall_edge_b,"upload fabric wall second nodes");
			wall_carrier_node_=upload(carrier_node,"upload fabric wall carrier nodes");
			wall_carrier_level_=upload(carrier_level,"upload fabric wall carrier levels");
			wall_carrier_index_=upload(carrier_index,"upload fabric wall carrier indices");
			wall_carrier_axis_=upload(carrier_axis,"upload fabric wall carrier axes");
			wall_carrier_mass_=upload(carrier_mass,"upload fabric wall carrier masses");
			wall_unique_carrier_level_=upload(unique_carrier_level,"upload unique fabric wall carrier levels");
			wall_unique_carrier_index_=upload(unique_carrier_index,"upload unique fabric wall carrier indices");
			wall_unique_carrier_axis_=upload(unique_carrier_axis,"upload unique fabric wall carrier axes");
			wall_unique_carrier_mass_=upload(unique_carrier_mass,"upload unique fabric wall carrier masses");
			wall_patch_node_=upload(wall_patch_node,"upload smooth-wall patch nodes");
			wall_patch_normal_=upload(wall_patch_normal,"upload smooth-wall patch normals");
			wall_patch_area_=upload(wall_patch_area,"upload smooth-wall patch areas");
			wall_patch_distance_=upload(wall_patch_distance,"upload smooth-wall patch distances");
			wall_patch_coefficient_=allocate_zero<Real>(fabric_wall_patch_count_,
				"allocate smooth-wall patch coefficients");
			wall_patch_force_per_density_=allocate_zero<Real>(
				static_cast<std::size_t>(fabric_wall_patch_count_)*3,"allocate smooth-wall patch forces");
			wall_momentum_scratch_=allocate_zero<double>(3,"allocate fabric wall momentum reduction");
			wall_node_rate_=upload(wall_rate,"upload fabric wall viscous rates");
			wall_node_axis_sum_=allocate_zero<Real>(
				static_cast<std::size_t>(fabric_wall_node_count_)*3,"allocate fabric wall velocity sums");
			wall_node_axis_weight_=allocate_zero<Real>(
				static_cast<std::size_t>(fabric_wall_node_count_)*3,"allocate fabric wall velocity weights");
			wall_node_velocity_delta_=allocate_zero<Real>(
				static_cast<std::size_t>(fabric_wall_node_count_)*3,"allocate fabric wall velocity deltas");
			wall_node_matrix_=allocate_zero<Real>(
				static_cast<std::size_t>(fabric_wall_node_count_)*6,"allocate implicit smooth-wall matrices");
			bytes_+=2*wall_edge_a.size()*sizeof(int)+2*carrier_node.size()*sizeof(int)+
				carrier_index.size()*sizeof(std::uint64_t)+carrier_axis.size()*sizeof(std::int8_t)+
				carrier_mass.size()*sizeof(Real)+wall_patch_node.size()*sizeof(int)+
				(9*wall_patch_node.size()+static_cast<std::size_t>(fabric_wall_node_count_)*16)*sizeof(Real)+
				3*sizeof(double)+unique_carrier_level.size()*sizeof(int)+
				unique_carrier_index.size()*sizeof(std::uint64_t)+
				unique_carrier_axis.size()*sizeof(std::int8_t)+unique_carrier_mass.size()*sizeof(Real);
		}
		integrated_=allocate_zero<Real>(storage_size_,"allocate composite integrated flux");divergence_=allocate_zero<Real>(storage_size_,"allocate composite divergence");rhs_=allocate_zero<Real>(storage_size_,"allocate composite RHS");pressure_=allocate_zero<Real>(storage_size_,"allocate composite pressure");
		bytes_+=solver_.bytes()+views.size()*sizeof(DeviceCompositeAmrFluxLevelView)+system.active.size()+system.cut_face_mask.size()*sizeof(std::uint8_t)+volume.size()*sizeof(Real)+static_cast<std::size_t>(special_count_)*(2*sizeof(int)+2*sizeof(std::int8_t)+4*sizeof(Real))+static_cast<std::size_t>(4)*storage_size_*sizeof(Real)+sizeof(Real);
	}

	DeviceCompositeAmrProjection::~DeviceCompositeAmrProjection()
	{
		if(eb_carrier_unmapped_mass_)cudaFree(eb_carrier_unmapped_mass_);
		for(void* pointer:{(void*)pressure_gradient_node_dof_,(void*)pressure_gradient_offset_,(void*)pressure_gradient_neighbour_,(void*)pressure_gradient_weight_,(void*)pressure_gradient_,(void*)special_lower_gradient_node_,(void*)special_upper_gradient_node_,(void*)special_nonorthogonal_correction_,(void*)special_upper_gradient_weight_,(void*)regular_correction_lower_,(void*)regular_correction_upper_,(void*)regular_correction_lower_node_,(void*)regular_correction_upper_node_,(void*)regular_correction_level_,(void*)regular_correction_index_,(void*)regular_correction_mirror_index_,(void*)regular_correction_axis_,(void*)regular_correction_two_point_delta_,(void*)regular_correction_vector_,(void*)regular_correction_upper_weight_})if(pointer)cudaFree(pointer);
		for(void* pointer:{(void*)levels_,(void*)active_,(void*)cut_face_mask_,(void*)volume_,(void*)integrated_,(void*)divergence_,(void*)rhs_,(void*)pressure_,(void*)first_dof_,(void*)second_dof_,(void*)special_incident_dof_,(void*)special_incident_offset_,(void*)special_incident_edge_,(void*)special_incident_sign_,(void*)direction_,(void*)axis_,(void*)open_area_,(void*)centre_distance_,(void*)pressure_gradient_factor_,(void*)special_velocity_,(void*)max_abs_scratch_,(void*)cf_fine_level_,(void*)cf_fine_index_,(void*)cf_group_,(void*)cf_group_level_,(void*)cf_group_index_,(void*)cf_group_axis_,(void*)cf_group_area_,(void*)cf_group_sum_,(void*)eb_node_a_,(void*)eb_node_b_,(void*)eb_negative_edge_,(void*)eb_positive_edge_,(void*)eb_carrier_node_,(void*)eb_carrier_level_,(void*)eb_carrier_index_,(void*)eb_carrier_axis_,(void*)eb_carrier_area_,(void*)eb_carrier_mass_,(void*)eb_regular_lower_node_,(void*)eb_regular_upper_node_,(void*)eb_regular_level_,(void*)eb_regular_index_,(void*)eb_regular_axis_,(void*)eb_regular_area_,(void*)eb_node_flux_balance_,(void*)eb_node_axis_sum_,(void*)eb_node_axis_weight_,(void*)eb_node_gradient_sum_,(void*)eb_node_gradient_inverse_,(void*)eb_edge_ls_displacement_,(void*)eb_edge_ls_weight_,(void*)eb_transport_length_,(void*)eb_transport_scratch_,(void*)eb_diffusion_rate_,(void*)eb_node_diffusion_sum_,(void*)eb_node_neighbor_count_,(void*)eb_diffusion_scratch_,(void*)wall_edge_node_a_,(void*)wall_edge_node_b_,(void*)wall_patch_node_,(void*)wall_patch_normal_,(void*)wall_patch_area_,(void*)wall_patch_distance_,(void*)wall_patch_coefficient_,(void*)wall_patch_force_per_density_,(void*)wall_carrier_node_,(void*)wall_carrier_level_,(void*)wall_carrier_index_,(void*)wall_carrier_axis_,(void*)wall_carrier_mass_,(void*)wall_unique_carrier_level_,(void*)wall_unique_carrier_index_,(void*)wall_unique_carrier_axis_,(void*)wall_unique_carrier_mass_,(void*)wall_momentum_scratch_,(void*)wall_node_rate_,(void*)wall_node_axis_sum_,(void*)wall_node_axis_weight_,(void*)wall_node_velocity_delta_,(void*)wall_node_matrix_})if(pointer)cudaFree(pointer);
	}

	void DeviceCompositeAmrProjection::clear_special_fluxes(){if(special_count_)check(cudaMemset(special_velocity_,0,static_cast<std::size_t>(special_count_)*sizeof(Real)),"clear composite special fluxes");}
	void DeviceCompositeAmrProjection::initialize_special_freestream(Real speed){if(special_count_)initialize_special_freestream_kernel<<<(special_count_+255)/256,256>>>(axis_,special_velocity_,special_count_,speed);check(cudaDeviceSynchronize(),"initialize composite special freestream");}
	void DeviceCompositeAmrProjection::sync_coarse_fine_from_fields()
	{
		if(coarse_fine_count_)gather_coarse_fine_velocity_kernel<<<(coarse_fine_count_+255)/256,256>>>(levels_,cf_fine_level_,cf_fine_index_,axis_,special_velocity_,coarse_fine_count_);check(cudaDeviceSynchronize(),"gather transported coarse/fine velocities");
	}
	void DeviceCompositeAmrProjection::transport_embedded_apertures(Real dt,Real molecular_nu,Real smagorinsky_cs,bool apply_molecular_fabric_wall)
	{
		if(!(dt>Real(0))||molecular_nu<Real(0)||smagorinsky_cs<Real(0))throw std::invalid_argument("embedded aperture transport timestep/viscosity");
		if(!embedded_count_&&!fabric_wall_node_count_)return;
		if(embedded_count_)
		{
		const std::size_t node_components=static_cast<std::size_t>(embedded_node_count_)*3;
		const std::size_t node_gradients=static_cast<std::size_t>(embedded_node_count_)*9;
		check(cudaMemset(eb_node_axis_sum_,0,node_components*sizeof(Real)),"clear EB transport component sums");
		check(cudaMemset(eb_node_axis_weight_,0,node_components*sizeof(Real)),"clear EB transport component weights");
		check(cudaMemset(eb_node_gradient_sum_,0,node_gradients*sizeof(Real)),"clear EB least-squares gradient RHS");
		if(embedded_carrier_count_)accumulate_eb_carrier_kernel<<<(embedded_carrier_count_+255)/256,256>>>(levels_,eb_carrier_node_,eb_carrier_level_,eb_carrier_index_,eb_carrier_axis_,eb_carrier_area_,eb_node_axis_sum_,eb_node_axis_weight_,embedded_carrier_count_);
		accumulate_eb_aperture_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,axis_,open_area_,special_velocity_,eb_node_axis_sum_,eb_node_axis_weight_,coarse_fine_count_,embedded_count_);
		accumulate_eb_gradient_rhs_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,eb_edge_ls_displacement_,eb_edge_ls_weight_,eb_node_axis_sum_,eb_node_axis_weight_,eb_node_gradient_sum_,embedded_count_);
		solve_eb_gradient_kernel<<<(embedded_node_count_+255)/256,256>>>(eb_node_gradient_sum_,eb_node_gradient_inverse_,embedded_node_count_);
		transport_eb_aperture_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,eb_negative_edge_,eb_positive_edge_,direction_,axis_,eb_transport_length_,special_velocity_,eb_node_axis_sum_,eb_node_axis_weight_,eb_node_gradient_sum_,dt,molecular_nu,smagorinsky_cs,eb_transport_scratch_,eb_diffusion_rate_,coarse_fine_count_,embedded_count_);

		// Reconstruct a compatible fragment-centred vector state from the advected
		// aperture/carrier velocities. The dual mass A/g uses the exact pressure-
		// gradient factor, so this remap conserves momentum and dissipates only the
		// pressure-null circulation that cannot represent a fragment velocity.
		check(cudaMemset(eb_node_axis_sum_,0,node_components*sizeof(Real)),"clear compatible EB momentum sums");
		check(cudaMemset(eb_node_axis_weight_,0,node_components*sizeof(Real)),"clear compatible EB momentum weights");
		accumulate_eb_dual_momentum_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,axis_,open_area_,pressure_gradient_factor_,eb_transport_scratch_,eb_node_axis_sum_,eb_node_axis_weight_,coarse_fine_count_,embedded_count_);
		if(embedded_carrier_count_)accumulate_conservative_carrier_kernel<<<(embedded_carrier_count_+255)/256,256>>>(levels_,eb_carrier_node_,eb_carrier_level_,eb_carrier_index_,eb_carrier_axis_,eb_carrier_mass_,eb_node_axis_sum_,eb_node_axis_weight_,embedded_carrier_count_);
		reconstruct_compatible_eb_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,axis_+coarse_fine_count_,eb_node_axis_sum_,eb_node_axis_weight_,eb_transport_scratch_,eb_diffusion_scratch_,embedded_count_);
		if(embedded_carrier_count_)reconstruct_compatible_carrier_kernel<<<(embedded_carrier_count_+255)/256,256>>>(levels_,eb_carrier_node_,eb_carrier_level_,eb_carrier_index_,eb_carrier_axis_,eb_node_axis_sum_,eb_node_axis_weight_,embedded_carrier_count_);

		// Apply molecular/LES diffusion as conservative pairwise fragment momentum
		// exchange after the compatible reconstruction.
		check(cudaMemset(eb_node_axis_sum_,0,node_components*sizeof(Real)),"clear conservative EB momentum sums");
		check(cudaMemset(eb_node_axis_weight_,0,node_components*sizeof(Real)),"clear conservative EB momentum weights");
		check(cudaMemset(eb_node_diffusion_sum_,0,node_components*sizeof(Real)),"clear conservative EB momentum deltas");
		check(cudaMemset(eb_node_neighbor_count_,0,node_components*sizeof(Real)),"clear conservative EB neighbour counts");
		accumulate_eb_dual_momentum_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,axis_,open_area_,pressure_gradient_factor_,eb_diffusion_scratch_,eb_node_axis_sum_,eb_node_axis_weight_,coarse_fine_count_,embedded_count_);
		if(embedded_carrier_count_)accumulate_conservative_carrier_kernel<<<(embedded_carrier_count_+255)/256,256>>>(levels_,eb_carrier_node_,eb_carrier_level_,eb_carrier_index_,eb_carrier_axis_,eb_carrier_mass_,eb_node_axis_sum_,eb_node_axis_weight_,embedded_carrier_count_);
		count_node_diffusion_neighbours_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,eb_node_axis_weight_,eb_node_neighbor_count_,embedded_count_);
		exchange_node_momentum_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,eb_diffusion_rate_,eb_node_axis_sum_,eb_node_axis_weight_,eb_node_neighbor_count_,dt,eb_node_diffusion_sum_,embedded_count_);
		apply_node_momentum_to_eb_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,axis_+coarse_fine_count_,eb_node_axis_weight_,eb_node_diffusion_sum_,eb_diffusion_scratch_,eb_transport_scratch_,embedded_count_);
		if(embedded_carrier_count_)apply_node_momentum_to_carrier_kernel<<<(embedded_carrier_count_+255)/256,256>>>(levels_,eb_carrier_node_,eb_carrier_level_,eb_carrier_index_,eb_carrier_axis_,eb_node_axis_weight_,eb_node_diffusion_sum_,embedded_carrier_count_);
		commit_eb_aperture_kernel<<<(embedded_count_+255)/256,256>>>(special_velocity_,eb_transport_scratch_,coarse_fine_count_,embedded_count_);
		}
		if(apply_molecular_fabric_wall&&fabric_wall_node_count_&&molecular_nu>Real(0))
		{
			const std::size_t components=static_cast<std::size_t>(fabric_wall_node_count_)*3;check(cudaMemset(wall_node_axis_sum_,0,components*sizeof(Real)),"clear fabric wall velocity sums");check(cudaMemset(wall_node_axis_weight_,0,components*sizeof(Real)),"clear fabric wall velocity weights");if(embedded_count_)accumulate_wall_aperture_velocity_kernel<<<(embedded_count_+255)/256,256>>>(wall_edge_node_a_,wall_edge_node_b_,axis_,open_area_,pressure_gradient_factor_,special_velocity_,wall_node_axis_sum_,wall_node_axis_weight_,coarse_fine_count_,embedded_count_);if(fabric_wall_carrier_count_)accumulate_wall_carrier_velocity_kernel<<<(fabric_wall_carrier_count_+255)/256,256>>>(levels_,wall_carrier_node_,wall_carrier_level_,wall_carrier_index_,wall_carrier_axis_,wall_carrier_mass_,wall_node_axis_sum_,wall_node_axis_weight_,fabric_wall_carrier_count_);compute_no_slip_wall_delta_kernel<<<(fabric_wall_node_count_+255)/256,256>>>(wall_node_axis_sum_,wall_node_axis_weight_,wall_node_rate_,dt,molecular_nu,wall_node_velocity_delta_,fabric_wall_node_count_);if(embedded_count_)apply_no_slip_wall_to_aperture_kernel<<<(embedded_count_+255)/256,256>>>(wall_edge_node_a_,wall_edge_node_b_,axis_,wall_node_velocity_delta_,special_velocity_,coarse_fine_count_,embedded_count_);if(fabric_wall_carrier_count_)apply_no_slip_wall_to_carrier_kernel<<<(fabric_wall_carrier_count_+255)/256,256>>>(levels_,wall_carrier_node_,wall_carrier_level_,wall_carrier_index_,wall_carrier_axis_,wall_node_velocity_delta_,fabric_wall_carrier_count_);
		}
		check(cudaDeviceSynchronize(),"transport compact EB and fabric wall velocities");
	}
	void DeviceCompositeAmrProjection::apply_smooth_fabric_wall_model(Real dt,Real molecular_nu)
	{
		if(!(dt>Real(0))||!(molecular_nu>Real(0)))
			throw std::invalid_argument("smooth fabric wall timestep/viscosity");
		if(!fabric_wall_node_count_||!fabric_wall_patch_count_)return;
		const std::size_t components=static_cast<std::size_t>(fabric_wall_node_count_)*3;
		check(cudaMemset(wall_node_axis_sum_,0,components*sizeof(Real)),
			"clear smooth-wall velocity sums");
		check(cudaMemset(wall_node_axis_weight_,0,components*sizeof(Real)),
			"clear smooth-wall velocity weights");
		check(cudaMemset(wall_node_velocity_delta_,0,components*sizeof(Real)),
			"clear smooth-wall velocity deltas");
		check(cudaMemset(wall_node_matrix_,0,
			static_cast<std::size_t>(fabric_wall_node_count_)*6*sizeof(Real)),
			"clear implicit smooth-wall matrices");
		if(embedded_count_)
			accumulate_wall_aperture_velocity_kernel<<<(embedded_count_+255)/256,256>>>(
				wall_edge_node_a_,wall_edge_node_b_,axis_,open_area_,pressure_gradient_factor_,
				special_velocity_,wall_node_axis_sum_,wall_node_axis_weight_,coarse_fine_count_,
				embedded_count_);
		if(fabric_wall_carrier_count_)
			accumulate_wall_carrier_velocity_kernel<<<(fabric_wall_carrier_count_+255)/256,256>>>(
				levels_,wall_carrier_node_,wall_carrier_level_,wall_carrier_index_,
				wall_carrier_axis_,wall_carrier_mass_,wall_node_axis_sum_,wall_node_axis_weight_,
				fabric_wall_carrier_count_);
		accumulate_smooth_wall_node_matrix_kernel<<<(fabric_wall_patch_count_+255)/256,256>>>(
			wall_patch_node_,wall_patch_normal_,wall_patch_area_,wall_patch_distance_,
			wall_node_axis_sum_,wall_node_axis_weight_,molecular_nu,wall_node_matrix_,
			wall_patch_coefficient_,fabric_wall_patch_count_);
		solve_implicit_smooth_wall_node_kernel<<<(fabric_wall_node_count_+255)/256,256>>>(
			wall_node_axis_sum_,wall_node_axis_weight_,wall_node_matrix_,dt,
			wall_node_velocity_delta_,fabric_wall_node_count_);
		finalize_smooth_wall_patch_force_kernel<<<(fabric_wall_patch_count_+255)/256,256>>>(
			wall_patch_node_,wall_patch_normal_,wall_node_axis_sum_,wall_node_axis_weight_,
			wall_node_velocity_delta_,wall_patch_coefficient_,wall_patch_force_per_density_,
			fabric_wall_patch_count_);
		if(embedded_count_)
			apply_no_slip_wall_to_aperture_kernel<<<(embedded_count_+255)/256,256>>>(
				wall_edge_node_a_,wall_edge_node_b_,axis_,wall_node_velocity_delta_,
				special_velocity_,coarse_fine_count_,embedded_count_);
		if(fabric_wall_carrier_count_)
			apply_no_slip_wall_to_carrier_kernel<<<(fabric_wall_carrier_count_+255)/256,256>>>(
				levels_,wall_carrier_node_,wall_carrier_level_,wall_carrier_index_,
				wall_carrier_axis_,wall_node_velocity_delta_,fabric_wall_carrier_count_);
		check(cudaDeviceSynchronize(),"apply implicit smooth fabric wall model");
	}
	void DeviceCompositeAmrProjection::download_smooth_fabric_wall_loads(Real rho,std::vector<SmoothFabricWallPatchLoad>& host) const
	{
		if(!(rho>Real(0)))throw std::invalid_argument("smooth fabric wall load density");host.resize(fabric_wall_patch_count_);std::vector<Real> force(static_cast<std::size_t>(fabric_wall_patch_count_)*3);if(!force.empty())check(cudaMemcpy(force.data(),wall_patch_force_per_density_,force.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download smooth-wall patch forces");for(int patch=0;patch<fabric_wall_patch_count_;++patch){host[patch].source_triangle_id=wall_patch_source_triangle_id_[patch];host[patch].source_face_id=wall_patch_source_face_id_[patch];host[patch].centroid=wall_patch_centroid_[patch];host[patch].force={static_cast<double>(rho*force[3*patch]),static_cast<double>(rho*force[3*patch+1]),static_cast<double>(rho*force[3*patch+2])};}
	}
	std::array<double,3> DeviceCompositeAmrProjection::fabric_wall_physical_momentum() const
	{
		std::array<double,3> device_sum{};if(!wall_momentum_scratch_)return {};check(cudaMemset(wall_momentum_scratch_,0,3*sizeof(double)),"clear fabric wall momentum reduction");if(embedded_count_)accumulate_wall_aperture_momentum_kernel<<<(embedded_count_+255)/256,256>>>(wall_edge_node_a_,wall_edge_node_b_,axis_,open_area_,pressure_gradient_factor_,special_velocity_,wall_momentum_scratch_,coarse_fine_count_,embedded_count_);if(fabric_wall_unique_carrier_count_)accumulate_wall_carrier_momentum_kernel<<<(fabric_wall_unique_carrier_count_+255)/256,256>>>(levels_,wall_unique_carrier_level_,wall_unique_carrier_index_,wall_unique_carrier_axis_,wall_unique_carrier_mass_,wall_momentum_scratch_,fabric_wall_unique_carrier_count_);check(cudaMemcpy(device_sum.data(),wall_momentum_scratch_,3*sizeof(double),cudaMemcpyDeviceToHost),"download fabric wall physical momentum");return device_sum;
	}
	void DeviceCompositeAmrProjection::transport_embedded_internal_momentum(Real dt)
	{
		if(!(dt>Real(0)))throw std::invalid_argument("conservative embedded aperture transport timestep");
		if(!embedded_count_)return;
		const std::size_t node_components=static_cast<std::size_t>(embedded_node_count_)*3;
		check(cudaMemset(eb_node_axis_sum_,0,node_components*sizeof(Real)),"clear conservative EB state momentum");
		check(cudaMemset(eb_node_axis_weight_,0,node_components*sizeof(Real)),"clear conservative EB state mass");
		check(cudaMemset(eb_node_diffusion_sum_,0,node_components*sizeof(Real)),"clear conservative EB advective increments");
		accumulate_eb_dual_momentum_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,axis_,open_area_,pressure_gradient_factor_,special_velocity_,eb_node_axis_sum_,eb_node_axis_weight_,coarse_fine_count_,embedded_count_);
		if(embedded_carrier_count_)accumulate_conservative_carrier_kernel<<<(embedded_carrier_count_+255)/256,256>>>(levels_,eb_carrier_node_,eb_carrier_level_,eb_carrier_index_,eb_carrier_axis_,eb_carrier_mass_,eb_node_axis_sum_,eb_node_axis_weight_,embedded_carrier_count_);
		advect_eb_node_momentum_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,direction_+coarse_fine_count_,open_area_+coarse_fine_count_,special_velocity_+coarse_fine_count_,eb_node_axis_sum_,eb_node_axis_weight_,dt,eb_node_diffusion_sum_,embedded_count_);
		apply_node_momentum_to_eb_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,axis_+coarse_fine_count_,eb_node_axis_weight_,eb_node_diffusion_sum_,special_velocity_+coarse_fine_count_,eb_transport_scratch_,embedded_count_);
		if(embedded_carrier_count_)apply_node_momentum_to_carrier_kernel<<<(embedded_carrier_count_+255)/256,256>>>(levels_,eb_carrier_node_,eb_carrier_level_,eb_carrier_index_,eb_carrier_axis_,eb_node_axis_weight_,eb_node_diffusion_sum_,embedded_carrier_count_);
		commit_eb_aperture_kernel<<<(embedded_count_+255)/256,256>>>(special_velocity_,eb_transport_scratch_,coarse_fine_count_,embedded_count_);
		check(cudaDeviceSynchronize(),"conservative compact EB momentum transport");
	}
	std::array<double,3> DeviceCompositeAmrProjection::embedded_transport_momentum()const
	{
		std::array<double,3> result{};if(!embedded_count_)return result;const std::size_t node_components=static_cast<std::size_t>(embedded_node_count_)*3;
		check(cudaMemset(eb_node_axis_sum_,0,node_components*sizeof(Real)),"clear EB momentum diagnostic sums");
		check(cudaMemset(eb_node_axis_weight_,0,node_components*sizeof(Real)),"clear EB momentum diagnostic masses");
		accumulate_eb_dual_momentum_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,axis_,open_area_,pressure_gradient_factor_,special_velocity_,eb_node_axis_sum_,eb_node_axis_weight_,coarse_fine_count_,embedded_count_);
		if(embedded_carrier_count_)accumulate_conservative_carrier_kernel<<<(embedded_carrier_count_+255)/256,256>>>(levels_,eb_carrier_node_,eb_carrier_level_,eb_carrier_index_,eb_carrier_axis_,eb_carrier_mass_,eb_node_axis_sum_,eb_node_axis_weight_,embedded_carrier_count_);
		if(embedded_carrier_count_)accumulate_unmapped_carrier_momentum_kernel<<<(embedded_carrier_count_+255)/256,256>>>(levels_,eb_carrier_node_,eb_carrier_level_,eb_carrier_index_,eb_carrier_axis_,eb_carrier_unmapped_mass_,eb_node_axis_sum_,embedded_carrier_count_);
		std::vector<Real> host(node_components);check(cudaMemcpy(host.data(),eb_node_axis_sum_,node_components*sizeof(Real),cudaMemcpyDeviceToHost),"download EB momentum diagnostic");for(int node=0;node<embedded_node_count_;++node)for(int component=0;component<3;++component)result[component]+=static_cast<double>(host[node*3+component]);return result;
	}
	double DeviceCompositeAmrProjection::max_embedded_transport_mass_imbalance()const
	{
		if(!embedded_count_)return 0;
		check(cudaMemset(eb_node_flux_balance_,0,static_cast<std::size_t>(embedded_node_count_)*sizeof(Real)),"clear compact-node flux balance");
		accumulate_embedded_node_flux_kernel<<<(embedded_count_+255)/256,256>>>(eb_node_a_,eb_node_b_,direction_+coarse_fine_count_,open_area_+coarse_fine_count_,special_velocity_+coarse_fine_count_,eb_node_flux_balance_,embedded_count_);
		if(embedded_regular_connection_count_)accumulate_regular_perimeter_flux_kernel<<<(embedded_regular_connection_count_+255)/256,256>>>(levels_,eb_regular_lower_node_,eb_regular_upper_node_,eb_regular_level_,eb_regular_index_,eb_regular_axis_,eb_regular_area_,eb_node_flux_balance_,embedded_regular_connection_count_);
		return max_abs_device_values(eb_node_flux_balance_,embedded_node_count_,max_abs_scratch_);
	}
	void DeviceCompositeAmrProjection::upload_special_fluxes(const CompositeAmrFluxes& host)
	{
		if(host.coarse_fine_velocity.size()!=static_cast<std::size_t>(coarse_fine_count_)||host.embedded_velocity.size()!=static_cast<std::size_t>(special_count_-coarse_fine_count_))throw std::invalid_argument("composite special flux upload size");std::vector<Real> values(special_count_);for(int q=0;q<coarse_fine_count_;++q)values[q]=static_cast<Real>(host.coarse_fine_velocity[q]);for(int q=coarse_fine_count_;q<special_count_;++q)values[q]=static_cast<Real>(host.embedded_velocity[q-coarse_fine_count_]);if(special_count_)check(cudaMemcpy(special_velocity_,values.data(),values.size()*sizeof(Real),cudaMemcpyHostToDevice),"upload composite special fluxes");
	}
	void DeviceCompositeAmrProjection::download_special_fluxes(CompositeAmrFluxes& host)const
	{
		std::vector<Real> values(special_count_);if(special_count_)check(cudaMemcpy(values.data(),special_velocity_,values.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download composite special fluxes");host.coarse_fine_velocity.resize(coarse_fine_count_);host.embedded_velocity.resize(special_count_-coarse_fine_count_);for(int q=0;q<coarse_fine_count_;++q)host.coarse_fine_velocity[q]=values[q];for(int q=coarse_fine_count_;q<special_count_;++q)host.embedded_velocity[q-coarse_fine_count_]=values[q];
	}
	void DeviceCompositeAmrProjection::download_embedded_node_gradients(std::vector<Real>& gradients,std::vector<std::uint8_t>* component_rank)const
	{
		gradients.resize(static_cast<std::size_t>(embedded_node_count_)*9);if(!gradients.empty())check(cudaMemcpy(gradients.data(),eb_node_gradient_sum_,gradients.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download EB least-squares gradients");if(component_rank)*component_rank=embedded_gradient_rank_;
	}
	double DeviceCompositeAmrProjection::max_abs_special_velocity()const{return max_abs_device_values(special_velocity_,special_count_,max_abs_scratch_);}
	double DeviceCompositeAmrProjection::max_embedded_cfl_rate()const
	{
		if(!embedded_count_)return 0;embedded_cfl_rate_kernel<<<(embedded_count_+255)/256,256>>>(special_velocity_,eb_transport_length_,eb_transport_scratch_,coarse_fine_count_,embedded_count_);return max_abs_device_values(eb_transport_scratch_,embedded_count_,max_abs_scratch_);
	}

	namespace
	{
		void launch_integrated_flux(const DeviceCompositeAmrFluxLevelView* levels,const std::vector<int>& brick_counts,int level_count,int brick_size,const unsigned char* active,const std::uint8_t* cut_face_mask,int storage_size,const int* special_incident_dof,const int* special_incident_offset,const int* special_incident_edge,const std::int8_t* special_incident_sign,int special_incident_dof_count,const Real* area,const Real* velocity,Real* integrated)
		{
			check(cudaMemset(integrated,0,static_cast<std::size_t>(storage_size)*sizeof(Real)),"clear composite integrated flux");for(int level=0;level<level_count;++level){const int work=brick_counts[level]*brick_size*brick_size*brick_size;if(work)structured_integrated_flux_kernel<<<(work+255)/256,256>>>(levels,level,brick_size,active,cut_face_mask,integrated);}if(special_incident_dof_count)special_integrated_flux_gather_kernel<<<(special_incident_dof_count+255)/256,256>>>(special_incident_dof,special_incident_offset,special_incident_edge,special_incident_sign,special_incident_dof_count,area,velocity,integrated);
		}
	}

	void DeviceCompositeAmrProjection::compute_divergence()
	{
		launch_integrated_flux(levels_,brick_counts_,level_count_,brick_size_,active_,cut_face_mask_,storage_size_,special_incident_dof_,special_incident_offset_,special_incident_edge_,special_incident_sign_,special_incident_dof_count_,open_area_,special_velocity_,integrated_);normalize_divergence_rhs_kernel<<<(storage_size_+255)/256,256>>>(integrated_,volume_,active_,Real(0),divergence_,nullptr,storage_size_);check(cudaDeviceSynchronize(),"compute composite AMR divergence");
	}
	void DeviceCompositeAmrProjection::build_projection_rhs(Real rho,Real dt)
	{
		if(!(rho>Real(0))||!(dt>Real(0)))throw std::invalid_argument("composite projection rho/dt");launch_integrated_flux(levels_,brick_counts_,level_count_,brick_size_,active_,cut_face_mask_,storage_size_,special_incident_dof_,special_incident_offset_,special_incident_edge_,special_incident_sign_,special_incident_dof_count_,open_area_,special_velocity_,integrated_);normalize_divergence_rhs_kernel<<<(storage_size_+255)/256,256>>>(integrated_,volume_,active_,rho/dt,divergence_,rhs_,storage_size_);check(cudaDeviceSynchronize(),"build composite AMR projection RHS");
	}
	void DeviceCompositeAmrProjection::correct_fluxes(Real rho,Real dt)
	{
		if(!(rho>Real(0))||!(dt>Real(0)))throw std::invalid_argument("composite projection correction rho/dt");const Real scale=dt/rho;for(int level=0;level<level_count_;++level){const int work=brick_counts_[level]*brick_size_*brick_size_*brick_size_;if(work)structured_flux_correction_kernel<<<(work+255)/256,256>>>(levels_,level,brick_size_,outlet_,active_,cut_face_mask_,pressure_,scale);}if(pressure_gradient_node_count_)reconstruct_pressure_gradient_kernel<<<(pressure_gradient_node_count_+255)/256,256>>>(pressure_gradient_node_dof_,pressure_gradient_offset_,pressure_gradient_neighbour_,pressure_gradient_weight_,pressure_,pressure_gradient_,pressure_gradient_node_count_);if(regular_pressure_correction_count_)regular_nonorthogonal_flux_correction_kernel<<<(regular_pressure_correction_count_+255)/256,256>>>(levels_,regular_correction_level_,regular_correction_index_,regular_correction_mirror_index_,regular_correction_axis_,regular_correction_lower_,regular_correction_upper_,regular_correction_lower_node_,regular_correction_upper_node_,regular_correction_two_point_delta_,regular_correction_vector_,regular_correction_upper_weight_,pressure_gradient_,pressure_,scale,regular_pressure_correction_count_);if(special_count_){special_flux_correction_kernel<<<(special_count_+255)/256,256>>>(first_dof_,second_dof_,direction_,pressure_gradient_factor_,special_velocity_,special_count_,active_,pressure_,scale);if(pressure_gradient_node_count_)special_nonorthogonal_flux_correction_kernel<<<(special_count_+255)/256,256>>>(first_dof_,second_dof_,direction_,special_lower_gradient_node_,special_upper_gradient_node_,special_nonorthogonal_correction_,special_upper_gradient_weight_,pressure_gradient_,special_velocity_,special_count_,active_,scale);}if(coarse_fine_count_){check(cudaMemset(cf_group_sum_,0,static_cast<std::size_t>(coarse_fine_group_count_)*sizeof(Real)),"clear coarse/fine group flux sums");scatter_coarse_fine_tiles_kernel<<<(coarse_fine_count_+255)/256,256>>>(levels_,cf_fine_level_,cf_fine_index_,axis_,cf_group_,open_area_,special_velocity_,cf_group_sum_,coarse_fine_count_);scatter_coarse_fine_groups_kernel<<<(coarse_fine_group_count_+255)/256,256>>>(levels_,cf_group_level_,cf_group_index_,cf_group_axis_,cf_group_area_,cf_group_sum_,coarse_fine_group_count_);}for(int level=0;level<level_count_;++level){const int work=brick_counts_[level]*brick_size_*brick_size_*brick_size_;if(work)scatter_regular_pressure_kernel<<<(work+255)/256,256>>>(levels_,level,brick_size_,active_,pressure_);}check(cudaDeviceSynchronize(),"correct composite AMR fluxes");
	}
	AmrGpuSolveResult DeviceCompositeAmrProjection::project(Real rho,Real dt,double tolerance,int max_iterations,bool warm_start)
	{
		build_projection_rhs(rho,dt);AmrGpuSolveResult result=solver_.solve(pressure_,rhs_,tolerance,max_iterations,warm_start);if(result.converged){correct_fluxes(rho,dt);compute_divergence();}return result;
	}
	void DeviceCompositeAmrProjection::download_divergence(std::vector<Real>& host)const{host.resize(storage_size_);check(cudaMemcpy(host.data(),divergence_,host.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download composite divergence");}
	void DeviceCompositeAmrProjection::download_pressure(std::vector<Real>& host)const{host.resize(storage_size_);check(cudaMemcpy(host.data(),pressure_,host.size()*sizeof(Real),cudaMemcpyDeviceToHost),"download composite pressure");}
}
