#include "core/fluid/eb_pressure.h"

#include <cuda_runtime.h>

#include <thrust/iterator/counting_iterator.h>
#include <thrust/execution_policy.h>
#include <thrust/transform_reduce.h>

#include <stdexcept>
#include <string>
#include <cmath>

namespace paracfd::core
{
	namespace
	{
		void ck(cudaError_t e,const char*w){if(e!=cudaSuccess)throw std::runtime_error(std::string(w)+": "+cudaGetErrorString(e));}
		template<class T>T* up(const std::vector<T>& h){if(h.empty())return nullptr;T*d=nullptr;ck(cudaMalloc(&d,h.size()*sizeof(T)),"cudaMalloc EB operator");ck(cudaMemcpy(d,h.data(),h.size()*sizeof(T),cudaMemcpyHostToDevice),"cudaMemcpy EB operator");return d;}
		__global__ void regular_apply(const Real*p,Real*out,UniformEbGrid g,const int*cell_dof,const std::uint8_t*cut_face_mask,bool outlet)
		{
			int cell=blockIdx.x*blockDim.x+threadIdx.x;if(cell>=g.cell_count())return;int a=cell_dof[cell];if(a<0)return;int i=cell%g.nx,j=(cell/g.nx)%g.ny,k=cell/(g.nx*g.ny);Real v=0;const Real c=(Real)g.h;
			auto edge=[&](int nb){int b=cell_dof[nb];if(b>=0)v+=c*(p[a]-p[b]);};if(i>0&&!(cut_face_mask[g.cell_index(i-1,j,k)]&1u))edge(g.cell_index(i-1,j,k));if(i+1<g.nx&&!(cut_face_mask[cell]&1u))edge(g.cell_index(i+1,j,k));if(j>0&&!(cut_face_mask[g.cell_index(i,j-1,k)]&2u))edge(g.cell_index(i,j-1,k));if(j+1<g.ny&&!(cut_face_mask[cell]&2u))edge(g.cell_index(i,j+1,k));if(k>0&&!(cut_face_mask[g.cell_index(i,j,k-1)]&4u))edge(g.cell_index(i,j,k-1));if(k+1<g.nz&&!(cut_face_mask[cell]&4u))edge(g.cell_index(i,j,k+1));if(outlet&&i==g.nx-1)v+=(Real)(2*g.h)*p[a];out[a]=v;
		}
		__global__ void aperture_apply(const Real*p,Real*out,const FragmentRef*ar,const FragmentRef*br,const Real*area,const Real*dist,int n,const int*cell_dof,const int*fragment_dof)
		{
			int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;auto dof=[&](FragmentRef r){return r<0?cell_dof[-r-1]:fragment_dof[r-1];};int a=dof(ar[i]),b=dof(br[i]);if(a<0||b<0||a==b)return;Real v=(area[i]/dist[i])*(p[a]-p[b]);atomicAdd(out+a,v);atomicAdd(out+b,-v);
		}
		__global__ void precondition_kernel(const Real*r,Real*z,Real*d,const Real*diag,const unsigned char*active,int n){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n){Real v=active[i]&&diag[i]>Real(0)?r[i]/diag[i]:Real(0);z[i]=v;d[i]=v;}}
		__global__ void subtract_kernel(Real*r,const Real*value,const unsigned char*active,int n){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)r[i]=active[i]?r[i]-value[i]:Real(0);}
		__global__ void update_xr_kernel(Real*x,Real*r,const Real*d,const Real*Ad,Real alpha,const unsigned char*active,int n){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n&&active[i]){x[i]+=alpha*d[i];r[i]-=alpha*Ad[i];}}
		__global__ void precondition_only_kernel(const Real*r,Real*z,const Real*diag,const unsigned char*active,int n){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)z[i]=active[i]&&diag[i]>Real(0)?r[i]/diag[i]:Real(0);}
		__global__ void update_direction_kernel(const Real*z,Real*d,Real beta,const unsigned char*active,int n){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)d[i]=active[i]?z[i]+beta*d[i]:Real(0);}
		struct DotActive
		{
			const Real*a,*b;const unsigned char*active;
			__host__ __device__ double operator()(int i)const{return active[i]?(double)a[i]*(double)b[i]:0.0;}
		};
		double dot_active(const Real*a,const Real*b,const unsigned char*active,int n){auto c=thrust::counting_iterator<int>(0);return thrust::transform_reduce(thrust::device,c,c+n,DotActive{a,b,active},0.0,thrust::plus<double>());}
	}

	DeviceEbPressureOperator::DeviceEbPressureOperator(const EbPressureSystem&s)
	{
		grid_=s.eb->grid;storage_size_=s.storage_size;outlet_=s.pressure_outlet_xmax;cell_dof_=up(s.cell_dof);fragment_dof_=up(s.fragment_dof);std::vector<std::uint8_t>mask=s.eb->cut_face_mask;if(mask.empty())mask.assign(grid_.cell_count(),0);cut_face_mask_=up(mask);std::vector<FragmentRef>a,b;std::vector<Real>area,dist;for(const auto&x:s.eb->apertures){a.push_back(x.fragment_a);b.push_back(x.fragment_b);area.push_back((Real)x.area);dist.push_back((Real)std::max(1e-12,std::sqrt(length2(s.eb->fragment_centroid(x.fragment_a)-s.eb->fragment_centroid(x.fragment_b)))));}aperture_count_=(int)a.size();aperture_a_=up(a);aperture_b_=up(b);aperture_area_=up(area);aperture_distance_=up(dist);bytes_=s.cell_dof.size()*sizeof(int)+s.fragment_dof.size()*sizeof(int)+mask.size()*sizeof(std::uint8_t)+a.size()*(2*sizeof(FragmentRef)+2*sizeof(Real));
	}
	DeviceEbPressureOperator::~DeviceEbPressureOperator(){for(void*p:{(void*)cell_dof_,(void*)fragment_dof_,(void*)cut_face_mask_,(void*)aperture_a_,(void*)aperture_b_,(void*)aperture_area_,(void*)aperture_distance_})if(p)cudaFree(p);}
	void DeviceEbPressureOperator::apply(const Real*p,Real*out)const{ck(cudaMemset(out,0,storage_size_*sizeof(Real)),"clear EB apply");regular_apply<<<(grid_.cell_count()+255)/256,256>>>(p,out,grid_,cell_dof_,cut_face_mask_,outlet_);if(aperture_count_)aperture_apply<<<(aperture_count_+255)/256,256>>>(p,out,aperture_a_,aperture_b_,aperture_area_,aperture_distance_,aperture_count_,cell_dof_,fragment_dof_);ck(cudaDeviceSynchronize(),"EB apply");}

	DeviceEbPressureSolver::DeviceEbPressureSolver(const EbPressureSystem&s):op_(s),n_(s.storage_size)
	{
		std::vector<double> dd;s.diagonal_cpu(dd);std::vector<Real>d(dd.size());for(std::size_t i=0;i<d.size();++i)d[i]=(Real)dd[i];
		r_=up(std::vector<Real>(n_));z_=up(std::vector<Real>(n_));direction_=up(std::vector<Real>(n_));Ad_=up(std::vector<Real>(n_));diagonal_=up(d);active_=up(s.active);bytes_=op_.bytes()+4ull*n_*sizeof(Real)+d.size()*sizeof(Real)+s.active.size();
	}
	DeviceEbPressureSolver::~DeviceEbPressureSolver(){for(void*p:{(void*)r_,(void*)z_,(void*)direction_,(void*)Ad_,(void*)diagonal_,(void*)active_})if(p)cudaFree(p);}
	EbGpuSolveResult DeviceEbPressureSolver::solve(Real*x,const Real*b,double tol,int maxit,bool warm)
	{
		EbGpuSolveResult out;if(!warm)ck(cudaMemset(x,0,n_*sizeof(Real)),"clear EB pressure");if(warm){op_.apply(x,Ad_);ck(cudaMemcpy(r_,b,n_*sizeof(Real),cudaMemcpyDeviceToDevice),"copy EB rhs");subtract_kernel<<<(n_+255)/256,256>>>(r_,Ad_,active_,n_);}else ck(cudaMemcpy(r_,b,n_*sizeof(Real),cudaMemcpyDeviceToDevice),"initial EB residual");
		const double b2=dot_active(b,b,active_,n_);if(!(b2>0)){out.converged=true;return out;}const double initial_r2=dot_active(r_,r_,active_,n_);out.relative_residual=std::sqrt(initial_r2/b2);if(out.relative_residual<=tol){out.converged=true;return out;}precondition_kernel<<<(n_+255)/256,256>>>(r_,z_,direction_,diagonal_,active_,n_);double rz=dot_active(r_,z_,active_,n_);
		for(int it=0;it<maxit;++it){op_.apply(direction_,Ad_);double dAd=dot_active(direction_,Ad_,active_,n_);if(!(dAd>0))break;Real alpha=(Real)(rz/dAd);update_xr_kernel<<<(n_+255)/256,256>>>(x,r_,direction_,Ad_,alpha,active_,n_);double r2=dot_active(r_,r_,active_,n_);out.iterations=it+1;out.relative_residual=std::sqrt(r2/b2);if(out.relative_residual<=tol){out.converged=true;break;}precondition_only_kernel<<<(n_+255)/256,256>>>(r_,z_,diagonal_,active_,n_);double rz2=dot_active(r_,z_,active_,n_);Real beta=(Real)(rz2/rz);update_direction_kernel<<<(n_+255)/256,256>>>(z_,direction_,beta,active_,n_);rz=rz2;}ck(cudaDeviceSynchronize(),"EB pressure solve");return out;
	}
}
