// GPU/CPU parity oracle for the reusable uniform MAC mathematics.
//
// The paraglider solver uses a different FP32 block-AMR/embedded-boundary path,
// but these original FP64 kernels remain useful reference implementations.  The
// snapshots are verification-only: this executable deliberately has no --bless
// mode, so a changed numerical result has to be reviewed rather than accepted by
// an accidental command-line flag.
#include "core/fluid/mac_grid.h"
#include "core/fluid/mac_ops.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <vector>

using namespace paracfd::core;

#define CU(expr) do { cudaError_t e_=(expr); if(e_!=cudaSuccess) { \
	std::fprintf(stderr,"CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); \
	std::exit(2); } } while(0)

namespace
{
	std::vector<void*> allocations;

	double* device_copy(const std::vector<double>& host)
	{
		double* ptr=nullptr;CU(cudaMalloc(&ptr,host.size()*sizeof(double)));
		CU(cudaMemcpy(ptr,host.data(),host.size()*sizeof(double),cudaMemcpyHostToDevice));
		allocations.push_back(ptr);return ptr;
	}

	double* device_zero(std::size_t count)
	{
		double* ptr=nullptr;CU(cudaMalloc(&ptr,count*sizeof(double)));
		CU(cudaMemset(ptr,0,count*sizeof(double)));allocations.push_back(ptr);return ptr;
	}

	std::vector<double> download(const double* ptr,std::size_t count)
	{
		std::vector<double> host(count);
		CU(cudaMemcpy(host.data(),ptr,count*sizeof(double),cudaMemcpyDeviceToHost));return host;
	}

	std::vector<double> filled(std::size_t count,std::uint32_t seed,double lo,double hi)
	{
		std::mt19937 rng(seed);std::uniform_real_distribution<double> dist(lo,hi);
		std::vector<double> values(count);for(double& value:values)value=dist(rng);return values;
	}

	double relative_max_norm(const std::vector<double>& actual,const std::vector<double>& reference)
	{
		double numerator=0.0,denominator=1e-30;
		const std::size_t count=std::min(actual.size(),reference.size());
		for(std::size_t i=0;i<count;++i)
		{
			numerator=std::max(numerator,std::abs(actual[i]-reference[i]));
			denominator=std::max(denominator,std::abs(reference[i]));
		}
		return numerator/denominator;
	}

	struct Reporter
	{
		std::string golden_dir="tests/golden";int total=0,failed=0;

		void report(const std::string& name,const std::vector<double>& gpu,const std::vector<double>& cpu)
		{
			++total;const double parity=relative_max_norm(gpu,cpu);bool ok=parity<=1e-5&&gpu.size()==cpu.size();
			std::vector<double> golden;std::ifstream input(golden_dir+"/"+name+".f64",std::ios::binary);
			if(input)
			{
				input.seekg(0,std::ios::end);const std::streamoff bytes=input.tellg();input.seekg(0,std::ios::beg);
				golden.resize(static_cast<std::size_t>(bytes)/sizeof(double));
				input.read(reinterpret_cast<char*>(golden.data()),bytes);
			}
			const double snapshot=golden.empty()?1.0:relative_max_norm(gpu,golden);
			ok=ok&&golden.size()==gpu.size()&&snapshot<=1e-5;
			if(!ok)++failed;
			std::printf("  %-22s parity=%.2e snapshot=%.2e %s\n",name.c_str(),parity,snapshot,ok?"PASS":"FAIL");
		}
	};
}

int main(int argc,char** argv)
{
	Reporter reporter;
	for(int i=1;i<argc;++i)
	{
		if(std::strcmp(argv[i],"--golden")==0&&i+1<argc)reporter.golden_dir=argv[++i];
		else { std::fprintf(stderr,"usage: parity_probe [--golden directory]\n");return 2; }
	}

	MacGrid grid;grid.nx=16;grid.ny=16;grid.nz=16;grid.h=0.05;
	const int np=grid.p_count(),nu_count=grid.u_count(),nv_count=grid.v_count(),nw_count=grid.w_count();
	const int scratch_count=std::max({nu_count,nv_count,nw_count});
	const double rho=1.225,viscosity=1e-3,dt=1e-3,cs=0.16,omega=2.0/3.0;
	BC bc;
	const auto u=filled(nu_count,101,-1,1),v=filled(nv_count,102,-1,1),w=filled(nw_count,103,-1,1);
	const auto p=filled(np,104,-1,1),rhs=filled(np,105,-1,1),nut=filled(np,106,0,0.05);

	std::printf("parity_probe: reusable FP64 MAC reference, %dx%dx%d, immutable snapshots in %s\n",
		grid.nx,grid.ny,grid.nz,reporter.golden_dir.c_str());
	{
		double* du=device_copy(u),*dv=device_copy(v),*dw=device_copy(w);
		double* out_u=device_zero(nu_count),*out_v=device_zero(nv_count),*out_w=device_zero(nw_count),*scratch=device_zero(scratch_count);
		mac_advect_gpu(du,dv,dw,out_u,out_v,out_w,scratch,grid,bc,dt,1);
		std::vector<double> cpu_u(nu_count),cpu_v(nv_count),cpu_w(nw_count);mac_advect_cpu(u,v,w,cpu_u,cpu_v,cpu_w,grid,bc,dt,1);
		reporter.report("mac_advect_u",download(out_u,nu_count),cpu_u);
		reporter.report("mac_advect_v",download(out_v,nv_count),cpu_v);
		reporter.report("mac_advect_w",download(out_w,nw_count),cpu_w);
	}
	{
		double* du=device_copy(u),*dv=device_copy(v),*dw=device_copy(w),*out=device_zero(np);
		smagorinsky_nut_gpu(du,dv,dw,out,grid,bc,cs);std::vector<double> cpu(np);
		smagorinsky_nut_cpu(u,v,w,cpu,grid,bc,cs);reporter.report("smagorinsky_nut",download(out,np),cpu);
	}
	{
		double* du=device_copy(u),*dv=device_copy(v),*dw=device_copy(w),*dn=device_copy(nut);
		double* out_u=device_zero(nu_count),*out_v=device_zero(nv_count),*out_w=device_zero(nw_count);
		mac_diffuse_gpu(du,dv,dw,out_u,out_v,out_w,dn,grid,bc,dt,viscosity);
		std::vector<double> cpu_u(nu_count),cpu_v(nv_count),cpu_w(nw_count);mac_diffuse_cpu(u,v,w,cpu_u,cpu_v,cpu_w,&nut,grid,bc,dt,viscosity);
		reporter.report("mac_diffuse_u",download(out_u,nu_count),cpu_u);
		reporter.report("mac_diffuse_v",download(out_v,nv_count),cpu_v);
		reporter.report("mac_diffuse_w",download(out_w,nw_count),cpu_w);
	}
	{
		double* du=device_copy(u),*dv=device_copy(v),*dw=device_copy(w),*out=device_zero(np);
		poisson_rhs_gpu(du,dv,dw,out,grid,rho,dt);std::vector<double> cpu(np);poisson_rhs_cpu(u,v,w,cpu,grid,rho,dt);
		reporter.report("poisson_rhs",download(out,np),cpu);
	}
	{
		double* dp=device_copy(p),*out=device_zero(np);poisson_apply_gpu(dp,out,grid);
		std::vector<double> cpu(np);poisson_apply_cpu(p,cpu,grid);reporter.report("poisson_apply",download(out,np),cpu);
	}
	{
		double* dp=device_copy(p),*dr=device_copy(rhs),*scratch=device_zero(np);jacobi_smooth_gpu(dp,dr,scratch,grid,omega,2);
		auto cpu=p;jacobi_smooth_cpu(cpu,rhs,grid,omega,2);reporter.report("jacobi_smooth",download(dp,np),cpu);
	}
	{
		double* dp=device_copy(p),*dr=device_copy(rhs);gs_band_gpu(dp,dr,grid,2,2,true);
		auto cpu=p;gs_band_cpu(cpu,rhs,grid,2,2,true);reporter.report("gs_band",download(dp,np),cpu);
	}
	{
		double* dp=device_copy(p),*dr=device_copy(rhs),*out=device_zero(np);poisson_residual_gpu(dp,dr,out,grid);
		std::vector<double> cpu(np);poisson_residual_cpu(p,rhs,cpu,grid);reporter.report("poisson_residual",download(out,np),cpu);
	}
	{
		MacGrid coarse=grid;coarse.nx/=2;coarse.ny/=2;coarse.nz/=2;coarse.h*=2;
		const auto fine=filled(np,107,-1,1);double* df=device_copy(fine),*dc=device_zero(coarse.p_count());restrict_gpu(df,dc,grid,coarse);
		std::vector<double> cpu(coarse.p_count());restrict_cpu(fine,cpu,grid,coarse);reporter.report("restrict",download(dc,coarse.p_count()),cpu);
	}
	{
		MacGrid coarse=grid;coarse.nx/=2;coarse.ny/=2;coarse.nz/=2;coarse.h*=2;
		const auto coarse_values=filled(coarse.p_count(),108,-1,1),fine=filled(np,109,-1,1);
		double* dc=device_copy(coarse_values),*df=device_copy(fine);prolong_add_gpu(dc,df,coarse,grid);
		auto cpu=fine;prolong_add_cpu(coarse_values,cpu,coarse,grid);reporter.report("prolong_add",download(df,np),cpu);
	}
	{
		double* du=device_copy(u),*dv=device_copy(v),*dw=device_copy(w),*dp=device_copy(p);subtract_gradient_gpu(du,dv,dw,dp,grid,rho,dt);
		auto cpu_u=u,cpu_v=v,cpu_w=w;subtract_gradient_cpu(cpu_u,cpu_v,cpu_w,p,grid,rho,dt);
		reporter.report("subtract_gradient_u",download(du,nu_count),cpu_u);
		reporter.report("subtract_gradient_v",download(dv,nv_count),cpu_v);
		reporter.report("subtract_gradient_w",download(dw,nw_count),cpu_w);
	}
	{
		double* du=device_copy(u),*dv=device_copy(v),*dw=device_copy(w),*scratch=device_zero(np);
		const double gpu=max_abs_divergence_gpu(du,dv,dw,scratch,grid),cpu=max_abs_divergence_cpu(u,v,w,grid);
		reporter.report("max_abs_divergence",{gpu},{cpu});
	}

	for(void* ptr:allocations)cudaFree(ptr);
	std::printf("parity_probe: %d/%d checks passed\n",reporter.total-reporter.failed,reporter.total);
	return reporter.failed==0?0:1;
}
