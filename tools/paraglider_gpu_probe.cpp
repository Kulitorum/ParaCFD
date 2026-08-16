#include "core/fluid/eb_pressure.h"

#include <cuda_runtime.h>

#include <chrono>
#include <cstdio>
#include <vector>

using namespace paracfd::core;

int main()
{
	EmbeddedBoundary eb;eb.grid={{0,0,0},64,64,64,0.05};eb.cells.resize(eb.grid.cell_count());EbPressureSystem sys=build_eb_pressure_system(eb,true);DeviceEbPressureOperator op(sys);std::vector<Real> host(sys.storage_size,Real(0.25));Real *p=nullptr,*out=nullptr;cudaMalloc(&p,host.size()*sizeof(Real));cudaMalloc(&out,host.size()*sizeof(Real));cudaMemcpy(p,host.data(),host.size()*sizeof(Real),cudaMemcpyHostToDevice);
	for(int i=0;i<5;++i)op.apply(p,out);constexpr int repeats=50;auto t0=std::chrono::steady_clock::now();for(int i=0;i<repeats;++i)op.apply(p,out);auto t1=std::chrono::steady_clock::now();std::vector<Real> result(host.size());cudaMemcpy(result.data(),out,result.size()*sizeof(Real),cudaMemcpyDeviceToHost);size_t free_b=0,total_b=0;cudaMemGetInfo(&free_b,&total_b);cudaFree(p);cudaFree(out);double ms=std::chrono::duration<double,std::milli>(t1-t0).count()/repeats;
	std::printf("[paraglider-gpu] scalar=%s grid=64^3 cells=%d operator=%.4f ms/apply (includes synchronization)\n",sizeof(Real)==4?"FP32":"FP64",eb.grid.cell_count(),ms);
	std::printf("[paraglider-gpu] operator metadata=%.3f MiB pressure+output=%.3f MiB GPU used=%.1f/%.1f MiB\n",op.bytes()/(1024.0*1024.0),2.0*host.size()*sizeof(Real)/(1024.0*1024.0),(total_b-free_b)/(1024.0*1024.0),total_b/(1024.0*1024.0));
	return ms>0&&result.back()>0?0:1;
}
