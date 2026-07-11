// slice_gl.cu — see slice_gl.h. The one translation unit that includes CUDA-GL interop
// (and therefore GL) headers; keeps them out of the Qt/MSVC widget compile.
//
// Each registered buffer carries its OWN non-blocking CUDA stream so the main-thread
// map/fill/unmap runs concurrently with the worker thread's simulation step (which uses
// the default stream) — this is what lets the viewer hit display-rate fps regardless of
// how expensive a physics step is.
#include "gui/slice_gl.h"

// cuda_gl_interop.h pulls in <GL/gl.h>, which on Windows needs windows.h first.
#if defined(_WIN32)
#	define WIN32_LEAN_AND_MEAN
#	define NOMINMAX
#	include <windows.h>
#endif
#include <cuda_gl_interop.h>
#include <cuda_runtime.h>

#include <cstdio>

namespace scour::gui
{
	namespace
	{
		struct SliceResource
		{
			cudaGraphicsResource* res = nullptr;
			cudaStream_t stream = nullptr;
			float* range_dev = nullptr; // persistent 3-float scratch for the auto-range reduction
		};
	}

	void* slice_gl_register(unsigned int vbo)
	{
		SliceResource* h = new SliceResource();
		cudaError_t e = cudaGraphicsGLRegisterBuffer(&h->res, vbo, cudaGraphicsRegisterFlagsWriteDiscard);
		if (e != cudaSuccess)
		{
			std::fprintf(stderr, "slice_gl_register: cudaGraphicsGLRegisterBuffer failed: %s\n", cudaGetErrorString(e));
			delete h;
			return nullptr;
		}
		// Dedicated non-blocking stream => render does not serialise behind the worker's
		// default-stream step kernels.
		cudaStreamCreateWithFlags(&h->stream, cudaStreamNonBlocking);
		cudaMalloc(&h->range_dev, 3 * sizeof(float)); // reused every auto-range reduction
		return h;
	}

	void slice_gl_unregister(void* handle)
	{
		if (!handle) return;
		SliceResource* h = reinterpret_cast<SliceResource*>(handle);
		if (h->range_dev) cudaFree(h->range_dev);
		if (h->stream) cudaStreamDestroy(h->stream);
		if (h->res) cudaGraphicsUnregisterResource(h->res);
		delete h;
	}

	bool slice_gl_fill(void* handle, const double* u, const double* v, const double* w, const double* p,
		const SliceParams& sp)
	{
		if (!handle) return false;
		SliceResource* h = reinterpret_cast<SliceResource*>(handle);
		cudaError_t e = cudaGraphicsMapResources(1, &h->res, h->stream);
		if (e != cudaSuccess)
		{
			std::fprintf(stderr, "slice_gl_fill: map failed: %s\n", cudaGetErrorString(e));
			return false;
		}
		float4* dptr = nullptr;
		size_t bytes = 0;
		e = cudaGraphicsResourceGetMappedPointer(reinterpret_cast<void**>(&dptr), &bytes, h->res);
		if (e != cudaSuccess)
		{
			std::fprintf(stderr, "slice_gl_fill: get pointer failed: %s\n", cudaGetErrorString(e));
			cudaGraphicsUnmapResources(1, &h->res, h->stream);
			return false;
		}
		size_t need = (size_t)sp.nu * (size_t)sp.nv * sizeof(float4);
		if (bytes >= need)
			slice_fill_gpu(u, v, w, p, sp, dptr, h->stream);
		cudaError_t k = cudaGetLastError();
		cudaGraphicsUnmapResources(1, &h->res, h->stream); // stream-ordered: writes visible to GL after this
		if (k != cudaSuccess)
		{
			std::fprintf(stderr, "slice_gl_fill: kernel error: %s\n", cudaGetErrorString(k));
			return false;
		}
		return bytes >= need;
	}

	bool slice_gl_reduce(void* handle, const double* u, const double* v, const double* w, const double* p,
		const unsigned char* solid, scour::core::MacGrid g, Field field, FieldRange* out)
	{
		if (!handle || !out) return false;
		SliceResource* h = reinterpret_cast<SliceResource*>(handle);
		if (!h->range_dev) return false;

		slice_reduce_gpu(u, v, w, p, solid, g, field, h->range_dev, h->stream);
		float host3[3] = { 0.0f, 0.0f, 0.0f };
		cudaError_t e = cudaMemcpyAsync(host3, h->range_dev, 3 * sizeof(float), cudaMemcpyDeviceToHost, h->stream);
		if (e == cudaSuccess) e = cudaStreamSynchronize(h->stream); // waits ONLY on the interop stream
		if (e != cudaSuccess)
		{
			std::fprintf(stderr, "slice_gl_reduce: %s\n", cudaGetErrorString(e));
			out->valid = false;
			return false;
		}
		out->field_min = host3[0];
		out->field_max = host3[1];
		out->speed_max = host3[2];
		out->valid = host3[1] >= host3[0]; // max>=min ⇒ at least one finite fluid cell contributed
		return out->valid;
	}
}
