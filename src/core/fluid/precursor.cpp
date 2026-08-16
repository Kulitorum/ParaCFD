// precursor.cpp — see precursor.h. Host library + linear-in-time looped replay; the
// replay uploads interpolated planes and applies them via the shared plane_ops launchers.
#include "core/fluid/precursor.h"
#include "core/fluid/plane_ops.h"

#include <cuda_runtime.h>

#include <cmath>

namespace paracfd::core
{
	PrecursorLibrary::PrecursorLibrary(MacGrid g, double record_hz) : g_(g)
	{
		record_dt_ = (record_hz > 0.0) ? 1.0 / record_hz : 0.1;
	}

	bool PrecursorLibrary::maybe_record(double t, const std::vector<double>& up, const std::vector<double>& vp, const std::vector<double>& wp)
	{
		if (t + 1e-12 < next_t_) return false;
		PrecursorFrame f; f.t = t; f.up = up; f.vp = vp; f.wp = wp;
		frames_.push_back(std::move(f));
		next_t_ += record_dt_;
		return true;
	}

	void PrecursorLibrary::sample(double t, std::vector<double>& up, std::vector<double>& vp, std::vector<double>& wp) const
	{
		int M = (int)frames_.size();
		size_t nu = (size_t)g_.ny * g_.nz, nv = (size_t)(g_.ny + 1) * g_.nz, nw = (size_t)g_.ny * (g_.nz + 1);
		up.assign(nu, 0.0); vp.assign(nv, 0.0); wp.assign(nw, 0.0);
		if (M == 0) return;
		if (M == 1) { up = frames_[0].up; vp = frames_[0].vp; wp = frames_[0].wp; return; }
		double P = M * record_dt_;
		double tm = std::fmod(t, P); if (tm < 0.0) tm += P;
		double g = tm / record_dt_;
		int i0 = (int)std::floor(g);
		double f = g - i0;
		int a = i0 % M, b = (i0 + 1) % M; // wrap last→first for a seamless loop
		const PrecursorFrame& A = frames_[a];
		const PrecursorFrame& B = frames_[b];
		for (size_t q = 0; q < nu; ++q) up[q] = A.up[q] * (1.0 - f) + B.up[q] * f;
		for (size_t q = 0; q < nv; ++q) vp[q] = A.vp[q] * (1.0 - f) + B.vp[q] * f;
		for (size_t q = 0; q < nw; ++q) wp[q] = A.wp[q] * (1.0 - f) + B.wp[q] * f;
	}

	PrecursorReplay::PrecursorReplay(const PrecursorLibrary& lib) : lib_(&lib), g_(lib.grid())
	{
		size_t nu = (size_t)g_.ny * g_.nz, nv = (size_t)(g_.ny + 1) * g_.nz, nw = (size_t)g_.ny * (g_.nz + 1);
		hup_.assign(nu, 0.0); hvp_.assign(nv, 0.0); hwp_.assign(nw, 0.0);
		cudaMalloc(&up_, sizeof(double) * nu);
		cudaMalloc(&vp_, sizeof(double) * nv);
		cudaMalloc(&wp_, sizeof(double) * nw);
	}
	PrecursorReplay::~PrecursorReplay()
	{
		if (up_) cudaFree(up_);
		if (vp_) cudaFree(vp_);
		if (wp_) cudaFree(wp_);
	}

	void PrecursorReplay::advance(double dt)
	{
		clock_ += dt;
		lib_->sample(clock_, hup_, hvp_, hwp_);
		cudaMemcpy(up_, hup_.data(), sizeof(double) * hup_.size(), cudaMemcpyHostToDevice);
		cudaMemcpy(vp_, hvp_.data(), sizeof(double) * hvp_.size(), cudaMemcpyHostToDevice);
		cudaMemcpy(wp_, hwp_.data(), sizeof(double) * hwp_.size(), cudaMemcpyHostToDevice);
	}
	void PrecursorReplay::add_u(double* u, MacGrid g) { plane_add_u_gpu(u, up_, g); }
	void PrecursorReplay::set_vw(double* v, double* w, MacGrid g) { plane_set_vw_gpu(v, w, vp_, wp_, g); }
}
