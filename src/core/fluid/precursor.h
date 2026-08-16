// precursor.h — inlet-plane recording/replay for validation & hero runs (RESEARCH §8,
// research/14 §4). A precursor periodic channel (or any source) records inlet-plane
// fluctuations at 5–10 Hz into a library; PrecursorReplay feeds them back at the obstacle-
// domain inlet with linear-in-time interpolation, looping the ~300 s library. This removes
// the ~6δ SEM adjustment-length bias for the M6/V7 gate (the full turbulent precursor run
// itself is the first task of M6; here we build and unit-test the record/replay machinery).
//
// Planes use the SEM/plane_ops layout: up[ny·nz], vp[(ny+1)·nz], wp[ny·(nz+1)], index
// k·(dimj)+j, matching the i=0 MAC faces.
#pragma once

#include "core/fluid/inlet_fluct.h"
#include "core/fluid/mac_grid.h"

#include <vector>

namespace paracfd::core
{
	struct PrecursorFrame
	{
		double t = 0.0;
		std::vector<double> up, vp, wp;
	};

	// Ordered library of inlet-plane frames sampled at a fixed rate; loops seamlessly with
	// period = nframes·record_dt (frame nframes−1 interpolates back to frame 0).
	class PrecursorLibrary
	{
	public:
		PrecursorLibrary(MacGrid g, double record_hz);

		// Record the source planes if the scheduled time has arrived; returns true if kept.
		bool maybe_record(double t, const std::vector<double>& up, const std::vector<double>& vp, const std::vector<double>& wp);

		int size() const { return (int)frames_.size(); }
		double record_dt() const { return record_dt_; }
		double period() const { return frames_.empty() ? 0.0 : frames_.size() * record_dt_; }
		MacGrid grid() const { return g_; }
		const PrecursorFrame& frame(int i) const { return frames_[i]; }

		// Interpolated planes at time t (looped). No-op fill (zeros) if empty.
		void sample(double t, std::vector<double>& up, std::vector<double>& vp, std::vector<double>& wp) const;

	private:
		MacGrid g_;
		double record_dt_ = 0.1;
		double next_t_ = 0.0;
		std::vector<PrecursorFrame> frames_;
	};

	class PrecursorReplay final : public InletFluct
	{
	public:
		explicit PrecursorReplay(const PrecursorLibrary& lib);
		~PrecursorReplay() override;
		PrecursorReplay(const PrecursorReplay&) = delete;
		PrecursorReplay& operator=(const PrecursorReplay&) = delete;

		void advance(double dt) override;
		void add_u(double* u, MacGrid g) override;
		void set_vw(double* v, double* w, MacGrid g) override;

		double clock() const { return clock_; }
		const std::vector<double>& u_plane() const { return hup_; }

	private:
		const PrecursorLibrary* lib_;
		MacGrid g_;
		double clock_ = 0.0;
		std::vector<double> hup_, hvp_, hwp_;
		double *up_ = nullptr, *vp_ = nullptr, *wp_ = nullptr;
	};
}
