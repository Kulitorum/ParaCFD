// flow_particles.cpp — see flow_particles.h. CPU advection of the arrow tracers.
#include "gui/flow_particles.h"

#include <algorithm>
#include <cmath>

namespace scour::gui
{
	using scour::core::MacGrid;

	namespace
	{
		inline float frand(std::mt19937& r)
		{
			return (float)((r() >> 8) & 0xFFFFFF) / (float)0x1000000; // [0,1)
		}
		inline int clampi(int x, int lo, int hi) { return x < lo ? lo : (x > hi ? hi : x); }
	}

	void FlowParticles::set_count(int n)
	{
		n = std::max(0, std::min(n, 20000));
		count_ = n;
		px_.resize(n); py_.resize(n); pz_.resize(n);
		age_.resize(n); life_.resize(n);
		dirx_.resize(n); diry_.resize(n); dirz_.resize(n);
		spd_.resize(n); alpha_.resize(n);
		inst_.resize((size_t)n * 8);
		needs_reset_ = true; // fresh particles spawn on the next advance() (needs the field)
	}

	bool FlowParticles::is_solid(const FlowField& f, float x, float y, float z) const
	{
		const MacGrid& g = f.grid;
		if (!f.solid) return false;
		int i = clampi((int)std::floor(x / (float)g.h), 0, g.nx - 1);
		int j = clampi((int)std::floor(y / (float)g.h), 0, g.ny - 1);
		int k = clampi((int)std::floor(z / (float)g.h), 0, g.nz - 1);
		return f.solid[(size_t)g.pidx(i, j, k)] != 0;
	}

	// Cell-centred velocity component helpers (average of the two bracketing MAC faces), with
	// solid cells treated as zero flow so tracers slow + respawn into obstacles rather than
	// tunnelling. Reads clamp to the valid index range.
	void FlowParticles::sample(const FlowField& f, float x, float y, float z, double& uu, double& vv, double& ww) const
	{
		const MacGrid& g = f.grid;
		const float h = (float)g.h;
		// Continuous cell-centre index space: centre of cell (i,j,k) is at ((i+0.5)h,...).
		float gx = x / h - 0.5f, gy = y / h - 0.5f, gz = z / h - 0.5f;
		int i0 = clampi((int)std::floor(gx), 0, g.nx - 1);
		int j0 = clampi((int)std::floor(gy), 0, g.ny - 1);
		int k0 = clampi((int)std::floor(gz), 0, g.nz - 1);
		int i1 = std::min(i0 + 1, g.nx - 1);
		int j1 = std::min(j0 + 1, g.ny - 1);
		int k1 = std::min(k0 + 1, g.nz - 1);
		float fx = std::clamp(gx - (float)i0, 0.0f, 1.0f);
		float fy = std::clamp(gy - (float)j0, 0.0f, 1.0f);
		float fz = std::clamp(gz - (float)k0, 0.0f, 1.0f);

		auto cc = [&](int i, int j, int k, double& cu, double& cv, double& cw)
		{
			if (f.solid && f.solid[(size_t)g.pidx(i, j, k)]) { cu = cv = cw = 0.0; return; }
			cu = 0.5 * (f.u[g.uidx(i, j, k)] + f.u[g.uidx(i + 1, j, k)]);
			cv = 0.5 * (f.v[g.vidx(i, j, k)] + f.v[g.vidx(i, j + 1, k)]);
			cw = 0.5 * (f.w[g.widx(i, j, k)] + f.w[g.widx(i, j, k + 1)]);
		};

		double U[2][2][2], V[2][2][2], W[2][2][2];
		const int ii[2] = { i0, i1 }, jj[2] = { j0, j1 }, kk[2] = { k0, k1 };
		for (int a = 0; a < 2; ++a)
			for (int b = 0; b < 2; ++b)
				for (int c = 0; c < 2; ++c)
					cc(ii[a], jj[b], kk[c], U[a][b][c], V[a][b][c], W[a][b][c]);

		auto trilerp = [&](double q[2][2][2]) -> double
		{
			double c00 = q[0][0][0] * (1 - fx) + q[1][0][0] * fx;
			double c10 = q[0][1][0] * (1 - fx) + q[1][1][0] * fx;
			double c01 = q[0][0][1] * (1 - fx) + q[1][0][1] * fx;
			double c11 = q[0][1][1] * (1 - fx) + q[1][1][1] * fx;
			double c0 = c00 * (1 - fy) + c10 * fy;
			double c1 = c01 * (1 - fy) + c11 * fy;
			return c0 * (1 - fz) + c1 * fz;
		};
		uu = trilerp(U); vv = trilerp(V); ww = trilerp(W);
	}

	void FlowParticles::spawn(int k, const ArrowView& view, const FlowField& f, bool initial)
	{
		const MacGrid& g = f.grid;
		float Lx = (float)(g.nx * g.h), Ly = (float)(g.ny * g.h), Lz = (float)(g.nz * g.h);
		float x = 0, y = 0, z = 0;
		for (int t = 0; t < 8; ++t) // reject solid spawn points (up to 8 tries)
		{
			if (view.three_d)
			{
				x = frand(rng_) * Lx; y = frand(rng_) * Ly; z = frand(rng_) * Lz;
			}
			else if (view.axis == 0) // X-normal: in-plane (y,z)
			{
				x = view.plane_pos; y = frand(rng_) * Ly; z = frand(rng_) * Lz;
			}
			else if (view.axis == 1) // Y-normal: (x,z)
			{
				x = frand(rng_) * Lx; y = view.plane_pos; z = frand(rng_) * Lz;
			}
			else // Z-normal: (x,y)
			{
				x = frand(rng_) * Lx; y = frand(rng_) * Ly; z = view.plane_pos;
			}
			if (!is_solid(f, x, y, z)) break;
		}
		px_[k] = x; py_[k] = y; pz_[k] = z;
		life_[k] = 2.5f + frand(rng_) * 2.5f;              // 2.5..5 s
		age_[k] = initial ? frand(rng_) * life_[k] : 0.0f; // stagger initial ages so it looks alive
		alpha_[k] = 0.0f; spd_[k] = 0.0f;
		dirx_[k] = diry_[k] = dirz_[k] = 0.0f;
	}

	void FlowParticles::advance(float dt, const ArrowView& view, const FlowField& f)
	{
		if (count_ == 0 || !f.u || f.grid.nx <= 0) return;
		const MacGrid& g = f.grid;
		float Lx = (float)(g.nx * g.h), Ly = (float)(g.ny * g.h), Lz = (float)(g.nz * g.h);
		const float eps = 1e-6f;

		if (needs_reset_)
		{
			for (int k = 0; k < count_; ++k) spawn(k, view, f, true);
			needs_reset_ = false;
		}

		dt = std::clamp(dt, 0.0f, 0.05f);
		const float move = dt * view.gain; // metres advanced per (m/s)
		// Sub-step so the fastest tracer never jumps more than ~½ cell (curved paths stay accurate).
		float maxstep = view.speed_scale * move / (0.5f * (float)g.h);
		int nsub = std::clamp((int)std::ceil(maxstep), 1, 8);
		const float hstep = move / (float)nsub;
		const float inv_scale = view.speed_scale > eps ? 1.0f / view.speed_scale : 1.0f;

		for (int k = 0; k < count_; ++k)
		{
			age_[k] += dt;
			if (age_[k] > life_[k]) { spawn(k, view, f, false); continue; }

			float x = px_[k], y = py_[k], z = pz_[k];
			// 2D: keep the tracer pinned to the (possibly moved) slice plane so dragging the plane
			// slider drags the arrows onto it, without a jarring full respawn.
			if (!view.three_d)
			{
				if (view.axis == 0) x = view.plane_pos;
				else if (view.axis == 1) y = view.plane_pos;
				else z = view.plane_pos;
			}
			double vx = 0, vy = 0, vz = 0;
			bool dead = false;
			for (int s = 0; s < nsub; ++s)
			{
				if (x < 0 || x > Lx || y < 0 || y > Ly || z < 0 || z > Lz || is_solid(f, x, y, z)) { dead = true; break; }
				double uu, vv, ww;
				sample(f, x, y, z, uu, vv, ww);
				if (view.three_d) { vx = uu; vy = vv; vz = ww; }
				else if (view.axis == 0) { vx = 0; vy = vv; vz = ww; }
				else if (view.axis == 1) { vx = uu; vy = 0; vz = ww; }
				else { vx = uu; vy = vv; vz = 0; }
				x += (float)vx * hstep; y += (float)vy * hstep; z += (float)vz * hstep;
			}
			if (dead) { spawn(k, view, f, false); continue; }

			double sp = std::sqrt(vx * vx + vy * vy + vz * vz);
			spd_[k] = std::clamp((float)sp * inv_scale, 0.0f, 1.0f);
			if (sp < 1e-5)
			{
				dirx_[k] = diry_[k] = dirz_[k] = 0.0f;
			}
			else
			{
				dirx_[k] = (float)(vx / sp); diry_[k] = (float)(vy / sp); dirz_[k] = (float)(vz / sp);
			}
			px_[k] = x; py_[k] = y; pz_[k] = z;
			float a = age_[k], life = life_[k];
			alpha_[k] = std::clamp(a / 0.3f, 0.0f, 1.0f) * std::clamp((life - a) / 0.5f, 0.0f, 1.0f);
		}

		// Pack the interleaved instance array (stride 8).
		for (int k = 0; k < count_; ++k)
		{
			float* o = &inst_[(size_t)k * 8];
			o[0] = px_[k]; o[1] = py_[k]; o[2] = pz_[k];
			o[3] = dirx_[k]; o[4] = diry_[k]; o[5] = dirz_[k];
			o[6] = spd_[k]; o[7] = alpha_[k];
		}
	}
}
