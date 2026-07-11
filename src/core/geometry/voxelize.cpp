// voxelize.cpp — see voxelize.h. Pure host C++ (no CUDA, no OCC): a one-time preprocessing
// step that runs on the CPU (the mask is then uploaded to the core like the cylinder mask).
#include "core/geometry/voxelize.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace windcfd::core
{
	namespace
	{
		// A triangle's three placed vertices (metres, double).
		struct Tri
		{
			double ax, ay, az, bx, by, bz, cx, cy, cz;
		};
	}

	double mesh_signed_volume_abs(const TriMesh& mesh)
	{
		// V = (1/6) Σ v0 · (v1 × v2) over triangles (signed). |V| is placement-invariant.
		double v6 = 0.0;
		for (std::size_t t = 0; t + 2 < mesh.indices.size(); t += 3)
		{
			const float* p0 = &mesh.positions[3 * mesh.indices[t + 0]];
			const float* p1 = &mesh.positions[3 * mesh.indices[t + 1]];
			const float* p2 = &mesh.positions[3 * mesh.indices[t + 2]];
			const double cx = (double)p1[1] * p2[2] - (double)p1[2] * p2[1];
			const double cy = (double)p1[2] * p2[0] - (double)p1[0] * p2[2];
			const double cz = (double)p1[0] * p2[1] - (double)p1[1] * p2[0];
			v6 += (double)p0[0] * cx + (double)p0[1] * cy + (double)p0[2] * cz;
		}
		return std::fabs(v6) / 6.0;
	}

	std::vector<unsigned char> voxelize_mesh(const TriMesh& mesh, MacGrid g, const ModelPlacement& place,
		double* out_mesh_volume, double* out_voxel_volume,
		std::vector<float>* out_fraction, int* out_thin_cells, int supersample)
	{
		const int ncell = g.p_count();
		std::vector<unsigned char> mask((std::size_t)ncell, 0);
		if (out_fraction) out_fraction->assign((std::size_t)ncell, 0.0f);
		if (out_thin_cells) *out_thin_cells = 0;

		// Place-correct the mesh volume by |det(M)| so the mesh-vs-voxel volume gate stays
		// meaningful when the placement scales the model (rotation/translation ⇒ |det| = 1,
		// so the byte-identical translation-only path is unchanged).
		const double mesh_vol = mesh_signed_volume_abs(mesh) * std::fabs(place.linear_det());
		if (out_mesh_volume) *out_mesh_volume = mesh_vol;
		if (mesh.empty() || g.h <= 0.0 || ncell <= 0)
		{
			if (out_voxel_volume) *out_voxel_volume = 0.0;
			return mask;
		}

		const int K = std::max(1, supersample);
		const double h = g.h;
		const double hs = h / (double)K; // sub-cell edge

		// --- Placed triangles + placed xy/z bbox ---------------------------------
		const std::size_t ntri = mesh.triangle_count();
		std::vector<Tri> tris;
		tris.reserve(ntri);
		double lo[3] = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
		double hi[3] = { -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max() };
		for (std::size_t t = 0; t + 2 < mesh.indices.size(); t += 3)
		{
			const float* p0 = &mesh.positions[3 * mesh.indices[t + 0]];
			const float* p1 = &mesh.positions[3 * mesh.indices[t + 1]];
			const float* p2 = &mesh.positions[3 * mesh.indices[t + 2]];
			Tri tr;
			place.apply(p0[0], p0[1], p0[2], tr.ax, tr.ay, tr.az);
			place.apply(p1[0], p1[1], p1[2], tr.bx, tr.by, tr.bz);
			place.apply(p2[0], p2[1], p2[2], tr.cx, tr.cy, tr.cz);
			tris.push_back(tr);
			for (const double* v : { &tr.ax, &tr.bx, &tr.cx })
			{
				lo[0] = std::min(lo[0], v[0]); hi[0] = std::max(hi[0], v[0]);
				lo[1] = std::min(lo[1], v[1]); hi[1] = std::max(hi[1], v[1]);
				lo[2] = std::min(lo[2], v[2]); hi[2] = std::max(hi[2], v[2]);
			}
		}

		// Cell index range overlapping the placed bbox (only these cells can be solid).
		auto clampi = [](int v, int a, int b) { return v < a ? a : (v > b ? b : v); };
		// Graded grid: the model must live in the UNIFORM fine core (h_fine cubes; windloads asserts
			// bbox ⊆ core). There xf(i)=i·h+Δ for a per-axis CONSTANT Δ, so this uniform voxelizer is exact
			// once the placed geometry is shifted by −Δ: then floor(x/h) equals the correct graded CORE cell
			// index. Δ is read from the cell containing the model bbox centre (a core cell). Uniform grid ⇒
			// xfa null ⇒ Δ=0 (no shift, byte-identical). A model straddling the graded transition would
			// mis-voxelize, but the load integration rejects that case (task 4.3).
			if (g.xfa)
			{
				const int ic = clampi((int)std::floor(grid_fx(g, 0.5 * (lo[0] + hi[0]))), 0, g.nx - 1);
				const int jc = clampi((int)std::floor(grid_fy(g, 0.5 * (lo[1] + hi[1]))), 0, g.ny - 1);
				const int kc = clampi((int)std::floor(grid_fz(g, 0.5 * (lo[2] + hi[2]))), 0, g.nz - 1);
				const double Dx = g.xf(ic) - ic * h, Dy = g.yf(jc) - jc * h, Dz = g.zf(kc) - kc * h;
				for (Tri& tr : tris)
				{
					tr.ax -= Dx; tr.bx -= Dx; tr.cx -= Dx;
					tr.ay -= Dy; tr.by -= Dy; tr.cy -= Dy;
					tr.az -= Dz; tr.bz -= Dz; tr.cz -= Dz;
				}
				lo[0] -= Dx; hi[0] -= Dx; lo[1] -= Dy; hi[1] -= Dy; lo[2] -= Dz; hi[2] -= Dz;
			}

			const int i0 = clampi((int)std::floor(lo[0] / h) - 1, 0, g.nx - 1);
		const int i1 = clampi((int)std::floor(hi[0] / h) + 1, 0, g.nx - 1);
		const int j0 = clampi((int)std::floor(lo[1] / h) - 1, 0, g.ny - 1);
		const int j1 = clampi((int)std::floor(hi[1] / h) + 1, 0, g.ny - 1);
		const int k0 = clampi((int)std::floor(lo[2] / h) - 1, 0, g.nz - 1);
		const int k1 = clampi((int)std::floor(hi[2] / h) + 1, 0, g.nz - 1);
		if (i1 < i0 || j1 < j0 || k1 < k0) // model entirely outside the domain
		{
			if (out_voxel_volume) *out_voxel_volume = 0.0;
			std::fprintf(stderr, "[voxelize] model bbox is outside the domain — empty mask.\n");
			return mask;
		}

		// --- Sub-column crossing lists over the bbox footprint -------------------
		// Sub-column grid: sub_i in [i0*K, (i1+1)*K), sub_j in [j0*K, (j1+1)*K). Stored dense
		// (a = sub_i - i0*K, b = sub_j - j0*K) so the memory is bounded by the model footprint.
		const int SA = (i1 - i0 + 1) * K; // sub-columns along x
		const int SB = (j1 - j0 + 1) * K; // sub-columns along y
		const int subi0 = i0 * K, subj0 = j0 * K;
		std::vector<std::vector<float>> cross((std::size_t)SA * SB);

		// Triangle-centric scatter: append each triangle's +z crossing z into the sub-columns
		// inside its xy projection (bbox-culled), tested with 2D barycentrics.
		//
		// BRUTE FORCE by design — NO BVH / k-d tree. The input is a watertight, LOW-POLY solid
		// (V000 ≈ 70 tris); column-parity is O(sub-columns × tris) and runs ONCE at load (the
		// obstacle is rigid — never re-voxelized per step), so a spatial tree would be pure
		// overhead. The only acceleration is the per-triangle xy-AABB cull below (each triangle
		// only touches the sub-columns it overlaps). If a future case is genuinely high-poly,
		// bucket triangles by xy-column here (or add a BVH) — not before.
		//
		// TIE-BREAK: axis-aligned input (e.g. a cube's flat face split into two coplanar
		// triangles along the x=y diagonal) puts sub-column sample points exactly on a shared
		// edge, where an inclusive test double-counts the crossing (parity flips the WRONG way)
		// and an exclusive test misses it. Both corrupt the even-odd count. We nudge the sample
		// point by a tiny ASYMMETRIC offset (≪ sub-cell) so it never lands on an axis-aligned
		// edge/vertex/diagonal — each shared edge is then owned by exactly one triangle ⇒ exactly
		// one crossing. The offset is far below feature size, so the inside-fraction is unchanged.
		const double eps = 1e-12;
		const double ox = 0.0013 * hs, oy = 0.0029 * hs; // asymmetric sub-cell nudge (ox ≠ oy)
		for (const Tri& tr : tris)
		{
			const double txlo = std::min({ tr.ax, tr.bx, tr.cx });
			const double txhi = std::max({ tr.ax, tr.bx, tr.cx });
			const double tylo = std::min({ tr.ay, tr.by, tr.cy });
			const double tyhi = std::max({ tr.ay, tr.by, tr.cy });
			int sa_lo = (int)std::floor(txlo / hs) - 1, sa_hi = (int)std::floor(txhi / hs) + 1;
			int sb_lo = (int)std::floor(tylo / hs) - 1, sb_hi = (int)std::floor(tyhi / hs) + 1;
			sa_lo = std::max(sa_lo, subi0); sa_hi = std::min(sa_hi, subi0 + SA - 1);
			sb_lo = std::max(sb_lo, subj0); sb_hi = std::min(sb_hi, subj0 + SB - 1);
			// 2D edge functions of the triangle projected to xy.
			const double denom = (tr.by - tr.cy) * (tr.ax - tr.cx) + (tr.cx - tr.bx) * (tr.ay - tr.cy);
			if (std::fabs(denom) < eps) continue; // triangle edge-on to +z: contributes no crossing
			const double invden = 1.0 / denom;
			for (int sa = sa_lo; sa <= sa_hi; ++sa)
			{
				const double px = (sa + 0.5) * hs + ox;
				for (int sb = sb_lo; sb <= sb_hi; ++sb)
				{
					const double py = (sb + 0.5) * hs + oy;
					const double l1 = ((tr.by - tr.cy) * (px - tr.cx) + (tr.cx - tr.bx) * (py - tr.cy)) * invden;
					const double l2 = ((tr.cy - tr.ay) * (px - tr.cx) + (tr.ax - tr.cx) * (py - tr.cy)) * invden;
					const double l3 = 1.0 - l1 - l2;
					if (l1 < 0.0 || l2 < 0.0 || l3 < 0.0) continue; // strictly inside exactly one triangle
					const double z = l1 * tr.az + l2 * tr.bz + l3 * tr.cz;
					cross[(std::size_t)(sa - subi0) * SB + (sb - subj0)].push_back((float)z);
				}
			}
		}

		// --- Supersampled majority fill -----------------------------------------
		// Sub-column-major sweep: for each sub-column walk its sorted crossings and toggle parity
		// across sub-z levels, accumulating inside-counts into the owning cells.
		std::vector<int> counts((std::size_t)ncell, 0);
		const int Zlo = k0 * K, Zhi = (k1 + 1) * K; // sub-z range
		for (int a = 0; a < SA; ++a)
		{
			const int cell_i = i0 + a / K;
			for (int b = 0; b < SB; ++b)
			{
				const int cell_j = j0 + b / K;
				std::vector<float>& lst = cross[(std::size_t)a * SB + b];
				if (lst.empty()) continue;
				std::sort(lst.begin(), lst.end());
				std::size_t ic = 0;
				int parity = 0;
				for (int Z = Zlo; Z < Zhi; ++Z)
				{
					const double subz = (Z + 0.5) * hs;
					while (ic < lst.size() && (double)lst[ic] < subz) { parity ^= 1; ++ic; }
					if (parity) ++counts[(std::size_t)g.pidx(cell_i, cell_j, Z / K)];
				}
			}
		}

		const double invK3 = 1.0 / (double)(K * K * K);
		int solid_count = 0;
		for (int k = k0; k <= k1; ++k)
			for (int j = j0; j <= j1; ++j)
				for (int i = i0; i <= i1; ++i)
				{
					const int idx = g.pidx(i, j, k);
					const float frac = (float)(counts[(std::size_t)idx] * invK3);
					if (out_fraction) (*out_fraction)[(std::size_t)idx] = frac;
					if (frac > 0.5f) { mask[(std::size_t)idx] = 1; ++solid_count; }
				}

		// --- Thin-wall safeguard: seal sub-cell sheets the majority rule dropped -
		// A partial cell (0<frac≤0.5) that the majority left fluid AND that has NO solid
		// face-neighbour is an isolated sub-cell sheet with no solid backing — i.e. a wall
		// thinner than ~h that would otherwise leak. Keep it solid (over-thicken beats leak).
		// A thick body's partial skin cells DO have a solid neighbour (the interior), so they
		// are correctly left fluid and not over-thickened. Neighbours are read from the
		// post-majority mask (single deterministic pass).
		int thin_kept = 0;
		auto solid_at = [&](int i, int j, int k) -> bool
		{
			if (i < 0 || i >= g.nx || j < 0 || j >= g.ny || k < 0 || k >= g.nz) return false;
			return mask[(std::size_t)g.pidx(i, j, k)] != 0;
		};
		std::vector<int> promote;
		for (int k = k0; k <= k1; ++k)
			for (int j = j0; j <= j1; ++j)
				for (int i = i0; i <= i1; ++i)
				{
					const int idx = g.pidx(i, j, k);
					if (mask[(std::size_t)idx]) continue;                  // already solid
					if (counts[(std::size_t)idx] <= 0) continue;          // surface never touches this cell
					if (solid_at(i - 1, j, k) || solid_at(i + 1, j, k) ||
						solid_at(i, j - 1, k) || solid_at(i, j + 1, k) ||
						solid_at(i, j, k - 1) || solid_at(i, j, k + 1)) continue; // has solid backing
					promote.push_back(idx);
				}
		for (int idx : promote)
		{
			mask[(std::size_t)idx] = 1;
			++solid_count;
			++thin_kept;
		}
		if (out_thin_cells) *out_thin_cells = thin_kept;

		// --- Volumes + log -------------------------------------------------------
		const double voxel_vol = (double)solid_count * h * h * h;
		if (out_voxel_volume) *out_voxel_volume = voxel_vol;
		const double pct = 100.0 * (double)solid_count / (double)ncell;
		const bool vol_gate = mesh_vol > 1e-12;
		if (vol_gate)
		{
			const double err_pct = 100.0 * std::fabs(mesh_vol - voxel_vol) / mesh_vol;
			std::fprintf(stderr,
				"[voxelize] %zu tris -> %d solid cells (%.3f%% of domain); mesh vol %.6g m^3 vs voxel vol %.6g m^3 (err %.3f%%)\n",
				ntri, solid_count, pct, mesh_vol, voxel_vol, err_pct);
		}
		else
		{
			std::fprintf(stderr,
				"[voxelize] %zu tris -> %d solid cells (%.3f%% of domain); voxel vol %.6g m^3 (volume gate N/A: open shell, |mesh vol|~0)\n",
				ntri, solid_count, pct, voxel_vol);
		}
		if (thin_kept > 0)
			std::fprintf(stderr, "[voxelize] WARNING: %d thin cells kept conservatively (feature < h) to seal the barrier.\n", thin_kept);
		return mask;
	}

	std::vector<unsigned char> voxelize_mesh_instances(const TriMesh& mesh, MacGrid g,
		const std::vector<ModelPlacement>& placements, int* out_solid_count, int supersample)
	{
		std::vector<unsigned char> mask((std::size_t)g.p_count(), 0);
		for (const ModelPlacement& p : placements)
		{
			std::vector<unsigned char> m = voxelize_mesh(mesh, g, p, nullptr, nullptr, nullptr, nullptr, supersample);
			const std::size_t nn = std::min(mask.size(), m.size());
			for (std::size_t i = 0; i < nn; ++i)
				if (m[i]) mask[i] = 1;
		}
		if (out_solid_count)
		{
			int c = 0;
			for (unsigned char v : mask) if (v) ++c;
			*out_solid_count = c;
		}
		return mask;
	}
}
