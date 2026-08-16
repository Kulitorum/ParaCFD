#include "core/fluid/amr_exchange.h"

#include <stdexcept>

namespace paracfd::core
{
	void restrict_fine_pressure_to_coarse(const AmrHierarchy& h, AmrHostFields& f, int fl)
	{
		if (fl <= 0 || fl >= static_cast<int>(h.levels().size())) throw std::invalid_argument("invalid fine AMR level");
		const int cl = fl - 1, bs = h.brick_size(); auto& coarse = f.levels()[cl]; auto& fine = f.levels()[fl];
		for (int pb = 0; pb < static_cast<int>(h.levels()[cl].bricks.size()); ++pb)
		{
			const auto& p = h.levels()[cl].bricks[pb]; if (p.active()) continue;
			for (int k = 0; k < bs; ++k) for (int j = 0; j < bs; ++j) for (int i = 0; i < bs; ++i)
			{
				double sum = 0;
				for (int dz = 0; dz < 2; ++dz) for (int dy = 0; dy < 2; ++dy) for (int dx = 0; dx < 2; ++dx)
				{
					int gx = 2 * i + dx, gy = 2 * j + dy, gz = 2 * k + dz;
					int cx = gx / bs, cy = gy / bs, cz = gz / bs; int child = p.children[(cz * 2 + cy) * 2 + cx];
					if (child < 0) throw std::runtime_error("covered parent missing AMR child");
					sum += fine.p[fine.layout.cell_index(child, gx % bs, gy % bs, gz % bs)];
				}
				coarse.p[coarse.layout.cell_index(pb, i, j, k)] = static_cast<Real>(sum / 8.0);
			}
		}
	}

	void prolong_coarse_pressure_to_fine(const AmrHierarchy& h, AmrHostFields& f, int fl)
	{
		if (fl <= 0 || fl >= static_cast<int>(h.levels().size())) throw std::invalid_argument("invalid fine AMR level");
		const int cl = fl - 1, bs = h.brick_size(); auto& coarse = f.levels()[cl]; auto& fine = f.levels()[fl];
		for (int pb = 0; pb < static_cast<int>(h.levels()[cl].bricks.size()); ++pb)
		{
			const auto& p = h.levels()[cl].bricks[pb]; if (p.active()) continue;
			for (int slot = 0; slot < 8; ++slot)
			{
				int child = p.children[slot]; if (child < 0) continue; int cx = slot % 2, cy = (slot / 2) % 2, cz = slot / 4;
				for (int k = 0; k < bs; ++k) for (int j = 0; j < bs; ++j) for (int i = 0; i < bs; ++i)
				{
					int pi = (cx * bs + i) / 2, pj = (cy * bs + j) / 2, pk = (cz * bs + k) / 2;
					fine.p[fine.layout.cell_index(child, i, j, k)] = coarse.p[coarse.layout.cell_index(pb, pi, pj, pk)];
				}
			}
		}
	}

	FluxMatchResult match_two_to_one_flux(double uc, double Ac, const std::array<double, 4>& uf, const std::array<double, 4>& Af)
	{
		FluxMatchResult r; r.original_coarse_flux = uc * Ac;
		for (int i = 0; i < 4; ++i) r.fine_flux_sum += uf[i] * Af[i];
		r.corrected_coarse_velocity = Ac > 0 ? r.fine_flux_sum / Ac : 0;
		r.reflux_correction = r.fine_flux_sum - r.original_coarse_flux;
		return r;
	}
}
