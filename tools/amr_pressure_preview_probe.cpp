// Regression for the explicitly qualitative pressure fallback.
//
// A deliberately rank-deficient skew aperture cannot represent its deferred
// non-orthogonal correction. The certified/default build must reject it; the
// labelled preview policy may retain only the conservative orthogonal A/d flux.

#include "core/fluid/amr_eb.h"
#include "core/fluid/amr_grid.h"
#include "core/fluid/amr_pressure.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

using namespace paracfd::core;

namespace
{
	int failures = 0;

	void check(bool condition, const char* description)
	{
		std::printf("[pressure-preview] %-66s %s\n", description,
			condition ? "PASS" : "FAIL");
		if (!condition) ++failures;
	}

	bool near(double a, double b, double tolerance = 1e-12)
	{
		return std::abs(a - b) <= tolerance;
	}

	AmrEmbeddedBoundaryAtlas unsupported_skew_aperture_fixture()
	{
		AmrEmbeddedBoundaryAtlas atlas;
		AmrEbLevelAtlas level;
		level.level = 0;
		EmbeddedBoundary& eb = level.topology;
		eb.grid = {{0, 0, 0}, 8, 4, 4, 1};
		eb.cells.resize(eb.grid.cell_count());
		eb.cut_face_mask.assign(eb.grid.cell_count(), 0);
		level.owned_cell.assign(eb.grid.cell_count(), 1);

		const int lower_cell = eb.grid.cell_index(0, 0, 0);
		const int upper_cell = eb.grid.cell_index(1, 0, 0);
		for (int fragment = 0; fragment < 2; ++fragment)
		{
			const int cell = fragment == 0 ? lower_cell : upper_cell;
			EbCellTopology& topology = eb.cells[cell];
			topology.state = EbCellState::split;
			topology.first_fragment = fragment;
			topology.fragment_count = 1;

			FluidFragment value;
			value.parent_cell = cell;
			value.volume = 1.0;
			value.centroid = fragment == 0 ? Vec3d{0.4, 0.25, 0.5} :
				Vec3d{1.6, 0.75, 0.5};
			value.merge_target = irregular_fragment(fragment);
			eb.fragments.push_back(value);
		}

		// The only same-fluid path is this skew X aperture, leaving both endpoint
		// WLS stencils rank one. Its face is nevertheless correctly bracketed.
		eb.apertures.push_back({lower_cell, 0, 0.5, {1.0, 0.5, 0.5},
			irregular_fragment(0), irregular_fragment(1)});
		atlas.levels.push_back(std::move(level));
		return atlas;
	}

	AmrEmbeddedBoundaryAtlas supported_skew_agglomeration_fixture()
	{
		AmrEmbeddedBoundaryAtlas atlas;
		AmrEbLevelAtlas level;
		level.level = 0;
		EmbeddedBoundary& eb = level.topology;
		eb.grid = {{0, 0, 0}, 8, 4, 4, 1};
		eb.cells.resize(eb.grid.cell_count());
		eb.cut_face_mask.assign(eb.grid.cell_count(), 0);
		level.owned_cell.assign(eb.grid.cell_count(), 1);

		const int merged_cell = eb.grid.cell_index(0, 0, 0);
		const int split_cell = eb.grid.cell_index(1, 0, 0);
		const int downstream_cell = eb.grid.cell_index(2, 0, 0);
		EbCellTopology& topology = eb.cells[split_cell];
		topology.state = EbCellState::split;
		topology.first_fragment = 0;
		topology.fragment_count = 1;

		FluidFragment fragment;
		fragment.parent_cell = split_cell;
		fragment.volume = 0.25;
		fragment.centroid = {1.25, 0.8, 0.7};
		fragment.merge_target = regular_fragment(merged_cell);
		eb.fragments.push_back(fragment);
		eb.apertures.push_back({split_cell, 0, 0.4, {2.0, 0.5, 0.5},
			irregular_fragment(0), regular_fragment(downstream_cell)});
		atlas.levels.push_back(std::move(level));
		return atlas;
	}
}

int main()
{
	const AmrHierarchy hierarchy = AmrHierarchy::uniform({{0, 0, 0}, {8, 4, 4}},
		1.0, 4, 1);
	const AmrEmbeddedBoundaryAtlas atlas = unsupported_skew_aperture_fixture();

	bool strict_rejected = false;
	std::string strict_message;
	try
	{
		(void)build_composite_amr_pressure_system(hierarchy, atlas, false);
	}
	catch (const std::runtime_error& error)
	{
		strict_message = error.what();
		strict_rejected = strict_message.find(
			"unsupported nonorthogonal pressure correction") != std::string::npos;
	}
	check(strict_rejected, "default/certified policy rejects unsupported WLS correction");

	CompositeAmrPressureBuildOptions preview_options;
	preview_options.pressure_outlet_xmax = false;
	preview_options.unsupported_nonorthogonal_correction =
		UnsupportedNonorthogonalCorrectionPolicy::qualitative_preview_orthogonal_fallback;
	const CompositeAmrPressureSystem preview =
		build_composite_amr_pressure_system(hierarchy, atlas, preview_options);

	check(preview.qualitative_preview_orthogonal_fallback_count() == 1 &&
		preview.qualitative_preview_orthogonal_fallback_embedded == 1 &&
		preview.qualitative_preview_orthogonal_fallback_regular == 0 &&
		preview.qualitative_preview_orthogonal_fallback_coarse_fine == 0,
		"preview records exactly the one dropped embedded correction");
	const bool has_one_preview_edge = preview.embedded.size() == 1;
	check(has_one_preview_edge &&
		length2(preview.embedded.front().nonorthogonal_correction) == 0.0 &&
		preview.embedded.front().lower_gradient_node == -1 &&
		preview.embedded.front().upper_gradient_node == -1,
		"preview drops only the unsupported deferred correction");
	if (!has_one_preview_edge) return 1;

	const CoarseFinePressureConnection& edge = preview.embedded.front();
	const int lower = edge.direction > 0 ? edge.coarse_dof : edge.fine_dof;
	const int upper = edge.direction > 0 ? edge.fine_dof : edge.coarse_dof;
	const double expected_conductance = edge.open_area / 1.2;
	check(near(pressure_gradient_factor(edge), 1.0 / 1.2) &&
		near(edge.open_area * pressure_gradient_factor(edge), expected_conductance),
		"preview retains the orthogonal A/d conductance");

	std::vector<double> pressure(preview.storage_size, 0.0), applied;
	pressure[lower] = 1.0;
	preview.apply_cpu(pressure, applied);
	double sum = 0.0;
	for (double value : applied) sum += value;
	check(near(applied[lower], expected_conductance) &&
		near(applied[upper], -expected_conductance) && near(sum, 0.0),
		"fallback two-point flux remains pairwise conservative");

	// Unlike the unsupported-only fallback above, the first-order policy drops
	// every deferred term, including valid rank-3 WLS corrections around an
	// agglomerated regular control volume.
	const AmrEmbeddedBoundaryAtlas supported_atlas =
		supported_skew_agglomeration_fixture();
	const CompositeAmrPressureSystem supported_strict =
		build_composite_amr_pressure_system(hierarchy, supported_atlas, false);
	check(!supported_strict.pressure_gradient_dof.empty() &&
		!supported_strict.regular_pressure_corrections.empty() &&
		length2(supported_strict.embedded.front().nonorthogonal_correction) > 0.0,
		"manufactured agglomeration has supported regular and embedded WLS terms");

	CompositeAmrPressureBuildOptions first_order_options;
	first_order_options.pressure_outlet_xmax = false;
	first_order_options.unsupported_nonorthogonal_correction =
		UnsupportedNonorthogonalCorrectionPolicy::qualitative_preview_first_order_orthogonal;
	const CompositeAmrPressureSystem first_order =
		build_composite_amr_pressure_system(hierarchy, supported_atlas, first_order_options);
	bool all_special_orthogonal = true;
	for (const CoarseFinePressureConnection& connection : first_order.coarse_fine)
		all_special_orthogonal = all_special_orthogonal &&
			length2(connection.nonorthogonal_correction) == 0.0;
	for (const CoarseFinePressureConnection& connection : first_order.embedded)
		all_special_orthogonal = all_special_orthogonal &&
			length2(connection.nonorthogonal_correction) == 0.0;
	bool all_regular_orthogonal = true;
	for (const RegularPressureCorrection& correction : first_order.regular_pressure_corrections)
		all_regular_orthogonal = all_regular_orthogonal &&
			length2(correction.nonorthogonal_correction) == 0.0;
	check(first_order.qualitative_preview_orthogonal_fallback_regular > 0 &&
		first_order.qualitative_preview_orthogonal_fallback_embedded > 0 &&
		all_special_orthogonal && all_regular_orthogonal,
		"first-order preview counts and drops every represented skew term");
	check(first_order.pressure_gradient_dof.empty() &&
		first_order.pressure_gradient_offset.empty() &&
		first_order.pressure_gradient_neighbour.empty() &&
		first_order.pressure_gradient_weight.empty(),
		"first-order preview allocates no unnecessary WLS topology");
	check(near(pressure_gradient_factor(first_order.embedded.front()),
		pressure_gradient_factor(supported_strict.embedded.front())) &&
		!first_order.regular_pressure_corrections.empty(),
		"first-order preview retains special A/d and regular two-point deltas");

	std::vector<double> first_pressure(first_order.storage_size, 0.0),
		first_applied, constant(first_order.storage_size, 0.0), constant_applied;
	for (int dof = 0; dof < first_order.storage_size; ++dof)
		if (first_order.active[dof])
		{
			first_pressure[dof] = std::sin(0.17 * dof) + 0.01 * dof;
			constant[dof] = 1.0;
		}
	first_order.apply_cpu(first_pressure, first_applied);
	first_order.apply_cpu(constant, constant_applied);
	double conservative_sum = 0.0, constant_maximum = 0.0;
	for (int dof = 0; dof < first_order.storage_size; ++dof)
	{
		conservative_sum += first_applied[dof];
		constant_maximum = std::max(constant_maximum, std::abs(constant_applied[dof]));
	}
	check(std::abs(conservative_sum) < 1e-10 && constant_maximum < 1e-12,
		"first-order orthogonal operator is conservative with constant nullspace");

	first_order_options.pressure_outlet_xmax = true;
	const CompositeAmrPressureSystem anchored =
		build_composite_amr_pressure_system(hierarchy, supported_atlas, first_order_options);
	std::vector<double> p(anchored.storage_size, 0.0), q(anchored.storage_size, 0.0), Ap, Aq;
	for (int dof = 0; dof < anchored.storage_size; ++dof)
		if (anchored.active[dof])
		{
			p[dof] = std::sin(0.11 * dof) + 0.003 * dof;
			q[dof] = std::cos(0.07 * dof) - 0.002 * dof;
		}
	anchored.apply_cpu(p, Ap);
	anchored.apply_cpu(q, Aq);
	double pAp = 0.0, pAq = 0.0, qAp = 0.0;
	for (int dof = 0; dof < anchored.storage_size; ++dof)
	{
		pAp += p[dof] * Ap[dof];
		pAq += p[dof] * Aq[dof];
		qAp += q[dof] * Ap[dof];
	}
	check(pAp > 0.0 && std::abs(pAq - qAp) < 1e-10,
		"outlet-anchored first-order pressure operator remains SPD");

	if (failures)
		std::fprintf(stderr, "[pressure-preview] %d regression(s) failed\n", failures);
	return failures ? 1 : 0;
}
