#include "core/fluid/amr_eb.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace paracfd::core
{
	std::size_t AmrEmbeddedBoundaryAtlas::owned_unresolved_count() const
	{
		std::size_t count=0;for(const AmrEbLevelAtlas& atlas:levels)for(const UnresolvedEbCell& problem:atlas.topology.unresolved)if(problem.parent_cell>=0&&problem.parent_cell<static_cast<int>(atlas.owned_cell.size())&&atlas.owned_cell[problem.parent_cell])++count;return count;
	}

	AmrEmbeddedBoundaryAtlas build_amr_embedded_boundary_atlas(const AmrHierarchy& hierarchy,
		const TriMesh& mesh, const TriangleBvh& bvh, const EmbeddedBoundaryBuildOptions& options)
	{
		AmrEmbeddedBoundaryAtlas result;const int bs=hierarchy.brick_size();
		for(int level=0;level<static_cast<int>(hierarchy.levels().size());++level)
		{
			const AmrLevel& source=hierarchy.levels()[level];Int3 minimum{std::numeric_limits<int>::max(),std::numeric_limits<int>::max(),std::numeric_limits<int>::max()},maximum{std::numeric_limits<int>::min(),std::numeric_limits<int>::min(),std::numeric_limits<int>::min()};bool found=false;
			for(const BrickMetadata& brick:source.bricks)if(brick.active())
			{
				const double width=bs*static_cast<double>(brick.h);const Aabb3d brick_box{brick.origin,brick.origin+Vec3d{width,width,width}};if(!brick.embedded_boundary()&&bvh.query_aabb(brick_box).empty())continue;found=true;minimum.x=std::min(minimum.x,brick.coord.x);minimum.y=std::min(minimum.y,brick.coord.y);minimum.z=std::min(minimum.z,brick.coord.z);maximum.x=std::max(maximum.x,brick.coord.x);maximum.y=std::max(maximum.y,brick.coord.y);maximum.z=std::max(maximum.z,brick.coord.z);
			}
			if(!found)continue;const double h=source.h;const Vec3d extent=hierarchy.domain().hi-hierarchy.domain().lo;const Int3 domain_cells{static_cast<int>(std::llround(extent.x/h)),static_cast<int>(std::llround(extent.y/h)),static_cast<int>(std::llround(extent.z/h))};const Int3 cell_lo{std::max(0,minimum.x*bs-1),std::max(0,minimum.y*bs-1),std::max(0,minimum.z*bs-1)},cell_hi{std::min(domain_cells.x,(maximum.x+1)*bs+1),std::min(domain_cells.y,(maximum.y+1)*bs+1),std::min(domain_cells.z,(maximum.z+1)*bs+1)};const Vec3d origin=hierarchy.domain().lo+Vec3d{cell_lo.x*h,cell_lo.y*h,cell_lo.z*h};const UniformEbGrid grid{origin,cell_hi.x-cell_lo.x,cell_hi.y-cell_lo.y,cell_hi.z-cell_lo.z,h};
			AmrEbLevelAtlas atlas;atlas.level=level;atlas.minimum_brick_coord=minimum;atlas.maximum_brick_coord=maximum;atlas.topology=build_embedded_boundary(mesh,bvh,grid,options);atlas.owned_cell.assign(grid.cell_count(),0);
			for(int cell=0;cell<grid.cell_count();++cell){const BrickLocation owner=hierarchy.locate_finest(grid.cell_centroid(cell));atlas.owned_cell[cell]=static_cast<unsigned char>(owner.found()&&owner.level==level);}
			result.levels.push_back(std::move(atlas));
		}
		return result;
	}

	OneLevelEmbeddedBoundary build_one_level_embedded_boundary(const AmrHierarchy& hierarchy,
		const TriMesh& mesh, const TriangleBvh& bvh, const EmbeddedBoundaryBuildOptions& options)
	{
		if (hierarchy.levels().size() != 1) throw std::invalid_argument("one-level EB bridge requires exactly one AMR level");
		const AmrLevel& level = hierarchy.levels().front();
		if (!(level.h > 0) || level.bricks.empty()) throw std::invalid_argument("one-level EB bridge requires a non-empty level");
		const Aabb3d& domain = hierarchy.domain(); const Vec3d extent = domain.hi - domain.lo;
		auto cell_count = [&](double length)
		{
			const double cells = length / level.h; const int rounded = static_cast<int>(std::llround(cells));
			if (rounded <= 0 || std::abs(cells - rounded) > 1e-5) throw std::invalid_argument("AMR domain is not cell-aligned");
			return rounded;
		};
		UniformEbGrid grid{domain.lo, cell_count(extent.x), cell_count(extent.y), cell_count(extent.z), level.h};
		OneLevelEmbeddedBoundary out; out.topology = build_embedded_boundary(mesh, bvh, grid, options);
		out.field_layout = BrickFieldLayout::make(hierarchy.brick_size(), hierarchy.ghost_cells());
		out.cell_field_index.resize(grid.cell_count(), std::numeric_limits<std::size_t>::max());
		const int bs = hierarchy.brick_size();
		for (int k = 0; k < grid.nz; ++k) for (int j = 0; j < grid.ny; ++j) for (int i = 0; i < grid.nx; ++i)
		{
			const Int3 brick_coord{i / bs, j / bs, k / bs}; const int brick = hierarchy.find_brick(0, brick_coord);
			if (brick < 0 || !level.bricks[brick].active()) throw std::runtime_error("one-level AMR tiling has a missing/covered brick");
			out.cell_field_index[grid.cell_index(i, j, k)] = out.field_layout.cell_index(brick, i % bs, j % bs, k % bs);
		}
		return out;
	}

	void gather_one_level_pressure(const OneLevelEmbeddedBoundary& composite, const EbPressureSystem& system,
		const AmrHostLevelFields& level, std::vector<Real>& pressure)
	{
		if (system.eb != &composite.topology) throw std::invalid_argument("pressure system/topology mismatch");
		if (pressure.size() != static_cast<std::size_t>(system.storage_size)) pressure.assign(system.storage_size, Real{});
		for (int cell = 0; cell < composite.topology.grid.cell_count(); ++cell)
		{
			const int dof = system.cell_dof[cell]; if (dof < 0) continue; const std::size_t field = composite.cell_field_index[cell];
			if (field >= level.p.size()) throw std::out_of_range("one-level pressure field mapping"); pressure[dof] = level.p[field];
		}
	}

	void scatter_one_level_pressure(const OneLevelEmbeddedBoundary& composite, const EbPressureSystem& system,
		const std::vector<Real>& pressure, AmrHostLevelFields& level)
	{
		if (system.eb != &composite.topology || pressure.size() != static_cast<std::size_t>(system.storage_size)) throw std::invalid_argument("pressure system/vector mismatch");
		for (int cell = 0; cell < composite.topology.grid.cell_count(); ++cell)
		{
			const int dof = system.cell_dof[cell]; if (dof < 0) continue; const std::size_t field = composite.cell_field_index[cell];
			if (field >= level.p.size()) throw std::out_of_range("one-level pressure field mapping"); level.p[field] = pressure[dof];
		}
	}
}
