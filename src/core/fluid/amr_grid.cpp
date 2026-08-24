#include "core/fluid/amr_grid.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unordered_set>

namespace paracfd::core
{
	namespace
	{
		std::uint64_t key_of(Int3 c)
		{
			return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(c.x)) << 42) ^
				(static_cast<std::uint64_t>(static_cast<std::uint32_t>(c.y)) << 21) ^ static_cast<std::uint32_t>(c.z);
		}
		Aabb3d expanded(Aabb3d b, double d) { b.lo = b.lo - Vec3d{d, d, d}; b.hi = b.hi + Vec3d{d, d, d}; return b; }
		Aabb3d brick_box(const BrickMetadata& b, int bs)
		{
			const double w = bs * static_cast<double>(b.h); return {b.origin, b.origin + Vec3d{w, w, w}};
		}
		bool box_overlap(const Aabb3d& a, const Aabb3d& b) { return a.overlaps(b); }
		template<class Work>
		void parallel_for_ranges(int count,Work&& work)
		{
			if(count<=0)return;unsigned workers=std::max(1u,std::thread::hardware_concurrency());
			if(const char* configured=std::getenv("PARACFD_GRID_THREADS"))
			{
				char* end=nullptr;const long parsed=std::strtol(configured,&end,10);
				if(end!=configured&&parsed>0)workers=static_cast<unsigned>(parsed);
			}
			workers=std::min(workers,static_cast<unsigned>(count));
			if(workers<=1){work(0,count);return;}
			std::atomic<bool> failed{false};std::exception_ptr error;std::mutex error_mutex;
			auto run=[&](unsigned worker)
			{
				try
				{
					const int begin=static_cast<int>((static_cast<long long>(count)*worker)/workers);
					const int end=static_cast<int>((static_cast<long long>(count)*(worker+1))/workers);
					if(!failed.load(std::memory_order_relaxed))work(begin,end);
				}
				catch(...)
				{
					failed.store(true,std::memory_order_relaxed);const std::lock_guard<std::mutex> lock(error_mutex);
					if(!error)error=std::current_exception();
				}
			};
			std::vector<std::thread> threads;threads.reserve(workers-1);
			for(unsigned worker=1;worker<workers;++worker)threads.emplace_back(run,worker);
			run(0);for(auto& thread:threads)thread.join();if(error)std::rethrow_exception(error);
		}
	}

	std::uint64_t AmrHierarchy::coord_hash(Int3 c)
	{
		std::uint64_t x = key_of(c) + 0x9e3779b97f4a7c15ull;
		x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull; x = (x ^ (x >> 27)) * 0x94d049bb133111ebull; return x ^ (x >> 31);
	}

	void AmrHierarchy::initialize_base(const Aabb3d& requested, double h, bool anchor_y_min)
	{
		if (!requested.valid() || !(h > 0.0) || brick_size_ <= 0) throw std::invalid_argument("invalid AMR domain/cell/brick size");
		const double bw = h * brick_size_;
		const Vec3d size = requested.hi - requested.lo;
		const Int3 dims{std::max(1, static_cast<int>(std::ceil(size.x / bw))), std::max(1, static_cast<int>(std::ceil(size.y / bw))), std::max(1, static_cast<int>(std::ceil(size.z / bw)))};
		const Vec3d padded_size{dims.x * bw, dims.y * bw, dims.z * bw};
		const Vec3d excess = padded_size - size;
		// +X is the prescribed freestream direction, so an indivisible base-brick
		// remainder is useful downstream. There is no analogous preferred side in Y
		// or Z: distribute padding symmetrically so the aerodynamic object does not
		// drift into a corner merely because the requested domain is not an exact
		// multiple of a base-brick width.
		// A half-wing run uses Y-min as its exact mirror plane. Put any indivisible
		// base-brick padding on the outer-span side so the symmetry plane cannot drift.
		domain_.lo = {requested.lo.x, anchor_y_min ? requested.lo.y : requested.lo.y - 0.5 * excess.y,
			requested.lo.z - 0.5 * excess.z};
		domain_.hi = domain_.lo + padded_size;
		levels_.assign(max_levels_, {});
		// Grid topology is constructed in world-space doubles. Keeping h as FP32 while
		// domain padding and base-brick origins used the requested FP64 value produced
		// two subtly different lattices for values such as 0.2. At a brick boundary,
		// point location could then select the brick on the wrong side and corrupt the
		// pressure graph. Retain one exact host lattice; GPU field values remain FP32.
		for (int l = 0; l < max_levels_; ++l) { levels_[l].level = l; levels_[l].h = h / static_cast<double>(1 << l); }
		AmrLevel& base = levels_[0]; base.bricks.reserve(static_cast<std::size_t>(dims.x) * dims.y * dims.z);
		for (int z = 0; z < dims.z; ++z) for (int y = 0; y < dims.y; ++y) for (int x = 0; x < dims.x; ++x)
		{
			BrickMetadata b; b.level = 0; b.coord = {x, y, z}; b.h = base.h; b.origin = domain_.lo + Vec3d{x * bw, y * bw, z * bw}; base.bricks.push_back(b);
		}
	}

	AmrHierarchy AmrHierarchy::uniform(const Aabb3d& domain, double h, int bs, int ghosts)
	{
		AmrHierarchy out; out.brick_size_ = bs; out.ghost_cells_ = ghosts; out.max_levels_ = 1; out.initialize_base(domain, h); out.rebuild_level_tables_and_metadata(nullptr); return out;
	}

	AmrHierarchy AmrHierarchy::build_static(const Aabb3d& domain, const TriMesh& wing, const TriangleBvh& bvh,
		const AmrConfig& c, bool anchor_y_min)
	{
		AmrHierarchy out; out.brick_size_ = c.brick_size; out.ghost_cells_ = c.ghost_cells; out.max_levels_ = std::max(1, c.max_levels + c.topology_refinement_levels); out.initialize_base(domain, c.base_cell_size, anchor_y_min);
		out.rebuild_level_tables_and_metadata(&bvh);
		Aabb3d wb; wb.lo = {wing.bbox_min[0], wing.bbox_min[1], wing.bbox_min[2]}; wb.hi = {wing.bbox_max[0], wing.bbox_max[1], wing.bbox_max[2]};
		const Vec3d wc = (wb.lo + wb.hi) * 0.5;
		Aabb3d wake{{wb.hi.x, wc.y - c.wake_radius, wc.z - c.wake_radius}, {wb.hi.x + c.wake_length, wc.y + c.wake_radius, wc.z + c.wake_radius}};
		for (int l = 0; l + 1 < out.max_levels_; ++l)
		{
			const std::size_t count = out.levels_[l].bricks.size();
			std::vector<unsigned char> refine(count,0);
			parallel_for_ranges(static_cast<int>(count),[&](int begin,int end)
			{
			for(int id=begin;id<end;++id)
			{
				const BrickMetadata& b = out.levels_[l].bricks[id]; if (!b.active()) continue;
				const Aabb3d bb = brick_box(b, out.brick_size_);
				const bool wing_region = l == 0 && box_overlap(bb, expanded(wb, c.wing_refinement_distance));
				// Normal flow levels retain the configured physical collar. Optional
				// topology-only levels add one parent-brick halo around intersected
				// parents, keeping every EB cell away from a 2:1 interface without
				// multiplying the whole physical collar at the extra resolution.
				const double surface_distance = l + 1 < c.max_levels
					? c.surface_refinement_distance
					: out.brick_size_ * static_cast<double>(b.h);
				const bool surface_region = !bvh.query_aabb(expanded(bb, surface_distance)).empty();
				// The finest level is reserved for fabric topology. Refining the entire
				// volumetric wake to that level multiplies memory while doing nothing to
				// resolve ribs, openings, or trailing edges. Keep the near wake one level
				// coarser; its length/radius remain independently configurable.
				const bool wake_region = l + 1 < c.max_levels - 1 && box_overlap(bb, wake);
				refine[id]=static_cast<unsigned char>(wing_region || surface_region || wake_region);
			}
			});
			for(std::size_t id=0;id<count;++id)if(refine[id])out.refine_brick(l,static_cast<int>(id));
			out.rebuild_level_tables_and_metadata(&bvh);
		}
		out.enforce_balance(); out.rebuild_level_tables_and_metadata(&bvh); return out;
	}

	void AmrHierarchy::refine_brick(int level, int id)
	{
		if (level < 0 || level + 1 >= max_levels_ || id < 0 || id >= static_cast<int>(levels_[level].bricks.size())) return;
		BrickMetadata parent_copy = levels_[level].bricks[id]; if (!parent_copy.active()) return;
		levels_[level].bricks[id].flags |= BRICK_COVERED;
		AmrLevel& fine = levels_[level + 1];
		for (int dz = 0; dz < 2; ++dz) for (int dy = 0; dy < 2; ++dy) for (int dx = 0; dx < 2; ++dx)
		{
			const int slot = (dz * 2 + dy) * 2 + dx; const Int3 cc{parent_copy.coord.x * 2 + dx, parent_copy.coord.y * 2 + dy, parent_copy.coord.z * 2 + dz};
			int child = find_brick(level + 1, cc);
			if (child < 0)
			{
				BrickMetadata b; b.level = level + 1; b.coord = cc; b.h = fine.h; b.parent = id;
				const double bw = brick_size_ * static_cast<double>(b.h); b.origin = domain_.lo + Vec3d{cc.x * bw, cc.y * bw, cc.z * bw};
				child = static_cast<int>(fine.bricks.size()); fine.bricks.push_back(b);
			}
			levels_[level].bricks[id].children[slot] = child;
		}
	}

	void AmrHierarchy::rebuild_level_tables_and_metadata(const TriangleBvh* bvh)
	{
		for (AmrLevel& lev : levels_)
		{
			std::size_t cap = 1; while (cap < std::max<std::size_t>(2, lev.bricks.size() * 2)) cap <<= 1;
			lev.lookup.assign(cap, {}); lev.lookup_mask = static_cast<std::uint32_t>(cap - 1); lev.active_bricks = 0;
			for (int id = 0; id < static_cast<int>(lev.bricks.size()); ++id)
			{
				BrickMetadata& b = lev.bricks[id]; b.same_level_neighbor.fill(-1); b.flags &= (BRICK_COVERED);
				std::uint32_t slot = static_cast<std::uint32_t>(coord_hash(b.coord)) & lev.lookup_mask;
				while (lev.lookup[slot].brick_id >= 0) slot = (slot + 1) & lev.lookup_mask;
				lev.lookup[slot] = {b.coord, id}; if (b.active()) ++lev.active_bricks;
			}
			parallel_for_ranges(static_cast<int>(lev.bricks.size()),[&](int begin,int end)
			{
			for(int id=begin;id<end;++id)
			{
				BrickMetadata& b = lev.bricks[id];
				const Int3 d[6] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
				for (int q = 0; q < 6; ++q) b.same_level_neighbor[q] = find_brick(lev.level, {b.coord.x + d[q].x, b.coord.y + d[q].y, b.coord.z + d[q].z});
				const double bw = brick_size_ * static_cast<double>(b.h); const Vec3d bh = b.origin + Vec3d{bw,bw,bw}; const double eps = 1e-9 * std::max(1.0, bw);
				if (std::abs(b.origin.x - domain_.lo.x) <= eps) b.flags |= BRICK_XMIN; if (std::abs(bh.x - domain_.hi.x) <= eps) b.flags |= BRICK_XMAX;
				if (std::abs(b.origin.y - domain_.lo.y) <= eps) b.flags |= BRICK_YMIN; if (std::abs(bh.y - domain_.hi.y) <= eps) b.flags |= BRICK_YMAX;
				if (std::abs(b.origin.z - domain_.lo.z) <= eps) b.flags |= BRICK_ZMIN; if (std::abs(bh.z - domain_.hi.z) <= eps) b.flags |= BRICK_ZMAX;
				if (bvh && !bvh->query_aabb(brick_box(b, brick_size_)).empty()) b.flags |= BRICK_EMBEDDED_BOUNDARY;
			}
			});
		}
	}

	int AmrHierarchy::find_brick(int level, Int3 c) const
	{
		if (level < 0 || level >= static_cast<int>(levels_.size())) return -1; const AmrLevel& lev = levels_[level]; if (lev.lookup.empty()) return -1;
		std::uint32_t slot = static_cast<std::uint32_t>(coord_hash(c)) & lev.lookup_mask;
		for (std::size_t probe = 0; probe < lev.lookup.size(); ++probe)
		{
			const BrickLookupEntry& e = lev.lookup[slot]; if (e.brick_id < 0) return -1; if (e.coord == c) return e.brick_id; slot = (slot + 1) & lev.lookup_mask;
		}
		return -1;
	}

	bool AmrHierarchy::point_in_domain(Vec3d p) const { return p.x >= domain_.lo.x && p.y >= domain_.lo.y && p.z >= domain_.lo.z && p.x < domain_.hi.x && p.y < domain_.hi.y && p.z < domain_.hi.z; }

	BrickLocation AmrHierarchy::locate_finest(Vec3d p) const
	{
		BrickLocation out; if (!point_in_domain(p)) return out;
		for (int l = static_cast<int>(levels_.size()) - 1; l >= 0; --l)
		{
			const double h = levels_[l].h, bw = brick_size_ * h; const Vec3d q = p - domain_.lo;
			const Int3 bc{static_cast<int>(std::floor(q.x / bw)), static_cast<int>(std::floor(q.y / bw)), static_cast<int>(std::floor(q.z / bw))}; const int id = find_brick(l, bc);
			if (id < 0 || !levels_[l].bricks[id].active()) continue;
			const Vec3d r = p - levels_[l].bricks[id].origin; out.level = l; out.brick = id;
			out.cell = {std::clamp(static_cast<int>(std::floor(r.x / h)), 0, brick_size_ - 1), std::clamp(static_cast<int>(std::floor(r.y / h)), 0, brick_size_ - 1), std::clamp(static_cast<int>(std::floor(r.z / h)), 0, brick_size_ - 1)}; return out;
		}
		return out;
	}

	void AmrHierarchy::enforce_balance()
	{
		for (int iteration = 0; iteration < max_levels_ * 4; ++iteration)
		{
			std::vector<std::pair<int,int>> refine;
			for (int l = 0; l < static_cast<int>(levels_.size()); ++l) for (int id = 0; id < static_cast<int>(levels_[l].bricks.size()); ++id)
			{
				const BrickMetadata& b = levels_[l].bricks[id]; if (!b.active()) continue; const double w = brick_size_ * static_cast<double>(b.h), eps = 1e-7 * w;
				for (int face = 0; face < 6; ++face) for (double a : {0.25, 0.75}) for (double c : {0.25, 0.75})
				{
					Vec3d p = b.origin + Vec3d{0.5*w,0.5*w,0.5*w};
					if (face < 2) { p.x = b.origin.x + (face ? w + eps : -eps); p.y = b.origin.y + a*w; p.z = b.origin.z + c*w; }
					else if (face < 4) { p.y = b.origin.y + ((face&1) ? w + eps : -eps); p.x = b.origin.x + a*w; p.z = b.origin.z + c*w; }
					else { p.z = b.origin.z + ((face&1) ? w + eps : -eps); p.x = b.origin.x + a*w; p.y = b.origin.y + c*w; }
					const BrickLocation n = locate_finest(p); if (!n.found() || std::abs(n.level - l) <= 1) continue;
					if (n.level < l) refine.push_back({n.level, n.brick}); else refine.push_back({l, id});
				}
			}
			std::sort(refine.begin(), refine.end()); refine.erase(std::unique(refine.begin(), refine.end()), refine.end()); if (refine.empty()) break;
			for (auto [l,id] : refine) refine_brick(l,id); rebuild_level_tables_and_metadata(nullptr);
		}
	}

	bool AmrHierarchy::is_two_to_one_balanced() const
	{
		for (int l = 0; l < static_cast<int>(levels_.size()); ++l) for (const BrickMetadata& b : levels_[l].bricks) if (b.active())
		{
			const double w = brick_size_ * static_cast<double>(b.h), eps=1e-7*w; const Vec3d c=b.origin+Vec3d{0.5*w,0.5*w,0.5*w};
			for (Vec3d d : {Vec3d{-0.5*w-eps,0,0},Vec3d{0.5*w+eps,0,0},Vec3d{0,-0.5*w-eps,0},Vec3d{0,0.5*w+eps,0},Vec3d{0,0,-0.5*w-eps},Vec3d{0,0,0.5*w+eps}})
			{ BrickLocation n=locate_finest(c+d); if(n.found() && std::abs(n.level-l)>1) return false; }
		}
		return true;
	}

	std::size_t AmrHierarchy::active_brick_count() const { std::size_t n=0; for(const auto& l:levels_) n+=l.active_bricks; return n; }
	std::size_t AmrHierarchy::active_cell_count() const { const std::size_t c=static_cast<std::size_t>(brick_size_)*brick_size_*brick_size_; return active_brick_count()*c; }
}
