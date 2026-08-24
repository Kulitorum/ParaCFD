#include "core/fluid/amr_eb.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>
#include <thread>
#include <utility>

namespace paracfd::core
{
	namespace
	{
		double elapsed_ms(std::chrono::steady_clock::time_point begin,
			std::chrono::steady_clock::time_point end)
		{
			return std::chrono::duration<double,std::milli>(end-begin).count();
		}

		template<class Work>
		unsigned parallel_for_cells(int count,Work&& work)
		{
			if(count<=0)return 0;
			unsigned workers=std::max(1u,std::thread::hardware_concurrency());
			if(const char* configured=std::getenv("PARACFD_GRID_THREADS"))
			{
				char* end=nullptr;const long parsed=std::strtol(configured,&end,10);
				if(end!=configured&&parsed>0)workers=static_cast<unsigned>(parsed);
			}
			workers=std::min(workers,static_cast<unsigned>(count));
			if(workers<=1)
			{
				for(int cell=0;cell<count;++cell)work(cell);
				return 1;
			}
			std::atomic<int> next{0};std::atomic<bool> failed{false};
			std::exception_ptr error;std::mutex error_mutex;
			constexpr int grain=4096;
			auto run=[&]
			{
				try
				{
					while(!failed.load(std::memory_order_relaxed))
					{
						const int begin=next.fetch_add(grain,std::memory_order_relaxed);
						if(begin>=count)break;const int end=std::min(count,begin+grain);
						for(int cell=begin;cell<end;++cell)work(cell);
					}
				}
				catch(...)
				{
					failed.store(true,std::memory_order_relaxed);
					std::lock_guard<std::mutex> lock(error_mutex);if(!error)error=std::current_exception();
				}
			};
			std::vector<std::thread> threads;threads.reserve(workers-1);
			for(unsigned worker=1;worker<workers;++worker)threads.emplace_back(run);
			run();for(auto& thread:threads)thread.join();if(error)std::rethrow_exception(error);
			return workers;
		}
	}

	std::size_t AmrEmbeddedBoundaryAtlas::unresolved_count() const
	{
		std::size_t count=0;for(const AmrEbLevelAtlas& atlas:levels)count+=atlas.topology.unresolved.size();return count;
	}

	std::size_t AmrEmbeddedBoundaryAtlas::owned_unresolved_count() const
	{
		std::size_t count=0;for(const AmrEbLevelAtlas& atlas:levels)for(const UnresolvedEbCell& problem:atlas.topology.unresolved)if(problem.parent_cell>=0&&problem.parent_cell<static_cast<int>(atlas.owned_cell.size())&&atlas.owned_cell[problem.parent_cell])++count;return count;
	}

	AmrEmbeddedBoundaryAtlas build_amr_embedded_boundary_atlas(const AmrHierarchy& hierarchy,
		const TriMesh& mesh, const TriangleBvh& bvh, const EmbeddedBoundaryBuildOptions& options)
	{
		if(!options.closed_solid)
			throw std::invalid_argument("AMR embedded-boundary construction requires a closed solid");
		AmrEmbeddedBoundaryAtlas result;const int bs=hierarchy.brick_size();
		for(int level=0;level<static_cast<int>(hierarchy.levels().size());++level)
		{
			const auto level_begin=std::chrono::steady_clock::now();
			const AmrLevel& source=hierarchy.levels()[level];Int3 minimum{std::numeric_limits<int>::max(),std::numeric_limits<int>::max(),std::numeric_limits<int>::max()},maximum{std::numeric_limits<int>::min(),std::numeric_limits<int>::min(),std::numeric_limits<int>::min()};bool found=false;
			for(const BrickMetadata& brick:source.bricks)if(brick.active())
			{
				const double width=bs*static_cast<double>(brick.h);const Aabb3d brick_box{brick.origin,brick.origin+Vec3d{width,width,width}};if(!brick.embedded_boundary()&&bvh.query_aabb(brick_box).empty())continue;found=true;minimum.x=std::min(minimum.x,brick.coord.x);minimum.y=std::min(minimum.y,brick.coord.y);minimum.z=std::min(minimum.z,brick.coord.z);maximum.x=std::max(maximum.x,brick.coord.x);maximum.y=std::max(maximum.y,brick.coord.y);maximum.z=std::max(maximum.z,brick.coord.z);
			}
			if(!found)continue;const double h=source.h;const Vec3d extent=hierarchy.domain().hi-hierarchy.domain().lo;const Int3 domain_cells{static_cast<int>(std::llround(extent.x/h)),static_cast<int>(std::llround(extent.y/h)),static_cast<int>(std::llround(extent.z/h))};const Int3 cell_lo{std::max(0,minimum.x*bs-1),std::max(0,minimum.y*bs-1),std::max(0,minimum.z*bs-1)},cell_hi{std::min(domain_cells.x,(maximum.x+1)*bs+1),std::min(domain_cells.y,(maximum.y+1)*bs+1),std::min(domain_cells.z,(maximum.z+1)*bs+1)};const Vec3d origin=hierarchy.domain().lo+Vec3d{cell_lo.x*h,cell_lo.y*h,cell_lo.z*h};const UniformEbGrid grid{origin,cell_hi.x-cell_lo.x,cell_hi.y-cell_lo.y,cell_hi.z-cell_lo.z,h};
			EmbeddedBoundaryBuildOptions level_options=options;
			AmrEbLevelAtlas atlas;atlas.level=level;atlas.minimum_brick_coord=minimum;atlas.maximum_brick_coord=maximum;atlas.topology=build_closed_solid_embedded_boundary(mesh,bvh,*level_options.closed_solid,grid,level_options);const auto topology_end=std::chrono::steady_clock::now();atlas.owned_cell.assign(grid.cell_count(),0);
			const unsigned ownership_workers=parallel_for_cells(grid.cell_count(),[&](int cell){const BrickLocation owner=hierarchy.locate_finest(grid.cell_centroid(cell));atlas.owned_cell[cell]=static_cast<unsigned char>(owner.found()&&owner.level==level);});const auto ownership_end=std::chrono::steady_clock::now();
			// Rebuild stabilization from OWNED fluid only. The geometry builder sees a
			// one-cell topology halo and cannot know which volumes survive composite AMR
			// ownership. Every accepted union below follows a real open aperture, remains
			// local, keeps every surviving aperture bracketed by its two centroids, and
			// preserves distinct raw plus/minus pressure sides of every fabric patch.
			EmbeddedBoundary& eb=atlas.topology;const int fragment_count=static_cast<int>(eb.fragments.size());const double cell_volume=h*h*h,minimum_volume=options.min_volume_fraction*cell_volume,maximum_merge_span=2*h,bracket_epsilon=1e-10*h;
			auto ref_cell=[&](FragmentRef ref){if(ref==invalid_fragment)return -1;return fragment_is_regular(ref)?regular_fragment_cell(ref):eb.fragments[irregular_fragment_index(ref)].parent_cell;};
			auto ref_owned=[&](FragmentRef ref){const int cell=ref_cell(ref);return cell>=0&&cell<static_cast<int>(atlas.owned_cell.size())&&atlas.owned_cell[cell]!=0;};
			std::vector<unsigned char> reported_cell(grid.cell_count(),0);auto report_unresolved=[&](int cell,std::string reason,std::vector<std::uint32_t> triangles=std::vector<std::uint32_t>{}){if(cell<0||cell>=grid.cell_count()||!atlas.owned_cell[cell]||reported_cell[cell])return;reported_cell[cell]=1;eb.cells[cell].state=EbCellState::unresolved;eb.unresolved.push_back({cell,std::move(triangles),std::move(reason)});};
			auto finite_vec=[](const Vec3d& value){return std::isfinite(value.x)&&std::isfinite(value.y)&&std::isfinite(value.z);};
			auto valid_aperture_geometry=[&](const FaceAperture& aperture){return aperture.axis>=0&&aperture.axis<=2&&std::isfinite(aperture.area)&&aperture.area>0&&finite_vec(aperture.centroid);};
			// Audit every record touching owned fluid before any centroid indexing or merge.
			// A mixed ownership edge needs an aperture-aware coarse/fine transaction; it
			// must not disappear merely because this sparse level atlas cannot own both ends.
			for(const FaceAperture& aperture:eb.apertures)
			{
				const bool a_owned=ref_owned(aperture.fragment_a),b_owned=ref_owned(aperture.fragment_b);
				if(a_owned!=b_owned){const FragmentRef owned=a_owned?aperture.fragment_a:aperture.fragment_b;report_unresolved(ref_cell(owned),"fluid aperture crosses an AMR ownership boundary; aperture-aware coarse/fine topology is required");continue;}
				if(a_owned&&(!valid_aperture_geometry(aperture)||aperture.fragment_a==aperture.fragment_b))report_unresolved(ref_cell(aperture.fragment_a),"owned fluid aperture has invalid endpoints or non-positive/non-finite geometry; refine and rebuild");
			}
			for(const SurfacePatch& patch:eb.patches)
			{
				if(patch.plus_fragment==invalid_fragment||patch.minus_fragment==invalid_fragment)
				{
					const FragmentRef fluid=patch.plus_fragment!=invalid_fragment?
						patch.plus_fragment:patch.minus_fragment;
					if(fluid==invalid_fragment)
						report_unresolved(-1,"wall patch has no adjacent fluid fragment");
					continue;
				}
				const bool plus_owned=ref_owned(patch.plus_fragment),minus_owned=ref_owned(patch.minus_fragment);
				if(plus_owned!=minus_owned){const FragmentRef owned=plus_owned?patch.plus_fragment:patch.minus_fragment;report_unresolved(ref_cell(owned),"fabric patch crosses an AMR ownership boundary; both pressure sides must be represented on one composite topology");}
			}
			struct Aggregate{double volume=0;Vec3d weighted_centroid{},lo{},hi{};std::vector<FragmentRef> members;};
			std::vector<FragmentRef> merge_parent(fragment_count,invalid_fragment);std::map<FragmentRef,Aggregate> aggregates;
			for(int fragment=0;fragment<fragment_count;++fragment)if(atlas.owned_cell[eb.fragments[fragment].parent_cell]){const FragmentRef ref=irregular_fragment(fragment);merge_parent[fragment]=ref;Aggregate value;value.volume=eb.fragments[fragment].volume;value.weighted_centroid=eb.fragments[fragment].centroid*value.volume;value.lo=value.hi=eb.fragments[fragment].centroid;value.members.push_back(ref);aggregates.emplace(ref,std::move(value));eb.fragments[fragment].merge_target=ref;eb.fragments[fragment].pressure_dof=eb.fragments[fragment].parent_cell;eb.fragments[fragment].pressure_static=false;}
			std::function<FragmentRef(FragmentRef)> find_root=[&](FragmentRef ref)->FragmentRef{if(ref==invalid_fragment||fragment_is_regular(ref))return ref;const int fragment=irregular_fragment_index(ref);if(fragment<0||fragment>=fragment_count||merge_parent[fragment]==invalid_fragment)return invalid_fragment;const FragmentRef parent=merge_parent[fragment];if(parent==ref)return ref;return merge_parent[fragment]=find_root(parent);};
			auto aggregate=[&](FragmentRef ref)->Aggregate&{ref=find_root(ref);auto found=aggregates.find(ref);if(found!=aggregates.end())return found->second;if(!fragment_is_regular(ref)||!ref_owned(ref))throw std::runtime_error("missing owned EB aggregate");const Vec3d centre=eb.grid.cell_centroid(regular_fragment_cell(ref));Aggregate value;value.volume=cell_volume;value.weighted_centroid=centre*cell_volume;value.lo=value.hi=centre;value.members.push_back(ref);return aggregates.emplace(ref,std::move(value)).first->second;};
			auto centroid=[&](FragmentRef ref){const Aggregate& value=aggregate(ref);return value.weighted_centroid/value.volume;};
			std::map<FragmentRef,std::vector<int>> incident_apertures,all_owned_incident_apertures;for(int q=0;q<static_cast<int>(eb.apertures.size());++q){const FaceAperture& aperture=eb.apertures[q];if(ref_owned(aperture.fragment_a))all_owned_incident_apertures[aperture.fragment_a].push_back(q);if(ref_owned(aperture.fragment_b))all_owned_incident_apertures[aperture.fragment_b].push_back(q);if(!ref_owned(aperture.fragment_a)||!ref_owned(aperture.fragment_b)||!valid_aperture_geometry(aperture)||aperture.fragment_a==aperture.fragment_b)continue;incident_apertures[aperture.fragment_a].push_back(q);incident_apertures[aperture.fragment_b].push_back(q);aggregate(aperture.fragment_a);aggregate(aperture.fragment_b);}
			std::map<FragmentRef,std::vector<FragmentRef>> forbidden_patch_pair;for(const SurfacePatch& patch:eb.patches)if(patch.plus_fragment!=invalid_fragment&&patch.minus_fragment!=invalid_fragment&&patch.plus_fragment!=patch.minus_fragment&&ref_owned(patch.plus_fragment)&&ref_owned(patch.minus_fragment)){forbidden_patch_pair[patch.plus_fragment].push_back(patch.minus_fragment);forbidden_patch_pair[patch.minus_fragment].push_back(patch.plus_fragment);aggregate(patch.plus_fragment);aggregate(patch.minus_fragment);}
			auto combined_span=[&](FragmentRef a,FragmentRef b){const Aggregate& aa=aggregate(a);const Aggregate& bb=aggregate(b);const Vec3d lo{std::min(aa.lo.x,bb.lo.x),std::min(aa.lo.y,bb.lo.y),std::min(aa.lo.z,bb.lo.z)},hi{std::max(aa.hi.x,bb.hi.x),std::max(aa.hi.y,bb.hi.y),std::max(aa.hi.z,bb.hi.z)};return std::max({hi.x-lo.x,hi.y-lo.y,hi.z-lo.z});};
			enum class UnionRejection
			{
				none,
				invalid_or_same,
				ordinary_pair,
				span,
				fabric_side,
				missing_neighbour,
				aperture_bracketing
			};
			auto union_rejection=[&](FragmentRef a,FragmentRef b)
			{
				a=find_root(a);b=find_root(b);
				if(a==invalid_fragment||b==invalid_fragment||a==b)return UnionRejection::invalid_or_same;
				if(fragment_is_regular(a)&&fragment_is_regular(b))return UnionRejection::ordinary_pair;
				Aggregate& aa=aggregate(a);Aggregate& bb=aggregate(b);
				if(combined_span(a,b)>maximum_merge_span+1e-12*h)return UnionRejection::span;
				auto crosses_patch=[&](const Aggregate& source,FragmentRef other){for(const FragmentRef member:source.members){const auto found=forbidden_patch_pair.find(member);if(found==forbidden_patch_pair.end())continue;for(const FragmentRef forbidden:found->second)if(find_root(forbidden)==other)return true;}return false;};
				if(crosses_patch(aa,b)||crosses_patch(bb,a))return UnionRejection::fabric_side;
				const double volume=aa.volume+bb.volume;const Vec3d merged_centroid=(aa.weighted_centroid+bb.weighted_centroid)/volume;std::vector<int> edges;for(const FragmentRef member:aa.members){const auto found=incident_apertures.find(member);if(found!=incident_apertures.end())edges.insert(edges.end(),found->second.begin(),found->second.end());}for(const FragmentRef member:bb.members){const auto found=incident_apertures.find(member);if(found!=incident_apertures.end())edges.insert(edges.end(),found->second.begin(),found->second.end());}std::sort(edges.begin(),edges.end());edges.erase(std::unique(edges.begin(),edges.end()),edges.end());for(const int edge:edges){const FaceAperture& aperture=eb.apertures[edge];FragmentRef lower=find_root(aperture.fragment_a),upper=find_root(aperture.fragment_b);const bool lower_merged=lower==a||lower==b,upper_merged=upper==a||upper==b;if(lower_merged&&upper_merged)continue;if(lower==invalid_fragment||upper==invalid_fragment)return UnionRejection::missing_neighbour;const Vec3d lower_centroid=lower_merged?merged_centroid:centroid(lower),upper_centroid=upper_merged?merged_centroid:centroid(upper);const double coordinate=aperture.centroid[aperture.axis];if(coordinate-lower_centroid[aperture.axis]<=bracket_epsilon||upper_centroid[aperture.axis]-coordinate<=bracket_epsilon)return UnionRejection::aperture_bracketing;}return UnionRejection::none;
			};
			auto legal_union=[&](FragmentRef a,FragmentRef b){return union_rejection(a,b)==UnionRejection::none;};
			auto unite=[&](FragmentRef a,FragmentRef b){a=find_root(a);b=find_root(b);if(a==b)return a;Aggregate& aa=aggregate(a);Aggregate& bb=aggregate(b);FragmentRef keep=a,drop=b;if(fragment_is_regular(b)||(!fragment_is_regular(a)&&(bb.volume>aa.volume||(bb.volume==aa.volume&&b<a)))){keep=b;drop=a;}if(fragment_is_regular(drop))throw std::runtime_error("attempted to agglomerate two ordinary EB control volumes");Aggregate& kept=aggregate(keep);Aggregate dropped=aggregate(drop);kept.volume+=dropped.volume;kept.weighted_centroid=kept.weighted_centroid+dropped.weighted_centroid;kept.lo={std::min(kept.lo.x,dropped.lo.x),std::min(kept.lo.y,dropped.lo.y),std::min(kept.lo.z,dropped.lo.z)};kept.hi={std::max(kept.hi.x,dropped.hi.x),std::max(kept.hi.y,dropped.hi.y),std::max(kept.hi.z,dropped.hi.z)};kept.members.insert(kept.members.end(),dropped.members.begin(),dropped.members.end());merge_parent[irregular_fragment_index(drop)]=keep;aggregates.erase(drop);return keep;};
			// Grow every sub-threshold aggregate by one legal aperture at a time. Prefer
			// irregular donors so an ordinary structured cell is not displaced unnecessarily.
			for(int pass=0;pass<=fragment_count;++pass)
			{
				bool found_small=false,progress=false;for(int fragment=0;fragment<fragment_count;++fragment){FragmentRef root=find_root(irregular_fragment(fragment));if(root!=irregular_fragment(fragment)||root==invalid_fragment||aggregate(root).volume+1e-12*cell_volume>=minimum_volume)continue;found_small=true;struct Choice{FragmentRef other=invalid_fragment;bool sufficient=false,irregular=false;double span=0,volume=0,area=0;}best;std::map<FragmentRef,double> neighbours;for(const FragmentRef member:aggregate(root).members){const auto incident=incident_apertures.find(member);if(incident==incident_apertures.end())continue;for(const int edge:incident->second){const FaceAperture& aperture=eb.apertures[edge];FragmentRef other=find_root(find_root(aperture.fragment_a)==root?aperture.fragment_b:aperture.fragment_a);if(other!=invalid_fragment&&other!=root)neighbours[other]=std::max(neighbours[other],aperture.area);}}for(const auto& candidate:neighbours){FragmentRef other=find_root(candidate.first);if(other==invalid_fragment||other==root||!legal_union(root,other))continue;Choice value;value.other=other;value.sufficient=aggregate(root).volume+aggregate(other).volume+1e-12*cell_volume>=minimum_volume;value.irregular=!fragment_is_regular(other);value.span=combined_span(root,other);value.volume=aggregate(other).volume;value.area=candidate.second;bool better=best.other==invalid_fragment;if(!better)better=value.sufficient!=best.sufficient?value.sufficient:((value.irregular!=best.irregular)?value.irregular:(value.span!=best.span?value.span<best.span:(value.volume!=best.volume?value.volume>best.volume:(value.area!=best.area?value.area>best.area:value.other<best.other))));if(better)best=value;}if(best.other!=invalid_fragment){unite(root,best.other);progress=true;}}
				if(!found_small||!progress)break;if(pass==fragment_count)throw std::runtime_error("constrained EB small-cell agglomeration did not terminate");
			}
			// Repair an initially unbracketed raw connection only if removing that edge by
			// union preserves every other topology invariant. Otherwise request refinement.
			for(int pass=0;pass<=fragment_count;++pass){bool changed=false,bad=false;for(const FaceAperture& aperture:eb.apertures){if(!ref_owned(aperture.fragment_a)||!ref_owned(aperture.fragment_b)||!valid_aperture_geometry(aperture))continue;FragmentRef a=find_root(aperture.fragment_a),b=find_root(aperture.fragment_b);if(a==invalid_fragment||b==invalid_fragment||a==b)continue;const double coordinate=aperture.centroid[aperture.axis],lower_distance=coordinate-centroid(a)[aperture.axis],upper_distance=centroid(b)[aperture.axis]-coordinate;if(lower_distance>bracket_epsilon&&upper_distance>bracket_epsilon)continue;bad=true;if(legal_union(a,b)){unite(a,b);++atlas.quality_agglomerations;changed=true;break;}report_unresolved(ref_cell(aperture.fragment_a),"no local topology-safe agglomeration brackets aperture: axis="+std::to_string(static_cast<int>(aperture.axis))+" lower_distance/h="+std::to_string(lower_distance/h)+" upper_distance/h="+std::to_string(upper_distance/h)+"; refine and rebuild");}if(!changed){(void)bad;break;}if(pass==fragment_count)throw std::runtime_error("constrained EB quality agglomeration did not terminate");}
			// Pressure activity follows represented fluid topology, not whether a component
			// happens to touch an ordinary Cartesian cell.  A closed sub-cell channel can be
			// made entirely of irregular control volumes and still carry conservative flux
			// through its positive-area apertures.  Such a component must remain dynamic;
			// the composite pressure builder supplies its otherwise-missing pressure gauge.
			// Only a single pressure root with no represented open edge is genuinely static.
			std::map<FragmentRef,std::vector<FragmentRef>> aperture_graph;
			for(const auto& item:aggregates)aperture_graph[find_root(item.first)];
			for(const FaceAperture& aperture:eb.apertures)
				if(valid_aperture_geometry(aperture)&&ref_owned(aperture.fragment_a)&&ref_owned(aperture.fragment_b))
				{
					const FragmentRef a=find_root(aperture.fragment_a),b=find_root(aperture.fragment_b);
					if(a==invalid_fragment||b==invalid_fragment||a==b)continue;
					aperture_graph[a].push_back(b);aperture_graph[b].push_back(a);
				}
			for(auto& item:aperture_graph)
			{
				auto& neighbours=item.second;std::sort(neighbours.begin(),neighbours.end());
				neighbours.erase(std::unique(neighbours.begin(),neighbours.end()),neighbours.end());
			}
			std::set<FragmentRef> visited;
			for(const auto& item:aperture_graph)if(!visited.count(item.first))
			{
				std::vector<FragmentRef> todo{item.first},component;bool regular_anchor=false,represented_edge=false;
				while(!todo.empty())
				{
					const FragmentRef root=find_root(todo.back());todo.pop_back();
					if(root==invalid_fragment||!visited.insert(root).second)continue;
					component.push_back(root);regular_anchor=regular_anchor||fragment_is_regular(root);
					const auto found=aperture_graph.find(root);if(found==aperture_graph.end())continue;
					represented_edge=represented_edge||!found->second.empty();
					for(const FragmentRef neighbour:found->second)todo.push_back(neighbour);
				}
				bool component_reported=false;
				for(const FragmentRef root:component)for(const FragmentRef member:aggregate(root).members){const int cell=ref_cell(member);component_reported=component_reported||(cell>=0&&cell<grid.cell_count()&&reported_cell[cell]);}
				if(regular_anchor||represented_edge||component_reported)continue;
				++atlas.static_subcell_components;
				for(const FragmentRef root:component)if(!fragment_is_regular(root))
				{
					++atlas.static_subcell_aggregates;atlas.static_subcell_volume+=aggregate(root).volume;
					eb.fragments[irregular_fragment_index(root)].pressure_static=true;
				}
			}
			// Materialize the deterministic parent forest for the composite builder.
			for(int fragment=0;fragment<fragment_count;++fragment)if(atlas.owned_cell[eb.fragments[fragment].parent_cell]){const FragmentRef root=find_root(irregular_fragment(fragment));if(root==invalid_fragment)throw std::runtime_error("owned constrained EB merge escaped its level");eb.fragments[fragment].merge_target=root;eb.fragments[fragment].pressure_dof=fragment_is_regular(root)?regular_fragment_cell(root):eb.fragments[irregular_fragment_index(root)].parent_cell;}
			for(const auto& item:aggregates){const Vec3d extent=item.second.hi-item.second.lo;atlas.maximum_aggregate_span_cells=std::max(atlas.maximum_aggregate_span_cells,std::max({extent.x,extent.y,extent.z})/h);}
			// Final invariants are intentionally redundant with legal_union: a future merge
			// policy change must fail preprocessing rather than silently short-circuit fabric.
			for(const SurfacePatch& patch:eb.patches)if(patch.plus_fragment!=invalid_fragment&&patch.minus_fragment!=invalid_fragment&&patch.plus_fragment!=patch.minus_fragment&&ref_owned(patch.plus_fragment)&&ref_owned(patch.minus_fragment)&&find_root(patch.plus_fragment)==find_root(patch.minus_fragment))report_unresolved(ref_cell(patch.plus_fragment),"agglomeration collapsed distinct plus/minus fabric pressure sides; refine and rebuild",{patch.source_triangle_id});
			for(const FaceAperture& aperture:eb.apertures)if(ref_owned(aperture.fragment_a)&&ref_owned(aperture.fragment_b)&&valid_aperture_geometry(aperture)){const FragmentRef a=find_root(aperture.fragment_a),b=find_root(aperture.fragment_b);if(a==invalid_fragment||b==invalid_fragment||a==b)continue;const double coordinate=aperture.centroid[aperture.axis],lower_distance=coordinate-centroid(a)[aperture.axis],upper_distance=centroid(b)[aperture.axis]-coordinate;if(!std::isfinite(lower_distance)||!std::isfinite(upper_distance)||lower_distance<=bracket_epsilon||upper_distance<=bracket_epsilon)report_unresolved(ref_cell(aperture.fragment_a),"final EB aggregate centroids do not bracket aperture; refine and rebuild");}
			for(const auto& item:aggregates)
			{
				if(fragment_is_regular(item.first)||find_root(item.first)!=item.first||
					item.second.volume+1e-12*cell_volume>=minimum_volume)continue;
				const int root_fragment=irregular_fragment_index(item.first);
				if(eb.fragments[root_fragment].pressure_static)continue;
				std::vector<int> edges;
				for(const FragmentRef member:item.second.members)
				{
					const auto incident=all_owned_incident_apertures.find(member);
					if(incident!=all_owned_incident_apertures.end())
						edges.insert(edges.end(),incident->second.begin(),incident->second.end());
				}
				std::sort(edges.begin(),edges.end());
				edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
				int open=0,invalid_incident=0;
				const Vec3d root_centroid=centroid(item.first);
				bool raw_topology_valid=std::isfinite(item.second.volume)&&item.second.volume>0&&
					finite_vec(root_centroid),raw_apertures_bracketed=true;
				std::map<FragmentRef,UnionRejection> donor_rejections;
				for(const int edge:edges)
				{
					const FaceAperture& aperture=eb.apertures[edge];
					if(!valid_aperture_geometry(aperture)||!ref_owned(aperture.fragment_a)||
						!ref_owned(aperture.fragment_b)||aperture.fragment_a==aperture.fragment_b)
					{
						raw_topology_valid=false;
						++invalid_incident;
						continue;
					}
					const FragmentRef a=find_root(aperture.fragment_a),b=find_root(aperture.fragment_b);
					if(a==b)continue;
					++open;
					if(a==invalid_fragment||b==invalid_fragment)
					{
						raw_topology_valid=false;
						raw_apertures_bracketed=false;
						++invalid_incident;
						continue;
					}
					const Vec3d lower_centroid=centroid(a),upper_centroid=centroid(b);
					const double coordinate=aperture.centroid[aperture.axis];
					const double lower_distance=coordinate-lower_centroid[aperture.axis];
					const double upper_distance=upper_centroid[aperture.axis]-coordinate;
					raw_topology_valid=raw_topology_valid&&std::isfinite(coordinate)&&
						finite_vec(lower_centroid)&&finite_vec(upper_centroid)&&
						std::isfinite(lower_distance)&&std::isfinite(upper_distance);
					raw_apertures_bracketed=raw_apertures_bracketed&&
						lower_distance>bracket_epsilon&&upper_distance>bracket_epsilon;
					const FragmentRef other=a==item.first?b:(b==item.first?a:invalid_fragment);
					if(other==invalid_fragment)raw_apertures_bracketed=false;
					else donor_rejections[other]=union_rejection(item.first,other);
				}
				if(open==0&&invalid_incident==0)
				{
					eb.fragments[root_fragment].pressure_static=true;
					continue;
				}

				int reject_span=0,reject_fabric=0,reject_bracketing=0,reject_other=0,
					regular_donors=0;
				bool legal_donor_remaining=false;
				for(const auto& rejection:donor_rejections)
				{
					regular_donors+=fragment_is_regular(rejection.first);
					switch(rejection.second)
					{
					case UnionRejection::span:++reject_span;break;
					case UnionRejection::aperture_bracketing:++reject_bracketing;break;
					case UnionRejection::fabric_side:
						++reject_fabric;break;
					case UnionRejection::none:
						legal_donor_remaining=true;break;
					default:++reject_other;break;
					}
				}
				const double volume_fraction=item.second.volume/cell_volume;
				const Real packed_volume=static_cast<Real>(item.second.volume);
				const bool production_volume_representable=
					std::isfinite(packed_volume)&&
					packed_volume>=std::numeric_limits<Real>::min();
				bool aggregate_already_unresolved=false;
				for(const FragmentRef member:item.second.members){const int cell=ref_cell(member);aggregate_already_unresolved=aggregate_already_unresolved||(cell>=0&&cell<grid.cell_count()&&reported_cell[cell]);}
				if(options.retain_signed_bracketed_small_roots_for_face_state&&
					!aggregate_already_unresolved&&open>0&&invalid_incident==0&&
					raw_topology_valid&&raw_apertures_bracketed&&!legal_donor_remaining&&
					production_volume_representable)
				{
					// This root remains an independent pressure state. The production MAC path
					// evolves aperture/dual momentum and projects integrated flux, so it has no
					// explicit update proportional to 1 / this fragment volume.
					eb.fragments[root_fragment].face_state_retained=true;
					++atlas.face_state_retained_small_roots;
					atlas.face_state_retained_small_volume+=item.second.volume;
					atlas.minimum_face_state_retained_volume_fraction=std::min(
						atlas.minimum_face_state_retained_volume_fraction,volume_fraction);
					continue;
				}

				const Vec3d extent=item.second.hi-item.second.lo;
				const double span=std::max({extent.x,extent.y,extent.z})/h;
				report_unresolved(eb.fragments[root_fragment].parent_cell,
					"owned fluid aggregate remains below min_volume_fraction after constrained agglomeration: level="+
					std::to_string(level)+" volume_fraction="+std::to_string(volume_fraction)+
					" span/h="+std::to_string(span)+" open_apertures="+std::to_string(open)+
					" invalid_incident="+std::to_string(invalid_incident)+
					" donor_roots="+std::to_string(donor_rejections.size())+
					" regular_donors="+std::to_string(regular_donors)+
					" raw_topology_valid="+std::to_string(raw_topology_valid?1:0)+
					" raw_bracketed="+std::to_string(raw_apertures_bracketed?1:0)+
					" production_volume_representable="+
						std::to_string(production_volume_representable?1:0)+
					" legal_donor_remaining="+std::to_string(legal_donor_remaining?1:0)+
					" rejected(span/fabric/bracketing/other)="+std::to_string(reject_span)+"/"+
					std::to_string(reject_fabric)+"/"+std::to_string(reject_bracketing)+"/"+
					std::to_string(reject_other)+"; refine and rebuild");
			}
			const auto level_end=std::chrono::steady_clock::now();
			std::fprintf(stderr,"[grid-profile] EB level %d: topology %.1f ms, ownership %.1f ms (%u CPU workers), stabilization %.1f ms, cells=%d\n",
				level,elapsed_ms(level_begin,topology_end),elapsed_ms(topology_end,ownership_end),
				ownership_workers,elapsed_ms(ownership_end,level_end),grid.cell_count());
			std::fflush(stderr);result.levels.push_back(std::move(atlas));
		}
		return result;
	}

	OneLevelEmbeddedBoundary build_one_level_embedded_boundary(const AmrHierarchy& hierarchy,
		const TriMesh& mesh, const TriangleBvh& bvh, const EmbeddedBoundaryBuildOptions& options)
	{
		if(!options.closed_solid)
			throw std::invalid_argument("one-level embedded-boundary construction requires a closed solid");
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
		OneLevelEmbeddedBoundary out; out.topology =
			build_closed_solid_embedded_boundary(mesh,bvh,*options.closed_solid,grid,options);
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
