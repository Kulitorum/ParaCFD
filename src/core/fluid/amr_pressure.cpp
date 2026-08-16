#include "core/fluid/amr_pressure.h"
#include "core/fluid/amr_eb.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace paracfd::core
{
	namespace
	{
		int local_index(int bs, int i, int j, int k) { return (k * bs + j) * bs + i; }
		void add_edge(std::vector<double>& output, const std::vector<double>& pressure, int a, int b, double coefficient)
		{
			if (a < 0 || b < 0 || a == b) return; const double flux = coefficient * (pressure[a] - pressure[b]); output[a] += flux; output[b] -= flux;
		}
		void finalize_component_gauges(CompositeAmrPressureSystem& system)
		{
			system.gauges.clear();if(!system.pressure_outlet_xmax)return;const int n=system.storage_size,bs=system.brick_size;std::vector<int> parent(n,-1);for(int q=0;q<n;++q)if(system.active[q])parent[q]=q;auto root=[&](int q){while(parent[q]!=q){parent[q]=parent[parent[q]];q=parent[q];}return q;};auto join=[&](int a,int b){if(a<0||b<0||!system.active[a]||!system.active[b])return;a=root(a);b=root(b);if(a!=b)parent[std::max(a,b)]=std::min(a,b);};
			for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level){const AmrLevel& source=system.hierarchy->levels()[level];for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick){const BrickMetadata& meta=source.bricks[brick];if(!meta.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=system.dof(level,brick,i,j,k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if(system.cut_face_mask[a]&(1u<<axis))continue;int c[3]={i,j,k},b=-1;if(++c[axis]<bs)b=system.dof(level,brick,c[0],c[1],c[2]);else{const int neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&source.bricks[neighbour].active()){c[axis]=0;b=system.dof(level,neighbour,c[0],c[1],c[2]);}}join(a,b);}}}}
			for(const CoarseFinePressureConnection& edge:system.coarse_fine)join(edge.coarse_dof,edge.fine_dof);for(const CoarseFinePressureConnection& edge:system.embedded)join(edge.coarse_dof,edge.fine_dof);
			std::vector<unsigned char> outlet_root(n,0);for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level){const AmrLevel& source=system.hierarchy->levels()[level];for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick){const BrickMetadata& meta=source.bricks[brick];if(!meta.active()||!(meta.flags&BRICK_XMAX))continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j){const int q=system.dof(level,brick,bs-1,j,k);if(system.active[q])outlet_root[root(q)]=1;}}}
			std::vector<unsigned char> emitted(n,0);for(int q=0;q<n;++q)if(system.active[q]){const int r=root(q);if(!outlet_root[r]&&!emitted[r]){emitted[r]=1;system.gauges.push_back({q,std::max(1e-12,std::cbrt(std::max(0.0,system.volume[q])))});}}
		}
		int base_aggregate_dof(const CompositeAmrPressureSystem& system,Vec3d point)
		{
			const AmrLevel& base=system.hierarchy->levels().front();const double brick_width=system.brick_size*base.h;const Aabb3d& domain=system.hierarchy->domain();Int3 coord{static_cast<int>(std::floor((point.x-domain.lo.x)/brick_width)),static_cast<int>(std::floor((point.y-domain.lo.y)/brick_width)),static_cast<int>(std::floor((point.z-domain.lo.z)/brick_width))};const int brick=system.hierarchy->find_brick(0,coord);if(brick<0)return -1;const BrickMetadata& metadata=base.bricks[brick];const int i=std::clamp(static_cast<int>(std::floor((point.x-metadata.origin.x)/base.h)),0,system.brick_size-1),j=std::clamp(static_cast<int>(std::floor((point.y-metadata.origin.y)/base.h)),0,system.brick_size-1),k=std::clamp(static_cast<int>(std::floor((point.z-metadata.origin.z)/base.h)),0,system.brick_size-1);return system.dof(0,brick,i,j,k);
		}
	}

	int CompositeAmrPressureSystem::dof(int level, int brick, int i, int j, int k) const
	{
		return level_offset[level] + brick * brick_size * brick_size * brick_size + local_index(brick_size, i, j, k);
	}

	CompositeAmrPressureSystem build_composite_amr_pressure_system(const AmrHierarchy& hierarchy, bool outlet)
	{
		CompositeAmrPressureSystem system; system.hierarchy = &hierarchy; system.brick_size = hierarchy.brick_size(); system.pressure_outlet_xmax = outlet;
		const int bs = system.brick_size, cells_per_brick = bs * bs * bs; system.level_offset.resize(hierarchy.levels().size());
		for (int level = 0; level < static_cast<int>(hierarchy.levels().size()); ++level)
		{
			system.level_offset[level] = system.storage_size; system.storage_size += static_cast<int>(hierarchy.levels()[level].bricks.size()) * cells_per_brick;
		}
		system.active.assign(system.storage_size, 0); system.volume.assign(system.storage_size, 0.0);system.centroid.assign(system.storage_size,{});system.preconditioner_aggregate.assign(system.storage_size,-1);
		system.cut_face_mask.assign(system.storage_size, 0);
		for (int level = 0; level < static_cast<int>(hierarchy.levels().size()); ++level)
		{
			const AmrLevel& source = hierarchy.levels()[level]; const double cell_volume = static_cast<double>(source.h) * source.h * source.h;
			for (int brick = 0; brick < static_cast<int>(source.bricks.size()); ++brick) if (source.bricks[brick].active())
				for (int k = 0; k < bs; ++k) for (int j = 0; j < bs; ++j) for (int i = 0; i < bs; ++i) { const int q = system.dof(level, brick, i, j, k); system.active[q] = 1; system.volume[q] = cell_volume;const Vec3d point=source.bricks[brick].origin+Vec3d{(i+0.5)*source.h,(j+0.5)*source.h,(k+0.5)*source.h};system.centroid[q]=point;system.preconditioner_aggregate[q]=base_aggregate_dof(system,point); }
		}

		// Each active coarse brick scans all six faces for a covered same-level
		// neighbour. Its four fine children per coarse face cell are the unique 2:1
		// flux tiles. Balancing guarantees those interface children are active.
		for (int level = 0; level + 1 < static_cast<int>(hierarchy.levels().size()); ++level)
		{
			const AmrLevel& coarse = hierarchy.levels()[level]; const AmrLevel& fine = hierarchy.levels()[level + 1]; const double hc = coarse.h, hf = fine.h;
			for (int brick = 0; brick < static_cast<int>(coarse.bricks.size()); ++brick)
			{
				const BrickMetadata& owner = coarse.bricks[brick]; if (!owner.active()) continue;
				for (int face = 0; face < 6; ++face)
				{
					const int covered_id = owner.same_level_neighbor[face]; if (covered_id < 0 || coarse.bricks[covered_id].active()) continue; const BrickMetadata& covered = coarse.bricks[covered_id];
					const int axis = face / 2; const bool positive = (face & 1) != 0;
					for (int t1 = 0; t1 < bs; ++t1) for (int t0 = 0; t0 < bs; ++t0) for (int s1 = 0; s1 < 2; ++s1) for (int s0 = 0; s0 < 2; ++s0)
					{
						int ci = t0, cj = t1, ck = 0; if (axis == 0) { ci = positive ? bs - 1 : 0; cj = t0; ck = t1; } else if (axis == 1) { ci = t0; cj = positive ? bs - 1 : 0; ck = t1; } else { ci = t0; cj = t1; ck = positive ? bs - 1 : 0; }
						int parent_fi = 0, parent_fj = 0, parent_fk = 0;
						if (axis == 0) { parent_fi = positive ? 0 : 2 * bs - 1; parent_fj = 2 * t0 + s0; parent_fk = 2 * t1 + s1; }
						else if (axis == 1) { parent_fi = 2 * t0 + s0; parent_fj = positive ? 0 : 2 * bs - 1; parent_fk = 2 * t1 + s1; }
						else { parent_fi = 2 * t0 + s0; parent_fj = 2 * t1 + s1; parent_fk = positive ? 0 : 2 * bs - 1; }
						const int child_x = parent_fi / bs, child_y = parent_fj / bs, child_z = parent_fk / bs; const int child_slot = (child_z * 2 + child_y) * 2 + child_x; const int child = covered.children[child_slot];
						if (child < 0 || !fine.bricks[child].active()) throw std::runtime_error("2:1 interface has missing or further-covered fine child");
						const int fi=parent_fi%bs,fj=parent_fj%bs,fk=parent_fk%bs,face_coordinate=axis==0?fi:(axis==1?fj:fk);Vec3d face_centroid=fine.bricks[child].origin+Vec3d{(fi+0.5)*hf,(fj+0.5)*hf,(fk+0.5)*hf};face_centroid[axis]=fine.bricks[child].origin[axis]+(face_coordinate+(positive?0:1))*hf;system.coarse_fine.push_back({system.dof(level, brick, ci, cj, ck), system.dof(level + 1, child, fi, fj, fk), hf * hf, 0.5 * hc + 0.5 * hf, static_cast<std::int8_t>(axis),static_cast<std::int8_t>(positive?1:-1),face_centroid});
					}
				}
			}
		}
		finalize_component_gauges(system);return system;
	}

	CompositeAmrPressureSystem build_composite_amr_pressure_system(const AmrHierarchy& hierarchy,
		const AmrEmbeddedBoundaryAtlas& atlas, bool outlet)
	{
		CompositeAmrPressureSystem system=build_composite_amr_pressure_system(hierarchy,outlet);system.gauges.clear();
		for(const AmrEbLevelAtlas& level_atlas:atlas.levels)
		{
			if(level_atlas.level<0||level_atlas.level>=static_cast<int>(hierarchy.levels().size()))throw std::invalid_argument("EB atlas level is outside AMR hierarchy");const EmbeddedBoundary& eb=level_atlas.topology;std::vector<int> cell_global(eb.grid.cell_count(),-1),fragment_global(eb.fragments.size(),-1);
			for(int cell=0;cell<eb.grid.cell_count();++cell)if(level_atlas.owned_cell[cell]){const BrickLocation owner=hierarchy.locate_finest(eb.grid.cell_centroid(cell));if(!owner.found()||owner.level!=level_atlas.level)throw std::runtime_error("owned EB atlas cell has no matching AMR owner");const int dof=system.dof(owner.level,owner.brick,owner.cell.x,owner.cell.y,owner.cell.z);cell_global[cell]=dof;if(eb.cells[cell].state==EbCellState::split){system.active[dof]=0;system.volume[dof]=0;}else if(eb.cells[cell].state==EbCellState::regular&&!eb.cut_face_mask.empty())system.cut_face_mask[dof]|=eb.cut_face_mask[cell];}
			// Allocate one appended DOF for every owned root fragment. Volume and centroid
			// are accumulated below from owned fragments only; using the atlas-local EB
			// aggregates here would incorrectly include covered halo fragments.
			for(int fragment=0;fragment<static_cast<int>(eb.fragments.size());++fragment)if(level_atlas.owned_cell[eb.fragments[fragment].parent_cell]&&eb.fragments[fragment].merge_target==irregular_fragment(fragment))
			{
				const int global_dof=system.storage_size++;fragment_global[fragment]=global_dof;system.active.push_back(eb.fragments[fragment].pressure_static?0:1);system.volume.push_back(0);system.centroid.push_back({});system.cut_face_mask.push_back(0);system.preconditioner_aggregate.push_back(-1);
			}
			std::vector<unsigned char> resolving(eb.fragments.size(),0);std::function<int(FragmentRef)> map_ref=[&](FragmentRef ref)->int
			{
				if(ref==invalid_fragment)return -1;if(fragment_is_regular(ref)){const int cell=regular_fragment_cell(ref);return cell>=0&&cell<static_cast<int>(cell_global.size())?cell_global[cell]:-1;}const int fragment=irregular_fragment_index(ref);if(fragment<0||fragment>=static_cast<int>(fragment_global.size()))return -1;if(fragment_global[fragment]>=0)return fragment_global[fragment];if(resolving[fragment])throw std::runtime_error("cyclic AMR EB fragment merge mapping");resolving[fragment]=1;fragment_global[fragment]=map_ref(eb.fragments[fragment].merge_target);resolving[fragment]=0;return fragment_global[fragment];
			};
			for(int fragment=0;fragment<static_cast<int>(eb.fragments.size());++fragment)if(level_atlas.owned_cell[eb.fragments[fragment].parent_cell]&&map_ref(irregular_fragment(fragment))<0)throw std::runtime_error("owned AMR EB fragment maps outside its level atlas");
			// Conservative merge transfer. A tiny fragment may resolve to an appended
			// irregular root or to an ordinary structured cell. In both cases its volume
			// and first moment belong to that same-side control volume in the composite
			// projection; omitting this transfer changes divergence normalization.
			for(int fragment=0;fragment<static_cast<int>(eb.fragments.size());++fragment)if(level_atlas.owned_cell[eb.fragments[fragment].parent_cell])
			{
				const int target=map_ref(irregular_fragment(fragment));const double add=eb.fragments[fragment].volume;
				if(!(add>=0)||target<0)throw std::runtime_error("invalid owned AMR EB merge volume");const double old=system.volume[target];
				if(old+add>0)system.centroid[target]=(system.centroid[target]*old+eb.fragments[fragment].centroid*add)/(old+add);system.volume[target]=old+add;
			}
			for(int fragment=0;fragment<static_cast<int>(eb.fragments.size());++fragment)if(level_atlas.owned_cell[eb.fragments[fragment].parent_cell])
			{
				const int target=map_ref(irregular_fragment(fragment));system.preconditioner_aggregate[target]=base_aggregate_dof(system,system.centroid[target]);
			}
			for(const FaceAperture& aperture:eb.apertures)
			{
				const int a=map_ref(aperture.fragment_a),b=map_ref(aperture.fragment_b);if(a<0||b<0||a==b)continue;const bool a_active=system.active[a]!=0,b_active=system.active[b]!=0;if(!a_active&&!b_active)continue;const double distance=std::max(1e-12,std::sqrt(length2(system.centroid[a]-system.centroid[b])));system.embedded.push_back({a,b,aperture.area,distance,aperture.axis,1,aperture.centroid});
			}
			for(const SurfacePatch& patch:eb.patches)
			{
				const BrickLocation owner=hierarchy.locate_finest(patch.centroid);if(!owner.found()||owner.level!=level_atlas.level)continue;const int plus=map_ref(patch.plus_fragment),minus=map_ref(patch.minus_fragment);if(plus<0||minus<0)throw std::runtime_error("owned AMR surface patch has no two-sided pressure mapping");system.surface_patches.push_back({patch.source_triangle_id,patch.source_face_id,patch.area,patch.centroid,patch.normal,plus,minus});
			}
		}
		for(const CoarseFinePressureConnection& connection:system.coarse_fine)if(!system.active[connection.coarse_dof]||!system.active[connection.fine_dof])throw std::runtime_error("embedded boundary intersects a 2:1 interface; aperture-aware cross-level topology is required");
		finalize_component_gauges(system);return system;
	}

	CompositeAmrFluxes make_zero_composite_fluxes(const CompositeAmrPressureSystem& system){CompositeAmrFluxes flux;flux.coarse_fine_velocity.assign(system.coarse_fine.size(),0);flux.embedded_velocity.assign(system.embedded.size(),0);return flux;}

	namespace
	{
		double regular_face_value(const AmrHostLevelFields& level,int brick,int axis,int i,int j,int k)
		{
			if(axis==0)return level.u[level.layout.u_index(brick,i+1,j,k)];if(axis==1)return level.v[level.layout.v_index(brick,i,j+1,k)];return level.w[level.layout.w_index(brick,i,j,k+1)];
		}
		void set_regular_face(AmrHostLevelFields& level,int brick,int axis,int i,int j,int k,double value)
		{
			if(axis==0)level.u[level.layout.u_index(brick,i+1,j,k)]=static_cast<Real>(value);else if(axis==1)level.v[level.layout.v_index(brick,i,j+1,k)]=static_cast<Real>(value);else level.w[level.layout.w_index(brick,i,j,k+1)]=static_cast<Real>(value);
		}
		void add_special_flux(const CoarseFinePressureConnection& connection,double velocity,std::vector<double>& integrated)
		{
			const int lower=connection.direction>0?connection.coarse_dof:connection.fine_dof,upper=connection.direction>0?connection.fine_dof:connection.coarse_dof;const double flux=velocity*connection.open_area;integrated[lower]+=flux;integrated[upper]-=flux;
		}
	}

	void composite_amr_divergence_cpu(const CompositeAmrPressureSystem& system,const AmrHostFields& fields,const CompositeAmrFluxes& special,std::vector<double>& divergence)
	{
		if(special.coarse_fine_velocity.size()!=system.coarse_fine.size()||special.embedded_velocity.size()!=system.embedded.size())throw std::invalid_argument("composite AMR special flux size");std::vector<double> integrated(system.storage_size,0);const int bs=system.brick_size;
		for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level){const AmrLevel& metadata=system.hierarchy->levels()[level];const AmrHostLevelFields& values=fields.levels()[level];const double area=metadata.h*metadata.h;for(int brick=0;brick<static_cast<int>(metadata.bricks.size());++brick){const BrickMetadata& meta=metadata.bricks[brick];if(!meta.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=system.dof(level,brick,i,j,k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if(system.cut_face_mask[a]&(1u<<axis))continue;const int local[3]={i,j,k};int coordinate[3]={i,j,k},b=-1;if(local[axis]+1<bs){++coordinate[axis];b=system.dof(level,brick,coordinate[0],coordinate[1],coordinate[2]);}else{const int neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&metadata.bricks[neighbour].active()){coordinate[axis]=0;b=system.dof(level,neighbour,coordinate[0],coordinate[1],coordinate[2]);}}if(b>=0&&system.active[b]){const double flux=regular_face_value(values,brick,axis,i,j,k)*area;integrated[a]+=flux;integrated[b]-=flux;}}if((meta.flags&BRICK_XMIN)&&i==0)integrated[a]-=values.u[values.layout.u_index(brick,0,j,k)]*area;if((meta.flags&BRICK_XMAX)&&i==bs-1)integrated[a]+=values.u[values.layout.u_index(brick,bs,j,k)]*area;}}}
		for(std::size_t edge=0;edge<system.coarse_fine.size();++edge)add_special_flux(system.coarse_fine[edge],special.coarse_fine_velocity[edge],integrated);for(std::size_t edge=0;edge<system.embedded.size();++edge)add_special_flux(system.embedded[edge],special.embedded_velocity[edge],integrated);divergence.assign(system.storage_size,0);for(int q=0;q<system.storage_size;++q)if(system.active[q]&&system.volume[q]>0)divergence[q]=integrated[q]/system.volume[q];
	}

	void composite_amr_projection_rhs_cpu(const CompositeAmrPressureSystem& system,const AmrHostFields& fields,const CompositeAmrFluxes& special,double rho,double dt,std::vector<double>& rhs)
	{
		std::vector<double> divergence;composite_amr_divergence_cpu(system,fields,special,divergence);rhs.assign(system.storage_size,0);for(int q=0;q<system.storage_size;++q)if(system.active[q])rhs[q]=-(rho/dt)*divergence[q]*system.volume[q];
	}

	void composite_amr_correct_fluxes_cpu(const CompositeAmrPressureSystem& system,const std::vector<double>& pressure,double rho,double dt,AmrHostFields& fields,CompositeAmrFluxes& special)
	{
		if(pressure.size()!=static_cast<std::size_t>(system.storage_size))throw std::invalid_argument("composite AMR correction pressure size");const int bs=system.brick_size;
		for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level){const AmrLevel& metadata=system.hierarchy->levels()[level];AmrHostLevelFields& values=fields.levels()[level];for(int brick=0;brick<static_cast<int>(metadata.bricks.size());++brick){const BrickMetadata& meta=metadata.bricks[brick];if(!meta.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=system.dof(level,brick,i,j,k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if(system.cut_face_mask[a]&(1u<<axis))continue;int coordinate[3]={i,j,k},b=-1,neighbour=-1;if(coordinate[axis]+1<bs){++coordinate[axis];b=system.dof(level,brick,coordinate[0],coordinate[1],coordinate[2]);}else{neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&metadata.bricks[neighbour].active()){coordinate[axis]=0;b=system.dof(level,neighbour,coordinate[0],coordinate[1],coordinate[2]);}}if(b>=0&&system.active[b]){double velocity=regular_face_value(values,brick,axis,i,j,k)-(dt/rho)*(pressure[b]-pressure[a])/metadata.h;set_regular_face(values,brick,axis,i,j,k,velocity);if(neighbour>=0){int opposite[3]={i,j,k};opposite[axis]=-1;set_regular_face(values,neighbour,axis,opposite[0],opposite[1],opposite[2],velocity);}}}if((meta.flags&BRICK_XMAX)&&i==bs-1&&system.pressure_outlet_xmax){double velocity=values.u[values.layout.u_index(brick,bs,j,k)]+(dt/rho)*pressure[a]/(0.5*metadata.h);values.u[values.layout.u_index(brick,bs,j,k)]=static_cast<Real>(velocity);}}}}
		auto correct_special=[&](const std::vector<CoarseFinePressureConnection>& connections,std::vector<double>& velocity){for(std::size_t edge=0;edge<connections.size();++edge){const auto& connection=connections[edge];const int lower=connection.direction>0?connection.coarse_dof:connection.fine_dof,upper=connection.direction>0?connection.fine_dof:connection.coarse_dof;if(system.active[lower]&&system.active[upper])velocity[edge]-=(dt/rho)*(pressure[upper]-pressure[lower])/connection.centre_distance;}};correct_special(system.coarse_fine,special.coarse_fine_velocity);correct_special(system.embedded,special.embedded_velocity);
	}

	void CompositeAmrPressureSystem::apply_cpu(const std::vector<double>& pressure, std::vector<double>& output) const
	{
		if (pressure.size() != static_cast<std::size_t>(storage_size)) throw std::invalid_argument("composite AMR pressure vector size"); output.assign(storage_size, 0.0); const int bs = brick_size;
		for (int level = 0; level < static_cast<int>(hierarchy->levels().size()); ++level)
		{
			const AmrLevel& source = hierarchy->levels()[level]; const double coefficient = source.h;
			for (int brick = 0; brick < static_cast<int>(source.bricks.size()); ++brick)
			{
				const BrickMetadata& meta = source.bricks[brick]; if (!meta.active()) continue;
				for (int k = 0; k < bs; ++k) for (int j = 0; j < bs; ++j) for (int i = 0; i < bs; ++i)
				{
					const int a = dof(level, brick, i, j, k);if(!active[a])continue; const int local[3] = {i,j,k};
					for (int axis = 0; axis < 3; ++axis)
					{
						if(cut_face_mask[a]&(1u<<axis))continue;if (local[axis] + 1 < bs) { int q[3] = {i,j,k}; ++q[axis];const int b=dof(level,brick,q[0],q[1],q[2]);if(active[b])add_edge(output,pressure,a,b,coefficient); }
						else { const int neighbour = meta.same_level_neighbor[2 * axis + 1]; if (neighbour >= 0 && source.bricks[neighbour].active()) { int q[3] = {i,j,k}; q[axis] = 0;const int b=dof(level,neighbour,q[0],q[1],q[2]);if(active[b])add_edge(output,pressure,a,b,coefficient); } }
					}
					if (pressure_outlet_xmax && (meta.flags & BRICK_XMAX) && i == bs - 1) output[a] += 2.0 * source.h * pressure[a];
				}
			}
		}
		for (const CoarseFinePressureConnection& edge : coarse_fine)if(active[edge.coarse_dof]&&active[edge.fine_dof])add_edge(output, pressure, edge.coarse_dof, edge.fine_dof, edge.open_area / edge.centre_distance);
		for (const CoarseFinePressureConnection& edge : embedded)if(active[edge.coarse_dof]&&active[edge.fine_dof])add_edge(output, pressure, edge.coarse_dof, edge.fine_dof, edge.open_area / edge.centre_distance);
		for(const CompositePressureGauge& gauge:gauges)if(gauge.dof>=0&&active[gauge.dof])output[gauge.dof]+=gauge.coefficient*pressure[gauge.dof];
	}

	void CompositeAmrPressureSystem::diagonal_cpu(std::vector<double>& diagonal) const
	{
		diagonal.assign(storage_size, 0.0); const int bs = brick_size;auto add_diagonal=[&](int a,int b,double coefficient){if(a<0||b<0||a==b||!active[a]||!active[b])return;diagonal[a]+=coefficient;diagonal[b]+=coefficient;};
		for (int level = 0; level < static_cast<int>(hierarchy->levels().size()); ++level)
		{
			const AmrLevel& source = hierarchy->levels()[level]; const double coefficient = source.h;
			for (int brick = 0; brick < static_cast<int>(source.bricks.size()); ++brick) { const BrickMetadata& meta = source.bricks[brick]; if (!meta.active()) continue; for (int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i) { const int a=dof(level,brick,i,j,k);if(!active[a])continue;for(int axis=0;axis<3;++axis){if(cut_face_mask[a]&(1u<<axis))continue;int c[3]={i,j,k};++c[axis];if(c[axis]<bs)add_diagonal(a,dof(level,brick,c[0],c[1],c[2]),coefficient);else{int neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&source.bricks[neighbour].active()){c[axis]=0;add_diagonal(a,dof(level,neighbour,c[0],c[1],c[2]),coefficient);}}}if(pressure_outlet_xmax&&(meta.flags&BRICK_XMAX)&&i==bs-1)diagonal[a]+=2.0*source.h;} }
		}
		for(const CoarseFinePressureConnection& edge:coarse_fine)add_diagonal(edge.coarse_dof,edge.fine_dof,edge.open_area/edge.centre_distance);for(const CoarseFinePressureConnection& edge:embedded)add_diagonal(edge.coarse_dof,edge.fine_dof,edge.open_area/edge.centre_distance);
		for(const CompositePressureGauge& gauge:gauges)if(gauge.dof>=0&&active[gauge.dof])diagonal[gauge.dof]+=gauge.coefficient;
	}
}
