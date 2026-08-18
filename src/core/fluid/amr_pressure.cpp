#include "core/fluid/amr_pressure.h"
#include "core/fluid/amr_eb.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

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
			system.gauges.clear();const int n=system.storage_size,bs=system.brick_size;system.freestream_connected.assign(n,0);if(!system.pressure_outlet_xmax){for(int q=0;q<n;++q)system.freestream_connected[q]=system.active[q];return;}std::vector<int> parent(n,-1);for(int q=0;q<n;++q)if(system.active[q])parent[q]=q;auto root=[&](int q){while(parent[q]!=q){parent[q]=parent[parent[q]];q=parent[q];}return q;};auto join=[&](int a,int b){if(a<0||b<0||!system.active[a]||!system.active[b])return;a=root(a);b=root(b);if(a!=b)parent[std::max(a,b)]=std::min(a,b);};
			for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level){const AmrLevel& source=system.hierarchy->levels()[level];for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick){const BrickMetadata& meta=source.bricks[brick];if(!meta.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=system.dof(level,brick,i,j,k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if(system.cut_face_mask[a]&(1u<<axis))continue;int c[3]={i,j,k},b=-1;if(++c[axis]<bs)b=system.dof(level,brick,c[0],c[1],c[2]);else{const int neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&source.bricks[neighbour].active()){c[axis]=0;b=system.dof(level,neighbour,c[0],c[1],c[2]);}}join(a,b);}}}}
			for(const CoarseFinePressureConnection& edge:system.coarse_fine)join(edge.coarse_dof,edge.fine_dof);for(const CoarseFinePressureConnection& edge:system.embedded)join(edge.coarse_dof,edge.fine_dof);
			std::vector<unsigned char> outlet_root(n,0),external_root(n,0);for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level){const AmrLevel& source=system.hierarchy->levels()[level];for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick){const BrickMetadata& meta=source.bricks[brick];if(!meta.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j){if(meta.flags&BRICK_XMIN){const int q=system.dof(level,brick,0,j,k);if(system.active[q])external_root[root(q)]=1;}if(meta.flags&BRICK_XMAX){const int q=system.dof(level,brick,bs-1,j,k);if(system.active[q])outlet_root[root(q)]=external_root[root(q)]=1;}}}}
			std::vector<unsigned char> emitted(n,0);for(int q=0;q<n;++q)if(system.active[q]){const int r=root(q);system.freestream_connected[q]=external_root[r];if(!outlet_root[r]&&!emitted[r]){emitted[r]=1;system.gauges.push_back({q,std::max(1e-12,std::cbrt(std::max(0.0,system.volume[q])))});}}
		}
		int base_aggregate_dof(const CompositeAmrPressureSystem& system,Vec3d point)
		{
			const AmrLevel& base=system.hierarchy->levels().front();const double brick_width=system.brick_size*base.h;const Aabb3d& domain=system.hierarchy->domain();Int3 coord{static_cast<int>(std::floor((point.x-domain.lo.x)/brick_width)),static_cast<int>(std::floor((point.y-domain.lo.y)/brick_width)),static_cast<int>(std::floor((point.z-domain.lo.z)/brick_width))};const int brick=system.hierarchy->find_brick(0,coord);if(brick<0)return -1;const BrickMetadata& metadata=base.bricks[brick];const int i=std::clamp(static_cast<int>(std::floor((point.x-metadata.origin.x)/base.h)),0,system.brick_size-1),j=std::clamp(static_cast<int>(std::floor((point.y-metadata.origin.y)/base.h)),0,system.brick_size-1),k=std::clamp(static_cast<int>(std::floor((point.z-metadata.origin.z)/base.h)),0,system.brick_size-1);return system.dof(0,brick,i,j,k);
		}

		int symmetric_pseudoinverse_3x3(const double* input,double* inverse)
		{
			double a[9],vectors[9]={1,0,0,0,1,0,0,0,1};for(int q=0;q<9;++q)a[q]=input[q];
			for(int sweep=0;sweep<20;++sweep)
			{
				bool changed=false;for(int pair=0;pair<3;++pair)
				{
					const int p=pair==0?0:(pair==1?0:1),q=pair==0?1:(pair==1?2:2);const double apq=a[p*3+q],scale=std::max({std::abs(a[p*3+p]),std::abs(a[q*3+q]),1e-300});if(std::abs(apq)<=1e-14*scale)continue;changed=true;const double phi=0.5*std::atan2(2*apq,a[q*3+q]-a[p*3+p]),c=std::cos(phi),s=std::sin(phi),app=a[p*3+p],aqq=a[q*3+q];for(int k=0;k<3;++k)if(k!=p&&k!=q){const double akp=a[k*3+p],akq=a[k*3+q];a[k*3+p]=a[p*3+k]=c*akp-s*akq;a[k*3+q]=a[q*3+k]=s*akp+c*akq;}a[p*3+p]=c*c*app-2*s*c*apq+s*s*aqq;a[q*3+q]=s*s*app+2*s*c*apq+c*c*aqq;a[p*3+q]=a[q*3+p]=0;for(int k=0;k<3;++k){const double vkp=vectors[k*3+p],vkq=vectors[k*3+q];vectors[k*3+p]=c*vkp-s*vkq;vectors[k*3+q]=s*vkp+c*vkq;}
				}
				if(!changed)break;
			}
			for(int q=0;q<9;++q)inverse[q]=0;const double maximum=std::max({std::abs(a[0]),std::abs(a[4]),std::abs(a[8])}),threshold=maximum*(sizeof(Real)==4?1e-4:1e-10);int rank=0;for(int mode=0;mode<3;++mode){const double eigenvalue=a[mode*3+mode];if(!(eigenvalue>threshold))continue;++rank;for(int row=0;row<3;++row)for(int column=0;column<3;++column)inverse[row*3+column]+=vectors[row*3+mode]*vectors[column*3+mode]/eigenvalue;}return rank;
		}

		struct PressureGradientSample { int neighbour=-1; double weight=0; };

		void finalize_nonorthogonal_pressure_topology(CompositeAmrPressureSystem& system,
			UnsupportedNonorthogonalCorrectionPolicy unsupported_policy)
		{
			if(unsupported_policy!=UnsupportedNonorthogonalCorrectionPolicy::reject&&
				unsupported_policy!=UnsupportedNonorthogonalCorrectionPolicy::qualitative_preview_orthogonal_fallback&&
				unsupported_policy!=UnsupportedNonorthogonalCorrectionPolicy::qualitative_preview_first_order_orthogonal)
				throw std::invalid_argument("unknown unsupported nonorthogonal pressure-correction policy");
			system.regular_pressure_corrections.clear();system.pressure_gradient_dof.clear();system.pressure_gradient_offset.clear();system.pressure_gradient_neighbour.clear();system.pressure_gradient_weight.clear();system.pressure_gradient_rank.clear();system.pressure_gradient_ring.clear();system.numerically_orthogonal_regular=system.numerically_orthogonal_coarse_fine=system.numerically_orthogonal_embedded=0;system.maximum_numerically_orthogonal_correction=0;
			system.qualitative_preview_orthogonal_fallback_regular=0;
			system.qualitative_preview_orthogonal_fallback_coarse_fine=0;
			system.qualitative_preview_orthogonal_fallback_embedded=0;
			auto prepare=[](CompositeAmrPressureSystem& owner,CoarseFinePressureConnection& edge,std::size_t& below_resolution)
			{
				if(edge.coarse_dof<0||edge.fine_dof<0||edge.coarse_dof>=owner.storage_size||edge.fine_dof>=owner.storage_size||edge.coarse_dof==edge.fine_dof||!(edge.open_area>0))throw std::runtime_error("invalid compact pressure connection");
				const int lower=edge.direction>0?edge.coarse_dof:edge.fine_dof,upper=edge.direction>0?edge.fine_dof:edge.coarse_dof;const Vec3d displacement=owner.centroid[upper]-owner.centroid[lower];const double distance2=length2(displacement),signed_normal_distance=displacement[edge.axis];if(!(distance2>0))throw std::runtime_error("invalid compact pressure centroid displacement");const double distance=std::sqrt(distance2),bracket_epsilon=1e-10*std::max(distance,1e-12),lower_face_distance=edge.face_centroid[edge.axis]-owner.centroid[lower][edge.axis],upper_face_distance=owner.centroid[upper][edge.axis]-edge.face_centroid[edge.axis];if(!(signed_normal_distance>2*bracket_epsilon&&lower_face_distance>bracket_epsilon&&upper_face_distance>bracket_epsilon))throw std::runtime_error("compact pressure face is not signed-bracketed by its control-volume centroids");const double normal_path=lower_face_distance+upper_face_distance;edge.centre_distance=distance;edge.normal_distance=signed_normal_distance;const double over_relaxed_alpha=1.0/signed_normal_distance;edge.orthogonal_gradient_factor=over_relaxed_alpha;Vec3d normal{};normal[edge.axis]=1;edge.nonorthogonal_correction=normal-displacement*edge.orthogonal_gradient_factor;if(nonorthogonal_correction_is_below_resolution(edge.nonorthogonal_correction)){owner.maximum_numerically_orthogonal_correction=std::max(owner.maximum_numerically_orthogonal_correction,std::sqrt(length2(edge.nonorthogonal_correction)));edge.nonorthogonal_correction={};++below_resolution;}edge.upper_gradient_weight=lower_face_distance/normal_path;edge.lower_gradient_node=edge.upper_gradient_node=-1;
			};
			for(auto& edge:system.coarse_fine)prepare(system,edge,system.numerically_orthogonal_coarse_fine);for(auto& edge:system.embedded)prepare(system,edge,system.numerically_orthogonal_embedded);
			const int bs=system.brick_size;
			for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level)
			{
				const AmrLevel& source=system.hierarchy->levels()[level];
				const double area=source.h*source.h,base_factor=1.0/source.h;
				for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick)
				{
					const BrickMetadata& meta=source.bricks[brick];if(!meta.active())continue;
					for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i)
					{
						const int lower=system.dof(level,brick,i,j,k);if(!system.active[lower])continue;
						for(int axis=0;axis<3;++axis)
						{
							if(system.cut_face_mask[lower]&(1u<<axis))continue;
							int coordinate[3]={i,j,k},upper_brick=brick;++coordinate[axis];
							if(coordinate[axis]>=bs){upper_brick=meta.same_level_neighbor[2*axis+1];if(upper_brick<0||!source.bricks[upper_brick].active())continue;coordinate[axis]=0;}
							const int upper=system.dof(level,upper_brick,coordinate[0],coordinate[1],coordinate[2]);if(!system.active[upper])continue;
							const Vec3d displacement=system.centroid[upper]-system.centroid[lower];const double distance2=length2(displacement);if(!(distance2>1e-24))throw std::runtime_error("implicit pressure face has coincident control-volume centroids");
							const double over_relaxed_alpha=1.0/displacement[axis];Vec3d face{meta.origin.x+(i+0.5)*source.h,meta.origin.y+(j+0.5)*source.h,meta.origin.z+(k+0.5)*source.h};face[axis]+=0.5*source.h;const double bracket_epsilon=1e-10*source.h,lower_distance=face[axis]-system.centroid[lower][axis],upper_distance=system.centroid[upper][axis]-face[axis];if(!(displacement[axis]>2*bracket_epsilon&&lower_distance>bracket_epsilon&&upper_distance>bracket_epsilon))throw std::runtime_error("implicit pressure face is not signed-bracketed by its control-volume centroids");const double normal_path=lower_distance+upper_distance,alpha=over_relaxed_alpha;Vec3d normal{};normal[axis]=1;Vec3d correction=normal-displacement*alpha;if(nonorthogonal_correction_is_below_resolution(correction)){system.maximum_numerically_orthogonal_correction=std::max(system.maximum_numerically_orthogonal_correction,std::sqrt(length2(correction)));correction={};++system.numerically_orthogonal_regular;}const double two_point_delta=alpha-base_factor;if(std::abs(two_point_delta)*source.h<=1e-12&&length2(correction)<=1e-24)continue;
							RegularPressureCorrection record;record.lower_dof=lower;record.upper_dof=upper;record.level=level;record.brick=brick;record.i=i;record.j=j;record.k=k;if(axis==0)++record.i;else if(axis==1)++record.j;else ++record.k;record.axis=static_cast<std::int8_t>(axis);record.open_area=area;record.two_point_delta=two_point_delta;record.nonorthogonal_correction=correction;record.upper_gradient_weight=lower_distance/normal_path;system.regular_pressure_corrections.push_back(record);
						}
					}
				}
			}
			if(unsupported_policy==UnsupportedNonorthogonalCorrectionPolicy::qualitative_preview_first_order_orthogonal)
			{
				auto drop_connection=[](CoarseFinePressureConnection& edge,std::size_t& count)
				{
					if(length2(edge.nonorthogonal_correction)>1e-24)++count;
					edge.nonorthogonal_correction={};edge.lower_gradient_node=edge.upper_gradient_node=-1;
				};
				for(auto& edge:system.coarse_fine)
					drop_connection(edge,system.qualitative_preview_orthogonal_fallback_coarse_fine);
				for(auto& edge:system.embedded)
					drop_connection(edge,system.qualitative_preview_orthogonal_fallback_embedded);
				for(auto& edge:system.regular_pressure_corrections)
				{
					if(length2(edge.nonorthogonal_correction)>1e-24)
						++system.qualitative_preview_orthogonal_fallback_regular;
					edge.nonorthogonal_correction={};edge.lower_gradient_node=edge.upper_gradient_node=-1;
				}
				// A skew same-level face whose signed-normal distance is still exactly h
				// now needs neither a two-point delta nor a deferred term. Removing that
				// empty record keeps the first-order preview's compact work list minimal.
				system.regular_pressure_corrections.erase(std::remove_if(
					system.regular_pressure_corrections.begin(),system.regular_pressure_corrections.end(),
					[&](const RegularPressureCorrection& edge)
					{
						const double h=system.hierarchy->levels()[edge.level].h;
						return std::abs(edge.two_point_delta)*h<=1e-12;
					}),system.regular_pressure_corrections.end());
				// No deferred term survives, so no WLS nodes, graph rings, response matrix,
				// or GPU gradient buffers are required for this explicitly first-order path.
				return;
			}

			std::unordered_map<int,int> node_map;auto register_connection=[&](const CoarseFinePressureConnection& edge)
			{
				if(length2(edge.nonorthogonal_correction)<=1e-24)return;const int lower=edge.direction>0?edge.coarse_dof:edge.fine_dof,upper=edge.direction>0?edge.fine_dof:edge.coarse_dof;for(int dof:{lower,upper})if(node_map.emplace(dof,static_cast<int>(node_map.size())).second)system.pressure_gradient_dof.push_back(dof);
			};
			for(const auto& edge:system.coarse_fine)register_connection(edge);for(const auto& edge:system.embedded)register_connection(edge);for(const auto& edge:system.regular_pressure_corrections)for(int dof:{edge.lower_dof,edge.upper_dof})if(node_map.emplace(dof,static_cast<int>(node_map.size())).second)system.pressure_gradient_dof.push_back(dof);
			std::vector<std::vector<PressureGradientSample>> samples(system.pressure_gradient_dof.size());auto add_sample=[&](int endpoint,int neighbour,double area)
			{
				auto found=node_map.find(endpoint);if(found==node_map.end()||neighbour<0||neighbour>=system.storage_size||!system.active[neighbour]||endpoint==neighbour)return;const Vec3d displacement=system.centroid[neighbour]-system.centroid[endpoint];const double distance2=length2(displacement);if(distance2>1e-24&&area>0)samples[found->second].push_back({neighbour,area/distance2});
			};
			for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level)
			{
				const AmrLevel& source=system.hierarchy->levels()[level];const double area=source.h*source.h;for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick){const BrickMetadata& meta=source.bricks[brick];if(!meta.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=system.dof(level,brick,i,j,k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if(system.cut_face_mask[a]&(1u<<axis))continue;int c[3]={i,j,k},b=-1;if(++c[axis]<bs)b=system.dof(level,brick,c[0],c[1],c[2]);else{const int neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&source.bricks[neighbour].active()){c[axis]=0;b=system.dof(level,neighbour,c[0],c[1],c[2]);}}if(b>=0&&system.active[b]){add_sample(a,b,area);add_sample(b,a,area);}}}}
			}
			auto add_special_samples=[&](const CoarseFinePressureConnection& edge){add_sample(edge.coarse_dof,edge.fine_dof,edge.open_area);add_sample(edge.fine_dof,edge.coarse_dof,edge.open_area);};for(const auto& edge:system.coarse_fine)add_special_samples(edge);for(const auto& edge:system.embedded)add_special_samples(edge);

			// A junction or thin local fragment can have a rank-deficient immediate WLS
			// neighbourhood even though its same-fluid graph is fully three-dimensional.
			// Extend only deficient correction nodes, one graph ring at a time. Every path
			// is composed of real open connections, so the stencil cannot cross fabric.
			auto sample_rank=[&](int node,double* inverse)
			{
				const int centre=system.pressure_gradient_dof[node];double normal_matrix[9]{};for(const auto& sample:samples[node]){const Vec3d d=system.centroid[sample.neighbour]-system.centroid[centre];for(int row=0;row<3;++row)for(int column=0;column<3;++column)normal_matrix[row*3+column]+=sample.weight*d[row]*d[column];}double local_inverse[9];return symmetric_pseudoinverse_3x3(normal_matrix,inverse?inverse:local_inverse);
			};
			std::vector<std::uint8_t> ring_used(samples.size(),1);
			auto extend_one_ring=[&](const std::vector<int>& deficient)
			{
				struct RingSource{int node=-1;double path_area=0;};std::unordered_map<int,std::vector<RingSource>> sources;
				for(const int node:deficient)
				{
					const int centre=system.pressure_gradient_dof[node];for(const auto& sample:samples[node]){const Vec3d d=system.centroid[sample.neighbour]-system.centroid[centre];sources[sample.neighbour].push_back({node,sample.weight*length2(d)});}
				}
				std::vector<std::vector<PressureGradientSample>> additions(samples.size());
				auto propose=[&](int middle,int other,double area)
				{
					const auto found=sources.find(middle);if(found==sources.end())return;for(const RingSource& source:found->second)
					{
						const int centre=system.pressure_gradient_dof[source.node];if(other==centre)continue;bool present=false;for(const auto& sample:samples[source.node])if(sample.neighbour==other){present=true;break;}if(present)continue;const double distance2=length2(system.centroid[other]-system.centroid[centre]);if(!(distance2>1e-24))continue;const double weight=std::min(source.path_area,area)/distance2;auto& candidates=additions[source.node];auto duplicate=std::find_if(candidates.begin(),candidates.end(),[&](const PressureGradientSample& value){return value.neighbour==other;});if(duplicate==candidates.end())candidates.push_back({other,weight});else duplicate->weight=std::max(duplicate->weight,weight);
					}
				};
				for(int level=0;level<static_cast<int>(system.hierarchy->levels().size());++level){const AmrLevel& source=system.hierarchy->levels()[level];const double area=source.h*source.h;for(int brick=0;brick<static_cast<int>(source.bricks.size());++brick){const BrickMetadata& meta=source.bricks[brick];if(!meta.active())continue;for(int k=0;k<bs;++k)for(int j=0;j<bs;++j)for(int i=0;i<bs;++i){const int a=system.dof(level,brick,i,j,k);if(!system.active[a])continue;for(int axis=0;axis<3;++axis){if(system.cut_face_mask[a]&(1u<<axis))continue;int c[3]={i,j,k},b=-1;if(++c[axis]<bs)b=system.dof(level,brick,c[0],c[1],c[2]);else{const int neighbour=meta.same_level_neighbor[2*axis+1];if(neighbour>=0&&source.bricks[neighbour].active()){c[axis]=0;b=system.dof(level,neighbour,c[0],c[1],c[2]);}}if(b>=0&&system.active[b]){propose(a,b,area);propose(b,a,area);}}}}}
				auto propose_special=[&](const CoarseFinePressureConnection& edge){propose(edge.coarse_dof,edge.fine_dof,edge.open_area);propose(edge.fine_dof,edge.coarse_dof,edge.open_area);};for(const auto& edge:system.coarse_fine)propose_special(edge);for(const auto& edge:system.embedded)propose_special(edge);
				for(const int node:deficient){auto& candidates=additions[node];std::sort(candidates.begin(),candidates.end(),[](const PressureGradientSample& a,const PressureGradientSample& b){return a.neighbour<b.neighbour;});samples[node].insert(samples[node].end(),candidates.begin(),candidates.end());}
			};
			std::vector<int> deficient;for(int node=0;node<static_cast<int>(samples.size());++node)if(sample_rank(node,nullptr)<3)deficient.push_back(node);
			// Three same-fluid graph rings resolve normal smooth-surface cells. Closed
			// all-irregular components remain active; exactly orthogonal narrow channels
			// need no WLS correction, while skew ones use the same bounded graph search.
			for(int ring=2;ring<=3&&!deficient.empty();++ring)
			{
				for(const int node:deficient)ring_used[node]=static_cast<std::uint8_t>(ring);
				extend_one_ring(deficient);
				if(ring<3){std::vector<int> remaining;for(const int node:deficient)if(sample_rank(node,nullptr)<3)remaining.push_back(node);deficient=std::move(remaining);}
			}
			// A rare active sheet-edge stencil can remain planar after three shells. One
			// final bounded extension is allowed only for rank-2 nodes; rank-0/1 topology
			// is not repaired by reaching farther and remains subject to the hard invariant.
			deficient.clear();for(int node=0;node<static_cast<int>(samples.size());++node)if(sample_rank(node,nullptr)==2){deficient.push_back(node);ring_used[node]=4;}extend_one_ring(deficient);

			system.pressure_gradient_offset.push_back(0);for(int node=0;node<static_cast<int>(system.pressure_gradient_dof.size());++node)
			{
				const int centre=system.pressure_gradient_dof[node];double inverse[9];system.pressure_gradient_rank.push_back(static_cast<std::uint8_t>(sample_rank(node,inverse)));system.pressure_gradient_ring.push_back(ring_used[node]);for(const auto& sample:samples[node]){const Vec3d d=system.centroid[sample.neighbour]-system.centroid[centre];Vec3d stencil{};for(int row=0;row<3;++row)for(int column=0;column<3;++column)stencil[row]+=inverse[row*3+column]*sample.weight*d[column];system.pressure_gradient_neighbour.push_back(sample.neighbour);system.pressure_gradient_weight.push_back(stencil);}system.pressure_gradient_offset.push_back(static_cast<int>(system.pressure_gradient_neighbour.size()));
			}
			auto bind_nodes=[&](CoarseFinePressureConnection& edge){if(length2(edge.nonorthogonal_correction)<=1e-24)return;const int lower=edge.direction>0?edge.coarse_dof:edge.fine_dof,upper=edge.direction>0?edge.fine_dof:edge.coarse_dof;edge.lower_gradient_node=node_map.at(lower);edge.upper_gradient_node=node_map.at(upper);};for(auto& edge:system.coarse_fine)bind_nodes(edge);for(auto& edge:system.embedded)bind_nodes(edge);for(auto& edge:system.regular_pressure_corrections){edge.lower_gradient_node=node_map.at(edge.lower_dof);edge.upper_gradient_node=node_map.at(edge.upper_dof);}

			// An unsupported correction silently drops part of the pressure derivative and
			// can create a non-conservative feedback loop. Audit the final affine response
			// here, once during preprocessing, and fail with the exact offending topology.
			std::vector<std::array<double,9>> response(system.pressure_gradient_dof.size());for(std::size_t node=0;node<response.size();++node){const int centre=system.pressure_gradient_dof[node];for(int q=system.pressure_gradient_offset[node];q<system.pressure_gradient_offset[node+1];++q){const Vec3d d=system.centroid[system.pressure_gradient_neighbour[q]]-system.centroid[centre],weight=system.pressure_gradient_weight[q];for(int row=0;row<3;++row)for(int column=0;column<3;++column)response[node][3*row+column]+=weight[row]*d[column];}}
			auto reject_or_drop=[&](std::string message,Vec3d& correction,int& lower_node,
				int& upper_node,std::size_t& fallback_count)
			{
				if(unsupported_policy==UnsupportedNonorthogonalCorrectionPolicy::reject)
					throw std::runtime_error(message);
				if(unsupported_policy!=UnsupportedNonorthogonalCorrectionPolicy::qualitative_preview_orthogonal_fallback)
					throw std::invalid_argument("unknown unsupported nonorthogonal pressure-correction policy");
				// Preserve the implicit conservative A/d term selected by prepare(). Only
				// the unsupported deferred correction is removed from this preview edge.
				correction={};lower_node=upper_node=-1;++fallback_count;
			};
			auto validate_correction=[&](const char* kind,std::size_t index,int lower_dof,
				int upper_dof,int& lower_node,int& upper_node,double upper_weight,
				double open_area,Vec3d face_centroid,Vec3d& correction,std::size_t& fallback_count)
			{
				const double magnitude=std::sqrt(length2(correction));if(magnitude<=1e-12)return;
				if(lower_node<0||upper_node<0||lower_node>=static_cast<int>(response.size())||
					upper_node>=static_cast<int>(response.size()))
				{
					reject_or_drop(std::string("nonorthogonal pressure correction has no WLS nodes: ")+
						kind+" edge "+std::to_string(index),correction,lower_node,upper_node,fallback_count);
					return;
				}
				Vec3d represented{};for(int column=0;column<3;++column)for(int row=0;row<3;++row)represented[column]+=((1-upper_weight)*response[lower_node][3*row+column]+upper_weight*response[upper_node][3*row+column])*correction[row];const double error=std::sqrt(length2(represented-correction));if(error<=1e-8+1e-3*magnitude)return;std::ostringstream message;message<<"unsupported nonorthogonal pressure correction: "<<kind<<" edge "<<index<<" dofs ["<<lower_dof<<','<<upper_dof<<"] nodes ["<<lower_node<<','<<upper_node<<"] ranks ["<<static_cast<int>(system.pressure_gradient_rank[lower_node])<<','<<static_cast<int>(system.pressure_gradient_rank[upper_node])<<"] rings ["<<static_cast<int>(system.pressure_gradient_ring[lower_node])<<','<<static_cast<int>(system.pressure_gradient_ring[upper_node])<<"] samples ["<<(system.pressure_gradient_offset[lower_node+1]-system.pressure_gradient_offset[lower_node])<<','<<(system.pressure_gradient_offset[upper_node+1]-system.pressure_gradient_offset[upper_node])<<"] volumes ["<<system.volume[lower_dof]<<','<<system.volume[upper_dof]<<"] centroids [["<<system.centroid[lower_dof].x<<','<<system.centroid[lower_dof].y<<','<<system.centroid[lower_dof].z<<"],["<<system.centroid[upper_dof].x<<','<<system.centroid[upper_dof].y<<','<<system.centroid[upper_dof].z<<"]] face=["<<face_centroid.x<<','<<face_centroid.y<<','<<face_centroid.z<<"] area="<<open_area<<" k=["<<correction.x<<','<<correction.y<<','<<correction.z<<"] relative error="<<error/magnitude;reject_or_drop(message.str(),correction,lower_node,upper_node,fallback_count);
			};
			for(std::size_t q=0;q<system.coarse_fine.size();++q){auto& edge=system.coarse_fine[q];const int lower=edge.direction>0?edge.coarse_dof:edge.fine_dof,upper=edge.direction>0?edge.fine_dof:edge.coarse_dof;validate_correction("coarse/fine",q,lower,upper,edge.lower_gradient_node,edge.upper_gradient_node,edge.upper_gradient_weight,edge.open_area,edge.face_centroid,edge.nonorthogonal_correction,system.qualitative_preview_orthogonal_fallback_coarse_fine);}for(std::size_t q=0;q<system.embedded.size();++q){auto& edge=system.embedded[q];const int lower=edge.direction>0?edge.coarse_dof:edge.fine_dof,upper=edge.direction>0?edge.fine_dof:edge.coarse_dof;validate_correction("embedded",q,lower,upper,edge.lower_gradient_node,edge.upper_gradient_node,edge.upper_gradient_weight,edge.open_area,edge.face_centroid,edge.nonorthogonal_correction,system.qualitative_preview_orthogonal_fallback_embedded);}for(std::size_t q=0;q<system.regular_pressure_corrections.size();++q){auto& edge=system.regular_pressure_corrections[q];validate_correction("regular",q,edge.lower_dof,edge.upper_dof,edge.lower_gradient_node,edge.upper_gradient_node,edge.upper_gradient_weight,edge.open_area,{},edge.nonorthogonal_correction,system.qualitative_preview_orthogonal_fallback_regular);}
		}

		std::vector<Vec3d> reconstruct_pressure_gradients(const CompositeAmrPressureSystem& system,const std::vector<double>& pressure)
		{
			std::vector<Vec3d> gradient(system.pressure_gradient_dof.size());for(int node=0;node<static_cast<int>(gradient.size());++node){const double centre=pressure[system.pressure_gradient_dof[node]];for(int q=system.pressure_gradient_offset[node];q<system.pressure_gradient_offset[node+1];++q)gradient[node]=gradient[node]+system.pressure_gradient_weight[q]*(pressure[system.pressure_gradient_neighbour[q]]-centre);}return gradient;
		}
		void add_nonorthogonal_pressure_flux(const CompositeAmrPressureSystem& system,const CoarseFinePressureConnection& edge,const std::vector<Vec3d>& gradient,std::vector<double>& output)
		{
			if(edge.lower_gradient_node<0||edge.upper_gradient_node<0)return;const int lower=edge.direction>0?edge.coarse_dof:edge.fine_dof,upper=edge.direction>0?edge.fine_dof:edge.coarse_dof;if(!system.active[lower]||!system.active[upper])return;const Vec3d face_gradient=gradient[edge.lower_gradient_node]*(1-edge.upper_gradient_weight)+gradient[edge.upper_gradient_node]*edge.upper_gradient_weight;const double correction=-edge.open_area*dot(edge.nonorthogonal_correction,face_gradient);output[lower]+=correction;output[upper]-=correction;
		}
		double regular_pressure_derivative_delta(const RegularPressureCorrection& edge,
			const std::vector<double>& pressure,const std::vector<Vec3d>& gradient)
		{
			double derivative=edge.two_point_delta*(pressure[edge.upper_dof]-pressure[edge.lower_dof]);
			if(edge.lower_gradient_node>=0&&edge.upper_gradient_node>=0)
			{
				const Vec3d face_gradient=gradient[edge.lower_gradient_node]*(1-edge.upper_gradient_weight)+
					gradient[edge.upper_gradient_node]*edge.upper_gradient_weight;
				derivative+=dot(edge.nonorthogonal_correction,face_gradient);
			}
			return derivative;
		}
	}

	int CompositeAmrPressureSystem::dof(int level, int brick, int i, int j, int k) const
	{
		return level_offset[level] + brick * brick_size * brick_size * brick_size + local_index(brick_size, i, j, k);
	}

	int CompositeAmrPressureSystem::pressure_dof_at_point(
		const AmrEmbeddedBoundaryAtlas& embedded_boundary, Vec3d point) const
	{
		if(!hierarchy)return -1;const BrickLocation owner=hierarchy->locate_finest(point);if(!owner.found())return -1;
		for(std::size_t atlas_index=0;atlas_index<embedded_boundary.levels.size();++atlas_index)
		{
			const AmrEbLevelAtlas& atlas=embedded_boundary.levels[atlas_index];if(atlas.level!=owner.level)continue;
			if(atlas_index>=eb_sampling_maps.size())break;const CompositeEbPressureSamplingMap& map=eb_sampling_maps[atlas_index];const EmbeddedBoundary& eb=atlas.topology;const UniformEbGrid& grid=eb.grid;
			const int i=static_cast<int>(std::floor((point.x-grid.origin.x)/grid.h)),j=static_cast<int>(std::floor((point.y-grid.origin.y)/grid.h)),k=static_cast<int>(std::floor((point.z-grid.origin.z)/grid.h));
			if(i<0||j<0||k<0||i>=grid.nx||j>=grid.ny||k>=grid.nz)break;const int cell=grid.cell_index(i,j,k);if(cell<0||cell>=static_cast<int>(eb.cells.size())||cell>=static_cast<int>(map.cell_dof.size()))break;
			const EbCellTopology& topology=eb.cells[cell];if(topology.state==EbCellState::regular)return map.cell_dof[cell];if(topology.state!=EbCellState::split)return -1;
			const FragmentRef ref=eb.fragment_containing_point(cell,point);
			if(ref==invalid_fragment)return -1;if(fragment_is_regular(ref)){const int regular=regular_fragment_cell(ref);return regular>=0&&regular<static_cast<int>(map.cell_dof.size())?map.cell_dof[regular]:-1;}const int fragment=irregular_fragment_index(ref);return fragment>=0&&fragment<static_cast<int>(map.fragment_dof.size())?map.fragment_dof[fragment]:-1;
		}
		return dof(owner.level,owner.brick,owner.cell.x,owner.cell.y,owner.cell.z);
	}

	CompositeAmrPressureSystem build_composite_amr_pressure_system(const AmrHierarchy& hierarchy, bool outlet)
	{
		CompositeAmrPressureBuildOptions options;options.pressure_outlet_xmax=outlet;
		return build_composite_amr_pressure_system(hierarchy,options);
	}

	CompositeAmrPressureSystem build_composite_amr_pressure_system(const AmrHierarchy& hierarchy,
		const CompositeAmrPressureBuildOptions& options)
	{
		CompositeAmrPressureSystem system; system.hierarchy = &hierarchy; system.brick_size = hierarchy.brick_size(); system.pressure_outlet_xmax = options.pressure_outlet_xmax;
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
						const int fi=parent_fi%bs,fj=parent_fj%bs,fk=parent_fk%bs,face_coordinate=axis==0?fi:(axis==1?fj:fk);Vec3d face_centroid=fine.bricks[child].origin+Vec3d{(fi+0.5)*hf,(fj+0.5)*hf,(fk+0.5)*hf};face_centroid[axis]=fine.bricks[child].origin[axis]+(face_coordinate+(positive?0:1))*hf;const double centre_distance=0.5*hc+0.5*hf;system.coarse_fine.push_back({system.dof(level, brick, ci, cj, ck), system.dof(level + 1, child, fi, fj, fk), hf * hf, centre_distance, centre_distance, static_cast<std::int8_t>(axis),static_cast<std::int8_t>(positive?1:-1),face_centroid});
					}
				}
			}
		}
		finalize_component_gauges(system);finalize_nonorthogonal_pressure_topology(system,
			options.unsupported_nonorthogonal_correction);return system;
	}

	CompositeAmrPressureSystem build_composite_amr_pressure_system(const AmrHierarchy& hierarchy,
		const AmrEmbeddedBoundaryAtlas& atlas, bool outlet)
	{
		CompositeAmrPressureBuildOptions options;options.pressure_outlet_xmax=outlet;
		return build_composite_amr_pressure_system(hierarchy,atlas,options);
	}

	CompositeAmrPressureSystem build_composite_amr_pressure_system(const AmrHierarchy& hierarchy,
		const AmrEmbeddedBoundaryAtlas& atlas,const CompositeAmrPressureBuildOptions& options)
	{
		if(!atlas.ready_for_flow())throw std::invalid_argument(
			"cannot build a composite pressure operator from unresolved AMR embedded-boundary topology");
		CompositeAmrPressureSystem system=build_composite_amr_pressure_system(hierarchy,options);system.gauges.clear();
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
			CompositeEbPressureSamplingMap sampling;sampling.level=level_atlas.level;sampling.cell_dof=cell_global;sampling.fragment_dof.resize(eb.fragments.size(),-1);for(int fragment=0;fragment<static_cast<int>(eb.fragments.size());++fragment)if(level_atlas.owned_cell[eb.fragments[fragment].parent_cell])sampling.fragment_dof[fragment]=map_ref(irregular_fragment(fragment));system.eb_sampling_maps.push_back(std::move(sampling));
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
				const int a=map_ref(aperture.fragment_a),b=map_ref(aperture.fragment_b);if(a<0||b<0||a==b)continue;const bool a_active=system.active[a]!=0,b_active=system.active[b]!=0;if(a_active!=b_active)throw std::runtime_error("embedded aperture has exactly one active pressure endpoint");if(!a_active)continue;const double distance=std::sqrt(length2(system.centroid[a]-system.centroid[b])),normal_distance=std::abs(system.centroid[a][aperture.axis]-system.centroid[b][aperture.axis]);system.embedded.push_back({a,b,aperture.area,distance,normal_distance,aperture.axis,1,aperture.centroid});
			}
			for(const SurfacePatch& patch:eb.patches)
			{
				const BrickLocation owner=hierarchy.locate_finest(patch.centroid);if(!owner.found()||owner.level!=level_atlas.level)continue;const int plus=map_ref(patch.plus_fragment),minus=map_ref(patch.minus_fragment);if(plus<0||minus<0)throw std::runtime_error("owned AMR surface patch has no two-sided pressure mapping");if(patch.plus_fragment!=patch.minus_fragment&&plus==minus)throw std::runtime_error("agglomeration collapsed distinct fabric pressure sides");system.surface_patches.push_back({patch.source_triangle_id,patch.source_face_id,patch.area,patch.centroid,patch.normal,plus,minus});
			}
		}
		for(const CoarseFinePressureConnection& connection:system.coarse_fine)if(!system.active[connection.coarse_dof]||!system.active[connection.fine_dof])throw std::runtime_error("embedded boundary intersects a 2:1 interface; aperture-aware cross-level topology is required");
		finalize_component_gauges(system);finalize_nonorthogonal_pressure_topology(system,
			options.unsupported_nonorthogonal_correction);return system;
	}

	CompositeAmrFluxes make_zero_composite_fluxes(const CompositeAmrPressureSystem& system){CompositeAmrFluxes flux;flux.coarse_fine_velocity.assign(system.coarse_fine.size(),0);flux.embedded_velocity.assign(system.embedded.size(),0);return flux;}

	double composite_mac_carrier_volume(const CompositeAmrPressureSystem& system, int dof_a, int dof_b)
	{
		if(dof_a<0||dof_b<0||dof_a==dof_b||dof_a>=system.storage_size||dof_b>=system.storage_size||
			!system.active[dof_a]||!system.active[dof_b]||system.volume[dof_a]<=0||system.volume[dof_b]<=0)
			throw std::invalid_argument("MAC carrier endpoints must be distinct active fluid control volumes: a="+
				std::to_string(dof_a)+" b="+std::to_string(dof_b)+" storage="+
				std::to_string(system.storage_size)+" active="+
				std::to_string(dof_a>=0&&dof_a<system.storage_size?system.active[dof_a]:0)+"/"+
				std::to_string(dof_b>=0&&dof_b<system.storage_size?system.active[dof_b]:0)+" volume="+
				std::to_string(dof_a>=0&&dof_a<system.storage_size?system.volume[dof_a]:-1)+"/"+
				std::to_string(dof_b>=0&&dof_b<system.storage_size?system.volume[dof_b]:-1));
		return 0.5*(system.volume[dof_a]+system.volume[dof_b]);
	}

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
		const std::vector<Vec3d> pressure_gradient=reconstruct_pressure_gradients(system,pressure);for(const RegularPressureCorrection& edge:system.regular_pressure_corrections){AmrHostLevelFields& values=fields.levels()[edge.level];const BrickFieldLayout& layout=values.layout;Real* face=edge.axis==0?&values.u[layout.u_index(edge.brick,edge.i,edge.j,edge.k)]:(edge.axis==1?&values.v[layout.v_index(edge.brick,edge.i,edge.j,edge.k)]:&values.w[layout.w_index(edge.brick,edge.i,edge.j,edge.k)]);*face=static_cast<Real>(static_cast<double>(*face)-(dt/rho)*regular_pressure_derivative_delta(edge,pressure,pressure_gradient));const int face_coordinate=edge.axis==0?edge.i:(edge.axis==1?edge.j:edge.k);if(face_coordinate==bs){const int neighbour=system.hierarchy->levels()[edge.level].bricks[edge.brick].same_level_neighbor[2*edge.axis+1];if(neighbour>=0){if(edge.axis==0)values.u[layout.u_index(neighbour,0,edge.j,edge.k)]=*face;else if(edge.axis==1)values.v[layout.v_index(neighbour,edge.i,0,edge.k)]=*face;else values.w[layout.w_index(neighbour,edge.i,edge.j,0)]=*face;}}}auto correct_special=[&](const std::vector<CoarseFinePressureConnection>& connections,std::vector<double>& velocity){for(std::size_t edge=0;edge<connections.size();++edge){const auto& connection=connections[edge];const int lower=connection.direction>0?connection.coarse_dof:connection.fine_dof,upper=connection.direction>0?connection.fine_dof:connection.coarse_dof;if(!system.active[lower]||!system.active[upper])continue;double derivative=(pressure[upper]-pressure[lower])*pressure_gradient_factor(connection);if(connection.lower_gradient_node>=0&&connection.upper_gradient_node>=0){const Vec3d face_gradient=pressure_gradient[connection.lower_gradient_node]*(1-connection.upper_gradient_weight)+pressure_gradient[connection.upper_gradient_node]*connection.upper_gradient_weight;derivative+=dot(connection.nonorthogonal_correction,face_gradient);}velocity[edge]-=(dt/rho)*derivative;}};correct_special(system.coarse_fine,special.coarse_fine_velocity);correct_special(system.embedded,special.embedded_velocity);
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
		for (const CoarseFinePressureConnection& edge : coarse_fine)if(active[edge.coarse_dof]&&active[edge.fine_dof])add_edge(output, pressure, edge.coarse_dof, edge.fine_dof, edge.open_area * pressure_gradient_factor(edge));
		for (const CoarseFinePressureConnection& edge : embedded)if(active[edge.coarse_dof]&&active[edge.fine_dof])add_edge(output, pressure, edge.coarse_dof, edge.fine_dof, edge.open_area * pressure_gradient_factor(edge));
		const std::vector<Vec3d> pressure_gradient=reconstruct_pressure_gradients(*this,pressure);for(const CoarseFinePressureConnection& edge:coarse_fine)add_nonorthogonal_pressure_flux(*this,edge,pressure_gradient,output);for(const CoarseFinePressureConnection& edge:embedded)add_nonorthogonal_pressure_flux(*this,edge,pressure_gradient,output);for(const RegularPressureCorrection& edge:regular_pressure_corrections){const double flux=-edge.open_area*regular_pressure_derivative_delta(edge,pressure,pressure_gradient);output[edge.lower_dof]+=flux;output[edge.upper_dof]-=flux;}
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
		for(const CoarseFinePressureConnection& edge:coarse_fine)add_diagonal(edge.coarse_dof,edge.fine_dof,edge.open_area*pressure_gradient_factor(edge));for(const CoarseFinePressureConnection& edge:embedded)add_diagonal(edge.coarse_dof,edge.fine_dof,edge.open_area*pressure_gradient_factor(edge));
		// Merged control-volume centroids move an otherwise implicit Cartesian
		// connection away from its h coefficient.  apply_orthogonal() adds this
		// two-point delta, so the Jacobi/PCG diagonal must contain the same term.
		// Omitting it makes the advertised SPD preconditioner inconsistent with
		// the operator it is solving (the non-orthogonal WLS term remains outside
		// this diagonal by design).
		for(const RegularPressureCorrection& edge:regular_pressure_corrections)
			add_diagonal(edge.lower_dof,edge.upper_dof,edge.open_area*edge.two_point_delta);
		for(const CompositePressureGauge& gauge:gauges)if(gauge.dof>=0&&active[gauge.dof])diagonal[gauge.dof]+=gauge.coefficient;
	}
}
