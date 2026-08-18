#include "core/geometry/local_surface_arrangement.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <utility>
#include <vector>

using namespace paracfd::core;

namespace
{
	int failures=0;
	void check(bool condition,const char* message)
	{
		std::printf("[%s] %s\n",condition?"PASS":"FAIL",message);if(!condition)++failures;
	}
	bool near(double a,double b,double tolerance=1e-10){return std::abs(a-b)<=tolerance;}

	void add_quad(std::vector<LocalSurfaceTriangle>& triangles,Vec3d a,Vec3d b,Vec3d c,Vec3d d,
		std::uint32_t face)
	{
		const std::uint32_t first=static_cast<std::uint32_t>(triangles.size());
		LocalSurfaceTriangle first_triangle{a,b,c,first,face};
		LocalSurfaceTriangle second_triangle{a,c,d,first+1,face};
		first_triangle.edge_kind={FabricEdgeKind::confirmed_free,FabricEdgeKind::confirmed_free,
			FabricEdgeKind::attached};
		second_triangle.edge_kind={FabricEdgeKind::attached,FabricEdgeKind::confirmed_free,
			FabricEdgeKind::confirmed_free};
		const std::uint64_t diagonal_attachment=(static_cast<std::uint64_t>(face)<<32)|first;
		first_triangle.edge_attachment_id[2]=diagonal_attachment;
		second_triangle.edge_attachment_id[0]=diagonal_attachment;
		triangles.push_back(first_triangle);triangles.push_back(second_triangle);
	}

	void reverse_winding(LocalSurfaceTriangle& triangle)
	{
		std::swap(triangle.b,triangle.c);
		std::swap(triangle.edge_kind[0],triangle.edge_kind[2]);
		std::swap(triangle.edge_attachment_id[0],triangle.edge_attachment_id[2]);
	}

	std::vector<LocalSurfaceTriangle> reversed(std::vector<LocalSurfaceTriangle> triangles)
	{
		for(auto& triangle : triangles) reverse_winding(triangle);
		return triangles;
	}

	double volume_sum(const LocalSurfaceArrangement& arrangement)
	{
		double value=0;for(const auto& fragment:arrangement.fragments)value+=fragment.volume;return value;
	}

	double boundary_area(const LocalSurfaceArrangement& arrangement)
	{
		double value=0;for(const auto& aperture:arrangement.boundary_apertures)value+=aperture.area;return value;
	}

	Vec3d manufactured_force(const LocalSurfaceArrangement& arrangement)
	{
		Vec3d force{};auto pressure=[&](int fragment){return arrangement.fragments[fragment].centroid.x<0.5?2.0:1.0;};
		for(const auto& patch:arrangement.surface_patches)
			force=force+patch.normal*((pressure(patch.minus_fragment)-pressure(patch.plus_fragment))*patch.area);
		return force;
	}

	void add_radial_panel(std::vector<LocalSurfaceTriangle>& triangles,double angle,std::uint32_t face)
	{
		const Vec3d centre0{0.5,0.5,0},centre1{0.5,0.5,1};const double radius=2.0;
		const Vec3d outer0{0.5+radius*std::cos(angle),0.5+radius*std::sin(angle),0};
		const Vec3d outer1{outer0.x,outer0.y,1};add_quad(triangles,centre0,outer0,outer1,centre1,face);
	}

	Vec3d map_to_box(const Aabb3d& box,Vec3d point)
	{
		const Vec3d extent=box.hi-box.lo;
		return {box.lo.x+extent.x*point.x,box.lo.y+extent.y*point.y,
			box.lo.z+extent.z*point.z};
	}

	std::vector<Vec3d> x_rectangle(double x,double y0,double y1,double z0,double z1)
	{
		return {{x,y0,z0},{x,y1,z0},{x,y1,z1},{x,y0,z1}};
	}

	std::array<Vec3d,4> x_face(double x,double y0=0,double y1=1,double z0=0,double z1=1)
	{
		return {{{x,y0,z0},{x,y1,z0},{x,y1,z1},{x,y0,z1}}};
	}

	LocalArrangementSharedFaceOwner shared_owner(std::vector<Vec3d> polygon,int fragment)
	{
		LocalArrangementSharedFaceOwner owner;owner.polygon=std::move(polygon);
		owner.fragment=fragment;return owner;
	}

	LocalArrangementSharedFaceBarrier shared_barrier(std::uint32_t triangle,
		std::uint32_t face,std::vector<Vec3d> polygon,Vec3d normal={1,0,0})
	{
		LocalArrangementSharedFaceBarrier barrier;barrier.source_triangle_id=triangle;
		barrier.source_face_id=face;barrier.normal=normal;
		barrier.polygon=std::move(polygon);return barrier;
	}

	LocalArrangementSharedFacePartitionPlane shared_plane(Vec3d normal,double offset)
	{
		return {normal,offset};
	}

	LocalArrangementSharedFaceRegion shared_region(std::vector<std::int8_t> signs,int fragment)
	{
		LocalArrangementSharedFaceRegion region;region.plane_side=std::move(signs);
		region.fragment=fragment;return region;
	}

	LocalArrangementSharedFaceSide regular_side(int fragment,std::int8_t inward)
	{
		LocalArrangementSharedFaceSide side;side.inward_axis_sign=inward;
		side.regions.push_back(shared_region({},fragment));return side;
	}

	std::array<double,4> shared_aperture_measure(
		const LocalArrangementSharedFaceAssembly& assembly)
	{
		double area=0;Vec3d moment{};
		for(const auto& aperture:assembly.apertures)
		{
			area+=aperture.area;moment=moment+aperture.centroid*aperture.area;
		}
		return {area,moment.x,moment.y,moment.z};
	}

	template<class Transform>
	std::vector<Vec3d> transformed_polygon(const std::vector<Vec3d>& polygon,Transform transform)
	{
		std::vector<Vec3d> result;result.reserve(polygon.size());
		for(Vec3d point:polygon)result.push_back(transform(point));return result;
	}
}

int main()
{
	std::setvbuf(stdout,nullptr,_IONBF,0);const Aabb3d cell{{0,0,0},{1,1,1}};
	LocalSurfaceArrangementOptions options;options.maximum_planes=96;options.maximum_atoms=8192;

	std::vector<LocalSurfaceTriangle> flat;
	add_quad(flat,{0.5,0,0},{0.5,1,0},{0.5,1,1},{0.5,0,1},10);
	const LocalSurfaceArrangement flat_result=build_local_surface_arrangement(cell,flat,options);
	double flat_area=0;for(const auto& patch:flat_result.surface_patches)flat_area+=patch.area;
	check(flat_result.valid&&flat_result.fragments.size()==2,"full membrane creates exactly two fluid fragments");
	check(near(volume_sum(flat_result),1)&&near(flat_area,1)&&near(boundary_area(flat_result),6),"full membrane conserves volume, fabric area, and boundary area");
	bool distinct_sides=true;for(const auto& patch:flat_result.surface_patches)distinct_sides=distinct_sides&&patch.plus_fragment!=patch.minus_fragment;
	check(distinct_sides,"full membrane has independent pressure states on its two sides");
	check(locate_local_arrangement_fragment(flat_result,{0.25,0.25,0.25})>=0&&locate_local_arrangement_fragment(flat_result,{0.75,0.25,0.25})>=0&&locate_local_arrangement_fragment(flat_result,{0.25,0.25,0.25})!=locate_local_arrangement_fragment(flat_result,{0.75,0.25,0.25}),"host point locator selects the two membrane sides");

	const LocalSurfaceArrangement reversed_result=build_local_surface_arrangement(cell,reversed(flat),options);
	const Vec3d forward_force=manufactured_force(flat_result),backward_force=manufactured_force(reversed_result);
	check(reversed_result.valid&&reversed_result.fragments.size()==flat_result.fragments.size()&&near(forward_force.x,1)&&near(forward_force.x,backward_force.x)&&near(forward_force.y,backward_force.y)&&near(forward_force.z,backward_force.z),"topology and pressure force are invariant to triangle winding");
	std::vector<LocalSurfaceTriangle> mixed_winding=flat;
	reverse_winding(mixed_winding.front());
	const LocalSurfaceArrangement mixed_winding_result=build_local_surface_arrangement(cell,mixed_winding,options);
	const Vec3d mixed_force=manufactured_force(mixed_winding_result);
	check(mixed_winding_result.valid&&mixed_winding_result.fragments.size()==2&&
		near(mixed_force.x,forward_force.x)&&near(mixed_force.y,forward_force.y)&&
		near(mixed_force.z,forward_force.z),
		"mixed triangle winding on one sheet preserves topology and physical force");

	std::vector<LocalSurfaceTriangle> terminating;
	add_quad(terminating,{0.5,0,0},{0.5,0.5,0},{0.5,0.5,1},{0.5,0,1},20);
	const LocalSurfaceArrangement terminating_result=build_local_surface_arrangement(cell,terminating,options);
	bool joined_sides=true;for(const auto& patch:terminating_result.surface_patches)joined_sides=joined_sides&&patch.plus_fragment==patch.minus_fragment;
	const bool terminating_surface_is_ambiguous=terminating_result.valid&&
		!terminating_result.surface_patches.empty()&&locate_local_arrangement_fragment(
			terminating_result,terminating_result.surface_patches.front().centroid)==-1;
	check(terminating_result.valid&&terminating_result.fragments.size()==1&&joined_sides&&
		terminating_surface_is_ambiguous,
		"finite terminating sheet reconnects only around its free edge and remains non-fluid on the patch");
	std::vector<LocalSurfaceTriangle> unknown_termination=terminating;
	unknown_termination[0].edge_kind[1]=FabricEdgeKind::unknown;
	const LocalSurfaceArrangement unknown_termination_result=build_local_surface_arrangement(
		cell,unknown_termination,options);
	check(!unknown_termination_result.valid,
		"an unknown interior fabric edge is fail-closed instead of treated as an opening");
	std::vector<LocalSurfaceTriangle> uncertified_attachment=terminating;
	uncertified_attachment[0].edge_kind[1]=FabricEdgeKind::attached;
	uncertified_attachment[0].edge_attachment_id[1]=8800;
	const LocalSurfaceArrangement uncertified_attachment_result=build_local_surface_arrangement(
		cell,uncertified_attachment,options);
	check(!uncertified_attachment_result.valid,
		"an attached interior edge without incident geometry is rejected rather than leaked");

	std::vector<LocalSurfaceTriangle> opening;
	add_quad(opening,{0.5,0,0},{0.5,0.4,0},{0.5,0.4,1},{0.5,0,1},21);
	add_quad(opening,{0.5,0.6,0},{0.5,1,0},{0.5,1,1},{0.5,0.6,1},22);
	const LocalSurfaceArrangement opening_result=build_local_surface_arrangement(cell,opening,options);
	bool opening_joins_sides=true;
	for(const auto& patch:opening_result.surface_patches)
		opening_joins_sides=opening_joins_sides&&patch.plus_fragment==patch.minus_fragment;
	check(opening_result.valid&&opening_result.fragments.size()==1&&opening_joins_sides,
		"two membrane pieces retain connectivity through a deliberate geometric opening");

	std::vector<LocalSurfaceTriangle> bounded_hole;
	add_quad(bounded_hole,{0.5,0,0},{0.5,0.35,0},{0.5,0.35,1},{0.5,0,1},24);
	add_quad(bounded_hole,{0.5,0.65,0},{0.5,1,0},{0.5,1,1},{0.5,0.65,1},25);
	add_quad(bounded_hole,{0.5,0.35,0},{0.5,0.65,0},{0.5,0.65,0.35},{0.5,0.35,0.35},26);
	add_quad(bounded_hole,{0.5,0.35,0.65},{0.5,0.65,0.65},{0.5,0.65,1},{0.5,0.35,1},27);
	const LocalSurfaceArrangement bounded_hole_result=build_local_surface_arrangement(cell,bounded_hole,options);
	check(bounded_hole_result.valid&&bounded_hole_result.fragments.size()==1,
		"a membrane with a bounded rectangular hole preserves one connected fluid region");

	LocalSurfaceArrangementOptions tolerant_options=options;
	tolerant_options.contact_tolerance=1e-5;
	tolerant_options.ambiguity_tolerance=5e-5;
	tolerant_options.angular_contact_tolerance=1e-10;
	tolerant_options.angular_ambiguity_tolerance=1e-8;
	std::vector<LocalSurfaceTriangle> sewn_seam;
	add_quad(sewn_seam,{0.5,0,0},{0.5,0.5,0},{0.5,0.5,1},{0.5,0,1},28);
	add_quad(sewn_seam,{0.5,0.5000002,0},{0.5,1,0},{0.5,1,1},{0.5,0.5000002,1},29);
	sewn_seam[0].edge_kind[1]=FabricEdgeKind::attached;
	sewn_seam[3].edge_kind[2]=FabricEdgeKind::attached;
	sewn_seam[0].edge_attachment_id[1]=9001;
	sewn_seam[3].edge_attachment_id[2]=9001;
	const LocalSurfaceArrangement sewn_result=build_local_surface_arrangement(cell,sewn_seam,tolerant_options);
	if(!sewn_result.valid)std::printf("[INFO] sewn seam rejection: %s (dV %.17g, dM %.17g)\n",
		sewn_result.error.c_str(),sewn_result.volume_conservation_error,
		sewn_result.first_moment_conservation_error);
	else std::printf("[INFO] sewn seam fragments: %zu\n",sewn_result.fragments.size());
	check(sewn_result.valid&&sewn_result.fragments.size()==2,
		"a contact-certified near-coincident sewn seam is closed rather than converted into an opening");
	std::vector<LocalSurfaceTriangle> free_micro_gap=sewn_seam;
	free_micro_gap[0].edge_kind[1]=FabricEdgeKind::confirmed_free;
	free_micro_gap[3].edge_kind[2]=FabricEdgeKind::confirmed_free;
	free_micro_gap[0].edge_attachment_id[1]=std::numeric_limits<std::uint64_t>::max();
	free_micro_gap[3].edge_attachment_id[2]=std::numeric_limits<std::uint64_t>::max();
	const LocalSurfaceArrangement free_micro_gap_result=build_local_surface_arrangement(
		cell,free_micro_gap,tolerant_options);
	if(!free_micro_gap_result.valid||free_micro_gap_result.fragments.size()!=1)
		std::printf("[INFO] free micro-gap: valid=%d fragments=%zu error='%s' represented=%.17g clipped=%.17g\n",
			free_micro_gap_result.valid?1:0,free_micro_gap_result.fragments.size(),
			free_micro_gap_result.error.c_str(),free_micro_gap_result.represented_surface_area,
			free_micro_gap_result.clipped_input_surface_area);
	check(free_micro_gap_result.valid&&free_micro_gap_result.fragments.size()==1,
		"the identical sub-resolution remainder stays open when its local edges are confirmed free");

	// A genuine fabric patch can be an arithmetic-scale fraction of its support
	// facet after later cusp planes subdivide it.  Its free perimeter keeps the fluid
	// connected, but the pressure patch itself must never be discarded by an area-
	// fraction shortcut.  The area is intentionally below the historical global
	// surface-conservation tolerance so this test exercises topology/provenance rather
	// than relying on the final aggregate audit to notice the deletion.
	std::vector<LocalSurfaceTriangle> tiny_fabric;
	add_quad(tiny_fabric,{0.5,0.499999,0.499999},{0.5,0.500001,0.499999},
		{0.5,0.500001,0.500001},{0.5,0.499999,0.500001},34);
	auto tiny_fabric_signature=[&](std::vector<LocalSurfaceTriangle> triangles)
	{
		const LocalSurfaceArrangement arrangement=build_local_surface_arrangement(
			cell,triangles,options);
		double area=0.0;bool same_fragment=true;
		for(const auto& patch:arrangement.surface_patches)
		{
			area+=patch.area;
			same_fragment=same_fragment&&patch.plus_fragment==patch.minus_fragment;
		}
		return std::array<double,4>{arrangement.valid?1.0:0.0,
			static_cast<double>(arrangement.fragments.size()),same_fragment?1.0:0.0,area};
	};
	const auto tiny_fabric_forward=tiny_fabric_signature(tiny_fabric);
	std::vector<LocalSurfaceTriangle> tiny_fabric_permuted=reversed(tiny_fabric);
	std::reverse(tiny_fabric_permuted.begin(),tiny_fabric_permuted.end());
	const auto tiny_fabric_reverse=tiny_fabric_signature(tiny_fabric_permuted);
	check(tiny_fabric_forward[0]==1&&tiny_fabric_forward[1]==1&&
		tiny_fabric_forward[2]==1&&near(tiny_fabric_forward[3],4e-12,1e-15)&&
		near(tiny_fabric_reverse[0],tiny_fabric_forward[0])&&
		near(tiny_fabric_reverse[1],tiny_fabric_forward[1])&&
		near(tiny_fabric_reverse[2],tiny_fabric_forward[2])&&
		near(tiny_fabric_reverse[3],tiny_fabric_forward[3],1e-15),
		"arithmetic-scale finite fabric is retained under winding and triangle permutation");

	// Four CAD faces tessellated independently form an open-ended square tube.  Each
	// longitudinal seam has reciprocal attachment provenance, but its two coordinate
	// copies differ inside the CAD contact tolerance.  The tube interior and exterior
	// must remain disconnected without relying on any profile-specific geometry.
	std::vector<LocalSurfaceTriangle> independently_tessellated_tube;
	const std::array<Vec3d,4> tube_corner{{{0.3,0.3,0},{0.7,0.3,0},
		{0.7,0.7,0},{0.3,0.7,0}}};
	const std::array<Vec3d,4> panel_offset{{{0,0,0},{2e-6,-1e-6,0},
		{-1.5e-6,1e-6,0},{1e-6,2e-6,0}}};
	for(int panel=0;panel<4;++panel)
	{
		const int next=(panel+1)%4;
		Vec3d a=tube_corner[panel]+panel_offset[panel];
		Vec3d d=tube_corner[next]+panel_offset[panel];
		Vec3d b=a,c=d;b.z=1;c.z=1;
		const int first=static_cast<int>(independently_tessellated_tube.size());
		add_quad(independently_tessellated_tube,a,b,c,d,100+panel);
		independently_tessellated_tube[first].edge_kind[0]=FabricEdgeKind::attached;
		independently_tessellated_tube[first].edge_attachment_id[0]=10000+panel;
		independently_tessellated_tube[first+1].edge_kind[1]=FabricEdgeKind::attached;
		independently_tessellated_tube[first+1].edge_attachment_id[1]=10000+next;
	}
	const LocalSurfaceArrangement tube_result=build_local_surface_arrangement(
		cell,independently_tessellated_tube,tolerant_options);
	if(!tube_result.valid)std::printf("[INFO] independently tessellated tube rejection: %s\n",
		tube_result.error.c_str());
	const int tube_inside=locate_local_arrangement_fragment(tube_result,{0.5,0.5,0.5});
	const int tube_outside=locate_local_arrangement_fragment(tube_result,{0.1,0.1,0.5});
	check(tube_result.valid&&tube_result.fragments.size()==2&&tube_inside>=0&&tube_outside>=0&&
		tube_inside!=tube_outside,
		"certified independently tessellated seams preserve a closed loop generically");

	std::vector<LocalSurfaceTriangle> resolved_gap;
	add_quad(resolved_gap,{0.5,0,0},{0.5,0.45,0},{0.5,0.45,1},{0.5,0,1},30);
	add_quad(resolved_gap,{0.5,0.55,0},{0.5,1,0},{0.5,1,1},{0.5,0.55,1},31);
	const LocalSurfaceArrangement resolved_gap_result=build_local_surface_arrangement(cell,resolved_gap,tolerant_options);
	check(resolved_gap_result.valid&&resolved_gap_result.fragments.size()==1,
		"a resolved gap bounded by confirmed-free edges remains an intentional opening");

	std::vector<LocalSurfaceTriangle> ambiguous_gap;
	add_quad(ambiguous_gap,{0.5,0,0},{0.5,0.5,0},{0.5,0.5,1},{0.5,0,1},32);
	add_quad(ambiguous_gap,{0.5,0.50002,0},{0.5,1,0},{0.5,1,1},{0.5,0.50002,1},33);
	ambiguous_gap[0].edge_kind[1]=FabricEdgeKind::unknown;
	ambiguous_gap[3].edge_kind[2]=FabricEdgeKind::unknown;
	const LocalSurfaceArrangement ambiguous_gap_result=build_local_surface_arrangement(cell,
		ambiguous_gap,tolerant_options);
	check(!ambiguous_gap_result.valid,
		"an ambiguity-band seam is rejected instead of silently opened or sewn");

	const Vec3d inclined_normal=normalized(Vec3d{1,2,3});
	const Vec3d inclined_u=normalized(cross(inclined_normal,Vec3d{0,0,1}));
	const Vec3d inclined_v=cross(inclined_normal,inclined_u);
	const Vec3d inclined_centre{0.5,0.5,0.5};
	std::vector<LocalSurfaceTriangle> inclined;
	add_quad(inclined,inclined_centre-inclined_u*2.0-inclined_v*2.0,
		inclined_centre+inclined_u*2.0-inclined_v*2.0,
		inclined_centre+inclined_u*2.0+inclined_v*2.0,
		inclined_centre-inclined_u*2.0+inclined_v*2.0,23);
	const LocalSurfaceArrangement inclined_result=build_local_surface_arrangement(cell,inclined,options);
	const int inclined_minus=locate_local_arrangement_fragment(inclined_result,
		inclined_centre-inclined_normal*0.1);
	const int inclined_plus=locate_local_arrangement_fragment(inclined_result,
		inclined_centre+inclined_normal*0.1);
	check(inclined_result.valid&&inclined_result.fragments.size()==2&&inclined_minus>=0&&
		inclined_plus>=0&&inclined_minus!=inclined_plus&&near(volume_sum(inclined_result),1),
		"arbitrarily inclined membrane divides the cell without axis-aligned assumptions");

	std::vector<LocalSurfaceTriangle> cross_sheets;
	add_quad(cross_sheets,{0.5,0,0},{0.5,1,0},{0.5,1,1},{0.5,0,1},30);
	add_quad(cross_sheets,{0,0.5,0},{1,0.5,0},{1,0.5,1},{0,0.5,1},31);
	const LocalSurfaceArrangement cross_result=build_local_surface_arrangement(cell,cross_sheets,options);
	check(cross_result.valid&&cross_result.fragments.size()==4,"two full intersecting sheets create four sectors without an at-most-two assumption");

	std::vector<LocalSurfaceTriangle> t_junction;
	add_quad(t_junction,{0,0,0.5},{1,0,0.5},{1,1,0.5},{0,1,0.5},40);
	add_quad(t_junction,{0.5,0,0},{0.5,1,0},{0.5,1,0.5},{0.5,0,0.5},41);
	const LocalSurfaceArrangement t_result=build_local_surface_arrangement(cell,t_junction,options);
	check(t_result.valid&&t_result.fragments.size()==3,"rib terminating on a skin creates the three geometric fluid sectors");

	std::vector<LocalSurfaceTriangle> skin_with_partial_rib;
	add_quad(skin_with_partial_rib,{0,0,0.5},{1,0,0.5},{1,1,0.5},{0,1,0.5},120);
	const int partial_rib_first=static_cast<int>(skin_with_partial_rib.size());
	add_quad(skin_with_partial_rib,{0.5,0.25,0.2},{0.5,0.75,0.2},
		{0.5,0.75,0.5},{0.5,0.25,0.5},121);
	// The rib top is CAD-certified against the interior of the uninterrupted skin;
	// its remaining perimeter is a real free edge.  That free edge must reconnect
	// the two rib sides without ever authorising a passage through the skin.
	skin_with_partial_rib[partial_rib_first+1].edge_kind[1]=FabricEdgeKind::attached;
	const LocalSurfaceArrangement partial_rib_result=build_local_surface_arrangement(
		cell,skin_with_partial_rib,options);
	bool skin_sides_distinct=partial_rib_result.valid;
	for(const auto& patch:partial_rib_result.surface_patches)
		if(patch.source_face_id==120)
			skin_sides_distinct=skin_sides_distinct&&patch.plus_fragment!=patch.minus_fragment;
	check(partial_rib_result.valid&&partial_rib_result.fragments.size()==2&&skin_sides_distinct,
		"a partial rib free edge cannot authorise reconnection through its full skin");

	std::vector<LocalSurfaceTriangle> wedge;add_radial_panel(wedge,-0.02,50);add_radial_panel(wedge,0.02,51);
	const LocalSurfaceArrangement wedge_result=build_local_surface_arrangement(cell,wedge,options);
	if(!wedge_result.valid) std::printf("[INFO] wedge rejection: %s\n",wedge_result.error.c_str());
	check(wedge_result.valid&&wedge_result.fragments.size()==2,"thin attached wedge retains its interior and exterior sectors independently of a subcell raster");
	bool shallow_attached_wedges_valid=true;
	for(double half_angle:{1e-3,1e-5})
	{
		std::vector<LocalSurfaceTriangle> shallow_wedge;
		add_radial_panel(shallow_wedge,-half_angle,54);
		add_radial_panel(shallow_wedge,half_angle,55);
		shallow_wedge[1].edge_kind[2]=FabricEdgeKind::attached;
		shallow_wedge[3].edge_kind[2]=FabricEdgeKind::attached;
		shallow_wedge[1].edge_attachment_id[2]=9300;
		shallow_wedge[3].edge_attachment_id[2]=9300;
		const LocalSurfaceArrangement shallow_result=build_local_surface_arrangement(
			cell,shallow_wedge,options);
		if(!shallow_result.valid||shallow_result.fragments.size()!=2)
			std::printf("[INFO] shallow attached wedge %.3g: valid=%d fragments=%zu error: %s\n",
				half_angle,shallow_result.valid?1:0,shallow_result.fragments.size(),
				shallow_result.error.c_str());
		shallow_attached_wedges_valid=shallow_attached_wedges_valid&&shallow_result.valid&&
			shallow_result.fragments.size()==2;
	}
	check(shallow_attached_wedges_valid,
		"certified shallow attached cusps keep their two fluid sectors");

	// The second panel's support plane only grazes an atom created by the first
	// panel's support plane. Both finite panels are far from that artificial sliver,
	// so discarding the BSP debris is safe; strict surface-area conservation must
	// nevertheless retain every positive-area piece of both panels.
	std::vector<LocalSurfaceTriangle> grazing_t;
	add_quad(grazing_t,{0.5,0.5,0},{0.5,1,0},{0.5,1,1},{0.5,0.5,1},52);
	const double seam_y=1e-10;
	add_quad(grazing_t,{0.5,seam_y,0},{1,0.05+seam_y,0},{1,0.05+seam_y,1},{0.5,seam_y,1},53);
	const LocalSurfaceArrangement grazing_t_result=build_local_surface_arrangement(cell,grazing_t,options);
	double grazing_area=0;for(const auto& patch:grazing_t_result.surface_patches)grazing_area+=patch.area;
	if(!grazing_t_result.valid||grazing_t_result.fragments.size()!=1||!near(volume_sum(grazing_t_result),1)||!near(grazing_area,0.5+std::sqrt(0.2525),1e-10))
		std::printf("[INFO] grazing support: valid=%d error='%s' fragments=%zu volume=%.17g area=%.17g expected-area=%.17g represented=%.17g clipped=%.17g dV=%.17g dM=%.17g\n",
			grazing_t_result.valid?1:0,grazing_t_result.error.c_str(),grazing_t_result.fragments.size(),
			volume_sum(grazing_t_result),grazing_area,0.5+std::sqrt(0.2525),
			grazing_t_result.represented_surface_area,grazing_t_result.clipped_input_surface_area,
			grazing_t_result.volume_conservation_error,grazing_t_result.first_moment_conservation_error);
	check(grazing_t_result.valid&&grazing_t_result.fragments.size()==1&&near(volume_sum(grazing_t_result),1)&&near(grazing_area,0.5+std::sqrt(0.2525),1e-10),"grazing support plane discards only a sub-tolerance BSP sliver while preserving fabric");
	std::vector<LocalSurfaceTriangle> grazing_permuted=reversed(grazing_t);
	std::reverse(grazing_permuted.begin(),grazing_permuted.end());
	const LocalSurfaceArrangement grazing_permuted_result=build_local_surface_arrangement(
		cell,grazing_permuted,options);
	double grazing_permuted_area=0.0;
	for(const auto& patch:grazing_permuted_result.surface_patches)
		grazing_permuted_area+=patch.area;
	check(grazing_permuted_result.valid&&grazing_permuted_result.fragments.size()==1&&
		near(volume_sum(grazing_permuted_result),1)&&
		near(grazing_permuted_area,0.5+std::sqrt(0.2525),1e-10),
		"grazing topology is invariant to support-plane insertion order and winding");

	bool fan_range_valid=true;
	for(int fan_count=3;fan_count<=8;++fan_count)
	{
		std::vector<LocalSurfaceTriangle> fan;
		for(int sector=0;sector<fan_count;++sector)
			add_radial_panel(fan,2.0*3.14159265358979323846*sector/fan_count,60+sector);
		const LocalSurfaceArrangement fan_result=build_local_surface_arrangement(cell,fan,options);
		if(!fan_result.valid) std::printf("[INFO] %d-sheet fan rejection: %s (represented %.17g, clipped %.17g)\n",
			fan_count,fan_result.error.c_str(),fan_result.represented_surface_area,
			fan_result.clipped_input_surface_area);
		if(!fan_result.valid) std::printf("[INFO] fan conservation dV %.17g, dM %.17g\n",
			fan_result.volume_conservation_error,fan_result.first_moment_conservation_error);
		if(!fan_result.valid)
		{
			for(int face=60;face<60+fan_count;++face)
			{
				double area=0.0;
				for(const auto& patch:fan_result.surface_patches)
					if(static_cast<int>(patch.source_face_id)==face) area+=patch.area;
				std::printf("[INFO] fan face %d represented area %.17g\n",face,area);
			}
		}
		fan_range_valid=fan_range_valid&&fan_result.valid&&
			static_cast<int>(fan_result.fragments.size())==fan_count;
	}
	check(fan_range_valid,"attached fans with three through eight sheets retain every fluid sector");

	std::vector<LocalSurfaceTriangle> near_pivot_pair;
	const double near_pivot_angle=0.25*3.14159265358979323846+1e-14;
	add_radial_panel(near_pivot_pair,near_pivot_angle,72);
	add_radial_panel(near_pivot_pair,near_pivot_angle+3.14159265358979323846,73);
	const LocalSurfaceArrangement near_pivot_result=build_local_surface_arrangement(cell,
		near_pivot_pair,tolerant_options);
	check(near_pivot_result.valid&&near_pivot_result.fragments.size()==2,
		"opposite near-45-degree half-sheets register as one unoriented support plane");

	const Aabb3d translated_cell{{1000000.0,-2000000.0,3000000.0},
		{1000002.0,-1999997.0,3000000.25}};
	std::vector<LocalSurfaceTriangle> translated_sheet;
	add_quad(translated_sheet,map_to_box(translated_cell,{0,0,0}),
		map_to_box(translated_cell,{1,1,0}),map_to_box(translated_cell,{1,1,1}),
		map_to_box(translated_cell,{0,0,1}),74);
	LocalSurfaceArrangementOptions translated_options=options;
	translated_options.contact_tolerance=1e-8;
	translated_options.ambiguity_tolerance=5e-8;
	translated_options.angular_contact_tolerance=1e-10;
	translated_options.angular_ambiguity_tolerance=1e-8;
	const LocalSurfaceArrangement translated_result=build_local_surface_arrangement(
		translated_cell,translated_sheet,translated_options);
	if(!translated_result.valid)std::printf("[INFO] translated rejection: %s (dV %.17g, dM %.17g)\n",
		translated_result.error.c_str(),translated_result.volume_conservation_error,
		translated_result.first_moment_conservation_error);
	else std::printf("[INFO] translated fragments/volume: %zu / %.17g\n",
		translated_result.fragments.size(),volume_sum(translated_result));
	check(translated_result.valid&&translated_result.fragments.size()==2&&
		near(volume_sum(translated_result),1.5,1e-8),
		"translated anisotropic cell preserves oblique topology with positive BVH-like tolerances");

	std::vector<LocalSurfaceTriangle> skipped_inputs;
	skipped_inputs.push_back({{0,0,0},{0,0,0},{0,0,0},800,80});
	skipped_inputs.push_back({{2,2,2},{3,2,2},{2,3,2},801,81});
	skipped_inputs.insert(skipped_inputs.end(),flat.begin(),flat.end());
	const LocalSurfaceArrangement skipped_result=build_local_surface_arrangement(cell,
		skipped_inputs,options);
	bool original_indices_preserved=false;
	for(const auto& plane:skipped_result.planes)
		if(plane.support_triangles==std::vector<int>{2,3})original_indices_preserved=true;
	check(skipped_result.valid&&original_indices_preserved,
		"public support-triangle provenance retains original indices after skipped inputs");

	std::vector<LocalSurfaceTriangle> touching_peer;
	add_quad(touching_peer,{0,0.5,0},{1,0.5,0},{1,0.5,1},{0,0.5,1},82);
	touching_peer[0].edge_kind[1]=FabricEdgeKind::attached;
	touching_peer[0].edge_attachment_id[1]=9200;
	LocalSurfaceTriangle zero_area_peer{{1,0.5,0},{1,0.5,1},{2,0.75,0.5},999,83};
	zero_area_peer.edge_kind={FabricEdgeKind::attached,FabricEdgeKind::confirmed_free,
		FabricEdgeKind::confirmed_free};
	zero_area_peer.edge_attachment_id[0]=9200;
	touching_peer.push_back(zero_area_peer);
	const LocalSurfaceArrangement touching_peer_result=build_local_surface_arrangement(cell,
		touching_peer,options);
	check(touching_peer_result.valid&&touching_peer_result.fragments.size()==2&&
		near(touching_peer_result.clipped_input_surface_area,1),
		"zero-area AABB-touching incident triangle certifies attachment without adding fabric area");

	std::vector<LocalSurfaceTriangle> boundary_sheet;
	add_quad(boundary_sheet,{1,0,0},{1,1,0},{1,1,1},{1,0,1},70);
	const LocalSurfaceArrangement boundary_result=build_local_surface_arrangement(cell,boundary_sheet,options);
	double boundary_fabric_area=0;for(const auto& patch:boundary_result.boundary_surface_patches)boundary_fabric_area+=patch.area;
	check(boundary_result.valid&&boundary_result.fragments.size()==1&&near(boundary_area(boundary_result),5)&&near(boundary_fabric_area,1),"fabric coincident with a Cartesian face removes only that face aperture");

	std::vector<LocalSurfaceTriangle> partial_boundary_sheet;
	add_quad(partial_boundary_sheet,{1,0.25,0},{1,0.75,0},{1,0.75,1},{1,0.25,1},71);
	const LocalSurfaceArrangement partial_boundary_result=build_local_surface_arrangement(
		cell,partial_boundary_sheet,options);
	double partial_boundary_open=0,partial_boundary_fabric=0;int partial_boundary_open_count=0;
	for(const auto& aperture:partial_boundary_result.boundary_apertures)
		if(aperture.axis==0&&aperture.upper)
		{
			partial_boundary_open+=aperture.area;++partial_boundary_open_count;
		}
	for(const auto& patch:partial_boundary_result.boundary_surface_patches)
		if(patch.axis==0&&patch.upper)partial_boundary_fabric+=patch.area;
	check(partial_boundary_result.valid&&partial_boundary_result.fragments.size()==1&&
		partial_boundary_open_count==2&&near(partial_boundary_open,0.5)&&
		near(partial_boundary_fabric,0.5)&&near(boundary_area(partial_boundary_result),5.5),
		"finite face-aligned fabric partitions exact open apertures without blocking the face");

	// CAD attachment tolerance can intentionally be much larger than the arithmetic
	// alignment tolerance used for Cartesian clipping.  A near-face interior membrane
	// must remain an interior barrier; otherwise one cell would emit a boundary owner
	// while the shared-face barrier scan correctly emits no barrier.
	std::vector<LocalSurfaceTriangle> near_boundary_sheet;
	add_quad(near_boundary_sheet,{0.999,0,0},{0.999,1,0},{0.999,1,1},{0.999,0,1},72);
	LocalSurfaceArrangementOptions separate_face_tolerance=options;
	separate_face_tolerance.contact_tolerance=0.01;
	separate_face_tolerance.ambiguity_tolerance=0.02;
	separate_face_tolerance.cartesian_face_tolerance=1e-8;
	separate_face_tolerance.angular_cartesian_face_tolerance=1e-10;
	const LocalSurfaceArrangement near_boundary_result=build_local_surface_arrangement(
		cell,near_boundary_sheet,separate_face_tolerance);
	double near_boundary_interior_area=0.0;
	for(const auto& patch:near_boundary_result.surface_patches)
		near_boundary_interior_area+=patch.area;
	check(near_boundary_result.valid&&near_boundary_result.fragments.size()==2&&
		near_boundary_result.boundary_surface_patches.empty()&&
		near(near_boundary_interior_area,1,1e-10)&&
		near(near_boundary_result.cartesian_face_tolerance,1e-8,1e-16),
		"Cartesian-face coplanarity is independent of a larger CAD contact tolerance");

	const LocalSurfaceArrangement empty_a=build_local_surface_arrangement(cell,{},options);
	const LocalSurfaceArrangement empty_b=build_local_surface_arrangement({{1,0,0},{2,1,1}}, {},options);
	const auto find_face=[](const LocalSurfaceArrangement& arrangement,int axis,bool upper)->const LocalArrangementBoundaryAperture*
	{
		for(const auto& aperture : arrangement.boundary_apertures)
			if(aperture.axis == axis && aperture.upper == upper) return &aperture;
		return nullptr;
	};
	const auto* a_face=find_face(empty_a,0,true);const auto* b_face=find_face(empty_b,0,false);LocalArrangementBoundaryOverlap overlap;
	check(empty_a.valid&&empty_b.valid&&a_face&&b_face&&intersect_local_boundary_apertures(*a_face,*b_face,overlap)&&near(overlap.area,1),"neighbour boundary pieces form an exact full Cartesian-face aperture");

	const Aabb3d large_a_cell{{1e8,2e8,-3e8},{1e8+1,2e8+1,-3e8+1}};
	const Aabb3d large_b_cell{{1e8+1,2e8,-3e8},{1e8+2,2e8+1,-3e8+1}};
	const LocalSurfaceArrangement large_a=build_local_surface_arrangement(large_a_cell,{},options);
	const LocalSurfaceArrangement large_b=build_local_surface_arrangement(large_b_cell,{},options);
	const auto* large_a_face=find_face(large_a,0,true);
	const auto* large_b_face=find_face(large_b,0,false);
	LocalArrangementBoundaryOverlap large_overlap;
	const bool large_intersects=large_a.valid&&large_b.valid&&large_a_face&&large_b_face&&
		intersect_local_boundary_apertures(*large_a_face,*large_b_face,large_overlap);
	if(!large_intersects)std::printf("[INFO] translated overlap failed; tolerances %.17g / %.17g\n",
		large_a.contact_tolerance,large_b.contact_tolerance);
	check(large_a.valid&&large_b.valid&&large_a_face&&large_b_face&&
		large_intersects&&
		near(large_overlap.area,1,1e-8),
		"default boundary-overlap tolerance scales across translated tangential coordinates");

	auto boundary_piece=[](bool upper,double y0,double y1,double z0,double z1,
		Vec3d normal,int fragment)
	{
		LocalArrangementBoundarySurfacePatch patch;
		patch.source_triangle_id=910;patch.source_face_id=91;patch.axis=0;patch.upper=upper;
		patch.area=(y1-y0)*(z1-z0);patch.centroid={1,0.5*(y0+y1),0.5*(z0+z1)};
		patch.normal=normal;patch.polygon={{1,y0,z0},{1,y1,z0},{1,y1,z1},{1,y0,z1}};
		patch.interior_fragment=fragment;
		const Vec3d inward=upper?Vec3d{-1,0,0}:Vec3d{1,0,0};
		patch.interior_is_plus=dot(normal,inward)>0;patch.contact_tolerance=1e-10;
		return patch;
	};
	std::vector<LocalArrangementBoundarySurfacePatch> subdivision_a{
		boundary_piece(true,0,0.5,0,1,{1,0,0},10),
		boundary_piece(true,0.5,1,0,1,{1,0,0},10)};
	std::vector<LocalArrangementBoundarySurfacePatch> subdivision_b{
		boundary_piece(false,0,1,0,0.25,{-1,0,0},20),
		boundary_piece(false,0,1,0.25,0.7,{-1,0,0},20),
		boundary_piece(false,0,1,0.7,1,{-1,0,0},20)};
	const LocalArrangementBoundarySurfaceAssembly assembly_ab=
		assemble_local_boundary_surface_patches(subdivision_a,subdivision_b);
	const LocalArrangementBoundarySurfaceAssembly assembly_ba=
		assemble_local_boundary_surface_patches(subdivision_b,subdivision_a);
	auto assembly_force=[](const LocalArrangementBoundarySurfaceAssembly& assembly)
	{
		Vec3d force{};
		auto pressure=[](int fragment){return fragment==10?2.0:1.0;};
		for(const auto& patch:assembly.patches)
			force=force+patch.normal*((pressure(patch.minus_fragment)-
				pressure(patch.plus_fragment))*patch.area);
		return force;
	};
	const Vec3d force_ab=assembly_force(assembly_ab),force_ba=assembly_force(assembly_ba);
	bool owner_ab=true,owner_ba=true;
	for(const auto& patch:assembly_ab.patches)owner_ab=owner_ab&&patch.owner_is_a;
	for(const auto& patch:assembly_ba.patches)owner_ba=owner_ba&&!patch.owner_is_a;
	check(assembly_ab.valid&&assembly_ba.valid&&assembly_ab.patches.size()==6&&
		assembly_ba.patches.size()==6&&near(assembly_ab.area,1)&&near(assembly_ba.area,1)&&
		owner_ab&&owner_ba&&near(force_ab.x,1)&&near(force_ba.x,1),
		"N-by-M boundary fabric assembly is conservative, winding-invariant, and owned once");
	std::vector<LocalArrangementBoundarySurfacePatch> incomplete_b=subdivision_b;
	incomplete_b.pop_back();const LocalArrangementBoundarySurfaceAssembly rejected_assembly=
		assemble_local_boundary_surface_patches(subdivision_a,incomplete_b);
	check(!rejected_assembly.valid&&rejected_assembly.patches.empty()&&near(rejected_assembly.area,0),
		"mismatched boundary subdivision fails transactionally without orphan patches");

	// Canonical shared-face assembly uses one common 2-D arrangement.  Side A and B
	// may subdivide an otherwise ordinary aperture differently; only their fragment
	// labels and total cover matter.
	const auto canonical_face=x_face(1.0);
	std::vector<LocalArrangementSharedFaceOwner> mismatched_owners_a{
		shared_owner(x_rectangle(1,0,0.5,0,1),10),
		shared_owner(x_rectangle(1,0.5,1,0,1),10)};
	std::vector<LocalArrangementSharedFaceOwner> mismatched_owners_b{
		shared_owner(x_rectangle(1,0,1,0,0.25),20),
		shared_owner(x_rectangle(1,0,1,0.25,0.7),20),
		shared_owner(x_rectangle(1,0,1,0.7,1),20)};
	const LocalArrangementSharedFaceAssembly canonical_mismatched=
		assemble_canonical_local_shared_face(canonical_face,0,{}, {},mismatched_owners_a,
			mismatched_owners_b);
	const auto canonical_mismatched_measure=shared_aperture_measure(canonical_mismatched);
	bool canonical_mismatched_labels=true;
	for(const auto& aperture:canonical_mismatched.apertures)
		canonical_mismatched_labels=canonical_mismatched_labels&&
			aperture.fragment_a==10&&aperture.fragment_b==20;
	check(canonical_mismatched.valid&&canonical_mismatched.apertures.size()==6&&
		canonical_mismatched.blocked_surfaces.empty()&&canonical_mismatched_labels&&
		near(canonical_mismatched_measure[0],1,1e-12)&&
		near(canonical_mismatched_measure[1],1,1e-12)&&
		near(canonical_mismatched_measure[2],0.5,1e-12)&&
		near(canonical_mismatched_measure[3],0.5,1e-12),
		"canonical shared face common-refines mismatched ownership subdivisions exactly once");

	// A transverse inclined sheet contributes only a one-dimensional trace on the
	// Cartesian face.  It selects the matching fluid fragments without blocking
	// positive face area.
	const LocalArrangementSharedFaceTrace inclined_trace{501,50,{1,0,0.2},{1,1,0.8}};
	const std::vector<Vec3d> inclined_below{{1,0,0},{1,1,0},{1,1,0.8},{1,0,0.2}};
	const std::vector<Vec3d> inclined_above{{1,0,0.2},{1,1,0.8},{1,1,1},{1,0,1}};
	std::vector<LocalArrangementSharedFaceOwner> inclined_a{
		shared_owner(inclined_below,30),shared_owner(inclined_above,31)};
	std::vector<LocalArrangementSharedFaceOwner> inclined_b{
		shared_owner(inclined_below,40),shared_owner(inclined_above,41)};
	const LocalArrangementSharedFaceAssembly canonical_inclined=
		assemble_canonical_local_shared_face(canonical_face,0,{inclined_trace},{},inclined_a,inclined_b);
	bool inclined_pairs=true;double inclined_area=0;
	for(const auto& aperture:canonical_inclined.apertures)
	{
		inclined_area+=aperture.area;
		inclined_pairs=inclined_pairs&&((aperture.fragment_a==30&&aperture.fragment_b==40)||
			(aperture.fragment_a==31&&aperture.fragment_b==41));
	}
	check(canonical_inclined.valid&&canonical_inclined.apertures.size()==2&&
		canonical_inclined.blocked_surfaces.empty()&&inclined_pairs&&near(inclined_area,1,1e-12),
		"inclined physical trace pairs same-region fragments without inventing blocked area");

	// An arithmetic-scale opening is still a physical positive-area tile.  The
	// canonical assembler has no CAD-contact or configured minimum-aperture cutoff.
	constexpr double microscopic_gap_target=3.4911e-16;
	const double microscopic_width=std::sqrt(microscopic_gap_target);
	const double microscopic_lo=0.5-0.5*microscopic_width;
	const double microscopic_hi=0.5+0.5*microscopic_width;
	std::vector<LocalArrangementSharedFaceBarrier> microscopic_barriers{
		shared_barrier(600,60,x_rectangle(1,0,microscopic_lo,0,1)),
		shared_barrier(601,60,x_rectangle(1,microscopic_hi,1,0,1)),
		shared_barrier(602,60,x_rectangle(1,microscopic_lo,microscopic_hi,0,microscopic_lo)),
		shared_barrier(603,60,x_rectangle(1,microscopic_lo,microscopic_hi,microscopic_hi,1))};
	std::vector<LocalArrangementSharedFaceOwner> whole_a{
		shared_owner(x_rectangle(1,0,1,0,1),50)};
	std::vector<LocalArrangementSharedFaceOwner> whole_b{
		shared_owner(x_rectangle(1,0,1,0,1),60)};
	const LocalArrangementSharedFaceAssembly canonical_microscopic=
		assemble_canonical_local_shared_face(canonical_face,0,{},microscopic_barriers,
			whole_a,whole_b);
	check(canonical_microscopic.valid&&canonical_microscopic.apertures.size()==1&&
		canonical_microscopic.apertures.front().fragment_a==50&&
		canonical_microscopic.apertures.front().fragment_b==60&&
		canonical_microscopic.open_area>0&&
		near(canonical_microscopic.open_area,microscopic_gap_target,
			0.35*microscopic_gap_target)&&
		near(canonical_microscopic.open_area+canonical_microscopic.blocked_area,1,1e-12),
		"canonical shared face retains an arithmetic-scale positive opening");

	const LocalArrangementSharedFaceBarrier full_barrier=
		shared_barrier(700,70,x_rectangle(1,0,1,0,1));
	const LocalArrangementSharedFaceAssembly canonical_full_barrier=
		assemble_canonical_local_shared_face(canonical_face,0,{}, {full_barrier},whole_a,whole_b);
	const LocalArrangementSharedFaceBarrier partial_barrier=
		shared_barrier(701,71,x_rectangle(1,0.25,0.75,0,1));
	const LocalArrangementSharedFaceAssembly canonical_partial_barrier=
		assemble_canonical_local_shared_face(canonical_face,0,{}, {partial_barrier},whole_a,whole_b);
	const bool full_barrier_sides=canonical_full_barrier.blocked_surfaces.size()==1&&
		canonical_full_barrier.blocked_surfaces.front().minus_fragment==50&&
		canonical_full_barrier.blocked_surfaces.front().plus_fragment==60;
	check(canonical_full_barrier.valid&&canonical_full_barrier.apertures.empty()&&
		near(canonical_full_barrier.blocked_area,1,1e-12)&&full_barrier_sides&&
		canonical_partial_barrier.valid&&near(canonical_partial_barrier.open_area,0.5,1e-12)&&
		near(canonical_partial_barrier.blocked_area,0.5,1e-12),
		"full and partial coplanar barriers block exactly their positive face area");

	// A labelled ownership cover is an interface contract.  Gaps and positive-area
	// overlaps are topology errors, not an invitation to probe a nearest surface.
	std::vector<LocalArrangementSharedFaceOwner> overlapping_a{
		shared_owner(x_rectangle(1,0,0.6,0,1),1),
		shared_owner(x_rectangle(1,0.4,1,0,1),2)};
	std::vector<LocalArrangementSharedFaceOwner> gapped_a{
		shared_owner(x_rectangle(1,0,0.4,0,1),1),
		shared_owner(x_rectangle(1,0.6,1,0,1),1)};
	const LocalArrangementSharedFaceAssembly rejected_overlap=
		assemble_canonical_local_shared_face(canonical_face,0,{}, {},overlapping_a,whole_b);
	const LocalArrangementSharedFaceAssembly rejected_gap=
		assemble_canonical_local_shared_face(canonical_face,0,{}, {},gapped_a,whole_b);
	check(!rejected_overlap.valid&&rejected_overlap.apertures.empty()&&
		rejected_overlap.blocked_surfaces.empty()&&!rejected_gap.valid&&
		rejected_gap.apertures.empty()&&rejected_gap.blocked_surfaces.empty(),
		"malformed overlapping/gapped ownership covers fail transactionally");

	// Reordering inputs, reversing barrier winding, and translating the entire face
	// cannot alter topology or measure.  Reversing winding swaps plus/minus and the
	// normal together, leaving the manufactured pressure force invariant.
	std::vector<LocalArrangementSharedFaceBarrier> invariance_barriers{
		shared_barrier(801,80,x_rectangle(1,0,0.25,0,1)),
		shared_barrier(802,80,x_rectangle(1,0.75,1,0,1))};
	std::vector<LocalArrangementSharedFaceOwner> invariance_a{
		shared_owner(x_rectangle(1,0,1,0,0.3),70),
		shared_owner(x_rectangle(1,0,1,0.3,1),70)};
	std::vector<LocalArrangementSharedFaceOwner> invariance_b{
		shared_owner(x_rectangle(1,0,1,0,0.7),80),
		shared_owner(x_rectangle(1,0,1,0.7,1),80)};
	std::vector<LocalArrangementSharedFaceTrace> invariance_traces{
		{900,90,{1,0,0.1},{1,1,0.9}}};
	const LocalArrangementSharedFaceAssembly invariance_reference=
		assemble_canonical_local_shared_face(canonical_face,0,invariance_traces,
			invariance_barriers,invariance_a,invariance_b);
	std::reverse(invariance_traces.begin(),invariance_traces.end());
	std::swap(invariance_traces.front().a,invariance_traces.front().b);
	std::reverse(invariance_barriers.begin(),invariance_barriers.end());
	std::reverse(invariance_a.begin(),invariance_a.end());
	std::reverse(invariance_b.begin(),invariance_b.end());
	const LocalArrangementSharedFaceAssembly invariance_reordered=
		assemble_canonical_local_shared_face(canonical_face,0,invariance_traces,
			invariance_barriers,invariance_a,invariance_b);
	auto translate=[](Vec3d point){return point+Vec3d{1.0e8,-2.0e8,3.0e8};};
	std::array<Vec3d,4> translated_face=canonical_face;
	for(Vec3d& point:translated_face)point=translate(point);
	for(auto& trace:invariance_traces)
	{
		trace.a=translate(trace.a);trace.b=translate(trace.b);
	}
	for(auto& barrier:invariance_barriers)
		barrier.polygon=transformed_polygon(barrier.polygon,translate);
	for(auto& owner:invariance_a)owner.polygon=transformed_polygon(owner.polygon,translate);
	for(auto& owner:invariance_b)owner.polygon=transformed_polygon(owner.polygon,translate);
	const LocalArrangementSharedFaceAssembly invariance_translated=
		assemble_canonical_local_shared_face(translated_face,0,invariance_traces,
			invariance_barriers,invariance_a,invariance_b);
	std::vector<LocalArrangementSharedFaceBarrier> reversed_barriers{
		shared_barrier(801,80,x_rectangle(1,0,0.25,0,1),{-1,0,0}),
		shared_barrier(802,80,x_rectangle(1,0.75,1,0,1),{-1,0,0})};
	const LocalArrangementSharedFaceAssembly invariance_winding=
		assemble_canonical_local_shared_face(canonical_face,0,
			{{900,90,{1,0,0.1},{1,1,0.9}}},reversed_barriers,
			{shared_owner(x_rectangle(1,0,1,0,0.3),70),
				shared_owner(x_rectangle(1,0,1,0.3,1),70)},
			{shared_owner(x_rectangle(1,0,1,0,0.7),80),
				shared_owner(x_rectangle(1,0,1,0.7,1),80)});
	auto blocked_force_x=[](const LocalArrangementSharedFaceAssembly& assembly)
	{
		double force=0;auto pressure=[](int fragment){return fragment==70?2.0:1.0;};
		for(const auto& surface:assembly.blocked_surfaces)
			force+=surface.normal.x*(pressure(surface.minus_fragment)-
				pressure(surface.plus_fragment))*surface.area;
		return force;
	};
	const bool invariant_measures=invariance_reference.valid&&invariance_reordered.valid&&
		invariance_translated.valid&&invariance_winding.valid&&
		near(invariance_reference.open_area,invariance_reordered.open_area,1e-12)&&
		near(invariance_reference.blocked_area,invariance_reordered.blocked_area,1e-12)&&
		near(invariance_reference.open_area,invariance_translated.open_area,1e-8)&&
		near(invariance_reference.blocked_area,invariance_translated.blocked_area,1e-8)&&
		near(blocked_force_x(invariance_reference),blocked_force_x(invariance_winding),1e-12);
	check(invariant_measures,
		"canonical shared face is invariant to input order, winding, and large translation");

	// Production shared-face assembly derives ownership directly from the two cells'
	// implicit plane/sign partitions.  A regular cell is the zero-plane special case.
	const LocalArrangementSharedFaceAssembly implicit_regular=
		assemble_canonical_local_shared_face(canonical_face,0,{},
			regular_side(100,-1),regular_side(200,1));
	check(implicit_regular.valid&&implicit_regular.apertures.size()==1&&
		implicit_regular.blocked_surfaces.empty()&&near(implicit_regular.open_area,1,1e-14)&&
		implicit_regular.apertures.front().fragment_a==100&&
		implicit_regular.apertures.front().fragment_b==200,
		"implicit shared face represents a regular full-face aperture with one empty sign row");

	LocalArrangementSharedFaceSide analytic_a;analytic_a.inward_axis_sign=-1;
	analytic_a.planes={shared_plane({0,1,0},0.4)};
	analytic_a.regions={shared_region({-1},101),shared_region({1},102)};
	const LocalArrangementSharedFaceAssembly implicit_analytic=
		assemble_canonical_local_shared_face(canonical_face,0,{},analytic_a,regular_side(200,1));
	bool analytic_labels=true;double analytic_minus_area=0,analytic_plus_area=0;
	for(const auto& aperture:implicit_analytic.apertures)
	{
		analytic_labels=analytic_labels&&aperture.fragment_b==200;
		if(aperture.fragment_a==101)analytic_minus_area+=aperture.area;
		else if(aperture.fragment_a==102)analytic_plus_area+=aperture.area;
		else analytic_labels=false;
	}
	check(implicit_analytic.valid&&implicit_analytic.apertures.size()==2&&analytic_labels&&
		near(analytic_minus_area,0.4,1e-14)&&near(analytic_plus_area,0.6,1e-14),
		"implicit shared face labels both sides of one analytic support plane");

	LocalArrangementSharedFaceSide implicit_n;implicit_n.inward_axis_sign=-1;
	implicit_n.planes={shared_plane({0,1,0},0.5)};
	implicit_n.regions={shared_region({-1},110),shared_region({1},111)};
	LocalArrangementSharedFaceSide implicit_m;implicit_m.inward_axis_sign=1;
	implicit_m.planes={shared_plane({0,0,1},0.25),shared_plane({0,0,1},0.7)};
	implicit_m.regions={shared_region({-1,-1},120),shared_region({1,-1},121),
		shared_region({1,1},122)};
	const LocalArrangementSharedFaceAssembly implicit_n_by_m=
		assemble_canonical_local_shared_face(canonical_face,0,{},implicit_n,implicit_m);
	bool n_by_m_labels=true;
	for(const auto& aperture:implicit_n_by_m.apertures)
		n_by_m_labels=n_by_m_labels&&(aperture.fragment_a==110||aperture.fragment_a==111)&&
			(aperture.fragment_b==120||aperture.fragment_b==121||aperture.fragment_b==122);
	check(implicit_n_by_m.valid&&implicit_n_by_m.apertures.size()==6&&n_by_m_labels&&
		near(implicit_n_by_m.open_area,1,1e-14),
		"implicit shared face common-refines differing N-by-M structural partitions");

	const LocalArrangementSharedFaceAssembly implicit_microscopic=
		assemble_canonical_local_shared_face(canonical_face,0,microscopic_barriers,
			regular_side(130,-1),regular_side(140,1));
	check(implicit_microscopic.valid&&implicit_microscopic.apertures.size()==1&&
		implicit_microscopic.open_area>0.0&&
		near(implicit_microscopic.open_area,microscopic_gap_target,
			0.35*microscopic_gap_target)&&
		near(implicit_microscopic.open_area+implicit_microscopic.blocked_area,1,1e-14),
		"implicit shared face preserves a 3.4911e-16 positive opening without an area cutoff");

	// A support plane coincident with the Cartesian face is classified from the
	// declared inward direction.  No epsilon-offset point is needed to choose the
	// adjacent fragment on either side of the barrier.
	LocalArrangementSharedFaceSide coplanar_a;coplanar_a.inward_axis_sign=-1;
	coplanar_a.planes={shared_plane({1,0,0},1)};
	coplanar_a.regions={shared_region({-1},150)};
	LocalArrangementSharedFaceSide coplanar_b;coplanar_b.inward_axis_sign=1;
	coplanar_b.planes={shared_plane({1,0,0},1)};
	coplanar_b.regions={shared_region({1},160)};
	const LocalArrangementSharedFaceAssembly implicit_coplanar_barrier=
		assemble_canonical_local_shared_face(canonical_face,0,{full_barrier},
			coplanar_a,coplanar_b);
	const bool coplanar_sides=implicit_coplanar_barrier.blocked_surfaces.size()==1&&
		implicit_coplanar_barrier.blocked_surfaces.front().minus_fragment==150&&
		implicit_coplanar_barrier.blocked_surfaces.front().plus_fragment==160;
	check(implicit_coplanar_barrier.valid&&implicit_coplanar_barrier.apertures.empty()&&
		near(implicit_coplanar_barrier.blocked_area,1,1e-14)&&coplanar_sides,
		"implicit face-coplanar plane uses inward-axis signs and preserves barrier plus/minus");

	// Reordering plane/region/barrier inputs only reorders the construction work.  A
	// translated plane updates its offset by n dot translation.  Reversing fabric
	// winding swaps plus/minus and normal together, leaving pressure force invariant.
	LocalArrangementSharedFaceSide invariant_side_a;invariant_side_a.inward_axis_sign=-1;
	invariant_side_a.planes={shared_plane({0,1,0},0.3),shared_plane({0,0,1},0.6)};
	invariant_side_a.regions={shared_region({-1,-1},170),shared_region({-1,1},171),
		shared_region({1,-1},172),shared_region({1,1},173)};
	LocalArrangementSharedFaceSide invariant_side_b;invariant_side_b.inward_axis_sign=1;
	invariant_side_b.planes={shared_plane({0,1,0},0.8)};
	invariant_side_b.regions={shared_region({-1},180),shared_region({1},181)};
	std::vector<LocalArrangementSharedFaceBarrier> structural_invariance_barriers{
		shared_barrier(801,80,x_rectangle(1,0,0.25,0,1)),
		shared_barrier(802,80,x_rectangle(1,0.75,1,0,1))};
	const LocalArrangementSharedFaceAssembly implicit_invariant_reference=
		assemble_canonical_local_shared_face(canonical_face,0,structural_invariance_barriers,
			invariant_side_a,invariant_side_b);
	LocalArrangementSharedFaceSide reordered_side_a=invariant_side_a;
	std::swap(reordered_side_a.planes[0],reordered_side_a.planes[1]);
	for(auto& region:reordered_side_a.regions)
		std::swap(region.plane_side[0],region.plane_side[1]);
	std::reverse(reordered_side_a.regions.begin(),reordered_side_a.regions.end());
	LocalArrangementSharedFaceSide reordered_side_b=invariant_side_b;
	std::reverse(reordered_side_b.regions.begin(),reordered_side_b.regions.end());
	std::vector<LocalArrangementSharedFaceBarrier> implicit_reordered_barriers=
		structural_invariance_barriers;
	std::reverse(implicit_reordered_barriers.begin(),implicit_reordered_barriers.end());
	const LocalArrangementSharedFaceAssembly implicit_invariant_reordered=
		assemble_canonical_local_shared_face(canonical_face,0,implicit_reordered_barriers,
			reordered_side_a,reordered_side_b);
	const Vec3d invariant_translation{1e8,-2e8,3e8};
	auto translate_partition=[&](LocalArrangementSharedFaceSide side)
	{
		for(auto& plane:side.planes)plane.offset+=dot(plane.normal,invariant_translation);
		return side;
	};
	std::array<Vec3d,4> implicit_translated_face=canonical_face;
	for(Vec3d& point:implicit_translated_face)point=point+invariant_translation;
	std::vector<LocalArrangementSharedFaceBarrier> implicit_translated_barriers=
		structural_invariance_barriers;
	for(auto& barrier:implicit_translated_barriers)
		for(Vec3d& point:barrier.polygon)point=point+invariant_translation;
	const LocalArrangementSharedFaceAssembly implicit_invariant_translated=
		assemble_canonical_local_shared_face(implicit_translated_face,0,
			implicit_translated_barriers,translate_partition(invariant_side_a),
			translate_partition(invariant_side_b));
	std::vector<LocalArrangementSharedFaceBarrier> implicit_reversed_barriers=
		structural_invariance_barriers;
	for(auto& barrier:implicit_reversed_barriers)barrier.normal=barrier.normal*-1.0;
	const LocalArrangementSharedFaceAssembly implicit_invariant_winding=
		assemble_canonical_local_shared_face(canonical_face,0,implicit_reversed_barriers,
			invariant_side_a,invariant_side_b);
	auto implicit_force=[](const LocalArrangementSharedFaceAssembly& assembly)
	{
		Vec3d force{};auto pressure=[](int fragment){return fragment<180?2.0:1.0;};
		for(const auto& surface:assembly.blocked_surfaces)
			force=force+surface.normal*((pressure(surface.minus_fragment)-
				pressure(surface.plus_fragment))*surface.area);
		return force;
	};
	const bool implicit_invariant=implicit_invariant_reference.valid&&
		implicit_invariant_reordered.valid&&implicit_invariant_translated.valid&&
		implicit_invariant_winding.valid&&
		near(implicit_invariant_reference.open_area,implicit_invariant_reordered.open_area,1e-14)&&
		near(implicit_invariant_reference.blocked_area,implicit_invariant_reordered.blocked_area,1e-14)&&
		near(implicit_invariant_reference.open_area,implicit_invariant_translated.open_area,1e-8)&&
		near(implicit_invariant_reference.blocked_area,implicit_invariant_translated.blocked_area,1e-8)&&
		near(implicit_force(implicit_invariant_reference).x,
			implicit_force(implicit_invariant_winding).x,1e-14);
	if(!implicit_invariant)
		std::printf("[INFO] implicit invariance ref=%d/%s %.17g/%.17g reorder=%d/%s %.17g/%.17g translated=%d/%s %.17g/%.17g winding=%d/%s force %.17g/%.17g\n",
			implicit_invariant_reference.valid,implicit_invariant_reference.error.c_str(),
			implicit_invariant_reference.open_area,implicit_invariant_reference.blocked_area,
			implicit_invariant_reordered.valid,implicit_invariant_reordered.error.c_str(),
			implicit_invariant_reordered.open_area,implicit_invariant_reordered.blocked_area,
			implicit_invariant_translated.valid,implicit_invariant_translated.error.c_str(),
			implicit_invariant_translated.open_area,implicit_invariant_translated.blocked_area,
			implicit_invariant_winding.valid,implicit_invariant_winding.error.c_str(),
			implicit_force(implicit_invariant_reference).x,
			implicit_force(implicit_invariant_winding).x);
	check(implicit_invariant,
		"implicit shared face is invariant to plane/region/barrier order, winding, and translation");

	LocalArrangementSharedFaceSide ambiguous_side=regular_side(190,-1);
	ambiguous_side.regions.push_back(shared_region({},191));
	LocalArrangementSharedFaceSide missing_side;missing_side.inward_axis_sign=-1;
	missing_side.planes={shared_plane({0,1,0},0.5)};
	missing_side.regions={shared_region({-1},192)};
	LocalArrangementSharedFaceSide wildcard_side;wildcard_side.inward_axis_sign=-1;
	wildcard_side.planes={shared_plane({0,1,0},0.5)};
	wildcard_side.regions={shared_region({0},193),shared_region({1},194)};
	const auto rejected_multiple=assemble_canonical_local_shared_face(canonical_face,0,{},
		ambiguous_side,regular_side(200,1));
	const auto rejected_missing=assemble_canonical_local_shared_face(canonical_face,0,{},
		missing_side,regular_side(200,1));
	const auto rejected_wildcard=assemble_canonical_local_shared_face(canonical_face,0,{},
		wildcard_side,regular_side(200,1));
	check(!rejected_multiple.valid&&!rejected_missing.valid&&!rejected_wildcard.valid&&
		rejected_multiple.apertures.empty()&&rejected_missing.apertures.empty()&&
		rejected_wildcard.apertures.empty(),
		"implicit shared face fails transactionally for missing, duplicate, or wildcard-ambiguous regions");

	std::printf("[local-arrangement] %s (%d failures)\n",failures?"FAIL":"PASS",failures);return failures?1:0;
}
