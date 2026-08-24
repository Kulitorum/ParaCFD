#include "core/geometry/embedded_boundary.h"

#include <algorithm>
#include <cmath>

namespace paracfd::core
{
	Aabb3d UniformEbGrid::cell_box(int i,int j,int k) const
	{
		const Vec3d lo=origin+Vec3d{i*h,j*h,k*h};
		return {lo,lo+Vec3d{h,h,h}};
	}

	Vec3d UniformEbGrid::cell_centroid(int cell) const
	{
		const auto coordinate=cell_coord(cell);
		return origin+Vec3d{(coordinate[0]+0.5)*h,(coordinate[1]+0.5)*h,
			(coordinate[2]+0.5)*h};
	}

	namespace
	{
		Vec3d closest_point_on_triangle(Vec3d point,
			const EbOrientedBoundaryTriangle& triangle)
		{
			const Vec3d a=triangle.a,b=triangle.b,c=triangle.c;
			const Vec3d ab=b-a,ac=c-a,ap=point-a;
			const double d1=dot(ab,ap),d2=dot(ac,ap);
			if(d1<=0.0&&d2<=0.0)return a;
			const Vec3d bp=point-b;
			const double d3=dot(ab,bp),d4=dot(ac,bp);
			if(d3>=0.0&&d4<=d3)return b;
			const double vc=d1*d4-d3*d2;
			if(vc<=0.0&&d1>=0.0&&d3<=0.0)return a+ab*(d1/(d1-d3));
			const Vec3d cp=point-c;
			const double d5=dot(ab,cp),d6=dot(ac,cp);
			if(d6>=0.0&&d5<=d6)return c;
			const double vb=d5*d2-d1*d6;
			if(vb<=0.0&&d2>=0.0&&d6<=0.0)return a+ac*(d2/(d2-d6));
			const double va=d3*d6-d5*d4;
			if(va<=0.0&&(d4-d3)>=0.0&&(d5-d6)>=0.0)
				return b+(c-b)*((d4-d3)/((d4-d3)+(d5-d6)));
			const double denominator=1.0/(va+vb+vc);
			return a+ab*(vb*denominator)+ac*(vc*denominator);
		}

		double oriented_solid_angle(const EbOrientedBoundaryTriangle& triangle,Vec3d query)
		{
			const Vec3d a=triangle.a-query,b=triangle.b-query,c=triangle.c-query;
			const double la=std::sqrt(length2(a)),lb=std::sqrt(length2(b)),
				lc=std::sqrt(length2(c));
			const double numerator=dot(a,cross(b,c));
			const double denominator=la*lb*lc+dot(a,b)*lc+dot(b,c)*la+dot(c,a)*lb;
			return 2.0*std::atan2(numerator,denominator);
		}

		int locate_exact_fragment(const ExactEbCellLocator& locator,Vec3d point,
			double requested_tolerance)
		{
			const double tolerance=std::max(requested_tolerance,locator.ambiguity_tolerance);
			bool have_convex_regions=false;int convex_containing=-1;
			for(int fragment=0;fragment<static_cast<int>(locator.fragments.size());++fragment)
			{
				for(const EbConvexRegion& region:locator.fragments[fragment].convex_regions)
				{
					have_convex_regions=true;bool inside=true;
					for(const EbConvexHalfspace& halfspace:region.halfspaces)
						if(dot(halfspace.outward_normal,point)>halfspace.offset+tolerance)
						{inside=false;break;}
					if(!inside)continue;
					if(convex_containing>=0&&convex_containing!=fragment)return -1;
					convex_containing=fragment;break;
				}
			}
			if(have_convex_regions)return convex_containing;
			const double tolerance2=tolerance*tolerance;
			const double winding_tolerance=std::max(1.0e-6,locator.winding_tolerance);
			int containing=-1;
			for(int fragment=0;fragment<static_cast<int>(locator.fragments.size());++fragment)
			{
				double sum=0.0,compensation=0.0;
				for(const EbOrientedBoundaryTriangle& triangle:
					locator.fragments[fragment].triangles)
				{
					if(length2(point-closest_point_on_triangle(point,triangle))<=tolerance2)
						return -1;
					const double value=oriented_solid_angle(triangle,point)-compensation;
					const double next=sum+value;
					compensation=(next-sum)-value;
					sum=next;
				}
				const double magnitude=std::abs(sum/(4.0*std::acos(-1.0)));
				if(magnitude<=winding_tolerance)continue;
				if(std::abs(magnitude-1.0)>winding_tolerance||containing>=0)return -1;
				containing=fragment;
			}
			return containing;
		}
	}

	FragmentRef EmbeddedBoundary::fragment_containing_point(int cell,Vec3d point,
		double tolerance) const
	{
		if(cell<0||cell>=static_cast<int>(cells.size()))return invalid_fragment;
		const EbCellTopology& topology=cells[cell];
		if(topology.state==EbCellState::regular)return regular_fragment(cell);
		if(topology.state!=EbCellState::split||topology.exact_locator_index<0||
			topology.exact_locator_index>=static_cast<int>(exact_locators.size()))
			return invalid_fragment;
		const int local=locate_exact_fragment(exact_locators[topology.exact_locator_index],
			point,tolerance);
		return local>=0&&local<topology.fragment_count
			?irregular_fragment(topology.first_fragment+local):invalid_fragment;
	}

	Vec3d EmbeddedBoundary::fragment_centroid(FragmentRef fragment) const
	{
		return fragment_is_regular(fragment)
			?grid.cell_centroid(regular_fragment_cell(fragment))
			:fragments[irregular_fragment_index(fragment)].centroid;
	}

	double EmbeddedBoundary::fragment_volume(FragmentRef fragment) const
	{
		return fragment_is_regular(fragment)?grid.h*grid.h*grid.h
			:fragments[irregular_fragment_index(fragment)].volume;
	}
}
