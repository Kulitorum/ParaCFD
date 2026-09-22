// Host-side composite AMR sampling for visualization and diagnostics.
//
// This mirrors the staggered trilinear reconstruction used by amr_advection.cu.
// In particular, interpolation stencils that cross a brick or 2:1 boundary are
// completed from the neighbouring active brick instead of clamping to the source
// brick.  A globally affine MAC field therefore remains affine across AMR levels.
#pragma once

#include "core/fluid/amr_fields.h"

#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace paracfd::core
{
	namespace amr_sampling_detail
	{
		inline bool valid_face_index(int component,int brick_size,int i,int j,int k)
		{
			return i>=0&&j>=0&&k>=0&&i<(component==0?brick_size+1:brick_size)&&
				j<(component==1?brick_size+1:brick_size)&&
				k<(component==2?brick_size+1:brick_size);
		}

		inline double value_at_clamped(const AmrHostLevelFields& level,int component,
			int brick,int i,int j,int k)
		{
			const int bs=level.layout.brick_size;
			if(component==0)
			{
				i=std::clamp(i,0,bs);j=std::clamp(j,0,bs-1);k=std::clamp(k,0,bs-1);
				return static_cast<double>(level.u[level.layout.u_index(brick,i,j,k)]);
			}
			if(component==1)
			{
				i=std::clamp(i,0,bs-1);j=std::clamp(j,0,bs);k=std::clamp(k,0,bs-1);
				return static_cast<double>(level.v[level.layout.v_index(brick,i,j,k)]);
			}
			i=std::clamp(i,0,bs-1);j=std::clamp(j,0,bs-1);k=std::clamp(k,0,bs);
			return static_cast<double>(level.w[level.layout.w_index(brick,i,j,k)]);
		}

		inline double value_at_extrapolated(const AmrHostLevelFields& level,int component,
			int brick,int i,int j,int k)
		{
			const int bs=level.layout.brick_size;
			const int n[3]={component==0?bs+1:bs,component==1?bs+1:bs,
				component==2?bs+1:bs};
			const int requested[3]={i,j,k};int index[3][2]{},count[3]{};
			double weight[3][2]{};
			for(int axis=0;axis<3;++axis)
			{
				if(requested[axis]<0)
				{
					index[axis][0]=0;index[axis][1]=1;
					weight[axis][0]=2;weight[axis][1]=-1;count[axis]=2;
				}
				else if(requested[axis]>=n[axis])
				{
					index[axis][0]=n[axis]-1;index[axis][1]=n[axis]-2;
					weight[axis][0]=2;weight[axis][1]=-1;count[axis]=2;
				}
				else
				{
					index[axis][0]=requested[axis];weight[axis][0]=1;count[axis]=1;
				}
			}
			double result=0;
			for(int az=0;az<count[2];++az)for(int ay=0;ay<count[1];++ay)
				for(int ax=0;ax<count[0];++ax)
					result+=weight[0][ax]*weight[1][ay]*weight[2][az]*
						value_at_clamped(level,component,brick,index[0][ax],index[1][ay],index[2][az]);
			return result;
		}

		inline Vec3d face_node_point(const BrickMetadata& record,int component,int i,int j,int k)
		{
			return {record.origin.x+(i+(component==0?0.0:0.5))*record.h,
				record.origin.y+(j+(component==1?0.0:0.5))*record.h,
				record.origin.z+(k+(component==2?0.0:0.5))*record.h};
		}

		inline double sample_local_extrapolated(const AmrHostLevelFields& level,
			const BrickMetadata& record,int component,int brick,Vec3d point)
		{
			double x=(point.x-record.origin.x)/record.h-(component==0?0.0:0.5);
			double y=(point.y-record.origin.y)/record.h-(component==1?0.0:0.5);
			double z=(point.z-record.origin.z)/record.h-(component==2?0.0:0.5);
			const int i0=static_cast<int>(std::floor(x));
			const int j0=static_cast<int>(std::floor(y));
			const int k0=static_cast<int>(std::floor(z));
			x-=i0;y-=j0;z-=k0;double result=0;
			for(int dz=0;dz<2;++dz)for(int dy=0;dy<2;++dy)for(int dx=0;dx<2;++dx)
			{
				const double wx=dx?x:1-x,wy=dy?y:1-y,wz=dz?z:1-z;
				result+=wx*wy*wz*value_at_extrapolated(level,component,brick,
					i0+dx,j0+dy,k0+dz);
			}
			return result;
		}

		inline int nearest_lattice_index(double coordinate,double origin,double h,double offset)
		{
			return static_cast<int>(std::floor((coordinate-origin)/h-offset+0.5));
		}

		inline bool refined_normal_face_value(const AmrHierarchy& hierarchy,
			const AmrHostFields& fields,int source_level,const BrickMetadata& source_record,
			int component,int source_brick,int i,int j,int k,double& value)
		{
			const int bs=fields.levels()[source_level].layout.brick_size;
			const int coordinate=component==0?i:(component==1?j:k);
			if(coordinate!=0&&coordinate!=bs)return false;
			const int face=2*component+(coordinate==bs);
			const int neighbour=source_record.same_level_neighbor[face];
			if(neighbour>=0&&hierarchy.levels()[source_level].bricks[neighbour].active())
				return false;
			const Vec3d node=face_node_point(source_record,component,i,j,k);
			Vec3d probe=node;
			probe[component]+=(coordinate==0?-1:1)*1e-5*source_record.h;
			const BrickLocation location=hierarchy.locate_finest(probe);
			if(!location.found()||location.level<=source_level)return false;
			const BrickMetadata& target_record=
				hierarchy.levels()[location.level].bricks[location.brick];
			value=sample_local_extrapolated(fields.levels()[location.level],target_record,
				component,location.brick,node);
			return true;
		}

		inline bool refined_interface_sample(const AmrHierarchy& hierarchy,
			const AmrHostFields& fields,int source_level,const BrickMetadata& source_record,
			int component,int face,Vec3d point,double& value)
		{
			const int neighbour=source_record.same_level_neighbor[face];
			if(neighbour>=0&&hierarchy.levels()[source_level].bricks[neighbour].active())
				return false;
			const int direction=(face&1)?1:-1;
			point[component]=source_record.origin[component]+
				((face&1)?hierarchy.brick_size()*source_record.h:0.0);
			Vec3d probe=point;probe[component]+=direction*1e-5*source_record.h;
			const BrickLocation location=hierarchy.locate_finest(probe);
			if(!location.found()||location.level<=source_level)return false;
			const BrickMetadata& target_record=
				hierarchy.levels()[location.level].bricks[location.brick];
			value=sample_local_extrapolated(fields.levels()[location.level],target_record,
				component,location.brick,point);
			return true;
		}

		inline double value_at_node(const AmrHierarchy& hierarchy,const AmrHostFields& fields,
			int source_level,const BrickMetadata& source_record,int component,int source_brick,
			int i,int j,int k)
		{
			const AmrHostLevelFields& source=fields.levels()[source_level];
			const int bs=source.layout.brick_size;
			if(valid_face_index(component,bs,i,j,k))
			{
				double refined=0;
				return refined_normal_face_value(hierarchy,fields,source_level,source_record,
					component,source_brick,i,j,k,refined)
					?refined:value_at_clamped(source,component,source_brick,i,j,k);
			}

			const int nx=component==0?bs+1:bs,ny=component==1?bs+1:bs,
				nz=component==2?bs+1:bs;
			int face=-1,outside=0;
			if(i<0){face=0;++outside;}else if(i>=nx){face=1;++outside;}
			if(j<0){face=2;++outside;}else if(j>=ny){face=3;++outside;}
			if(k<0){face=4;++outside;}else if(k>=nz){face=5;++outside;}
			if(outside==1)
			{
				const int neighbour=source_record.same_level_neighbor[face];
				if(neighbour>=0&&hierarchy.levels()[source_level].bricks[neighbour].active())
				{
					if(face==0)i+=bs;else if(face==1)i-=bs;
					else if(face==2)j+=bs;else if(face==3)j-=bs;
					else if(face==4)k+=bs;else k-=bs;
					return value_at_clamped(source,component,neighbour,i,j,k);
				}
			}

			const Vec3d node=face_node_point(source_record,component,i,j,k);
			const Aabb3d& domain=hierarchy.domain();
			if(node.x<domain.lo.x||node.y<domain.lo.y||node.z<domain.lo.z||
				node.x>domain.hi.x||node.y>domain.hi.y||node.z>domain.hi.z)
				return value_at_clamped(source,component,source_brick,i,j,k);
			const double epsilon=1e-5*source_record.h;
			Vec3d probe{std::clamp(node.x,domain.lo.x+epsilon,domain.hi.x-epsilon),
				std::clamp(node.y,domain.lo.y+epsilon,domain.hi.y-epsilon),
				std::clamp(node.z,domain.lo.z+epsilon,domain.hi.z-epsilon)};
			const BrickLocation location=hierarchy.locate_finest(probe);
			if(!location.found())return value_at_clamped(source,component,source_brick,i,j,k);
			const BrickMetadata& record=hierarchy.levels()[location.level].bricks[location.brick];
			const AmrHostLevelFields& target=fields.levels()[location.level];
			if(location.level!=source_level)
				return sample_local_extrapolated(target,record,component,location.brick,node);
			const int ni=nearest_lattice_index(node.x,record.origin.x,record.h,component==0?0.0:0.5);
			const int nj=nearest_lattice_index(node.y,record.origin.y,record.h,component==1?0.0:0.5);
			const int nk=nearest_lattice_index(node.z,record.origin.z,record.h,component==2?0.0:0.5);
			return value_at_clamped(target,component,location.brick,ni,nj,nk);
		}

		inline double sample_local(const AmrHierarchy& hierarchy,const AmrHostFields& fields,
			int level,const BrickMetadata& record,int component,int brick,Vec3d point)
		{
			double x=(point.x-record.origin.x)/record.h-(component==0?0.0:0.5);
			double y=(point.y-record.origin.y)/record.h-(component==1?0.0:0.5);
			double z=(point.z-record.origin.z)/record.h-(component==2?0.0:0.5);
			const int i0=static_cast<int>(std::floor(x));
			const int j0=static_cast<int>(std::floor(y));
			const int k0=static_cast<int>(std::floor(z));
			x-=i0;y-=j0;z-=k0;double result=0,boundary_contribution=0;
			const int bs=fields.levels()[level].layout.brick_size;
			const double normal=(point[component]-record.origin[component])/record.h;
			int transition_face=-1,boundary_index=-1;
			if(normal>=0&&normal<1){transition_face=2*component;boundary_index=0;}
			else if(normal>bs-1&&normal<=bs){transition_face=2*component+1;boundary_index=bs;}
			for(int dz=0;dz<2;++dz)for(int dy=0;dy<2;++dy)for(int dx=0;dx<2;++dx)
			{
				const int ni=i0+dx,nj=j0+dy,nk=k0+dz;
				const double wx=dx?x:1-x,wy=dy?y:1-y,wz=dz?z:1-z;
				double value=0;
				if(valid_face_index(component,bs,ni,nj,nk))
				{
					if(!refined_normal_face_value(hierarchy,fields,level,record,component,
						brick,ni,nj,nk,value))
						value=value_at_clamped(fields.levels()[level],component,brick,ni,nj,nk);
				}
				else value=value_at_node(hierarchy,fields,level,record,component,brick,ni,nj,nk);
				result+=wx*wy*wz*value;
				const int node_normal=component==0?ni:(component==1?nj:nk);
				if(node_normal==boundary_index)boundary_contribution+=wx*wy*wz*value;
			}
			if(transition_face>=0)
			{
				const double boundary_weight=transition_face&1?normal-(bs-1):1-normal;
				double exact_interface=0;
				if(boundary_weight>1e-14&&refined_interface_sample(hierarchy,fields,level,
					record,component,transition_face,point,exact_interface))
					result+=boundary_weight*(exact_interface-boundary_contribution/boundary_weight);
			}
			return result;
		}
	}

	inline double sample_amr_face_component(const AmrHierarchy& hierarchy,
		const AmrHostFields& fields,Vec3d point,int component,double fallback=0.0)
	{
		if(component<0||component>2||fields.levels().size()!=hierarchy.levels().size())
			return fallback;
		const Aabb3d& domain=hierarchy.domain();
		point.x=std::clamp(point.x,domain.lo.x,std::nextafter(domain.hi.x,domain.lo.x));
		point.y=std::clamp(point.y,domain.lo.y,std::nextafter(domain.hi.y,domain.lo.y));
		point.z=std::clamp(point.z,domain.lo.z,std::nextafter(domain.hi.z,domain.lo.z));
		const BrickLocation location=hierarchy.locate_finest(point);
		if(!location.found())return fallback;
		return amr_sampling_detail::sample_local(hierarchy,fields,location.level,
			hierarchy.levels()[location.level].bricks[location.brick],component,
			location.brick,point);
	}

	inline Vec3d sample_amr_velocity(const AmrHierarchy& hierarchy,
		const AmrHostFields& fields,Vec3d point,Vec3d fallback={})
	{
		if(fields.levels().size()!=hierarchy.levels().size())return fallback;
		const Aabb3d& domain=hierarchy.domain();
		point.x=std::clamp(point.x,domain.lo.x,std::nextafter(domain.hi.x,domain.lo.x));
		point.y=std::clamp(point.y,domain.lo.y,std::nextafter(domain.hi.y,domain.lo.y));
		point.z=std::clamp(point.z,domain.lo.z,std::nextafter(domain.hi.z,domain.lo.z));
		const BrickLocation location=hierarchy.locate_finest(point);
		if(!location.found())return fallback;
		const BrickMetadata& record=hierarchy.levels()[location.level].bricks[location.brick];
		return {amr_sampling_detail::sample_local(hierarchy,fields,location.level,record,0,location.brick,point),
			amr_sampling_detail::sample_local(hierarchy,fields,location.level,record,1,location.brick,point),
			amr_sampling_detail::sample_local(hierarchy,fields,location.level,record,2,location.brick,point)};
	}
    // A geometrically closed contour must be sampled at its actual positions.
    // Nearest-cell values move its four sides independently on an AMR grid and
    // create spurious circulation even for a curl-free affine velocity field.
    inline double rectangular_circulation_xz(const AmrHierarchy& hierarchy,
        const AmrHostFields& fields,double x0,double x1,double y,double z0,double z1,
        int samples=512)
    {
        if(!(x1>x0&&z1>z0)||samples<1)
            throw std::invalid_argument("invalid circulation contour");
        const double dx=(x1-x0)/samples,dz=(z1-z0)/samples;
        double circulation=0;
        for(int q=0;q<samples;++q)
        {
            const double x=x0+(q+0.5)*dx,z=z0+(q+0.5)*dz;
            circulation+=(sample_amr_velocity(hierarchy,fields,{x,y,z0}).x-
                sample_amr_velocity(hierarchy,fields,{x,y,z1}).x)*dx;
            circulation+=(sample_amr_velocity(hierarchy,fields,{x1,y,z}).z-
                sample_amr_velocity(hierarchy,fields,{x0,y,z}).z)*dz;
        }
        return circulation;
    }

}
