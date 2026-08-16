#include "core/fluid/eb_pressure.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace paracfd::core
{
	namespace
	{
		void add_edge(std::vector<double>& out,const std::vector<double>& p,int a,int b,double c)
		{if(a<0||b<0||a==b)return;const double v=c*(p[a]-p[b]);out[a]+=v;out[b]-=v;}
		double aperture_distance(const EmbeddedBoundary& eb,const FaceAperture& a){return std::max(1e-12,std::sqrt(length2(eb.fragment_centroid(a.fragment_a)-eb.fragment_centroid(a.fragment_b))));}
		std::uint8_t face_mask(const EmbeddedBoundary& eb,int cell){return eb.cut_face_mask.empty()?0:eb.cut_face_mask[cell];}
	}

	int EbPressureSystem::dof(FragmentRef r)const{return fragment_is_regular(r)?cell_dof[regular_fragment_cell(r)]:fragment_dof[irregular_fragment_index(r)];}

	EbPressureSystem build_eb_pressure_system(const EmbeddedBoundary& eb,bool outlet)
	{
		EbPressureSystem s;s.eb=&eb;s.storage_size=eb.grid.cell_count()+static_cast<int>(eb.fragments.size());s.cell_dof.assign(eb.grid.cell_count(),-1);s.fragment_dof.assign(eb.fragments.size(),-1);s.volume.assign(s.storage_size,0);s.centroid.assign(s.storage_size,{});s.active.assign(s.storage_size,0);s.pressure_outlet_xmax=outlet;
		for(int c=0;c<eb.grid.cell_count();++c)if(eb.cells[c].state==EbCellState::regular){s.cell_dof[c]=c;s.active[c]=1;s.volume[c]=eb.grid.h*eb.grid.h*eb.grid.h;s.centroid[c]=eb.grid.cell_centroid(c);}
		for(int i=0;i<static_cast<int>(eb.fragments.size());++i)if(eb.fragments[i].merge_target==irregular_fragment(i)){int q=eb.grid.cell_count()+i;s.fragment_dof[i]=q;s.active[q]=eb.fragments[i].pressure_static?0:1;s.volume[q]=eb.fragments[i].volume;s.centroid[q]=eb.fragments[i].centroid;}
		// Resolve arbitrary acyclic merge chains. Every link was selected through an open
		// same-side aperture to a strictly larger fragment, so a cycle indicates corrupt
		// preprocessing rather than a case that may be silently joined.
		std::vector<unsigned char> resolving(eb.fragments.size(),0);
		std::function<int(int)> resolve=[&](int fragment)->int
		{
			if(s.fragment_dof[fragment]>=0)return s.fragment_dof[fragment];
			if(resolving[fragment])throw std::runtime_error("cyclic EB small-fragment merge chain");
			resolving[fragment]=1;FragmentRef target_ref=eb.fragments[fragment].merge_target;int target=-1;
			if(fragment_is_regular(target_ref))target=s.cell_dof[regular_fragment_cell(target_ref)];
			else{int next=irregular_fragment_index(target_ref);if(next<0||next>=static_cast<int>(eb.fragments.size()))throw std::runtime_error("invalid EB small-fragment merge target");target=resolve(next);}
			if(target<0)throw std::runtime_error("unresolved EB small-fragment merge target");resolving[fragment]=0;s.fragment_dof[fragment]=target;return target;
		};
		for(int i=0;i<static_cast<int>(eb.fragments.size());++i)resolve(i);
		for(int i=0;i<static_cast<int>(eb.fragments.size());++i)if(eb.fragments[i].merge_target!=irregular_fragment(i))
		{
			const int target=s.fragment_dof[i];const double old=s.volume[target],add=eb.fragments[i].volume;
			s.centroid[target]=(s.centroid[target]*old+eb.fragments[i].centroid*add)/(old+add);s.volume[target]+=add;
		}
		return s;
	}

	void EbPressureSystem::apply_cpu(const std::vector<double>& p,std::vector<double>& out)const
	{
		if(static_cast<int>(p.size())!=storage_size)throw std::invalid_argument("EB pressure vector size");out.assign(storage_size,0);const auto& e=*eb;const double regular_c=e.grid.h; // h^2 / h
		for(int k=0;k<e.grid.nz;++k)for(int j=0;j<e.grid.ny;++j)for(int i=0;i<e.grid.nx;++i)
		{
			int a=e.grid.cell_index(i,j,k);if(cell_dof[a]<0)continue;
			if(i+1<e.grid.nx&&!(face_mask(e,a)&1u)){int b=e.grid.cell_index(i+1,j,k);if(cell_dof[b]>=0)add_edge(out,p,cell_dof[a],cell_dof[b],regular_c);}
			if(j+1<e.grid.ny&&!(face_mask(e,a)&2u)){int b=e.grid.cell_index(i,j+1,k);if(cell_dof[b]>=0)add_edge(out,p,cell_dof[a],cell_dof[b],regular_c);}
			if(k+1<e.grid.nz&&!(face_mask(e,a)&4u)){int b=e.grid.cell_index(i,j,k+1);if(cell_dof[b]>=0)add_edge(out,p,cell_dof[a],cell_dof[b],regular_c);}
			if(pressure_outlet_xmax&&i==e.grid.nx-1)out[cell_dof[a]]+=(e.grid.h*e.grid.h/(0.5*e.grid.h))*p[cell_dof[a]];
		}
		for(const FaceAperture& a:e.apertures)add_edge(out,p,dof(a.fragment_a),dof(a.fragment_b),a.area/aperture_distance(e,a));
		// A split fragment touching xmax also receives pressure-reference coupling over its boundary
		// aperture in the full solver. Uniform preprocessing currently emits internal apertures only;
		// such outlet-cut geometries are rejected by domain-margin validation rather than guessed here.
	}

	void EbPressureSystem::diagonal_cpu(std::vector<double>& d)const
	{
		std::vector<double> z(storage_size,0),a;d.assign(storage_size,0);const auto& e=*eb;const double c=e.grid.h;
		for(int k=0;k<e.grid.nz;++k)for(int j=0;j<e.grid.ny;++j)for(int i=0;i<e.grid.nx;++i){int ca=e.grid.cell_index(i,j,k),da=cell_dof[ca];if(da<0)continue;auto edge=[&](int cb){int db=cell_dof[cb];if(db>=0&&db!=da){d[da]+=c;d[db]+=c;}};if(i+1<e.grid.nx&&!(face_mask(e,ca)&1u))edge(e.grid.cell_index(i+1,j,k));if(j+1<e.grid.ny&&!(face_mask(e,ca)&2u))edge(e.grid.cell_index(i,j+1,k));if(k+1<e.grid.nz&&!(face_mask(e,ca)&4u))edge(e.grid.cell_index(i,j,k+1));if(pressure_outlet_xmax&&i==e.grid.nx-1)d[da]+=2*e.grid.h;}
		for(const auto& ap:e.apertures){int da=dof(ap.fragment_a),db=dof(ap.fragment_b);if(da<0||db<0||da==db)continue;double q=ap.area/aperture_distance(e,ap);d[da]+=q;d[db]+=q;}
	}

	double EbPressureSystem::volume_weighted_mean(const std::vector<double>& x)const{double s=0,v=0;for(int i=0;i<storage_size;++i)if(active[i]){s+=x[i]*volume[i];v+=volume[i];}return v>0?s/v:0;}

	EbFaceFluxes make_zero_fluxes(const EmbeddedBoundary& e){EbFaceFluxes f;f.x.assign((e.grid.nx+1)*e.grid.ny*e.grid.nz,0);f.y.assign(e.grid.nx*(e.grid.ny+1)*e.grid.nz,0);f.z.assign(e.grid.nx*e.grid.ny*(e.grid.nz+1),0);f.aperture_velocity.assign(e.apertures.size(),0);return f;}

	void eb_divergence_cpu(const EbPressureSystem& s,const EbFaceFluxes& f,std::vector<double>& div)
	{
		const auto& e=*s.eb;div.assign(s.storage_size,0);auto add=[&](int a,int b,double q){if(a>=0)div[a]+=q;if(b>=0)div[b]-=q;};auto xi=[&](int i,int j,int k){return(k*e.grid.ny+j)*(e.grid.nx+1)+i;};auto yi=[&](int i,int j,int k){return(k*(e.grid.ny+1)+j)*e.grid.nx+i;};auto zi=[&](int i,int j,int k){return(k*e.grid.ny+j)*e.grid.nx+i;};const double A=e.grid.h*e.grid.h;
		for(int k=0;k<e.grid.nz;++k)for(int j=0;j<e.grid.ny;++j)for(int i=0;i<e.grid.nx;++i){int ca=e.grid.cell_index(i,j,k),da=s.cell_dof[ca];if(da<0)continue;if(i+1<e.grid.nx&&!(face_mask(e,ca)&1u)){int db=s.cell_dof[e.grid.cell_index(i+1,j,k)];if(db>=0)add(da,db,f.x[xi(i+1,j,k)]*A);}if(j+1<e.grid.ny&&!(face_mask(e,ca)&2u)){int db=s.cell_dof[e.grid.cell_index(i,j+1,k)];if(db>=0)add(da,db,f.y[yi(i,j+1,k)]*A);}if(k+1<e.grid.nz&&!(face_mask(e,ca)&4u)){int db=s.cell_dof[e.grid.cell_index(i,j,k+1)];if(db>=0)add(da,db,f.z[zi(i,j,k+1)]*A);}}
		for(std::size_t i=0;i<e.apertures.size();++i){const auto& a=e.apertures[i];add(s.dof(a.fragment_a),s.dof(a.fragment_b),f.aperture_velocity[i]*a.area);}
		for(int i=0;i<s.storage_size;++i)if(s.active[i]&&s.volume[i]>0)div[i]/=s.volume[i];
	}

	void eb_projection_rhs_cpu(const EbPressureSystem& s,const EbFaceFluxes& f,double rho,double dt,std::vector<double>& rhs){std::vector<double>d;eb_divergence_cpu(s,f,d);rhs.assign(s.storage_size,0);for(int i=0;i<s.storage_size;++i)if(s.active[i])rhs[i]=-(rho/dt)*d[i]*s.volume[i];}

	void eb_correct_fluxes_cpu(const EbPressureSystem& s,const std::vector<double>& p,double rho,double dt,EbFaceFluxes& f)
	{
		const auto& e=*s.eb;auto xi=[&](int i,int j,int k){return(k*e.grid.ny+j)*(e.grid.nx+1)+i;};auto yi=[&](int i,int j,int k){return(k*(e.grid.ny+1)+j)*e.grid.nx+i;};auto zi=[&](int i,int j,int k){return(k*e.grid.ny+j)*e.grid.nx+i;};const double fac=dt/rho/e.grid.h;
		for(int k=0;k<e.grid.nz;++k)for(int j=0;j<e.grid.ny;++j)for(int i=0;i<e.grid.nx;++i){int cell=e.grid.cell_index(i,j,k),a=s.cell_dof[cell];if(a<0)continue;if(i+1<e.grid.nx&&!(face_mask(e,cell)&1u)){int b=s.cell_dof[e.grid.cell_index(i+1,j,k)];if(b>=0)f.x[xi(i+1,j,k)]-=fac*(p[b]-p[a]);}if(j+1<e.grid.ny&&!(face_mask(e,cell)&2u)){int b=s.cell_dof[e.grid.cell_index(i,j+1,k)];if(b>=0)f.y[yi(i,j+1,k)]-=fac*(p[b]-p[a]);}if(k+1<e.grid.nz&&!(face_mask(e,cell)&4u)){int b=s.cell_dof[e.grid.cell_index(i,j,k+1)];if(b>=0)f.z[zi(i,j,k+1)]-=fac*(p[b]-p[a]);}}
		for(std::size_t i=0;i<e.apertures.size();++i){const auto&a=e.apertures[i];int da=s.dof(a.fragment_a),db=s.dof(a.fragment_b);if(da!=db)f.aperture_velocity[i]-=(dt/rho)*(p[db]-p[da])/aperture_distance(e,a);}
	}

	EbCpuSolveResult solve_eb_pressure_cpu(const EbPressureSystem& s,const std::vector<double>& rhs,std::vector<double>& p,double tol,int maxit)
	{
		const int n=s.storage_size;if((int)rhs.size()!=n)throw std::invalid_argument("EB rhs size");if((int)p.size()!=n)p.assign(n,0);std::vector<double> b=rhs,r,z,d,Ad,diag;s.diagonal_cpu(diag);
		if(!s.pressure_outlet_xmax)
		{
			// Integrated FV rows have constant nullspace, so compatibility is sum(rhs)=0.
			// Remove any roundoff incompatibility as a uniform divergence density: rhs -= c*V.
			double sum=0,vol=0;for(int i=0;i<n;++i)if(s.active[i]){sum+=b[i];vol+=s.volume[i];}
			const double c=vol>0?sum/vol:0;for(int i=0;i<n;++i)if(s.active[i])b[i]-=c*s.volume[i];
		}
		s.apply_cpu(p,Ad);r.resize(n);z.resize(n);d.resize(n);double b2=0;for(int i=0;i<n;++i)if(s.active[i]){r[i]=b[i]-Ad[i];z[i]=diag[i]>0?r[i]/diag[i]:0;d[i]=z[i];b2+=b[i]*b[i];}if(!(b2>0))return {0,0,true};double rz=0;for(int i=0;i<n;++i)rz+=r[i]*z[i];EbCpuSolveResult out;
		for(int it=0;it<maxit;++it){s.apply_cpu(d,Ad);double dAd=0;for(int i=0;i<n;++i)dAd+=d[i]*Ad[i];if(!(dAd>0))break;double a=rz/dAd,r2=0;for(int i=0;i<n;++i)if(s.active[i]){p[i]+=a*d[i];r[i]-=a*Ad[i];r2+=r[i]*r[i];}out.iterations=it+1;out.relative_residual=std::sqrt(r2/b2);if(out.relative_residual<=tol){out.converged=true;break;}double rz2=0;for(int i=0;i<n;++i)if(s.active[i]){z[i]=diag[i]>0?r[i]/diag[i]:0;rz2+=r[i]*z[i];}double beta=rz2/rz;for(int i=0;i<n;++i)d[i]=z[i]+beta*d[i];rz=rz2;}
		if(!s.pressure_outlet_xmax){double mean=s.volume_weighted_mean(p);for(int i=0;i<n;++i)if(s.active[i])p[i]-=mean;}return out;
	}
}
