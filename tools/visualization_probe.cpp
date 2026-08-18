#include "gui/slice_field.h"

#include <cmath>
#include <cstdio>
#include <vector>

namespace
{
	bool near(float actual,float expected,float tolerance,const char* label)
	{
		if(std::abs(actual-expected)<=tolerance)return true;
		std::fprintf(stderr,"visualization probe: %s expected %.9g, got %.9g\n",label,expected,actual);return false;
	}
}

int main()
{
	using paracfd::core::MacGrid;using paracfd::gui::Field;
	MacGrid grid;grid.nx=grid.ny=grid.nz=4;grid.h=0.5;
	std::vector<double> u(grid.u_count()),v(grid.v_count()),w(grid.w_count()),p(grid.p_count());
	// Rigid rotation u=(-y,x,0) has |curl u|=2 and Q=1 everywhere. Pressure is a
	// separate affine field whose gradient magnitude is sqrt(1+4+9)=sqrt(14).
	for(int k=0;k<grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<=grid.nx;++i)
		u[grid.uidx(i,j,k)]=-(j+0.5)*grid.h;
	for(int k=0;k<grid.nz;++k)for(int j=0;j<=grid.ny;++j)for(int i=0;i<grid.nx;++i)
		v[grid.vidx(i,j,k)]=(i+0.5)*grid.h;
	for(int k=0;k<=grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<grid.nx;++i)
		w[grid.widx(i,j,k)]=0.0;
	for(int k=0;k<grid.nz;++k)for(int j=0;j<grid.ny;++j)for(int i=0;i<grid.nx;++i)
	{
		const double x=(i+0.5)*grid.h,y=(j+0.5)*grid.h,z=(k+0.5)*grid.h;
		p[grid.pidx(i,j,k)]=x+2*y+3*z;
	}
	float range[3]{};bool ok=true;
	paracfd::gui::slice_reduce_cpu(u.data(),v.data(),w.data(),p.data(),nullptr,grid,Field::VorticityMagnitude,range);
	ok&=near(range[0],2.0f,1e-5f,"vorticity min")&&near(range[1],2.0f,1e-5f,"vorticity max");
	paracfd::gui::slice_reduce_cpu(u.data(),v.data(),w.data(),p.data(),nullptr,grid,Field::QCriterion,range);
	ok&=near(range[0],1.0f,1e-5f,"Q min")&&near(range[1],1.0f,1e-5f,"Q max");
	paracfd::gui::slice_reduce_cpu(u.data(),v.data(),w.data(),p.data(),nullptr,grid,Field::PressureGradient,range);
	ok&=near(range[0],std::sqrt(14.0f),1e-5f,"pressure-gradient min")&&near(range[1],std::sqrt(14.0f),1e-5f,"pressure-gradient max");
	paracfd::gui::slice_reduce_cpu(u.data(),v.data(),w.data(),p.data(),nullptr,grid,Field::PressureCoefficient,range,2.0f,4.0f);
	ok&=near(range[0],1.5f/16.0f,1e-6f,"Cp min")&&near(range[1],10.5f/16.0f,1e-6f,"Cp max");
	if(ok)std::printf("visualization probe: derived scalar fields pass\n");return ok?0:1;
}
