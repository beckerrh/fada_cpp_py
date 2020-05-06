#include <math.h>
#include "operator.hpp"

/*-------------------------------------------------*/
void Operator::dot(vector& out, const vector& in, double d) const
{
  // Laplacien   elements finis q1  (9-point-stencil)
  int dim = out.dim();
  if(dim==2)
  {
    double d0 = 8.0/3.0 * d;
    double d1 = -1.0/3.0 * d;
    int nx = out.n(0), ny = out.n(1);
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        out(ix,iy) += d0 * in(ix,iy)
                         + d1* in(ix-1,iy  ) + d1* in(ix+1,iy)
                         + d1* in(ix  ,iy-1) + d1* in(ix  ,iy+1)
                         + d1* in(ix-1,iy-1) + d1* in(ix-1,iy+1)
                         + d1* in(ix+1,iy-1) + d1* in(ix+1,iy+1);
      }
    }
  }
  else if(dim==3)
  {
    int nx = out.n(0), ny = out.n(1), nz = out.n(2);
    double e = d/arma::mean(out.n());
    double d0 = 8.0/3.0 * e;
    double d1 = -0.0/3.0 * e;
    double d2 = -1.0/6.0 * e;
    double d3 = -1.0/12.0 * e;
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        for(int iz=1;iz<nz-1;iz++)
        {
        out(ix,iy,iz) +=
                      d0* in(ix,iy,iz)
          
                      +d1* in(ix-1,iy  ,iz  )
                      +d1* in(ix+1,iy  ,iz  )
                      +d1* in(ix  ,iy-1,iz  )
                      +d1* in(ix  ,iy+1,iz  )
                      +d1* in(ix  ,iy  ,iz-1)
                      +d1* in(ix+1,iy  ,iz+1)
          
                      +d2* in(ix-1,iy-1,iz  )
                      +d2* in(ix-1,iy+1,iz  )
                      +d2* in(ix+1,iy-1,iz  )
                      +d2* in(ix+1,iy+1,iz  )
                      +d2* in(ix-1,iy  ,iz-1)
                      +d2* in(ix-1,iy  ,iz+1)
                      +d2* in(ix+1,iy  ,iz-1)
                      +d2* in(ix+1,iy  ,iz+1)
                      +d2* in(ix  ,iy-1,iz-1)
                      +d2* in(ix  ,iy-1,iz+1)
                      +d2* in(ix  ,iy+1,iz-1)
                      +d2* in(ix  ,iy+1,iz+1)
          
                      +d3* in(ix-1,iy-1,iz-1)
                      +d3* in(ix-1,iy-1,iz+1)
                      +d3* in(ix-1,iy+1,iz-1)
                      +d3* in(ix-1,iy+1,iz+1)
                      +d3* in(ix+1,iy-1,iz-1)
                      +d3* in(ix+1,iy-1,iz+1)
                      +d3* in(ix+1,iy+1,iz-1)
                      +d3* in(ix+1,iy+1,iz+1)
                          ;
        }
      }
    }
  }
  // Conditions aux limites ( Dirichlet)
  
  //  out.boundary(in, d);
}


/*-------------------------------------------------*/
void Operator::jacobi(vector& out, const vector& in) const
{
  int dim = out.dim();
  if(dim==2)
  {
    double omega = 0.8;
    double d0inv = 3.0/8.0 * omega;
    int nx = out.n(0), ny = out.n(1);
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        out(ix,iy) = d0inv * in(ix,iy);
      }
    }
  }
  else if(dim==3)
  {
    int nx = out.n(0), ny = out.n(1), nz = out.n(2);
    double omega = 0.8;
    double d0inv = 3.0/8.0 * arma::mean(out.n()) * omega;
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        for(int iz=1;iz<nz-1;iz++)
        {
          out(ix,iy,iz) = d0inv * in(ix,iy,iz);
        }
      }
    }
  }
  else
  {
    assert(0);
  }
}

/*-------------------------------------------------*/
void Operator::gauss_seidel1(vector& out, const vector& in) const
{
  int dim = out.dim();
  if(dim==2)
  {
    double omega = 0.8;
    double d0inv = 3.0/8.0;
    double d1 = -1.0/3.0;
/*
     (ix+p)*ny + iy+q < ix*ny + iy
     p*ny +q < 0
     p=-1 q=-1,0,1
     p= 0 q=-1
     */
    int nx = out.n(0), ny = out.n(1);
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        out(ix,iy) = d0inv * (
                              in(ix,iy)
                              -d1* out(ix-1,iy-1)
                              -d1* out(ix-1,iy  )
                              -d1* out(ix-1,iy+1)
                              -d1* out(ix  ,iy-1)
                              );
      }
    }
  }
  else
  {
    /*
     (ix+p)*ny*nz + (iy+q)*nz + iz+r < ix*ny*nz + iy*nz + iz
     p*ny*nz +q*nz +r < 0
     p=-1 q=-1,0,1  r=-1,0,1
     p= 0 q=-1 r=-1,0,1 q=0 r=-1
     */
    int nx = out.n(0), ny = out.n(1), nz = out.n(2);
    double e = 1.0/arma::mean(out.n());
    double d0 = 8.0/3.0 * e;
    double d1 = -0.0/3.0 * e;
    double d2 = -1.0/6.0 * e;
    double d3 = -1.0/12.0 * e;
    double d0inv = 1.0/d0;
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        for(int iz=1;iz<nz-1;iz++)
        {
          out(ix,iy,iz) = d0inv*(
                             in(ix,iy,iz)
                           -d1* out(ix-1,iy  ,iz  )
                           -d1* out(ix  ,iy-1,iz  )
                           -d1* out(ix  ,iy  ,iz-1)
                           -d2* out(ix-1,iy-1,iz  )
                           -d2* out(ix-1,iy+1,iz  )
                           -d2* out(ix-1,iy  ,iz-1)
                           -d2* out(ix-1,iy  ,iz+1)
                           -d2* out(ix  ,iy-1,iz-1)
                           -d2* out(ix  ,iy-1,iz+1)
                           -d3* out(ix-1,iy-1,iz-1)
                           -d3* out(ix-1,iy-1,iz+1)
                           -d3* out(ix-1,iy+1,iz-1)
                           -d3* out(ix-1,iy+1,iz+1)
                              );
        }
      }
    }
  }
}

/*-------------------------------------------------*/
void Operator::gauss_seidel2(vector& out, const vector& in) const
{
  int dim = out.dim();
  if(dim==2)
  {
    int nx = out.n(0), ny = out.n(1);
    double omega = 0.8;
    double d0inv = 3.0/8.0;
    double d1 = -1.0/3.0;
    for(int ix=nx-2;ix>=1;ix--)
    {
      for(int iy=ny-2;iy>=1;iy--)
      {
        out(ix,iy) = d0inv * (
                            in(ix,iy)
                          -d1* out(ix+1,iy-1)
                          -d1* out(ix+1,iy  )
                          -d1* out(ix+1,iy+1)
                          -d1* out(ix  ,iy+1)
                              );
      }
    }
  }
  else
  {
    int nx = out.n(0), ny = out.n(1), nz = out.n(2);
    double e = 1.0/arma::mean(out.n());
    double d0 = 8.0/3.0 * e;
    double d1 = -0.0/3.0 * e;
    double d2 = -1.0/6.0 * e;
    double d3 = -1.0/12.0 * e;
    double d0inv = 1.0/d0;
    for(int ix=nx-2;ix>=1;ix--)
    {
      for(int iy=ny-2;iy>=1;iy--)
      {
        for(int iz=nz-2;iz>=1;iz--)
        {
          out(ix,iy,iz) = d0inv*(
                             in(ix,iy,iz)+
                           -d1* out(ix+1,iy  ,iz  )
                           -d1* out(ix  ,iy+1,iz  )
                           -d2* out(ix+1,iy  ,iz+1)
                           -d2* out(ix+1,iy-1,iz  )
                           -d2* out(ix+1,iy+1,iz  )
                           -d2* out(ix+1,iy  ,iz-1)
                           -d2* out(ix+1,iy  ,iz+1)
                           -d2* out(ix  ,iy+1,iz-1)
                           -d2* out(ix  ,iy+1,iz+1)
                           -d3* out(ix+1,iy-1,iz-1)
                           -d3* out(ix+1,iy-1,iz+1)
                           -d3* out(ix+1,iy+1,iz-1)
                           -d3* out(ix+1,iy+1,iz+1)
                          );
        }
      }
    }
  }
}
