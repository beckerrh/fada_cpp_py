#include <math.h>
#include "operator.hpp"

/*-------------------------------------------------*/
void Operator::_boundary(int l, Vector& out) const
{
  if(dim()==2)
  {
    int nx = _n(0,l), ny = _n(1,l);
    for(int ix=0;ix<nx;ix++)
    {
      out.atp(ix,0)    = 0.0;
      out.atp(ix,ny-1) = 0.0;
    }
    for(int iy=0;iy<ny;iy++)
    {
      out.atp(0,iy)    = 0.0;
      out.atp(nx-1,iy) = 0.0;
    }
  }
  else if(dim()==3)
  {
    int nx = _n(0,l), ny = _n(1,l), nz = _n(2,l);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        out.atp(ix,iy,0)    = 0.0;
        out.atp(ix,iy,nz-1) = 0.0;
      }
    }
    for(int ix=0;ix<nx;ix++)
    {
      for(int iz=0;iz<nz;iz++)
      {
        out.atp(ix,0,   iz) = 0.0;
        out.atp(ix,ny-1,iz) = 0.0;
      }
    }
    for(int iy=0;iy<ny;iy++)
    {
      for(int iz=0;iz<nz;iz++)
      {
        out.atp(0,   iy,iz) = 0.0;
        out.atp(nx-1,iy,iz) = 0.0;
      }
    }
  }
}

/*-------------------------------------------------*/
void Operator::dot(int l, Vector& out, const Vector& in, double d) const
{
  // Laplacien   elements finis q1  (9-point-stencil)
  if(dim()==2)
  {
    double d0 = 8.0/3.0 * d;
    double d1 = -1.0/3.0 * d;
    int nx = _n(0,l), ny = _n(1,l);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        out.atp(ix,iy) += d0 * in.atp(ix,iy)
        + d1* in.atp(ix-1,iy  ) + d1* in.atp(ix+1,iy)
        + d1* in.atp(ix  ,iy-1) + d1* in.atp(ix  ,iy+1)
        + d1* in.atp(ix-1,iy-1) + d1* in.atp(ix-1,iy+1)
        + d1* in.atp(ix+1,iy-1) + d1* in.atp(ix+1,iy+1);
      }
    }
  }
  else if(dim()==3)
  {
    int nx = _n(0,l), ny = _n(1,l), nz = _n(2,l);
    double e = d/arma::mean(out.n());
    double d0 = 8.0/3.0 * e;
    double d1 = -0.0/3.0 * e;
    double d2 = -1.0/6.0 * e;
    double d3 = -1.0/12.0 * e;
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<nz;iz++)
        {
          out.atp(ix,iy,iz) +=
          d0* in.atp(ix,iy,iz)
          
          +d1* in.atp(ix-1,iy  ,iz  )
          +d1* in.atp(ix+1,iy  ,iz  )
          +d1* in.atp(ix  ,iy-1,iz  )
          +d1* in.atp(ix  ,iy+1,iz  )
          +d1* in.atp(ix  ,iy  ,iz-1)
          +d1* in.atp(ix+1,iy  ,iz+1)
          
          +d2* in.atp(ix-1,iy-1,iz  )
          +d2* in.atp(ix-1,iy+1,iz  )
          +d2* in.atp(ix+1,iy-1,iz  )
          +d2* in.atp(ix+1,iy+1,iz  )
          +d2* in.atp(ix-1,iy  ,iz-1)
          +d2* in.atp(ix-1,iy  ,iz+1)
          +d2* in.atp(ix+1,iy  ,iz-1)
          +d2* in.atp(ix+1,iy  ,iz+1)
          +d2* in.atp(ix  ,iy-1,iz-1)
          +d2* in.atp(ix  ,iy-1,iz+1)
          +d2* in.atp(ix  ,iy+1,iz-1)
          +d2* in.atp(ix  ,iy+1,iz+1)
          
          +d3* in.atp(ix-1,iy-1,iz-1)
          +d3* in.atp(ix-1,iy-1,iz+1)
          +d3* in.atp(ix-1,iy+1,iz-1)
          +d3* in.atp(ix-1,iy+1,iz+1)
          +d3* in.atp(ix+1,iy-1,iz-1)
          +d3* in.atp(ix+1,iy-1,iz+1)
          +d3* in.atp(ix+1,iy+1,iz-1)
          +d3* in.atp(ix+1,iy+1,iz+1)
          ;
        }
      }
    }
  }
  // Conditions aux limites ( Dirichlet)
  _boundary(l, out);
}


/*-------------------------------------------------*/
void Operator::jacobi(int l, Vector& out, const Vector& in) const
{
  if(dim()==2)
  {
    double omega = 0.8;
    double d0inv = 3.0/8.0 * omega;
    int nx = _n(0,l), ny = _n(1,l);
    d0inv = 0.1;
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        out.atp(ix,iy) = d0inv * in.atp(ix,iy);
      }
    }
  }
  else if(dim()==3)
  {
    int nx = _n(0,l), ny = _n(1,l), nz = _n(2,l);
    double omega = 0.8;
    double d0inv = 3.0/8.0 * arma::mean(out.n()) * omega;
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<nz;iz++)
        {
          out.atp(ix,iy,iz) = d0inv * in.atp(ix,iy,iz);
        }
      }
    }
  }
  else
  {
    assert(0);
  }
  _boundary(l, out);
}

/*-------------------------------------------------*/
void Operator::gauss_seidel1(int l, Vector& out, const Vector& in) const
{
  if(dim()==2)
  {
    /*
     (ix+p)*ny + iy+q < ix*ny + iy
     p*ny +q < 0
     p=-1 q=-1,0,1
     p= 0 q=-1
     */
    int nx = _n(0,l), ny = _n(1,l);
    double omega = 0.8;
    double d0inv = 3.0/8.0;
    double d1 = -1.0/3.0;
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        out.atp(ix,iy) = d0inv * (
                                  in.atp(ix,iy)
                                  -d1* out.atp(ix-1,iy-1)
                                  -d1* out.atp(ix-1,iy  )
                                  -d1* out.atp(ix-1,iy+1)
                                  -d1* out.atp(ix  ,iy-1)
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
    int nx = _n(0,l), ny = _n(1,l), nz = _n(2,l);
    double e = 1.0/arma::mean(out.n());
    double d0 = 8.0/3.0 * e;
    double d1 = -0.0/3.0 * e;
    double d2 = -1.0/6.0 * e;
    double d3 = -1.0/12.0 * e;
    double d0inv = 1.0/d0;
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<nz;iz++)
        {
          out.atp(ix,iy,iz) = d0inv*(
                                     in.atp(ix,iy,iz)
                                     -d1* out.atp(ix-1,iy  ,iz  )
                                     -d1* out.atp(ix  ,iy-1,iz  )
                                     -d1* out.atp(ix  ,iy  ,iz-1)
                                     -d2* out.atp(ix-1,iy-1,iz  )
                                     -d2* out.atp(ix-1,iy+1,iz  )
                                     -d2* out.atp(ix-1,iy  ,iz-1)
                                     -d2* out.atp(ix-1,iy  ,iz+1)
                                     -d2* out.atp(ix  ,iy-1,iz-1)
                                     -d2* out.atp(ix  ,iy-1,iz+1)
                                     -d3* out.atp(ix-1,iy-1,iz-1)
                                     -d3* out.atp(ix-1,iy-1,iz+1)
                                     -d3* out.atp(ix-1,iy+1,iz-1)
                                     -d3* out.atp(ix-1,iy+1,iz+1)
                                     );
        }
      }
    }
  }
  _boundary(l, out);
}

/*-------------------------------------------------*/
void Operator::gauss_seidel2(int l, Vector& out, const Vector& in) const
{
  if(dim()==2)
  {
    int nx = _n(0,l), ny = _n(1,l);
    double omega = 0.8;
    double d0inv = 3.0/8.0;
    double d1 = -1.0/3.0;
    for(int ix=nx-1;ix>=0;ix--)
    {
      for(int iy=ny-1;iy>=0;iy--)
      {
        out.atp(ix,iy) = d0inv * (
                                  in.atp(ix,iy)
                                  -d1* out.atp(ix+1,iy-1)
                                  -d1* out.atp(ix+1,iy  )
                                  -d1* out.atp(ix+1,iy+1)
                                  -d1* out.atp(ix  ,iy+1)
                                  );
      }
    }
  }
  else
  {
    int nx = _n(0,l), ny = _n(1,l), nz = _n(2,l);
    double e = 1.0/arma::mean(out.n());
    double d0 = 8.0/3.0 * e;
    double d1 = -0.0/3.0 * e;
    double d2 = -1.0/6.0 * e;
    double d3 = -1.0/12.0 * e;
    double d0inv = 1.0/d0;
    for(int ix=nx-1;ix>=0;ix--)
    {
      for(int iy=ny-1;iy>=0;iy--)
      {
        for(int iz=nz-1;iz>=0;iz--)
        {
          out.atp(ix,iy,iz) = d0inv*(
                                     in.atp(ix,iy,iz)+
                                     -d1* out.atp(ix+1,iy  ,iz  )
                                     -d1* out.atp(ix  ,iy+1,iz  )
                                     -d2* out.atp(ix+1,iy  ,iz+1)
                                     -d2* out.atp(ix+1,iy-1,iz  )
                                     -d2* out.atp(ix+1,iy+1,iz  )
                                     -d2* out.atp(ix+1,iy  ,iz-1)
                                     -d2* out.atp(ix+1,iy  ,iz+1)
                                     -d2* out.atp(ix  ,iy+1,iz-1)
                                     -d2* out.atp(ix  ,iy+1,iz+1)
                                     -d3* out.atp(ix+1,iy-1,iz-1)
                                     -d3* out.atp(ix+1,iy-1,iz+1)
                                     -d3* out.atp(ix+1,iy+1,iz-1)
                                     -d3* out.atp(ix+1,iy+1,iz+1)
                                     );
        }
      }
    }
  }
  _boundary(l, out);
}
