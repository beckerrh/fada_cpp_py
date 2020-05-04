#include <math.h>
#include "operateur.hpp"

/*-------------------------------------------------*/
void Operateur::dot(vector& out, const vector& in, double d) const
{
  // Laplacien   elements finis q1  (9-point-stencil)
  int dim = out.dim();
  if(dim==2)
  {
    int nx = out.n(0), ny = out.n(1);
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        out(ix,iy) += d*(8. * in(ix,iy)
                         - in(ix-1,iy  ) - in(ix+1,iy)
                         - in(ix  ,iy-1) - in(ix  ,iy+1)
                         - in(ix-1,iy-1) - in(ix-1,iy+1)
                         - in(ix+1,iy-1) - in(ix+1,iy+1));
      }
    }
  }
  else if(dim==3)
  {
    int nx = out.n(0), ny = out.n(1), nz = out.n(2);
    double e = d/arma::mean(out.n());
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        for(int iz=1;iz<nz-1;iz++)
        {
        out(ix,iy,iz) += e*(
                      26. * in(ix,iy,iz)
                         - in(ix-1,iy  ,iz  ) - in(ix+1,iy  ,iz  )
                         - in(ix  ,iy-1,iz  ) - in(ix  ,iy+1,iz  )
                         - in(ix  ,iy  ,iz-1) - in(ix+1,iy  ,iz+1)
                         - in(ix-1,iy-1,iz  ) - in(ix-1,iy+1,iz  )
                         - in(ix+1,iy-1,iz  ) - in(ix+1,iy+1,iz  )
                         - in(ix-1,iy  ,iz-1) - in(ix-1,iy  ,iz+1)
                         - in(ix+1,iy  ,iz-1) - in(ix+1,iy  ,iz+1)
                         - in(ix  ,iy-1,iz-1) - in(ix  ,iy-1,iz+1)
                         - in(ix  ,iy+1,iz-1) - in(ix  ,iy+1,iz+1)
                         - in(ix-1,iy-1,iz-1) - in(ix-1,iy-1,iz+1)
                         - in(ix-1,iy+1,iz-1) - in(ix-1,iy+1,iz+1)
                         - in(ix+1,iy-1,iz-1) - in(ix+1,iy-1,iz+1)
                         - in(ix+1,iy+1,iz-1) - in(ix+1,iy+1,iz+1)
                            );
        }
      }
    }
  }
  // Conditions aux limites ( Dirichlet)
  
  //  out.boundary(in, d);
}


/*-------------------------------------------------*/
void Operateur::jacobi(vector& out, const vector& in) const
{
  int dim = out.dim();
  if(dim==2)
  {
    double omega = 0.1;
    int nx = out.n(0), ny = out.n(1);
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        out(ix,iy) = omega * in(ix,iy);
      }
    }
  }
  else if(dim==3)
  {
    int nx = out.n(0), ny = out.n(1), nz = out.n(2);
    double e = 1.0/arma::mean(out.n());
    double omega = 1./30/e;
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        for(int iz=1;iz<nz-1;iz++)
        {
          out(ix,iy,iz) = omega * in(ix,iy,iz);
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
void Operateur::gauss_seidel1(vector& out, const vector& in) const
{
  int dim = out.dim();
  if(dim==2)
  {
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
        out(ix,iy) = 0.125 * ( in(ix,iy)+ out(ix-1,iy-1) + out(ix-1,iy  )+ out(ix-1,iy+1) + out(ix  ,iy-1) );
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
    double omega = 1.0/26.0/e;
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        for(int iz=1;iz<nz-1;iz++)
        {
          out(ix,iy,iz) = omega*(
                             in(ix,iy,iz)
                           + out(ix-1,iy  ,iz  )
//                           + out(ix+1,iy  ,iz  )
                           + out(ix  ,iy-1,iz  )
//                           + out(ix  ,iy+1,iz  )
                           + out(ix  ,iy  ,iz-1)
//                           + out(ix+1,iy  ,iz+1)
                           + out(ix-1,iy-1,iz  )
                           + out(ix-1,iy+1,iz  )
//                           + out(ix+1,iy-1,iz  )
//                           + out(ix+1,iy+1,iz  )
                           + out(ix-1,iy  ,iz-1)
                           + out(ix-1,iy  ,iz+1)
//                           + out(ix+1,iy  ,iz-1)
//                           + out(ix+1,iy  ,iz+1)
                           + out(ix  ,iy-1,iz-1)
                           + out(ix  ,iy-1,iz+1)
//                           + out(ix  ,iy+1,iz-1)
//                           + out(ix  ,iy+1,iz+1)
                           + out(ix-1,iy-1,iz-1)
                           + out(ix-1,iy-1,iz+1)
                           + out(ix-1,iy+1,iz-1)
                           + out(ix-1,iy+1,iz+1)
//                           + out(ix+1,iy-1,iz-1)
//                           + out(ix+1,iy-1,iz+1)
//                           + out(ix+1,iy+1,iz-1)
//                           + out(ix+1,iy+1,iz+1)
                              );
        }
      }
    }
  }
}

/*-------------------------------------------------*/
void Operateur::gauss_seidel2(vector& out, const vector& in) const
{
  int dim = out.dim();
  if(dim==2)
  {
    int nx = out.n(0), ny = out.n(1);
    double d = 1.0/8.0;
    for(int ix=nx-2;ix>=1;ix--)
    {
      for(int iy=ny-2;iy>=1;iy--)
      {
        out(ix,iy) = d * ( in(ix,iy)
                          + out(ix+1,iy-1) + out(ix+1,iy  )
                          + out(ix+1,iy+1) + out(ix  ,iy+1) );
      }
    }
  }
  else
  {
    int nx = out.n(0), ny = out.n(1), nz = out.n(2);
    double e = 1.0/arma::mean(out.n());
    double d = 1.0/26.0/e;
    double omega = 0.99;
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        for(int iz=1;iz<nz-1;iz++)
        {
          out(ix,iy,iz) = (1-omega)*out(ix,iy,iz)
                        + d*omega*(
                             in(ix,iy,iz)+
                                   e*(
//                           + out(ix-1,iy  ,iz  )
                           + out(ix+1,iy  ,iz  )
//                           + out(ix  ,iy-1,iz  )
                           + out(ix  ,iy+1,iz  )
//                           + out(ix  ,iy  ,iz-1)
                           + out(ix+1,iy  ,iz+1)
//                           + out(ix-1,iy-1,iz  )
//                           + out(ix-1,iy+1,iz  )
                           + out(ix+1,iy-1,iz  )
                           + out(ix+1,iy+1,iz  )
//                           + out(ix-1,iy  ,iz-1)
//                           + out(ix-1,iy  ,iz+1)
                           + out(ix+1,iy  ,iz-1)
                           + out(ix+1,iy  ,iz+1)
//                           + out(ix  ,iy-1,iz-1)
//                           + out(ix  ,iy-1,iz+1)
                           + out(ix  ,iy+1,iz-1)
                           + out(ix  ,iy+1,iz+1)
//                           + out(ix-1,iy-1,iz-1)
//                           + out(ix-1,iy-1,iz+1)
//                           + out(ix-1,iy+1,iz-1)
//                           + out(ix-1,iy+1,iz+1)
                           + out(ix+1,iy-1,iz-1)
                           + out(ix+1,iy-1,iz+1)
                           + out(ix+1,iy+1,iz-1)
                           + out(ix+1,iy+1,iz+1)
                              ));
        }
      }
    }
  }
}
