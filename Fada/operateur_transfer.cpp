#include <math.h>
#include "operateur.hpp"

/*-------------------------------------------------*/
void Operateur::restrict(int l, VecteurMG& out, const VecteurMG& in) const
//   in_{l+1} --> out_l
{
  int dim = out(l).dim();
  if(dim==2)
  {
    int nx = out(l).n(0), ny = out(l).n(1);
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        out(l)(ix,iy) = in(l+1)(2*ix  ,2*iy  )
        +  0.5 *( in(l+1)(2*ix-1,2*iy  ) + in(l+1)(2*ix+1,2*iy  )
                 + in(l+1)(2*ix  ,2*iy-1) + in(l+1)(2*ix  ,2*iy+1) )
        +  0.25*( in(l+1)(2*ix-1,2*iy-1) + in(l+1)(2*ix-1,2*iy+1)
                 + in(l+1)(2*ix+1,2*iy-1) + in(l+1)(2*ix+1,2*iy+1) );
      }
    }
  }
  else if(dim==3)
  {
    int nx = out(l).n(0), ny = out(l).n(1), nz = out(l).n(2);
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        for(int iz=1;iz<nz-1;iz++)
        {
          out(l)(ix,iy, iz) = in(l+1)(2*ix  ,2*iy,  2*iz  )
          +  0.5 *(
                   in(l+1)(2*ix-1,2*iy  ,2*iz  )
                   + in(l+1)(2*ix+1,2*iy  ,2*iz  )
                   + in(l+1)(2*ix  ,2*iy-1,2*iz  )
                   + in(l+1)(2*ix  ,2*iy+1,2*iz  )
                   + in(l+1)(2*ix  ,2*iy  ,2*iz+1)
                   + in(l+1)(2*ix  ,2*iy  ,2*iz-1)
                   )
          +  0.25*(
                   in(l+1)(2*ix-1,2*iy-1,2*iz  )
                   + in(l+1)(2*ix-1,2*iy+1,2*iz  )
                   + in(l+1)(2*ix+1,2*iy-1,2*iz  )
                   + in(l+1)(2*ix+1,2*iy+1,2*iz  )
                   + in(l+1)(2*ix-1,2*iy  ,2*iz-1)
                   + in(l+1)(2*ix-1,2*iy  ,2*iz+1)
                   + in(l+1)(2*ix+1,2*iy  ,2*iz-1)
                   + in(l+1)(2*ix+1,2*iy  ,2*iz+1)
                   + in(l+1)(2*ix  ,2*iy-1,2*iz-1)
                   + in(l+1)(2*ix  ,2*iy-1,2*iz+1)
                   + in(l+1)(2*ix  ,2*iy+1,2*iz-1)
                   + in(l+1)(2*ix  ,2*iy+1,2*iz+1)
                   )
          +  0.125*(
                    in(l+1)(2*ix-1,2*iy-1,2*iz-1)
                    + in(l+1)(2*ix-1,2*iy-1,2*iz+1)
                    + in(l+1)(2*ix-1,2*iy+1,2*iz-1)
                    + in(l+1)(2*ix-1,2*iy+1,2*iz+1)
                    + in(l+1)(2*ix+1,2*iy-1,2*iz-1)
                    + in(l+1)(2*ix+1,2*iy-1,2*iz+1)
                    + in(l+1)(2*ix+1,2*iy+1,2*iz-1)
                    + in(l+1)(2*ix+1,2*iy+1,2*iz+1)
                    );
        }
      }
    }
  }
  //  boundary(out(l));
}

/*-------------------------------------------------*/
void Operateur::prolongate(int l, VecteurMG& out, const VecteurMG& in) const
//   in_{l-1} --> out_l
{
  int dim = out(l).dim();
  //  out(l) = 0.0;
  out(l).fill(0.0);
  if(dim==2)
  {
    int nx = in(l-1).n(0), ny = in(l-1).n(1);
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        out(l)(2*ix  ,2*iy  ) +=        in(l-1)(ix,iy);
        out(l)(2*ix-1,2*iy  ) += 0.5  * in(l-1)(ix,iy);
        out(l)(2*ix+1,2*iy  ) += 0.5  * in(l-1)(ix,iy);
        out(l)(2*ix  ,2*iy-1) += 0.5  * in(l-1)(ix,iy);
        out(l)(2*ix  ,2*iy+1) += 0.5  * in(l-1)(ix,iy);
        out(l)(2*ix-1,2*iy-1) += 0.25 * in(l-1)(ix,iy);
        out(l)(2*ix-1,2*iy+1) += 0.25 * in(l-1)(ix,iy);
        out(l)(2*ix+1,2*iy-1) += 0.25 * in(l-1)(ix,iy);
        out(l)(2*ix+1,2*iy+1) += 0.25 * in(l-1)(ix,iy);
      }
    }
  }
  else if(dim==3)
  {
    int nx = in(l-1).n(0), ny = in(l-1).n(1), nz = in(l-1).n(2);
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        for(int iz=1;iz<nz-1;iz++)
        {
          out(l)(2*ix  ,2*iy  ,2*iz  ) +=        in(l-1)(ix,iy,iz);
          out(l)(2*ix-1,2*iy  ,2*iz  ) += 0.5  * in(l-1)(ix,iy,iz);
          out(l)(2*ix+1,2*iy  ,2*iz  ) += 0.5  * in(l-1)(ix,iy,iz);
          out(l)(2*ix  ,2*iy-1,2*iz  ) += 0.5  * in(l-1)(ix,iy,iz);
          out(l)(2*ix  ,2*iy+1,2*iz  ) += 0.5  * in(l-1)(ix,iy,iz);
          out(l)(2*ix  ,2*iy  ,2*iz-1) += 0.5  * in(l-1)(ix,iy,iz);
          out(l)(2*ix  ,2*iy  ,2*iz+1) += 0.5  * in(l-1)(ix,iy,iz);

          out(l)(2*ix-1,2*iy-1,2*iz  ) += 0.25 * in(l-1)(ix,iy,iz);
          out(l)(2*ix-1,2*iy+1,2*iz  ) += 0.25 * in(l-1)(ix,iy,iz);
          out(l)(2*ix+1,2*iy-1,2*iz  ) += 0.25 * in(l-1)(ix,iy,iz);
          out(l)(2*ix+1,2*iy+1,2*iz  ) += 0.25 * in(l-1)(ix,iy,iz);
          
          out(l)(2*ix-1,2*iy  ,2*iz-1) += 0.25 * in(l-1)(ix,iy,iz);
          out(l)(2*ix-1,2*iy  ,2*iz+1) += 0.25 * in(l-1)(ix,iy,iz);
          out(l)(2*ix+1,2*iy  ,2*iz-1) += 0.25 * in(l-1)(ix,iy,iz);
          out(l)(2*ix+1,2*iy  ,2*iz+1) += 0.25 * in(l-1)(ix,iy,iz);
          
          out(l)(2*ix  ,2*iy-1,2*iz-1) += 0.25 * in(l-1)(ix,iy,iz);
          out(l)(2*ix  ,2*iy-1,2*iz+1) += 0.25 * in(l-1)(ix,iy,iz);
          out(l)(2*ix  ,2*iy+1,2*iz-1) += 0.25 * in(l-1)(ix,iy,iz);
          out(l)(2*ix  ,2*iy+1,2*iz+1) += 0.25 * in(l-1)(ix,iy,iz);

          
          out(l)(2*ix-1,2*iy-1,2*iz-1) += 0.125 * in(l-1)(ix,iy,iz);
          out(l)(2*ix-1,2*iy-1,2*iz+1) += 0.125 * in(l-1)(ix,iy,iz);
          out(l)(2*ix-1,2*iy+1,2*iz-1) += 0.125 * in(l-1)(ix,iy,iz);
          out(l)(2*ix-1,2*iy+1,2*iz+1) += 0.125 * in(l-1)(ix,iy,iz);
          
          out(l)(2*ix+1,2*iy-1,2*iz-1) += 0.125 * in(l-1)(ix,iy,iz);
          out(l)(2*ix+1,2*iy-1,2*iz+1) += 0.125 * in(l-1)(ix,iy,iz);
          out(l)(2*ix+1,2*iy+1,2*iz-1) += 0.125 * in(l-1)(ix,iy,iz);
          out(l)(2*ix+1,2*iy+1,2*iz+1) += 0.125 * in(l-1)(ix,iy,iz);
        }
      }
    }
  }
}
