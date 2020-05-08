#include <math.h>
#include "operator.hpp"

/*-------------------------------------------------*/
void Operator::restrict(int l, VectorMG& out, const VectorMG& in) const
//   in_{l+1} --> out_l
{
  out(l).fill(0.0);
  if(dim()==2)
  {
    int nx = _n(0,l), ny = _n(1,l);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        out(l).atp(ix,iy) = in(l+1).atp(2*ix  ,2*iy  )
        +  0.5 *( in(l+1).atp(2*ix-1,2*iy  ) + in(l+1).atp(2*ix+1,2*iy  )
                 + in(l+1).atp(2*ix  ,2*iy-1) + in(l+1).atp(2*ix  ,2*iy+1) )
        +  0.25*( in(l+1).atp(2*ix-1,2*iy-1) + in(l+1).atp(2*ix-1,2*iy+1)
                 + in(l+1).atp(2*ix+1,2*iy-1) + in(l+1).atp(2*ix+1,2*iy+1) );
      }
    }
  }
  else if(dim()==3)
  {
    int nx = _n(0,l), ny = _n(1,l), nz = _n(2,l);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<nz;iz++)
        {
          out(l).atp(ix,iy, iz) = in(l+1).atp(2*ix  ,2*iy,  2*iz  )
          +  0.5 *(
                   in(l+1).atp(2*ix-1,2*iy  ,2*iz  )
                   + in(l+1).atp(2*ix+1,2*iy  ,2*iz  )
                   + in(l+1).atp(2*ix  ,2*iy-1,2*iz  )
                   + in(l+1).atp(2*ix  ,2*iy+1,2*iz  )
                   + in(l+1).atp(2*ix  ,2*iy  ,2*iz+1)
                   + in(l+1).atp(2*ix  ,2*iy  ,2*iz-1)
                   )
          +  0.25*(
                   in(l+1).atp(2*ix-1,2*iy-1,2*iz  )
                   + in(l+1).atp(2*ix-1,2*iy+1,2*iz  )
                   + in(l+1).atp(2*ix+1,2*iy-1,2*iz  )
                   + in(l+1).atp(2*ix+1,2*iy+1,2*iz  )
                   + in(l+1).atp(2*ix-1,2*iy  ,2*iz-1)
                   + in(l+1).atp(2*ix-1,2*iy  ,2*iz+1)
                   + in(l+1).atp(2*ix+1,2*iy  ,2*iz-1)
                   + in(l+1).atp(2*ix+1,2*iy  ,2*iz+1)
                   + in(l+1).atp(2*ix  ,2*iy-1,2*iz-1)
                   + in(l+1).atp(2*ix  ,2*iy-1,2*iz+1)
                   + in(l+1).atp(2*ix  ,2*iy+1,2*iz-1)
                   + in(l+1).atp(2*ix  ,2*iy+1,2*iz+1)
                   )
          +  0.125*(
                    in(l+1).atp(2*ix-1,2*iy-1,2*iz-1)
                    + in(l+1).atp(2*ix-1,2*iy-1,2*iz+1)
                    + in(l+1).atp(2*ix-1,2*iy+1,2*iz-1)
                    + in(l+1).atp(2*ix-1,2*iy+1,2*iz+1)
                    + in(l+1).atp(2*ix+1,2*iy-1,2*iz-1)
                    + in(l+1).atp(2*ix+1,2*iy-1,2*iz+1)
                    + in(l+1).atp(2*ix+1,2*iy+1,2*iz-1)
                    + in(l+1).atp(2*ix+1,2*iy+1,2*iz+1)
                    );
        }
      }
    }
  }
  _boundary(l, out(l));
  //  boundary(out(l));
}

/*-------------------------------------------------*/
void Operator::prolongate(int l, VectorMG& out, const VectorMG& in) const
//   in_{l-1} --> out_l
{
  out(l).fill(0.0);
  if(dim()==2)
  {
    int nx = _n(0,l-1), ny = _n(1,l-1);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        out(l).atp(2*ix  ,2*iy  ) +=        in(l-1).atp(ix,iy);
        out(l).atp(2*ix-1,2*iy  ) += 0.5  * in(l-1).atp(ix,iy);
        out(l).atp(2*ix+1,2*iy  ) += 0.5  * in(l-1).atp(ix,iy);
        out(l).atp(2*ix  ,2*iy-1) += 0.5  * in(l-1).atp(ix,iy);
        out(l).atp(2*ix  ,2*iy+1) += 0.5  * in(l-1).atp(ix,iy);
        out(l).atp(2*ix-1,2*iy-1) += 0.25 * in(l-1).atp(ix,iy);
        out(l).atp(2*ix-1,2*iy+1) += 0.25 * in(l-1).atp(ix,iy);
        out(l).atp(2*ix+1,2*iy-1) += 0.25 * in(l-1).atp(ix,iy);
        out(l).atp(2*ix+1,2*iy+1) += 0.25 * in(l-1).atp(ix,iy);
      }
    }
  }
  else if(dim()==3)
  {
    int nx = _n(0,l-1), ny = _n(1,l-1), nz = _n(2,l-1);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<nz;iz++)
        {
          out(l).atp(2*ix  ,2*iy  ,2*iz  ) +=        in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix-1,2*iy  ,2*iz  ) += 0.5  * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix+1,2*iy  ,2*iz  ) += 0.5  * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix  ,2*iy-1,2*iz  ) += 0.5  * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix  ,2*iy+1,2*iz  ) += 0.5  * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix  ,2*iy  ,2*iz-1) += 0.5  * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix  ,2*iy  ,2*iz+1) += 0.5  * in(l-1).atp(ix,iy,iz);

          out(l).atp(2*ix-1,2*iy-1,2*iz  ) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix-1,2*iy+1,2*iz  ) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix+1,2*iy-1,2*iz  ) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix+1,2*iy+1,2*iz  ) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix-1,2*iy  ,2*iz-1) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix-1,2*iy  ,2*iz+1) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix+1,2*iy  ,2*iz-1) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix+1,2*iy  ,2*iz+1) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix  ,2*iy-1,2*iz-1) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix  ,2*iy-1,2*iz+1) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix  ,2*iy+1,2*iz-1) += 0.25 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix  ,2*iy+1,2*iz+1) += 0.25 * in(l-1).atp(ix,iy,iz);
          
          out(l).atp(2*ix-1,2*iy-1,2*iz-1) += 0.125 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix-1,2*iy-1,2*iz+1) += 0.125 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix-1,2*iy+1,2*iz-1) += 0.125 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix-1,2*iy+1,2*iz+1) += 0.125 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix+1,2*iy-1,2*iz-1) += 0.125 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix+1,2*iy-1,2*iz+1) += 0.125 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix+1,2*iy+1,2*iz-1) += 0.125 * in(l-1).atp(ix,iy,iz);
          out(l).atp(2*ix+1,2*iy+1,2*iz+1) += 0.125 * in(l-1).atp(ix,iy,iz);
        }
      }
    }
  }
}
