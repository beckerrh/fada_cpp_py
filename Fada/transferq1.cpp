//
//  transferq1.cpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "transferq1.hpp"
#include  "uniformgrid.hpp"
#include  "vector.hpp"

/*-------------------------------------------------*/
void TransferQ12d::set_grid(const GridInterface& grid)
{
  const UniformGrid* ug = dynamic_cast<const UniformGrid*>(&grid);
  assert(ug);
  assert(ug->dim()==2);
  _nx = ug->nx();
  _ny = ug->ny();
}
/*-------------------------------------------------*/
void TransferQ13d::set_grid(const GridInterface& grid)
{
  const UniformGrid* ug = dynamic_cast<const UniformGrid*>(&grid);
  assert(ug);
  assert(ug->dim()==3);
  _nx = ug->nx();
  _ny = ug->ny();
  _nz = ug->nz();
}
/*-------------------------------------------------*/
void TransferQ12d::restrict(Vector& out, const Vector& in) const
{
  out.fill(0.0);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      out.atp(ix,iy) = in.atp(2*ix  ,2*iy  )
      +  0.5 *( in.atp(2*ix-1,2*iy  ) + in.atp(2*ix+1,2*iy  )
               + in.atp(2*ix  ,2*iy-1) + in.atp(2*ix  ,2*iy+1) )
      +  0.25*( in.atp(2*ix-1,2*iy-1) + in.atp(2*ix-1,2*iy+1)
               + in.atp(2*ix+1,2*iy-1) + in.atp(2*ix+1,2*iy+1) );
    }
  }
}
/*-------------------------------------------------*/
void TransferQ12d::prolongate(Vector& out, const Vector& in) const
{
  out.fill(0.0);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      out.atp(2*ix  ,2*iy  ) +=        in.atp(ix,iy);
      out.atp(2*ix-1,2*iy  ) += 0.5  * in.atp(ix,iy);
      out.atp(2*ix+1,2*iy  ) += 0.5  * in.atp(ix,iy);
      out.atp(2*ix  ,2*iy-1) += 0.5  * in.atp(ix,iy);
      out.atp(2*ix  ,2*iy+1) += 0.5  * in.atp(ix,iy);
      out.atp(2*ix-1,2*iy-1) += 0.25 * in.atp(ix,iy);
      out.atp(2*ix-1,2*iy+1) += 0.25 * in.atp(ix,iy);
      out.atp(2*ix+1,2*iy-1) += 0.25 * in.atp(ix,iy);
      out.atp(2*ix+1,2*iy+1) += 0.25 * in.atp(ix,iy);
    }
  }
}
/*-------------------------------------------------*/
void TransferQ13d::restrict(Vector& out, const Vector& in) const
{
  out.fill(0.0);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        out.atp(ix,iy, iz) = in.atp(2*ix  ,2*iy,  2*iz  )
        +  0.5 *(
                 in.atp(2*ix-1,2*iy  ,2*iz  )
                 + in.atp(2*ix+1,2*iy  ,2*iz  )
                 + in.atp(2*ix  ,2*iy-1,2*iz  )
                 + in.atp(2*ix  ,2*iy+1,2*iz  )
                 + in.atp(2*ix  ,2*iy  ,2*iz+1)
                 + in.atp(2*ix  ,2*iy  ,2*iz-1)
                 )
        +  0.25*(
                 in.atp(2*ix-1,2*iy-1,2*iz  )
                 + in.atp(2*ix-1,2*iy+1,2*iz  )
                 + in.atp(2*ix+1,2*iy-1,2*iz  )
                 + in.atp(2*ix+1,2*iy+1,2*iz  )
                 + in.atp(2*ix-1,2*iy  ,2*iz-1)
                 + in.atp(2*ix-1,2*iy  ,2*iz+1)
                 + in.atp(2*ix+1,2*iy  ,2*iz-1)
                 + in.atp(2*ix+1,2*iy  ,2*iz+1)
                 + in.atp(2*ix  ,2*iy-1,2*iz-1)
                 + in.atp(2*ix  ,2*iy-1,2*iz+1)
                 + in.atp(2*ix  ,2*iy+1,2*iz-1)
                 + in.atp(2*ix  ,2*iy+1,2*iz+1)
                 )
        +  0.125*(
                  in.atp(2*ix-1,2*iy-1,2*iz-1)
                  + in.atp(2*ix-1,2*iy-1,2*iz+1)
                  + in.atp(2*ix-1,2*iy+1,2*iz-1)
                  + in.atp(2*ix-1,2*iy+1,2*iz+1)
                  + in.atp(2*ix+1,2*iy-1,2*iz-1)
                  + in.atp(2*ix+1,2*iy-1,2*iz+1)
                  + in.atp(2*ix+1,2*iy+1,2*iz-1)
                  + in.atp(2*ix+1,2*iy+1,2*iz+1)
                  );
      }
    }
  }
}
/*-------------------------------------------------*/
void TransferQ13d::prolongate(Vector& out, const Vector& in) const
{
  out.fill(0.0);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        out.atp(2*ix  ,2*iy  ,2*iz  ) +=        in.atp(ix,iy,iz);
        out.atp(2*ix-1,2*iy  ,2*iz  ) += 0.5  * in.atp(ix,iy,iz);
        out.atp(2*ix+1,2*iy  ,2*iz  ) += 0.5  * in.atp(ix,iy,iz);
        out.atp(2*ix  ,2*iy-1,2*iz  ) += 0.5  * in.atp(ix,iy,iz);
        out.atp(2*ix  ,2*iy+1,2*iz  ) += 0.5  * in.atp(ix,iy,iz);
        out.atp(2*ix  ,2*iy  ,2*iz-1) += 0.5  * in.atp(ix,iy,iz);
        out.atp(2*ix  ,2*iy  ,2*iz+1) += 0.5  * in.atp(ix,iy,iz);
        
        out.atp(2*ix-1,2*iy-1,2*iz  ) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix-1,2*iy+1,2*iz  ) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix+1,2*iy-1,2*iz  ) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix+1,2*iy+1,2*iz  ) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix-1,2*iy  ,2*iz-1) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix-1,2*iy  ,2*iz+1) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix+1,2*iy  ,2*iz-1) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix+1,2*iy  ,2*iz+1) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix  ,2*iy-1,2*iz-1) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix  ,2*iy-1,2*iz+1) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix  ,2*iy+1,2*iz-1) += 0.25 * in.atp(ix,iy,iz);
        out.atp(2*ix  ,2*iy+1,2*iz+1) += 0.25 * in.atp(ix,iy,iz);
        
        out.atp(2*ix-1,2*iy-1,2*iz-1) += 0.125 * in.atp(ix,iy,iz);
        out.atp(2*ix-1,2*iy-1,2*iz+1) += 0.125 * in.atp(ix,iy,iz);
        out.atp(2*ix-1,2*iy+1,2*iz-1) += 0.125 * in.atp(ix,iy,iz);
        out.atp(2*ix-1,2*iy+1,2*iz+1) += 0.125 * in.atp(ix,iy,iz);
        out.atp(2*ix+1,2*iy-1,2*iz-1) += 0.125 * in.atp(ix,iy,iz);
        out.atp(2*ix+1,2*iy-1,2*iz+1) += 0.125 * in.atp(ix,iy,iz);
        out.atp(2*ix+1,2*iy+1,2*iz-1) += 0.125 * in.atp(ix,iy,iz);
        out.atp(2*ix+1,2*iy+1,2*iz+1) += 0.125 * in.atp(ix,iy,iz);
      }
    }
  }
}
