//
//  transferq1.cpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "transferq1.hpp"

TransferQ12d::~TransferQ12d() {}
TransferQ13d::~TransferQ13d() {}

/*-------------------------------------------------*/
void TransferQ12d::set_grid(const armaicvec& n, const armavec& dx)
{
//  std::shared_ptr<UniformGrid> ug = std::dynamic_pointer_cast<UniformGrid>(grid);
//  assert(ug);
  assert(n.n_elem==2);
  _nx = n[0];
  _ny = n[1];
}
/*-------------------------------------------------*/
void TransferQ13d::set_grid(const armaicvec& n, const armavec& dx)
{
  assert(n.n_elem==3);
  _nx = n[0];
  _ny = n[1];
  _nz = n[2];
}
/*-------------------------------------------------*/
void TransferQ12d::_boundary(NodeVector& v) const
{
  for(int ix=0;ix<_nx;ix++)
  {
    v.atp(ix,0)    = 0;
    v.atp(ix,_ny-1) = 0;
  }
  for(int iy=0;iy<_ny;iy++)
  {
    v.atp(0,iy)    = 0;
    v.atp(_nx-1,iy) = 0;
  }

}
/*-------------------------------------------------*/
void TransferQ13d::_boundary(NodeVector& v) const
{
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      v.atp(ix,iy,0)    = 0;
      v.atp(ix,iy,_nz-1) = 0;
    }
  }
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iz=0;iz<_nz;iz++)
    {
      v.atp(ix,0,   iz) = 0;
      v.atp(ix,_ny-1,iz) = 0;
    }
  }
  for(int iy=0;iy<_ny;iy++)
  {
    for(int iz=0;iz<_nz;iz++)
    {
      v.atp(0,   iy,iz) = 0;
      v.atp(_nx-1,iy,iz) = 0;
    }
  }
}

/*-------------------------------------------------*/
void TransferQ12d::restrict(NodeVector& out, const NodeVector& in) const
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
  _boundary(out);
}
/*-------------------------------------------------*/
void TransferQ12d::prolongate(NodeVector& out, const NodeVector& in) const
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
void TransferQ13d::restrict(NodeVector& out, const NodeVector& in) const
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
  _boundary(out);
}
/*-------------------------------------------------*/
void TransferQ13d::prolongate(NodeVector& out, const NodeVector& in) const
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
