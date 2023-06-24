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
  _seam.set_size(2*n+2);
}
/*-------------------------------------------------*/
void TransferQ13d::set_grid(const armaicvec& n, const armavec& dx)
{
  assert(n.n_elem==3);
  _nx = n[0];
  _ny = n[1];
  _nz = n[2];
  _seam.set_size(2*n+2);
}
/*-------------------------------------------------*/
void TransferQ12d::_boundary(NodeVector& v) const
{
  for(int ix=0;ix<_nx;ix++)
  {
    v.at(ix,0)    = 0;
    v.at(ix,_ny-1) = 0;
  }
  for(int iy=0;iy<_ny;iy++)
  {
    v.at(0,iy)    = 0;
    v.at(_nx-1,iy) = 0;
  }

}
/*-------------------------------------------------*/
void TransferQ13d::_boundary(NodeVector& v) const
{
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      v.at(ix,iy,0)    = 0;
      v.at(ix,iy,_nz-1) = 0;
    }
  }
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iz=0;iz<_nz;iz++)
    {
      v.at(ix,0,   iz) = 0;
      v.at(ix,_ny-1,iz) = 0;
    }
  }
  for(int iy=0;iy<_ny;iy++)
  {
    for(int iz=0;iz<_nz;iz++)
    {
      v.at(0,   iy,iz) = 0;
      v.at(_nx-1,iy,iz) = 0;
    }
  }
}

/*-------------------------------------------------*/
void TransferQ12d::restrict(NodeVector& out, const NodeVector& in) const
{
  _seam.fromvector(in);
  out.fill(0.0);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      out.at(ix,iy) = _seam.atp(2*ix  ,2*iy  )
      +  0.5 *( _seam.atp(2*ix-1,2*iy  ) + _seam.atp(2*ix+1,2*iy  )
               + _seam.atp(2*ix  ,2*iy-1) + _seam.atp(2*ix  ,2*iy+1) )
      +  0.25*( _seam.atp(2*ix-1,2*iy-1) + _seam.atp(2*ix-1,2*iy+1)
               + _seam.atp(2*ix+1,2*iy-1) + _seam.atp(2*ix+1,2*iy+1) );
    }
  }
  _boundary(out);
}
/*-------------------------------------------------*/
void TransferQ12d::prolongate(NodeVector& out, const NodeVector& in) const
{
  _seam.fill(0.0);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      _seam.atp(2*ix  ,2*iy  ) +=        in.at(ix,iy);
      _seam.atp(2*ix-1,2*iy  ) += 0.5  * in.at(ix,iy);
      _seam.atp(2*ix+1,2*iy  ) += 0.5  * in.at(ix,iy);
      _seam.atp(2*ix  ,2*iy-1) += 0.5  * in.at(ix,iy);
      _seam.atp(2*ix  ,2*iy+1) += 0.5  * in.at(ix,iy);
      _seam.atp(2*ix-1,2*iy-1) += 0.25 * in.at(ix,iy);
      _seam.atp(2*ix-1,2*iy+1) += 0.25 * in.at(ix,iy);
      _seam.atp(2*ix+1,2*iy-1) += 0.25 * in.at(ix,iy);
      _seam.atp(2*ix+1,2*iy+1) += 0.25 * in.at(ix,iy);
    }
  }
  _seam.tovector(out);
}
/*-------------------------------------------------*/
void TransferQ13d::restrict(NodeVector& out, const NodeVector& in) const
{
  _seam.fromvector(in);
  out.fill(0.0);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        out.at(ix,iy, iz) = _seam.atp(2*ix  ,2*iy,  2*iz  )
        +  0.5 *(
                   _seam.atp(2*ix-1,2*iy  ,2*iz  )
                 + _seam.atp(2*ix+1,2*iy  ,2*iz  )
                 + _seam.atp(2*ix  ,2*iy-1,2*iz  )
                 + _seam.atp(2*ix  ,2*iy+1,2*iz  )
                 + _seam.atp(2*ix  ,2*iy  ,2*iz+1)
                 + _seam.atp(2*ix  ,2*iy  ,2*iz-1)
                 )
        +  0.25*(
                   _seam.atp(2*ix-1,2*iy-1,2*iz  )
                 + _seam.atp(2*ix-1,2*iy+1,2*iz  )
                 + _seam.atp(2*ix+1,2*iy-1,2*iz  )
                 + _seam.atp(2*ix+1,2*iy+1,2*iz  )
                 + _seam.atp(2*ix-1,2*iy  ,2*iz-1)
                 + _seam.atp(2*ix-1,2*iy  ,2*iz+1)
                 + _seam.atp(2*ix+1,2*iy  ,2*iz-1)
                 + _seam.atp(2*ix+1,2*iy  ,2*iz+1)
                 + _seam.atp(2*ix  ,2*iy-1,2*iz-1)
                 + _seam.atp(2*ix  ,2*iy-1,2*iz+1)
                 + _seam.atp(2*ix  ,2*iy+1,2*iz-1)
                 + _seam.atp(2*ix  ,2*iy+1,2*iz+1)
                 )
        +  0.125*(
                    _seam.atp(2*ix-1,2*iy-1,2*iz-1)
                  + _seam.atp(2*ix-1,2*iy-1,2*iz+1)
                  + _seam.atp(2*ix-1,2*iy+1,2*iz-1)
                  + _seam.atp(2*ix-1,2*iy+1,2*iz+1)
                  + _seam.atp(2*ix+1,2*iy-1,2*iz-1)
                  + _seam.atp(2*ix+1,2*iy-1,2*iz+1)
                  + _seam.atp(2*ix+1,2*iy+1,2*iz-1)
                  + _seam.atp(2*ix+1,2*iy+1,2*iz+1)
                  );
      }
    }
  }
  _boundary(out);
}
/*-------------------------------------------------*/
void TransferQ13d::prolongate(NodeVector& out, const NodeVector& in) const
{
  _seam.fill(0.0);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        _seam.atp(2*ix  ,2*iy  ,2*iz  ) +=        in.at(ix,iy,iz);
        _seam.atp(2*ix-1,2*iy  ,2*iz  ) += 0.5  * in.at(ix,iy,iz);
        _seam.atp(2*ix+1,2*iy  ,2*iz  ) += 0.5  * in.at(ix,iy,iz);
        _seam.atp(2*ix  ,2*iy-1,2*iz  ) += 0.5  * in.at(ix,iy,iz);
        _seam.atp(2*ix  ,2*iy+1,2*iz  ) += 0.5  * in.at(ix,iy,iz);
        _seam.atp(2*ix  ,2*iy  ,2*iz-1) += 0.5  * in.at(ix,iy,iz);
        _seam.atp(2*ix  ,2*iy  ,2*iz+1) += 0.5  * in.at(ix,iy,iz);

        _seam.atp(2*ix-1,2*iy-1,2*iz  ) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix-1,2*iy+1,2*iz  ) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix+1,2*iy-1,2*iz  ) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix+1,2*iy+1,2*iz  ) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix-1,2*iy  ,2*iz-1) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix-1,2*iy  ,2*iz+1) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix+1,2*iy  ,2*iz-1) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix+1,2*iy  ,2*iz+1) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix  ,2*iy-1,2*iz-1) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix  ,2*iy-1,2*iz+1) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix  ,2*iy+1,2*iz-1) += 0.25 * in.at(ix,iy,iz);
        _seam.atp(2*ix  ,2*iy+1,2*iz+1) += 0.25 * in.at(ix,iy,iz);

        _seam.atp(2*ix-1,2*iy-1,2*iz-1) += 0.125 * in.at(ix,iy,iz);
        _seam.atp(2*ix-1,2*iy-1,2*iz+1) += 0.125 * in.at(ix,iy,iz);
        _seam.atp(2*ix-1,2*iy+1,2*iz-1) += 0.125 * in.at(ix,iy,iz);
        _seam.atp(2*ix-1,2*iy+1,2*iz+1) += 0.125 * in.at(ix,iy,iz);
        _seam.atp(2*ix+1,2*iy-1,2*iz-1) += 0.125 * in.at(ix,iy,iz);
        _seam.atp(2*ix+1,2*iy-1,2*iz+1) += 0.125 * in.at(ix,iy,iz);
        _seam.atp(2*ix+1,2*iy+1,2*iz-1) += 0.125 * in.at(ix,iy,iz);
        _seam.atp(2*ix+1,2*iy+1,2*iz+1) += 0.125 * in.at(ix,iy,iz);
      }
    }
  }
  _seam.tovector(out);
}
