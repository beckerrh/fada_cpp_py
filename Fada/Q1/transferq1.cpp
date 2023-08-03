//
//  transferq1.cpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "transferq1.hpp"

// /*-------------------------------------------------*/
// void TransferBase::_boundary(GridVector& v) const
// {
//   // return;
//   // for(const auto p: *_boundaryconditions) for( auto q: p) std::cerr << q ;
//   if(_dim==2)
//   {
//       if ((*_boundaryconditions)[0][0] == "dir")
//       {
//         for (int iy = 0; iy < _ny; iy++)
//         {
//           v.at(0, iy) = 0;
//         }
//       }
//       if ((*_boundaryconditions)[0][1] == "dir")
//       {
//         for (int iy = 0; iy < _ny; iy++)
//         {
//           v.at(_nx - 1, iy) = 0;
//         }
//       }
//       if ((*_boundaryconditions)[1][0] == "dir")
//       {
//         for (int ix = 0; ix < _nx; ix++)
//         {
//           v.at(ix, 0) = 0;
//         }
//       }
//       if ((*_boundaryconditions)[1][1] == "dir")
//       {
//         for (int ix = 0; ix < _nx; ix++)
//         {
//           v.at(ix, _ny - 1) = 0;
//         }
//       }
//     }
//     else
//     {
//   if ((*_boundaryconditions)[0][0] == "dir")
//   {
//     for (int iy = 0; iy < _ny; iy++)
//     {
//       for (int iz = 0; iz < _nz; iz++)
//       {
//         v.at(0, iy, iz)       = 0;
//       }
//     }
//   }
//   if ((*_boundaryconditions)[0][1] == "dir")
//   {
//     for (int iy = 0; iy < _ny; iy++)
//     {
//       for (int iz = 0; iz < _nz; iz++)
//       {
//         v.at(_nx - 1, iy, iz) = 0;
//       }
//     }
//   }
//   if ((*_boundaryconditions)[1][0] == "dir")
//   {
//     for (int ix = 0; ix < _nx; ix++)
//     {
//       for (int iz = 0; iz < _nz; iz++)
//       {
//         v.at(ix, 0, iz)       = 0;
//       }
//     }
//   }
//   if ((*_boundaryconditions)[1][1] == "dir")
//   {
//     for (int ix = 0; ix < _nx; ix++)
//     {
//       for (int iz = 0; iz < _nz; iz++)
//       {
//         v.at(ix, _ny - 1, iz) = 0;
//       }
//     }
//   }
//   if ((*_boundaryconditions)[2][0] == "dir")
//   {
//     for (int ix = 0; ix < _nx; ix++)
//     {
//       for (int iy = 0; iy < _ny; iy++)
//       {
//         v.at(ix, iy, 0)       = 0;
//       }
//     }
//   }
//   if ((*_boundaryconditions)[2][1] == "dir")
//   {
//     for (int ix = 0; ix < _nx; ix++)
//     {
//       for (int iy = 0; iy < _ny; iy++)
//       {
//         v.at(ix, iy, _nz - 1) = 0;
//       }
//     }
//   }
// }
// }

/*-------------------------------------------------*/
void TransferBase::set_grid(const armaicvec& n, const armavec& dx)
{
  _dim = n.n_elem;
  _nx = n[0];
  _ny = n[1];
  if(_dim==3)
  {
    _nz = n[2];
  }
  _seam.set_size(2 * n + 2);
}
/*-------------------------------------------------*/
void TransferQ12d::restrict (GridVector & out, const GridVector& in) const
{
  // _seam.fromvector(in);
  // out.fill(0.0);
  // for(int ix=0;ix<_nx;ix++)
  // {
  //   for(int iy=0;iy<_ny;iy++)
  //   {
  //     out.at(ix,iy) = _seam.atp(2*ix  ,2*iy  )
  //     +  0.5 * _seam.atp(2*ix-1,2*iy  )
  //     +  0.5 * _seam.atp(2*ix+1,2*iy  )
  //     +  0.5 * _seam.atp(2*ix  ,2*iy-1)
  //     +  0.5 * _seam.atp(2*ix  ,2*iy+1)
  //     +  0.25* _seam.atp(2*ix-1,2*iy-1)
  //     +  0.25* _seam.atp(2*ix-1,2*iy+1)
  //     +  0.25* _seam.atp(2*ix+1,2*iy-1)
  //     +  0.25* _seam.atp(2*ix+1,2*iy+1)
  //   );
  //   }
  // }
  out.fill(0.0);
  out.at(0, 0) = in.at(0, 0)
                 + 0.5 * in.at(1, 0)
                 + 0.5 * in.at(0, 1)
                 + 0.25 * in.at(1, 1)
  ;
  for (int iy = 1; iy < _ny - 1; iy++)
  {
    out.at(0, iy) = in.at(0, 2 * iy)
                    + 0.5 * in.at(1, 2 * iy)
                    + 0.5 * in.at(0, 2 * iy - 1)
                    + 0.5 * in.at(0, 2 * iy + 1)
                    + 0.25 * in.at(1, 2 * iy - 1)
                    + 0.25 * in.at(1, 2 * iy + 1)
    ;
  }
  out.at(0, _ny - 1) = in.at(0, 2 * _ny - 2)
                       + 0.5 * in.at(1, 2 * _ny - 2)
                       + 0.5 * in.at(0, 2 * _ny - 3)
                       + 0.25 * in.at(1, 2 * _ny - 3)
  ;
  for (int ix = 1; ix < _nx - 1; ix++)
  {
    out.at(ix, 0) = in.at(2 * ix, 0)
                    + 0.5 * in.at(2 * ix - 1, 0)
                    + 0.5 * in.at(2 * ix + 1, 0)
                    + 0.5 * in.at(2 * ix, 1)
                    + 0.25 * in.at(2 * ix - 1, 1)
                    + 0.25 * in.at(2 * ix + 1, 1)
    ;
    for (int iy = 1; iy < _ny - 1; iy++)
    {
      out.at(ix, iy) = in.at(2 * ix, 2 * iy)
                       + 0.5 * in.at(2 * ix - 1, 2 * iy)
                       + 0.5 * in.at(2 * ix + 1, 2 * iy)
                       + 0.5 * in.at(2 * ix, 2 * iy - 1)
                       + 0.5 * in.at(2 * ix, 2 * iy + 1)
                       + 0.25 * in.at(2 * ix - 1, 2 * iy - 1)
                       + 0.25 * in.at(2 * ix - 1, 2 * iy + 1)
                       + 0.25 * in.at(2 * ix + 1, 2 * iy - 1)
                       + 0.25 * in.at(2 * ix + 1, 2 * iy + 1)
      ;
    }
    out.at(ix, _ny - 1) = in.at(2 * ix, 2 * _ny - 2)
                          + 0.5 * in.at(2 * ix - 1, 2 * _ny - 2)
                          + 0.5 * in.at(2 * ix + 1, 2 * _ny - 2)
                          + 0.5 * in.at(2 * ix, 2 * _ny - 3)
                          + 0.25 * in.at(2 * ix - 1, 2 * _ny - 3)
                          + 0.25 * in.at(2 * ix + 1, 2 * _ny - 3)
    ;
  }
  out.at(_nx - 1, 0) = in.at(2 * _nx - 2, 0)
                       + 0.5 * in.at(2 * _nx - 3, 0)
                       + 0.5 * in.at(2 * _nx - 2, 1)
                       + 0.25 * in.at(2 * _nx - 3, 1)
  ;
  for (int iy = 1; iy < _ny - 1; iy++)
  {
    out.at(_nx - 1, iy) = in.at(2 * _nx - 2, 2 * iy)
                          + 0.5 * in.at(2 * _nx - 3, 2 * iy)
                          + 0.5 * in.at(2 * _nx - 2, 2 * iy - 1)
                          + 0.5 * in.at(2 * _nx - 2, 2 * iy + 1)
                          + 0.25 * in.at(2 * _nx - 3, 2 * iy - 1)
                          + 0.25 * in.at(2 * _nx - 3, 2 * iy + 1)
    ;
  }
  out.at(_nx - 1, _ny - 1) = in.at(2 * _nx - 2, 2 * _ny - 2)
                             + 0.5 * in.at(2 * _nx - 3, 2 * _ny - 2)
                             + 0.5 * in.at(2 * _nx - 2, 2 * _ny - 3)
                             + 0.25 * in.at(2 * _nx - 3, 2 * _ny - 3)
  ;
  // _boundary(out);
  // out.boundary_zero();
}
/*-------------------------------------------------*/
void TransferQ12d::prolongate(GridVector& out, const GridVector& in) const
{
  _seam.fill(0.0);
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      _seam.atp(2 * ix, 2 * iy)         += in.at(ix, iy);
      _seam.atp(2 * ix - 1, 2 * iy)     += 0.5 * in.at(ix, iy);
      _seam.atp(2 * ix + 1, 2 * iy)     += 0.5 * in.at(ix, iy);
      _seam.atp(2 * ix, 2 * iy - 1)     += 0.5 * in.at(ix, iy);
      _seam.atp(2 * ix, 2 * iy + 1)     += 0.5 * in.at(ix, iy);
      _seam.atp(2 * ix - 1, 2 * iy - 1) += 0.25 * in.at(ix, iy);
      _seam.atp(2 * ix - 1, 2 * iy + 1) += 0.25 * in.at(ix, iy);
      _seam.atp(2 * ix + 1, 2 * iy - 1) += 0.25 * in.at(ix, iy);
      _seam.atp(2 * ix + 1, 2 * iy + 1) += 0.25 * in.at(ix, iy);
    }
  }
  _seam.tovector(out);
}

/*-------------------------------------------------*/
void TransferQ13d::restrict (GridVector & out, const GridVector& in) const
{
  _seam.fromvector(in);
  out.fill(0.0);
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      for (int iz = 0; iz < _nz; iz++)
      {
        out.at(ix, iy, iz) = _seam.atp(2 * ix, 2 * iy, 2 * iz)
                             + 0.5 * (
          _seam.atp(2 * ix - 1, 2 * iy, 2 * iz)
          + _seam.atp(2 * ix + 1, 2 * iy, 2 * iz)
          + _seam.atp(2 * ix, 2 * iy - 1, 2 * iz)
          + _seam.atp(2 * ix, 2 * iy + 1, 2 * iz)
          + _seam.atp(2 * ix, 2 * iy, 2 * iz + 1)
          + _seam.atp(2 * ix, 2 * iy, 2 * iz - 1)
          )
                             + 0.25 * (
          _seam.atp(2 * ix - 1, 2 * iy - 1, 2 * iz)
          + _seam.atp(2 * ix - 1, 2 * iy + 1, 2 * iz)
          + _seam.atp(2 * ix + 1, 2 * iy - 1, 2 * iz)
          + _seam.atp(2 * ix + 1, 2 * iy + 1, 2 * iz)
          + _seam.atp(2 * ix - 1, 2 * iy, 2 * iz - 1)
          + _seam.atp(2 * ix - 1, 2 * iy, 2 * iz + 1)
          + _seam.atp(2 * ix + 1, 2 * iy, 2 * iz - 1)
          + _seam.atp(2 * ix + 1, 2 * iy, 2 * iz + 1)
          + _seam.atp(2 * ix, 2 * iy - 1, 2 * iz - 1)
          + _seam.atp(2 * ix, 2 * iy - 1, 2 * iz + 1)
          + _seam.atp(2 * ix, 2 * iy + 1, 2 * iz - 1)
          + _seam.atp(2 * ix, 2 * iy + 1, 2 * iz + 1)
          )
                             + 0.125 * (
          _seam.atp(2 * ix - 1, 2 * iy - 1, 2 * iz - 1)
          + _seam.atp(2 * ix - 1, 2 * iy - 1, 2 * iz + 1)
          + _seam.atp(2 * ix - 1, 2 * iy + 1, 2 * iz - 1)
          + _seam.atp(2 * ix - 1, 2 * iy + 1, 2 * iz + 1)
          + _seam.atp(2 * ix + 1, 2 * iy - 1, 2 * iz - 1)
          + _seam.atp(2 * ix + 1, 2 * iy - 1, 2 * iz + 1)
          + _seam.atp(2 * ix + 1, 2 * iy + 1, 2 * iz - 1)
          + _seam.atp(2 * ix + 1, 2 * iy + 1, 2 * iz + 1)
          );
      }
    }
  }
  // out.boundary_zero();
  // _boundary(out);
}
/*-------------------------------------------------*/
void TransferQ13d::prolongate(GridVector& out, const GridVector& in) const
{
  _seam.fill(0.0);
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      for (int iz = 0; iz < _nz; iz++)
      {
        _seam.atp(2 * ix, 2 * iy, 2 * iz)     += in.at(ix, iy, iz);
        _seam.atp(2 * ix - 1, 2 * iy, 2 * iz) += 0.5 * in.at(ix, iy, iz);
        _seam.atp(2 * ix + 1, 2 * iy, 2 * iz) += 0.5 * in.at(ix, iy, iz);
        _seam.atp(2 * ix, 2 * iy - 1, 2 * iz) += 0.5 * in.at(ix, iy, iz);
        _seam.atp(2 * ix, 2 * iy + 1, 2 * iz) += 0.5 * in.at(ix, iy, iz);
        _seam.atp(2 * ix, 2 * iy, 2 * iz - 1) += 0.5 * in.at(ix, iy, iz);
        _seam.atp(2 * ix, 2 * iy, 2 * iz + 1) += 0.5 * in.at(ix, iy, iz);

        _seam.atp(2 * ix - 1, 2 * iy - 1, 2 * iz) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix - 1, 2 * iy + 1, 2 * iz) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix + 1, 2 * iy - 1, 2 * iz) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix + 1, 2 * iy + 1, 2 * iz) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix - 1, 2 * iy, 2 * iz - 1) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix - 1, 2 * iy, 2 * iz + 1) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix + 1, 2 * iy, 2 * iz - 1) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix + 1, 2 * iy, 2 * iz + 1) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix, 2 * iy - 1, 2 * iz - 1) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix, 2 * iy - 1, 2 * iz + 1) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix, 2 * iy + 1, 2 * iz - 1) += 0.25 * in.at(ix, iy, iz);
        _seam.atp(2 * ix, 2 * iy + 1, 2 * iz + 1) += 0.25 * in.at(ix, iy, iz);

        _seam.atp(2 * ix - 1, 2 * iy - 1, 2 * iz - 1) += 0.125 * in.at(ix, iy, iz);
        _seam.atp(2 * ix - 1, 2 * iy - 1, 2 * iz + 1) += 0.125 * in.at(ix, iy, iz);
        _seam.atp(2 * ix - 1, 2 * iy + 1, 2 * iz - 1) += 0.125 * in.at(ix, iy, iz);
        _seam.atp(2 * ix - 1, 2 * iy + 1, 2 * iz + 1) += 0.125 * in.at(ix, iy, iz);
        _seam.atp(2 * ix + 1, 2 * iy - 1, 2 * iz - 1) += 0.125 * in.at(ix, iy, iz);
        _seam.atp(2 * ix + 1, 2 * iy - 1, 2 * iz + 1) += 0.125 * in.at(ix, iy, iz);
        _seam.atp(2 * ix + 1, 2 * iy + 1, 2 * iz - 1) += 0.125 * in.at(ix, iy, iz);
        _seam.atp(2 * ix + 1, 2 * iy + 1, 2 * iz + 1) += 0.125 * in.at(ix, iy, iz);
      }
    }
  }
  _seam.tovector(out);
}
