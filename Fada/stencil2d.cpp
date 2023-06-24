//
//  stencil2d.cpp
//  Fada
//
//  Created by Roland Becker on 01/06/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include "stencil2d.hpp"
#include "sparsematrix.hpp"

/*-------------------------------------------------*/
void Stencil2d::_boundary(NodeVector& out) const
{
  // for(int ix=0;ix<_nx;ix++)
  // {
  //   out.atp(ix,0)    = 0.0;
  //   out.atp(ix,_ny-1) = 0.0;
  // }
  // for(int iy=0;iy<_ny;iy++)
  // {
  //   out.atp(0,iy)    = 0.0;
  //   out.atp(_nx-1,iy) = 0.0;
  // }
  for (int ix = 0; ix < _nx; ix++)
  {
    out.at(ix, 0)       = 0.0;
    out.at(ix, _ny - 1) = 0.0;
  }
  for (int iy = 0; iy < _ny; iy++)
  {
    out.at(0, iy)       = 0.0;
    out.at(_nx - 1, iy) = 0.0;
  }
}

/*-------------------------------------------------*/
void Stencil2d9::set_grid(const armaicvec& n, const armavec& coef)
{
  _seam.set_size(n + 2);
  assert(n.n_elem == 2);
  _nx   = n[0];
  _ny   = n[1];
  _coef = coef;
//  std::cerr << "Stencil2d9() _coef="<<_coef.t();
}

/*-------------------------------------------------*/
void Stencil2d9::dot(NodeVector& out, const NodeVector& in, double d) const
{
  arma::vec::fixed <9> coef = d * _coef;

  // for(int ix=0;ix<_nx;ix++)
  // {
  //   for(int iy=0;iy<_ny;iy++)
  //   {
  //     out.atp(ix,iy) +=
  //       coef[0]* in.atp(ix-1,iy-1)
  //     + coef[1]* in.atp(ix-1,iy+1)
  //     + coef[2]* in.atp(ix-1,iy  )
  //     + coef[3]* in.atp(ix  ,iy-1)
  //     + coef[4]* in.atp(ix  ,iy  )
  //     + coef[5]* in.atp(ix  ,iy+1)
  //     + coef[6]* in.atp(ix+1,iy-1)
  //     + coef[7]* in.atp(ix+1,iy)
  //     + coef[8]* in.atp(ix+1,iy+1);
  //   }
  // }
  _seam.fromvector(in);
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      out.at(ix, iy) +=
        coef[0] * _seam.atp(ix - 1, iy - 1)
        + coef[1] * _seam.atp(ix - 1, iy + 1)
        + coef[2] * _seam.atp(ix - 1, iy)
        + coef[3] * _seam.atp(ix, iy - 1)
        + coef[4] * _seam.atp(ix, iy)
        + coef[5] * _seam.atp(ix, iy + 1)
        + coef[6] * _seam.atp(ix + 1, iy - 1)
        + coef[7] * _seam.atp(ix + 1, iy)
        + coef[8] * _seam.atp(ix + 1, iy + 1);
    }
  }
  // Conditions aux limites ( Dirichlet)
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil2d9::jacobi(NodeVector& out, const NodeVector& in) const
{
  double d0inv = 1.0 / _coef[4];

  // for(int ix=0;ix<_nx;ix++)
  // {
  //   for(int iy=0;iy<_ny;iy++)
  //   {
  //     out.atp(ix,iy) = d0inv * in.atp(ix,iy);
  //   }
  // }
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      out.at(ix, iy) = d0inv * in.at(ix, iy);
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil2d9::gauss_seidel1(NodeVector& out, const NodeVector& in) const
{
  /*
   * (ix+p)*ny + iy+q < ix*ny + iy
   * p*ny +q < 0
   * p=-1 q=-1,0,1
   * p= 0 q=-1
   */
  double d0inv = 1.0 / _coef[4];

  // for(int ix=0;ix<_nx;ix++)
  // {
  //   for(int iy=0;iy<_ny;iy++)
  //   {
  //     out.atp(ix,iy) = d0inv * (
  //                               in.atp(ix,iy)
  //                               - _coef[0]* out.atp(ix-1,iy-1)
  //                               - _coef[1]* out.atp(ix-1,iy+1)
  //                               - _coef[2]* out.atp(ix-1,iy  )
  //                               - _coef[3]* out.atp(ix  ,iy-1)
  //                               );
  //   }
  // }
  _seam.fromvector(out);
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      out.at(ix, iy) = d0inv * (
        in.at(ix, iy)
        - _coef[0] * _seam.atp(ix - 1, iy - 1)
        - _coef[1] * _seam.atp(ix - 1, iy + 1)
        - _coef[2] * _seam.atp(ix - 1, iy)
        - _coef[3] * _seam.atp(ix, iy - 1)
        );
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil2d9::gauss_seidel2(NodeVector& out, const NodeVector& in) const
{
//  double omega = 0.8;
  double d0inv = 1.0 / _coef[4];

  // for(int ix=_nx-1;ix>=0;ix--)
  // {
  //   for(int iy=_ny-1;iy>=0;iy--)
  //   {
  //     out.atp(ix,iy) = d0inv * (
  //                               in.atp(ix,iy)
  //                               - _coef[5]* out.atp(ix  ,iy+1)
  //                               - _coef[6]* out.atp(ix+1,iy-1)
  //                               - _coef[7]* out.atp(ix+1,iy)
  //                               - _coef[8]* out.atp(ix+1,iy+1)
  //                               );
  //   }
  // }
  _seam.fromvector(out);
  for (int ix = _nx - 1; ix >= 0; ix--)
  {
    for (int iy = _ny - 1; iy >= 0; iy--)
    {
      out.at(ix, iy) = d0inv * (
        in.at(ix, iy)
        - _coef[5] * _seam.atp(ix, iy + 1)
        - _coef[6] * _seam.atp(ix + 1, iy - 1)
        - _coef[7] * _seam.atp(ix + 1, iy)
        - _coef[8] * _seam.atp(ix + 1, iy + 1)
        );
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil2d9::get_sparse_matrix(SparseMatrix& sp) const
{
  // (nx+2)*(ny+2) - nx*ny + nx*ny - (nx-2)*(ny-2) = 4*(nx+ny)
//   int size = 9*(_nx-2)*(_ny-2) + 4*(_nx+_ny);
//   int ofsx = _ny + 2;
//   int ofsp = ofsx+1;
//   int i, j;
//   Construct_Elements ce(size);
//   for(int ix=1;ix<_nx-1;ix++)
//   {
//     for(int iy=1;iy<_ny-1;iy++)
//     {
//       i = ofsx*ix + iy + ofsp;
//
//       j = ofsx*(ix-1) + iy-1 + ofsp;
//       ce.add(i, j, _coef[0]);
//       j = ofsx*(ix-1) + iy + ofsp;
//       ce.add(i, j, _coef[1]);
//       j = ofsx*(ix-1) + iy+1 + ofsp;
//       ce.add(i, j, _coef[2]);
//       j = ofsx*ix + iy-1 + ofsp;
//       ce.add(i, j, _coef[3]);
//       j = ofsx*ix + iy + ofsp;
//       ce.add(i, j, _coef[4]);
//       j = ofsx*ix + iy+1 + ofsp;
//       ce.add(i, j, _coef[5]);
//       j = ofsx*(ix+1) + iy-1 + ofsp;
//       ce.add(i, j, _coef[6]);
//       j = ofsx*(ix+1) + iy   + ofsp;
//       ce.add(i, j, _coef[7]);
//       j = ofsx*(ix+1) + iy+1 + ofsp;
//       ce.add(i, j, _coef[8]);
//     }
//   }
//   //bdry
//   for(int ix=0;ix<_nx;ix++)
//   {
//     i = ofsx*ix + 0 + ofsp;
//     ce.add(i, i, 1);
//     i = ofsx*ix + _ny-1 + ofsp;
//     ce.add(i, i, 1);
//   }
//   for(int iy=1;iy<_ny-1;iy++)
//   {
//     i = ofsx*0 + iy + ofsp;
//     ce.add(i, i, 1);
//     i = ofsx*(_nx-1) + iy + ofsp;
//     ce.add(i, i, 1);
//   }
//   //aux
//   for(int ix=0;ix<_nx+2;ix++)
//   {
//     i = ofsx*ix + 0;
//     ce.add(i, i, 1);
//     i = ofsx*ix + _ny+1;
//     ce.add(i, i, 1);
//   }
//   for(int iy=1;iy<_ny+1;iy++)
//   {
//     i = ofsx*0 + iy;
//     ce.add(i, i, 1);
//     i = ofsx*(_nx+1) + iy;
//     ce.add(i, i, 1);
//   }
// //  std::cerr << "locations i " << locations.row(0);
// //  std::cerr << "locations j " << locations.row(1);
// //  std::cerr << "values " << values.t();
//   sp.set_elements(ce.locations(), ce.values());
  // (nx+2)*(ny+2) - nx*ny + nx*ny - (nx-2)*(ny-2) = 4*(nx+ny)
  int size = 9 * (_nx - 2) * (_ny - 2) + 4 * (_nx + _ny);
  int ofsx = _ny;
  int i, j;
  Construct_Elements ce(size);

  for (int ix = 1; ix < _nx - 1; ix++)
  {
    for (int iy = 1; iy < _ny - 1; iy++)
    {
      i = ofsx * ix + iy;

      j = ofsx * (ix - 1) + iy - 1;
      ce.add(i, j, _coef[0]);
      j = ofsx * (ix - 1) + iy;
      ce.add(i, j, _coef[1]);
      j = ofsx * (ix - 1) + iy + 1;
      ce.add(i, j, _coef[2]);
      j = ofsx * ix + iy - 1;
      ce.add(i, j, _coef[3]);
      j = ofsx * ix + iy;
      ce.add(i, j, _coef[4]);
      j = ofsx * ix + iy + 1;
      ce.add(i, j, _coef[5]);
      j = ofsx * (ix + 1) + iy - 1;
      ce.add(i, j, _coef[6]);
      j = ofsx * (ix + 1) + iy;
      ce.add(i, j, _coef[7]);
      j = ofsx * (ix + 1) + iy + 1;
      ce.add(i, j, _coef[8]);
    }
  }
  //bdry
  for (int ix = 0; ix < _nx; ix++)
  {
    i = ofsx * ix + 0;
    ce.add(i, i, 1);
    i = ofsx * ix + _ny - 1;
    ce.add(i, i, 1);
  }
  for (int iy = 1; iy < _ny - 1; iy++)
  {
    i = ofsx * 0 + iy;
    ce.add(i, i, 1);
    i = ofsx * (_nx - 1) + iy;
    ce.add(i, i, 1);
  }
  //aux
  for (int ix = 0; ix < _nx + 2; ix++)
  {
    i = ofsx * ix + 0;
    ce.add(i, i, 1);
    i = ofsx * ix + _ny + 1;
    ce.add(i, i, 1);
  }
  for (int iy = 1; iy < _ny + 1; iy++)
  {
    i = ofsx * 0 + iy;
    ce.add(i, i, 1);
    i = ofsx * (_nx + 1) + iy;
    ce.add(i, i, 1);
  }
//  std::cerr << "locations i " << locations.row(0);
//  std::cerr << "locations j " << locations.row(1);
//  std::cerr << "values " << values.t();
  sp.set_elements(ce.locations(), ce.values());
}

/*-------------------------------------------------*/
void Stencil2d5::set_grid(const armaicvec& n, const armavec& coef)
{
  _seam.set_size(n + 2);
  assert(n.n_elem == 2);
  _nx   = n[0];
  _ny   = n[1];
  _coef = coef;
//  std::cerr << "Stencil2d9() _coef="<<_coef.t();
}

/*-------------------------------------------------*/
void Stencil2d5::dot(NodeVector& out, const NodeVector& in, double d) const
{
  arma::vec::fixed <5> coef = d * _coef;

  // for(int ix=0;ix<_nx;ix++)
  // {
  //   for(int iy=0;iy<_ny;iy++)
  //   {
  //     out.atp(ix,iy) +=
  //       coef[0]* in.atp(ix-1,iy  )
  //     + coef[1]* in.atp(ix  ,iy-1)
  //     + coef[2]* in.atp(ix  ,iy  )
  //     + coef[3]* in.atp(ix  ,iy+1)
  //     + coef[4]* in.atp(ix+1,iy  );
  //   }
  // }
  _seam.fromvector(in);
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      out.at(ix, iy) +=
        coef[0] * _seam.atp(ix - 1, iy)
        + coef[1] * _seam.atp(ix, iy - 1)
        + coef[2] * _seam.atp(ix, iy)
        + coef[3] * _seam.atp(ix, iy + 1)
        + coef[4] * _seam.atp(ix + 1, iy);
    }
  }
  // Conditions aux limites ( Dirichlet)
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil2d5::jacobi(NodeVector& out, const NodeVector& in) const
{
  double d0inv = 1.0 / _coef[2];

  // for(int ix=0;ix<_nx;ix++)
  // {
  //   for(int iy=0;iy<_ny;iy++)
  //   {
  //     out.atp(ix,iy) = d0inv * in.atp(ix,iy);
  //   }
  // }
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      out.at(ix, iy) = d0inv * in.at(ix, iy);
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil2d5::gauss_seidel1(NodeVector& out, const NodeVector& in) const
{
  /*
   * (ix+p)*ny + iy+q < ix*ny + iy
   * p*ny +q < 0
   * p=-1 q=-1,0,1
   * p= 0 q=-1
   */
  double d0inv = 1.0 / _coef[2];
  double d1    = -1.0;

  // for(int ix=0;ix<_nx;ix++)
  // {
  //   for(int iy=0;iy<_ny;iy++)
  //   {
  //     out.atp(ix,iy) = d0inv * (
  //                               in.atp(ix,iy)
  //                               -_coef[0]* out.atp(ix-1,iy  )
  //                               -_coef[1]* out.atp(ix  ,iy-1)
  //                               );
  //   }
  // }
  out.at(0, 0) = d0inv *in.at(0, 0);
  for (int iy = 1; iy < _ny; iy++)
  {
    out.at(0, iy) = d0inv * (
      in.at(0, iy)
      - _coef[1] * out.at(0, iy - 1)
      );
  }
  for (int ix = 1; ix < _nx; ix++)
  {
    out.at(ix, 0) = d0inv * (
      in.at(ix, 0)
      - _coef[0] * out.at(ix - 1, 0)
      );
    for (int iy = 1; iy < _ny; iy++)
    {
      out.at(ix, iy) = d0inv * (
        in.at(ix, iy)
        - _coef[0] * out.at(ix - 1, iy)
        - _coef[1] * out.at(ix, iy - 1)
        );
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil2d5::gauss_seidel2(NodeVector& out, const NodeVector& in) const
{
//  double omega = 0.8;
  double omega = 1.0;
  double d0inv = 1.0 / _coef[2] * omega;

  // for(int ix=_nx-1;ix>=0;ix--)
  // {
  //   for(int iy=_ny-1;iy>=0;iy--)
  //   {
  //     out.atp(ix,iy) = d0inv * (
  //                               in.atp(ix,iy)
  //                               - _coef[3]* out.atp(ix  ,iy+1)
  //                               - _coef[4]* out.atp(ix+1,iy  )
  //                               );
  //   }
  // }
  out.at(_nx - 1, _ny - 1) = d0inv * in.at(_nx - 1, _ny - 1);
  for (int iy = _ny - 2; iy >= 0; iy--)
  {
    out.at(_nx - 1, iy) = d0inv * (
      in.at(_nx - 1, iy)
      - _coef[3] * out.at(_nx - 1, iy + 1)
      );
  }
  for (int ix = _nx - 2; ix >= 0; ix--)
  {
    out.at(ix, _ny - 1) = d0inv * (
      in.at(ix, _ny - 1)
      - _coef[4] * out.at(ix + 1, _ny - 1)
      );
    for (int iy = _ny - 2; iy >= 0; iy--)
    {
      out.at(ix, iy) = d0inv * (
        in.at(ix, iy)
        - _coef[3] * out.at(ix, iy + 1)
        - _coef[4] * out.at(ix + 1, iy)
        );
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil2d5::get_sparse_matrix(SparseMatrix& sp) const
{
//   // (nx+2)*(ny+2) - nx*ny + nx*ny - (nx-2)*(ny-2) = 4*(nx+ny)
//   int size = 5*(_nx-2)*(_ny-2) + 4*(_nx+_ny);
//   int ofsx = _ny + 2;
//   int ofsp = ofsx+1;
// //  std::cerr << "_nx _ny " << _nx << " " << _ny << " ofsx " << ofsx << " ofsp " << ofsp << "\n";
//   int i, j;
//   Construct_Elements ce(size);
//   for(int ix=1;ix<_nx-1;ix++)
//   {
//     for(int iy=1;iy<_ny-1;iy++)
//     {
//       i = ofsx*ix + iy + ofsp;
//
//       j = ofsx*(ix-1) + iy + ofsp;
//       ce.add(i, j, _coef[0]);
//       j = ofsx*ix + iy-1 + ofsp;
//       ce.add(i, j, _coef[1]);
//       j = ofsx*ix + iy + ofsp;
//       ce.add(i, j, _coef[2]);
//       j = ofsx*ix + iy+1 + ofsp;
//       ce.add(i, j, _coef[3]);
//       j = ofsx*(ix+1) + iy   + ofsp;
//       ce.add(i, j, _coef[4]);
//     }
//   }
//   //bdry
//   for(int ix=0;ix<_nx;ix++)
//   {
//     i = ofsx*ix + 0 + ofsp;
//     ce.add(i, i, 1);
//     i = ofsx*ix + _ny-1 + ofsp;
//     ce.add(i, i, 1);
//   }
//   for(int iy=1;iy<_ny-1;iy++)
//   {
//     i = ofsx*0 + iy + ofsp;
//     ce.add(i, i, 1);
//     i = ofsx*(_nx-1) + iy + ofsp;
//     ce.add(i, i, 1);
//   }
//   //aux
//   for(int ix=0;ix<_nx+2;ix++)
//   {
//     i = ofsx*ix + 0;
//     ce.add(i, i, 1);
//     i = ofsx*ix + _ny+1;
//     ce.add(i, i, 1);
//   }
//   for(int iy=1;iy<_ny+1;iy++)
//   {
//     i = ofsx*0 + iy;
//     ce.add(i, i, 1);
//     i = ofsx*(_nx+1) + iy;
//     ce.add(i, i, 1);
//   }
// //  std::cerr << "locations i " << locations.row(0);
// //  std::cerr << "locations j " << locations.row(1);
// //  std::cerr << "values " << values.t();
//   sp.set_elements(ce.locations(), ce.values());



  int size = 5 * (_nx - 2) * (_ny - 2) + 2 * _nx + 2 * (_ny - 2);
  int ofsx = _ny;
  // std::cerr << "_nx _ny " << _nx << " " << _ny << " ofsx " << ofsx  << "\n";
  int i, j;
  Construct_Elements ce(size);

  for (int ix = 1; ix < _nx - 1; ix++)
  {
    for (int iy = 1; iy < _ny - 1; iy++)
    {
      i = ofsx * ix + iy;

      j = ofsx * (ix - 1) + iy;
      ce.add(i, j, _coef[0]);
      j = ofsx * ix + iy - 1;
      ce.add(i, j, _coef[1]);
      j = ofsx * ix + iy;
      ce.add(i, j, _coef[2]);
      j = ofsx * ix + iy + 1;
      ce.add(i, j, _coef[3]);
      j = ofsx * (ix + 1) + iy;
      ce.add(i, j, _coef[4]);
    }
  }
//bdry
  for (int ix = 0; ix < _nx; ix++)
  {
    i = ofsx * ix + 0;
    ce.add(i, i, 1);
    i = ofsx * ix + _ny - 1;
    ce.add(i, i, 1);
  }
  for (int iy = 1; iy < _ny - 1; iy++)
  {
    i = ofsx * 0 + iy;
    ce.add(i, i, 1);
    i = ofsx * (_nx - 1) + iy;
    ce.add(i, i, 1);
  }
//  std::cerr << "locations i " << locations.row(0);
//  std::cerr << "locations j " << locations.row(1);
//  std::cerr << "values " << values.t();
  sp.set_elements(ce.locations(), ce.values());
}
