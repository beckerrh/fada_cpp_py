//
//  stencil3d.cpp
//  Fada
//
//  Created by Roland Becker on 02/06/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "stencil.hpp"
#include  "../construct_elements_matrix.hpp"

// /*-------------------------------------------------*/
// template<int N>
// void Stencil3d<N>::_boundary(GridVector& out) const
// {
//   // for(int ix=0;ix<_nx;ix++)
//   // {
//   //   for(int iy=0;iy<_ny;iy++)
//   //   {
//   //     out.atp(ix,iy,0)    = 0.0;
//   //     out.atp(ix,iy,_nz-1) = 0.0;
//   //   }
//   // }
//   // for(int ix=0;ix<_nx;ix++)
//   // {
//   //   for(int iz=0;iz<_nz;iz++)
//   //   {
//   //     out.atp(ix,0,   iz) = 0.0;
//   //     out.atp(ix,_ny-1,iz) = 0.0;
//   //   }
//   // }
//   // for(int iy=0;iy<_ny;iy++)
//   // {
//   //   for(int iz=0;iz<_nz;iz++)
//   //   {
//   //     out.atp(0,   iy,iz) = 0.0;
//   //     out.atp(_nx-1,iy,iz) = 0.0;
//   //   }
//   // }
//   for (int ix = 0; ix < _nx; ix++)
//   {
//     for (int iy = 0; iy < _ny; iy++)
//     {
//       out.at(ix, iy, 0)       = 0.0;
//       out.at(ix, iy, _nz - 1) = 0.0;
//     }
//   }
//   for (int ix = 0; ix < _nx; ix++)
//   {
//     for (int iz = 0; iz < _nz; iz++)
//     {
//       out.at(ix, 0, iz)       = 0.0;
//       out.at(ix, _ny - 1, iz) = 0.0;
//     }
//   }
//   for (int iy = 0; iy < _ny; iy++)
//   {
//     for (int iz = 0; iz < _nz; iz++)
//     {
//       out.at(0, iy, iz)       = 0.0;
//       out.at(_nx - 1, iy, iz) = 0.0;
//     }
//   }
// }

/*-------------------------------------------------*/
void Stencil3d27::dot(GridVector& out, const GridVector& in, double d) const
{
  _seam.fromvector(in);
  arma::vec::fixed <27> coef = d * _coef;
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      for (int iz = 0; iz < _nz; iz++)
      {
        out.at(ix, iy, iz) +=
          coef[0] * _seam.atp(ix - 1, iy - 1, iz - 1)
          + coef[1] * _seam.atp(ix - 1, iy - 1, iz)
          + coef[2] * _seam.atp(ix - 1, iy - 1, iz + 1)
          + coef[3] * _seam.atp(ix - 1, iy, iz - 1)
          + coef[4] * _seam.atp(ix - 1, iy, iz)
          + coef[5] * _seam.atp(ix - 1, iy, iz + 1)
          + coef[6] * _seam.atp(ix - 1, iy + 1, iz - 1)
          + coef[7] * _seam.atp(ix - 1, iy + 1, iz)
          + coef[8] * _seam.atp(ix - 1, iy + 1, iz + 1)
          + coef[9] * _seam.atp(ix, iy - 1, iz - 1)
          + coef[10] * _seam.atp(ix, iy - 1, iz)
          + coef[11] * _seam.atp(ix, iy - 1, iz + 1)
          + coef[12] * _seam.atp(ix, iy, iz - 1)
          + coef[13] * _seam.atp(ix, iy, iz)
          + coef[14] * _seam.atp(ix, iy, iz + 1)
          + coef[15] * _seam.atp(ix, iy + 1, iz - 1)
          + coef[16] * _seam.atp(ix, iy + 1, iz)
          + coef[17] * _seam.atp(ix, iy + 1, iz + 1)
          + coef[18] * _seam.atp(ix + 1, iy - 1, iz - 1)
          + coef[19] * _seam.atp(ix + 1, iy - 1, iz)
          + coef[20] * _seam.atp(ix + 1, iy - 1, iz + 1)
          + coef[21] * _seam.atp(ix + 1, iy, iz - 1)
          + coef[22] * _seam.atp(ix + 1, iy, iz)
          + coef[23] * _seam.atp(ix + 1, iy, iz + 1)
          + coef[24] * _seam.atp(ix + 1, iy + 1, iz - 1)
          + coef[25] * _seam.atp(ix + 1, iy + 1, iz)
          + coef[26] * _seam.atp(ix + 1, iy + 1, iz + 1)
        ;
      }
    }
  }
  // Conditions aux limites ( Dirichlet)
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil3d27::jacobi(GridVector& out, const GridVector& in) const
{
  double d0inv = 1.0 / _coef[13];

  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      for (int iz = 0; iz < _nz; iz++)
      {
        out.at(ix, iy, iz) = d0inv * in.at(ix, iy, iz);
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil3d27::gauss_seidel1(GridVector& out, const GridVector& in) const
{
  /*
   * (ix+p)*ny*nz + (iy+q)*nz + iz+r < ix*ny*nz + iy*nz + iz
   * p*ny*nz +q*nz +r < 0
   * p=-1 q=-1,0,1  r=-1,0,1
   * p= 0 q=-1 r=-1,0,1 q=0 r=-1
   */
  _seam.fromvector(out);
  double d0inv = 1.0 / _coef[13];
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      for (int iz = 0; iz < _nz; iz++)
      {
        out.at(ix, iy, iz) = d0inv * (
          in.at(ix, iy, iz)
          - _coef[0] * _seam.atp(ix - 1, iy - 1, iz - 1)
          - _coef[1] * _seam.atp(ix - 1, iy - 1, iz)
          - _coef[2] * _seam.atp(ix - 1, iy - 1, iz + 1)
          - _coef[3] * _seam.atp(ix - 1, iy, iz - 1)
          - _coef[4] * _seam.atp(ix - 1, iy, iz)
          - _coef[5] * _seam.atp(ix - 1, iy, iz + 1)
          - _coef[6] * _seam.atp(ix - 1, iy + 1, iz - 1)
          - _coef[7] * _seam.atp(ix - 1, iy + 1, iz)
          - _coef[8] * _seam.atp(ix - 1, iy + 1, iz + 1)
          - _coef[9] * _seam.atp(ix, iy - 1, iz - 1)
          - _coef[10] * _seam.atp(ix, iy - 1, iz)
          - _coef[11] * _seam.atp(ix, iy - 1, iz + 1)
          - _coef[12] * _seam.atp(ix, iy, iz - 1)
          );
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil3d27::gauss_seidel2(GridVector& out, const GridVector& in) const
{
  _seam.fromvector(out);
  double d0inv = 1.0 / _coef[13];
  for (int ix = _nx - 1; ix >= 0; ix--)
  {
    for (int iy = _ny - 1; iy >= 0; iy--)
    {
      for (int iz = _nz - 1; iz >= 0; iz--)
      {
        out.at(ix, iy, iz) = d0inv * (
          in.at(ix, iy, iz)
          - _coef[14] * _seam.atp(ix, iy, iz + 1)
          - _coef[15] * _seam.atp(ix, iy + 1, iz - 1)
          - _coef[16] * _seam.atp(ix, iy + 1, iz)
          - _coef[17] * _seam.atp(ix, iy + 1, iz + 1)
          - _coef[18] * _seam.atp(ix + 1, iy - 1, iz - 1)
          - _coef[19] * _seam.atp(ix + 1, iy - 1, iz)
          - _coef[20] * _seam.atp(ix + 1, iy - 1, iz + 1)
          - _coef[21] * _seam.atp(ix + 1, iy, iz - 1)
          - _coef[22] * _seam.atp(ix + 1, iy, iz)
          - _coef[23] * _seam.atp(ix + 1, iy, iz + 1)
          - _coef[24] * _seam.atp(ix + 1, iy + 1, iz - 1)
          - _coef[25] * _seam.atp(ix + 1, iy + 1, iz)
          - _coef[26] * _seam.atp(ix + 1, iy + 1, iz + 1)
          );
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
// std::shared_ptr<MatrixInterface> Stencil3d27::get_sparse_matrix() const
void Stencil3d27::get_locations_values(arma::umat& locations, armavec& values) const
{
  int size = 27 * (_nx - 2) * (_ny - 2) * (_nz - 2) + 16 + 4 * (_nx * _ny + _nx * _nz + _ny * _nz);
  int ofsy = _nz;
  int ofsx = ofsy * _ny;
//  std::cerr << "ofsx " << ofsx << " ofsy " << ofsy << " ofsp " << ofsp<< " size " << size << "\n";
  int i, j;
  Construct_Elements_Matrix ce(locations, values, size);
  for (int ix = 1; ix < _nx - 1; ix++)
  {
    for (int iy = 1; iy < _ny - 1; iy++)
    {
      for (int iz = 1; iz < _nz - 1; iz++)
      {
        i = ofsx * ix + ofsy * iy + iz;
        j = ofsx * (ix - 1) + ofsy * (iy - 1) + iz - 1;
        ce.add(i, j, _coef[0]);
        j = ofsx * (ix - 1) + ofsy * (iy - 1) + iz;
        ce.add(i, j, _coef[1]);
        j = ofsx * (ix - 1) + ofsy * (iy - 1) + iz + 1;
        ce.add(i, j, _coef[2]);
        j = ofsx * (ix - 1) + ofsy * (iy) + iz - 1;
        ce.add(i, j, _coef[3]);
        j = ofsx * (ix - 1) + ofsy * (iy) + iz;
        ce.add(i, j, _coef[4]);
        j = ofsx * (ix - 1) + ofsy * (iy) + iz + 1;
        ce.add(i, j, _coef[5]);
        j = ofsx * (ix - 1) + ofsy * (iy + 1) + iz - 1;
        ce.add(i, j, _coef[6]);
        j = ofsx * (ix - 1) + ofsy * (iy + 1) + iz;
        ce.add(i, j, _coef[7]);
        j = ofsx * (ix - 1) + ofsy * (iy + 1) + iz + 1;
        ce.add(i, j, _coef[8]);

        j = ofsx * (ix) + ofsy * (iy - 1) + iz - 1;
        ce.add(i, j, _coef[9]);
        j = ofsx * (ix) + ofsy * (iy - 1) + iz;
        ce.add(i, j, _coef[10]);
        j = ofsx * (ix) + ofsy * (iy - 1) + iz + 1;
        ce.add(i, j, _coef[11]);
        j = ofsx * (ix) + ofsy * (iy) + iz - 1;
        ce.add(i, j, _coef[12]);
        j = ofsx * (ix) + ofsy * (iy) + iz;
        ce.add(i, j, _coef[13]);
        j = ofsx * (ix) + ofsy * (iy) + iz + 1;
        ce.add(i, j, _coef[14]);
        j = ofsx * (ix) + ofsy * (iy + 1) + iz - 1;
        ce.add(i, j, _coef[15]);
        j = ofsx * (ix) + ofsy * (iy + 1) + iz;
        ce.add(i, j, _coef[16]);
        j = ofsx * (ix) + ofsy * (iy + 1) + iz + 1;
        ce.add(i, j, _coef[17]);

        j = ofsx * (ix + 1) + ofsy * (iy - 1) + iz - 1;
        ce.add(i, j, _coef[18]);
        j = ofsx * (ix + 1) + ofsy * (iy - 1) + iz;
        ce.add(i, j, _coef[19]);
        j = ofsx * (ix + 1) + ofsy * (iy - 1) + iz + 1;
        ce.add(i, j, _coef[20]);
        j = ofsx * (ix + 1) + ofsy * (iy) + iz - 1;
        ce.add(i, j, _coef[21]);
        j = ofsx * (ix + 1) + ofsy * (iy) + iz;
        ce.add(i, j, _coef[22]);
        j = ofsx * (ix + 1) + ofsy * (iy) + iz + 1;
        ce.add(i, j, _coef[23]);
        j = ofsx * (ix + 1) + ofsy * (iy + 1) + iz - 1;
        ce.add(i, j, _coef[24]);
        j = ofsx * (ix + 1) + ofsy * (iy + 1) + iz;
        ce.add(i, j, _coef[25]);
        j = ofsx * (ix + 1) + ofsy * (iy + 1) + iz + 1;
        ce.add(i, j, _coef[26]);
      }
    }
  }
  //bdry
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      i = ofsx * ix + ofsy * iy + 0;
      ce.add(i, i, 1);
      i = ofsx * ix + ofsy * iy + _nz - 1;
      ce.add(i, i, 1);
    }
  }
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iz = 1; iz < _nz - 1; iz++)
    {
      i = ofsx * ix + ofsy * 0 + iz;
      ce.add(i, i, 1);
      i = ofsx * ix + ofsy * (_ny - 1) + iz;
      ce.add(i, i, 1);
    }
  }
  for (int iy = 1; iy < _ny - 1; iy++)
  {
    for (int iz = 1; iz < _nz - 1; iz++)
    {
      i = ofsx * 0 + ofsy * iy + iz;
      ce.add(i, i, 1);
      i = ofsx * (_nx - 1) + ofsy * iy + iz;
      ce.add(i, i, 1);
    }
  }
  //aux
  for (int ix = 0; ix < _nx + 2; ix++)
  {
    for (int iy = 0; iy < _ny + 2; iy++)
    {
      i = ofsx * ix + ofsy * iy + 0;
      ce.add(i, i, 1);
      i = ofsx * ix + ofsy * iy + _nz + 1;
      ce.add(i, i, 1);
    }
  }
  for (int ix = 0; ix < _nx + 2; ix++)
  {
    for (int iz = 1; iz < _nz + 1; iz++)
    {
      i = ofsx * ix + ofsy * 0 + iz;
      ce.add(i, i, 1);
      i = ofsx * ix + ofsy * (_ny + 1) + iz;
      ce.add(i, i, 1);
    }
  }
  for (int iy = 1; iy < _ny + 1; iy++)
  {
    for (int iz = 1; iz < _nz + 1; iz++)
    {
      i = ofsx * 0 + ofsy * iy + iz;
      ce.add(i, i, 1);
      i = ofsx * (_nx + 1) + ofsy * iy + iz;
      ce.add(i, i, 1);
    }
  }
}

/*-------------------------------------------------*/
void Stencil3d7::dot(GridVector& out, const GridVector& in, double d) const
{
  _seam.fromvector(in);
// std::cerr << "**_seam" << _seam.data().t() << "\n";
  arma::vec::fixed <7> coef = d * _coef;
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      for (int iz = 0; iz < _nz; iz++)
      {
        out.at(ix, iy, iz) +=
          coef[0] * _seam.atp(ix - 1, iy, iz)
          + coef[1] * _seam.atp(ix, iy - 1, iz)
          + coef[2] * _seam.atp(ix, iy, iz - 1)
          + coef[3] * _seam.atp(ix, iy, iz)
          + coef[4] * _seam.atp(ix, iy, iz + 1)
          + coef[5] * _seam.atp(ix, iy + 1, iz)
          + coef[6] * _seam.atp(ix + 1, iy, iz)
        ;
      }
    }
  }
  // Conditions aux limites ( Dirichlet)
  // std::cerr << "***out" << out.data().t() << "\n";
  _boundary(out);
  // std::cerr << "***out" << out.data().t() << "\n";
}

/*-------------------------------------------------*/
void Stencil3d7::jacobi(GridVector& out, const GridVector& in) const
{
  double d0inv = 1.0 / _coef[3];

  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      for (int iz = 0; iz < _nz; iz++)
      {
        out.at(ix, iy, iz) = d0inv * in.at(ix, iy, iz);
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil3d7::gauss_seidel1(GridVector& out, const GridVector& in) const
{
  double d0inv = 1.0 / _coef[3];
  out.at(0, 0, 0) = d0inv * in.at(0, 0, 0);
  for (int iz = 1; iz < _nz; iz++)
  {
    out.at(0, 0, iz) = d0inv * (
      in.at(0, 0, iz)
      - _coef[2] * out.at(0, 0, iz - 1)
      );
  }
  for (int iy = 1; iy < _ny; iy++)
  {
    out.at(0, iy, 0) = d0inv * (
      in.at(0, iy, 0)
      - _coef[1] * out.at(0, iy - 1, 0)
      );
    for (int iz = 1; iz < _nz; iz++)
    {
      out.at(0, iy, iz) = d0inv * (
        in.at(0, iy, iz)
        - _coef[1] * out.at(0, iy - 1, iz)
        - _coef[2] * out.at(0, iy, iz - 1)
        );
    }
  }
  for (int ix = 1; ix < _nx; ix++)
  {
    out.at(ix, 0, 0) = d0inv * (
      in.at(ix, 0, 0)
      - _coef[0] * out.at(ix - 1, 0, 0)
      );
    for (int iz = 1; iz < _nz; iz++)
    {
      out.at(ix, 0, iz) = d0inv * (
        in.at(ix, 0, iz)
        - _coef[0] * out.at(ix - 1, 0, iz)
        - _coef[2] * out.at(ix, 0, iz - 1)
        );
    }
    for (int iy = 1; iy < _ny; iy++)
    {
      out.at(ix, iy, 0) = d0inv * (
        in.at(ix, iy, 0)
        - _coef[0] * out.at(ix - 1, iy, 0)
        - _coef[1] * out.at(ix, iy - 1, 0)
        );
      for (int iz = 1; iz < _nz; iz++)
      {
        out.at(ix, iy, iz) = d0inv * (
          in.at(ix, iy, iz)
          - _coef[0] * out.at(ix - 1, iy, iz)
          - _coef[1] * out.at(ix, iy - 1, iz)
          - _coef[2] * out.at(ix, iy, iz - 1)
          );
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil3d7::gauss_seidel2(GridVector& out, const GridVector& in) const
{
  double d0inv = 1.0/_coef[3];
  out.at(_nx - 1, _ny - 1, _nz - 1) = d0inv * in.at(_nx - 1, _ny - 1, _nz - 1);
  for (int iz = _nz - 2; iz >= 0; iz--)
  {
    out.at(_nx - 1, _ny - 1, iz) = d0inv * (
      in.at(_nx - 1, _ny - 1, iz)
      - _coef[4] * out.at(_nx - 1, _ny - 1, iz + 1)
      );
  }
  for (int iy = _ny - 1; iy >= 0; iy--)
  {
    out.at(_nx - 1, iy, _nz - 1) = d0inv * (
      in.at(_nx - 1, iy, _nz - 1)
      - _coef[5] * out.at(_nx - 1, iy + 1, _nz - 1)
      );
    for (int iz = _nz - 2; iz >= 0; iz--)
    {
      out.at(_nx - 1, iy, iz) = d0inv * (
        in.at(_nx - 1, iy, iz)
        - _coef[4] * out.at(_nx - 1, iy, iz + 1)
        - _coef[5] * out.at(_nx - 1, iy + 1, iz)
        );
    }
  }
  for (int ix = _nx - 2; ix >= 0; ix--)
  {
    for (int iy = _ny - 1; iy >= 0; iy--)
    {
      for (int iz = _nz - 1; iz >= 0; iz--)
      {
        out.at(ix, iy, iz) = d0inv * (
          in.at(ix, iy, iz)
          - _coef[4] * out.at(ix, iy, iz + 1)
          - _coef[5] * out.at(ix, iy + 1, iz)
          - _coef[6] * out.at(ix + 1, iy, iz)
          );
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil3d7::get_locations_values(arma::umat& locations, armavec& values) const
{
  int size = 7 * (_nx - 2) * (_ny - 2) * (_nz - 2) + 2 * _nx * _ny + 2 * _nx * (_nz - 2) + 2 * (_ny - 2) * (_nz - 2);
  int ofsy = _nz;
  int ofsx = ofsy * _ny;
//  std::cerr << "ofsx " << ofsx << " ofsy " << ofsy << " ofsp " << ofsp<< " size " << size << "\n";
  int i, j;
  Construct_Elements_Matrix ce(locations, values, size);
  for (int ix = 1; ix < _nx - 1; ix++)
  {
    for (int iy = 1; iy < _ny - 1; iy++)
    {
      for (int iz = 1; iz < _nz - 1; iz++)
      {
        i = ofsx * ix + ofsy * iy + iz;
        j = ofsx * (ix - 1) + ofsy * (iy) + iz;
        ce.add(i, j, _coef[0]);

        j = ofsx * (ix) + ofsy * (iy - 1) + iz;
        ce.add(i, j, _coef[1]);
        j = ofsx * (ix) + ofsy * (iy) + iz - 1;
        ce.add(i, j, _coef[2]);
        j = ofsx * (ix) + ofsy * (iy) + iz;
        ce.add(i, j, _coef[3]);
        j = ofsx * (ix) + ofsy * (iy) + iz + 1;
        ce.add(i, j, _coef[4]);
        j = ofsx * (ix) + ofsy * (iy + 1) + iz;
        ce.add(i, j, _coef[5]);
        j = ofsx * (ix + 1) + ofsy * (iy) + iz;
        ce.add(i, j, _coef[6]);
      }
    }
  }
  //bdry
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iy = 0; iy < _ny; iy++)
    {
      i = ofsx * ix + ofsy * iy + 0;
      ce.add(i, i, 1);
      i = ofsx * ix + ofsy * iy + _nz - 1;
      ce.add(i, i, 1);
    }
  }
  for (int ix = 0; ix < _nx; ix++)
  {
    for (int iz = 1; iz < _nz - 1; iz++)
    {
      i = ofsx * ix + ofsy * 0 + iz;
      ce.add(i, i, 1);
      i = ofsx * ix + ofsy * (_ny - 1) + iz;
      ce.add(i, i, 1);
    }
  }
  for (int iy = 1; iy < _ny - 1; iy++)
  {
    for (int iz = 1; iz < _nz - 1; iz++)
    {
      i = ofsx * 0 + ofsy * iy + iz;
      ce.add(i, i, 1);
      i = ofsx * (_nx - 1) + ofsy * iy + iz;
      ce.add(i, i, 1);
    }
  }
  // std::cerr << "locations i " << locations.row(0);
  // std::cerr << "locations j " << locations.row(1);
  // std::cerr << "values " << values.t();
}
