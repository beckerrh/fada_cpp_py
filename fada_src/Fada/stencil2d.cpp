//
//  stencil2d.cpp
//  Fada
//
//  Created by Roland Becker on 01/06/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "stencil.hpp"
#include  "construct_elements_matrix.hpp"

/*-------------------------------------------------*/
void Stencil2d9::dot(GridVector& out, const GridVector& in, double d) const
{
  arma::vec::fixed <9> coef = d * _coef;
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
void Stencil2d9::jacobi(GridVector& out, const GridVector& in) const
{
  double d0inv = 1.0 / _coef[4];
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
void Stencil2d9::gauss_seidel1(GridVector& out, const GridVector& in) const
{
  /*
   * (ix+p)*ny + iy+q < ix*ny + iy
   * p*ny +q < 0
   * p=-1 q=-1,0,1
   * p= 0 q=-1
   */
  double d0inv = 1.0 / _coef[4];
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
void Stencil2d9::gauss_seidel2(GridVector& out, const GridVector& in) const
{
  double d0inv = 1.0 / _coef[4];
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
// std::shared_ptr<MatrixInterface> Stencil2d9::get_sparse_matrix() const
void Stencil2d9::get_locations_values(arma::umat& locations, armavec& values) const
{
  // (nx+2)*(ny+2) - nx*ny + nx*ny - (nx-2)*(ny-2) = 4*(nx+ny)
  int size = 9 * (_nx - 2) * (_ny - 2) + 4 * (_nx + _ny);
  int ofsx = _ny;
  int i, j;
  Construct_Elements_Matrix ce(locations, values, size);
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
}

/*-------------------------------------------------*/
void Stencil2d5::dot(GridVector& out, const GridVector& in, double d) const
{
  arma::vec::fixed <5> coef = d * _coef;
  // std::cerr << "Stencil2d5::dot() _nx="<<_nx << "_ny=" << _ny <<"\n";
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
void Stencil2d5::jacobi(GridVector& out, const GridVector& in) const
{
  double d0inv = 1.0 / _coef[2];
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
void Stencil2d5::gauss_seidel1(GridVector& out, const GridVector& in) const
{
  /*
   * (ix+p)*ny + iy+q < ix*ny + iy
   * p*ny +q < 0
   * p=-1 q=-1,0,1
   * p= 0 q=-1
   */
  double d0inv = 1.0 / _coef[2];
  double d1    = -1.0;
  out.at(0, 0) = d0inv * in.at(0, 0);
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
void Stencil2d5::gauss_seidel2(GridVector& out, const GridVector& in) const
{
  double omega = 1.0;
  double d0inv = 1.0 / _coef[2] * omega;
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
      // std::cerr << "outtt " << out.at(ix, iy) << " in " << in.at(ix, iy);
      out.at(ix, iy) = d0inv * (
        in.at(ix, iy)
        - _coef[3] * out.at(ix, iy + 1)
        - _coef[4] * out.at(ix + 1, iy)
        );
        // std::cerr << "  outtt " << out.at(ix, iy) << "\n";
    }
  }
  // std::cerr << "out="<<out<<"\n";
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil2d5::get_locations_values(arma::umat& locations, armavec& values) const
{
  int size = 5 * (_nx - 2) * (_ny - 2) + 2 * _nx + 2 * (_ny - 2);
  int ofsx = _ny;
  // std::cerr << "_nx _ny " << _nx << " " << _ny << " ofsx " << ofsx  << "\n";
  int i, j;
  Construct_Elements_Matrix ce(locations, values, size);
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
}
