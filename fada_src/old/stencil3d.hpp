//
//  stencil3d.hpp
//  Fada
//
//  Created by Roland Becker on 02/06/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef stencil3d_hpp
#define stencil3d_hpp

#include  "gridvector.hpp"
#include  "seamvector.hpp"

class MatrixInterface;
/*-------------------------------------------------*/
template<int N>
class Stencil3d
{
protected:
  arma::vec::fixed<N> _coef;
  mutable SeamVector _seam;
  int _nx, _ny, _nz;
  void _boundary(GridVector& out) const;
public:
  void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const{out<<_coef;}
};

/*-------------------------------------------------*/
class Stencil3d27 : public Stencil3d<27>
{
public:
  Stencil3d27() : Stencil3d<27>() {}
  Stencil3d27(const Stencil3d27& stencil) : Stencil3d<27>(stencil) {}
  Stencil3d27(const armaicvec& n, const armavec& coef) : Stencil3d<27>()
  {
    set_grid(n, coef);
  }

  void set_grid(const armaicvec& n, const armavec& coef);
  void jacobi       (GridVector& out, const GridVector& in) const;
  void gauss_seidel1(GridVector& out, const GridVector& in) const;
  void gauss_seidel2(GridVector& out, const GridVector& in) const;
  void dot(GridVector& out, const GridVector& in, double d) const;
  void get_locations_values(arma::umat& locations, armavec& values) const;
};

/*-------------------------------------------------*/
class Stencil3d7 : public Stencil3d<7>
{
public:
  Stencil3d7() : Stencil3d<7>() {}
  Stencil3d7(const Stencil3d7& stencil) : Stencil3d<7>(stencil) {}
  Stencil3d7(const armaicvec& n, const armavec& coef) : Stencil3d<7>()
  {
    set_grid(n, coef);
  }

  void set_grid(const armaicvec& n, const armavec& coef);
  void jacobi       (GridVector& out, const GridVector& in) const;
  void gauss_seidel1(GridVector& out, const GridVector& in) const;
  void gauss_seidel2(GridVector& out, const GridVector& in) const;
  void dot(GridVector& out, const GridVector& in, double d) const;
  void get_locations_values(arma::umat& locations, armavec& values) const;
};

#endif /* stencil3d_hpp */
