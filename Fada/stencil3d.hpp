//
//  stencil3d.hpp
//  Fada
//
//  Created by Roland Becker on 02/06/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef stencil3d_hpp
#define stencil3d_hpp

#include  "nodevector.hpp"
class SparseMatrix;

/*-------------------------------------------------*/
class Stencil3d
{
protected:
  int _nx, _ny, _nz;
  void _boundary(NodeVector& out) const;
};

/*-------------------------------------------------*/
class Stencil3d27 : public Stencil3d
{
protected:
  arma::vec::fixed<27> _coef;

public:
  Stencil3d27() {}
  Stencil3d27(const Stencil3d27& stencil) {}
  Stencil3d27(const armaicvec& n, const armavec& coef)
  {
    set_grid(n, coef);
  }

  void set_grid(const armaicvec& n, const armavec& coef);
  void jacobi       (NodeVector& out, const NodeVector& in) const;
  void gauss_seidel1(NodeVector& out, const NodeVector& in) const;
  void gauss_seidel2(NodeVector& out, const NodeVector& in) const;
  void dot(NodeVector& out, const NodeVector& in, double d) const;
  void get_sparse_matrix(SparseMatrix& sp) const;
};

#endif /* stencil3d_hpp */
