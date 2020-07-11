//
//  stencil2d.hpp
//  Fada
//
//  Created by Roland Becker on 01/06/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef stencil2d_hpp
#define stencil2d_hpp

#include  "nodevector.hpp"
class SparseMatrix;

/*-------------------------------------------------*/
class Stencil2d
{
protected:
  int _nx, _ny;
  void _boundary(NodeVector& out) const;
};

/*-------------------------------------------------*/
class Stencil2d9 : public Stencil2d
{
protected:
  arma::vec::fixed<9> _coef;

public:
  Stencil2d9() {}
  Stencil2d9(const Stencil2d9& stencil) {}
  Stencil2d9(const armaicvec& n, const armavec& coef)
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

/*-------------------------------------------------*/
class Stencil2d5 : public Stencil2d
{
protected:
  arma::vec::fixed<5> _coef;

public:
  Stencil2d5() {}
  Stencil2d5(const Stencil2d5& stencil) {}
  Stencil2d5(const armaicvec& n, const armavec& coef)
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

#endif /* stencil2d_hpp */
