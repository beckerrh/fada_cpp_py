
//
//  fullmatrix.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef fullmatrix_h
#define fullmatrix_h

#include  "gridvector.hpp"

class armasp;
class SparseMatrix;
/*-------------------------------------------------*/
class FullMatrix2d
{
protected:
  int _nx, _ny;
  double _vol, _dx, _dy;
  void _boundary(GridVector& out) const;

public:
  ~FullMatrix2d();
  FullMatrix2d() {}
  FullMatrix2d(const FullMatrix2d& fullmatrix) {}
  FullMatrix2d(const armaicvec& n, const armavec& dx)
  {
    set_grid(n, dx);
  }

  void set_grid(const armaicvec& n, const armavec& dx);
  void jacobi       (GridVector& out, const GridVector& in) const;
  void gauss_seidel1(GridVector& out, const GridVector& in) const;
  void gauss_seidel2(GridVector& out, const GridVector& in) const;
  void dot(GridVector& out, const GridVector& in, double d) const;
  void get_sparse_matrix(SparseMatrix& sp) const;
};

/*-------------------------------------------------*/
class FullMatrix3d
{
protected:
  int _nx, _ny, _nz;
  double _vol, _dx, _dy, _dz;
  void _boundary(GridVector& out) const;

public:
  ~FullMatrix3d();
  FullMatrix3d() {}
  FullMatrix3d(const FullMatrix3d& fullmatrix) {}
  FullMatrix3d(const armaicvec& n, const armavec& dx)
  {
    set_grid(n, dx);
  }

  void set_grid(const armaicvec& n, const armavec& dx);
//  void set_grid(std::shared_ptr<GridInterface> grid);
  void jacobi       (GridVector& out, const GridVector& in) const;
  void gauss_seidel1(GridVector& out, const GridVector& in) const;
  void gauss_seidel2(GridVector& out, const GridVector& in) const;
  void dot(GridVector& out, const GridVector& in, double d) const;
  void get_sparse_matrix(SparseMatrix& sp) const;
};

#endif /* fullmatrix_h */
