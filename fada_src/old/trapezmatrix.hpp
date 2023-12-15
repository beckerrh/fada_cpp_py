
//
//  trapezmatrix.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef trapezmatrix_h
#define trapezmatrix_h

#include  "matrixinterface.hpp"
#include  "uniformmultigrid.hpp"
#include  "uniformgrid.hpp"
#include  "gridvector.hpp"

/*-------------------------------------------------*/
class TrapezMatrix2d
{
protected:
  int _nx, _ny;
  double _vol, _dx, _dy;
  void _boundary(GridVector& out) const;

public:
  ~TrapezMatrix2d();
  TrapezMatrix2d() {}
  TrapezMatrix2d(const TrapezMatrix2d& trapezmatrix)  {}
  TrapezMatrix2d(const armaicvec& n, const armavec& dx)
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

/*-------------------------------------------------*/
class TrapezMatrix3d
{
protected:
  int _nx, _ny, _nz;
  double _vol, _dx, _dy, _dz;
  void _boundary(GridVector& out) const;

public:
  ~TrapezMatrix3d();
  TrapezMatrix3d() {}
  TrapezMatrix3d(const TrapezMatrix3d& trapezmatrix) {}
  TrapezMatrix3d(const armaicvec& n, const armavec& dx)
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

#endif /* trapezmatrix_h */
