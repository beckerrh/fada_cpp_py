
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
#include  "vector.hpp"

/*-------------------------------------------------*/
class TrapezMatrix2d : public MatrixInterface
{
protected:
  int _nx, _ny;
  void _boundary(Vector& out) const;

public:
  ~TrapezMatrix2d();
  TrapezMatrix2d() : MatrixInterface() {}
  TrapezMatrix2d(const TrapezMatrix2d& trapezmatrix) : MatrixInterface(trapezmatrix) {}

  void set_grid(std::shared_ptr<GridInterface> grid);
  void jacobi       (Vector& out, const Vector& in) const;
  void gauss_seidel1(Vector& out, const Vector& in) const;
  void gauss_seidel2(Vector& out, const Vector& in) const;
  void dot(Vector& out, const Vector& in, double d) const;
  void get_sparse_matrix(SparseMatrix& sp) const;
};

/*-------------------------------------------------*/
class TrapezMatrix3d : public MatrixInterface
{
protected:
  int _nx, _ny, _nz;
  void _boundary(Vector& out) const;

public:
  ~TrapezMatrix3d();
  TrapezMatrix3d() : MatrixInterface() {}
  TrapezMatrix3d(const TrapezMatrix3d& trapezmatrix) : MatrixInterface(trapezmatrix) {}

  void set_grid(std::shared_ptr<GridInterface> grid);
  void jacobi       (Vector& out, const Vector& in) const;
  void gauss_seidel1(Vector& out, const Vector& in) const;
  void gauss_seidel2(Vector& out, const Vector& in) const;
  void dot(Vector& out, const Vector& in, double d) const;
  void get_sparse_matrix(SparseMatrix& sp) const;
};

#endif /* trapezmatrix_h */
