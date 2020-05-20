
//
//  fullmatrix.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef fullmatrix_h
#define fullmatrix_h

#include  "matrixinterface.hpp"
#include  "vector.hpp"

class armasp;
/*-------------------------------------------------*/
class FullMatrix2d : public MatrixInterface
{
protected:
  int _nx, _ny;
  double _vol, _dx, _dy;
  void _boundary(Vector& out) const;

public:
  FullMatrix2d() : MatrixInterface() {}
  FullMatrix2d(const FullMatrix2d& fullmatrix) : MatrixInterface(fullmatrix) {}

  void set_grid(const GridInterface& grid);

  void jacobi       (Vector& out, const Vector& in) const;
  void gauss_seidel1(Vector& out, const Vector& in) const;
  void gauss_seidel2(Vector& out, const Vector& in) const;
  void dot(Vector& out, const Vector& in, double d) const;
  void get_sparse_matrix(SparseMatrix& sp) const;
};

/*-------------------------------------------------*/
class FullMatrix3d : public MatrixInterface
{
protected:
  int _nx, _ny, _nz;
  double _vol, _dx, _dy, _dz;
  void _boundary(Vector& out) const;

public:
  FullMatrix3d() : MatrixInterface() {}
  FullMatrix3d(const FullMatrix3d& fullmatrix) : MatrixInterface(fullmatrix) {}

  void set_grid(const GridInterface& grid);

  void jacobi       (Vector& out, const Vector& in) const;
  void gauss_seidel1(Vector& out, const Vector& in) const;
  void gauss_seidel2(Vector& out, const Vector& in) const;
  void dot(Vector& out, const Vector& in, double d) const;
  void get_sparse_matrix(SparseMatrix& sp) const;
};

#endif /* fullmatrix_h */
