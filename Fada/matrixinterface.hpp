
//
//  matrixinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef matrixinterface_h
#define matrixinterface_h

#include  <armadillo>

class GridInterface;
class Vector;
/*-------------------------------------------------*/
class MatrixInterface
{
public:
  virtual ~MatrixInterface() {}
  MatrixInterface() {}
  MatrixInterface(const MatrixInterface& updater) {}

  virtual void set_grid(const GridInterface& grid)=0;
  virtual void jacobi       (Vector& out, const Vector& in) const=0;
  virtual void gauss_seidel1(Vector& out, const Vector& in) const=0;
  virtual void gauss_seidel2(Vector& out, const Vector& in) const=0;
  virtual void dot(Vector& out, const Vector& in, double d) const=0;
  virtual arma::sp_mat set_sparse() const=0;
};


#endif /* matrixinterface_h */
