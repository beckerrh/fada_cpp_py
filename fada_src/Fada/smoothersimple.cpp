//
//  smoothersimple.cpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include "smoothersimple.hpp"
#include "matrixinterface.hpp"
#include "vectorinterface.hpp"
#include "sparsematrix.hpp"
#include "sparsematrix_arma.hpp"

/*-------------------------------------------------*/
template<typename T>
void SmootherSimple<T>::set_matrix(std::shared_ptr<MatrixInterface const> matrix)
{
  _matrix = std::dynamic_pointer_cast<T const>(matrix);
  assert(_matrix);
  //   auto stencil = std::dynamic_pointer_cast<FemAndMatrixAndSmootherInterface const>(matrix);
  // _matrix = matrix;
}
/*-------------------------------------------------*/
template<typename T>
void SmootherSimple<T>::presmooth(armavec& out, const armavec& in) const
{
  out.fill(0.0);
  if(_type=="Jac")
  {
    _matrix->jacobi(out, in);
  }
  else if(_type=="GS")
  {
    _matrix->gauss_seidel2(out, in);
  }
  else if(_type=="GS1")
  {
    _matrix->gauss_seidel1(out, in);
  }
  else if(_type=="GS2")
  {
    _matrix->gauss_seidel2(out, in);
  }
  else
  {
    throw std::runtime_error("unknown smoother " + _type);
  }
}
/*-------------------------------------------------*/
template<typename T>
void SmootherSimple<T>::postsmooth(armavec& out, const armavec& in) const
{
  out.fill(0.0);
  if(_type=="Jac")
  {
    _matrix->jacobi(out, in);
  }
  else if(_type=="GS")
  {
    _matrix->gauss_seidel1(out, in);
  }
  else if(_type=="GS1")
  {
    _matrix->gauss_seidel1(out, in);
  }
  else if(_type=="GS2")
  {
    _matrix->gauss_seidel2(out, in);
  }
  else
  {
    throw std::runtime_error("unknown smoother " + _type);
  }
}

template class SmootherSimple<SparseMatrix>;
template class SmootherSimple<SparseMatrix_arma>;
