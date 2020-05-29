//
//  smoothersimple.cpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include "smoothersimple.hpp"
#include "matrixinterface.hpp"
#include "vector.hpp"

/*-------------------------------------------------*/
void SmootherSimple::set_matrix(std::shared_ptr<MatrixInterface> matrix)
{
  _matrix = matrix;
}
/*-------------------------------------------------*/
void SmootherSimple::pre(Vector& out, const Vector& in) const
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
/*-------------------------------------------------*/
void SmootherSimple::post(Vector& out, const Vector& in) const
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
  }}
