//
//  smoothersimple.cpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include "../smoothersimple.hpp"
#include "../matrixinterface.hpp"
#include "../vectorinterface.hpp"

/*-------------------------------------------------*/
void SmootherSimple::set_matrix(std::shared_ptr<MatrixInterface const> matrix)
{
  _matrix = matrix;
}
/*-------------------------------------------------*/
void SmootherSimple::pre(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const
{
  out->fill(0.0);
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
void SmootherSimple::post(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const
{
  out->fill(0.0);
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
