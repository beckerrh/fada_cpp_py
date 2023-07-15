//
//  q1.cpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <cassert>
#include  "modelbase.hpp"

/*-------------------------------------------------*/
ModelBase::ModelBase(const std::map <std::string, std::string>& parameters, std::shared_ptr<BoundaryConditions const> boundaryconditions) : _ug(nullptr), _boundaryconditions(boundaryconditions)
{
  _stenciltype  = "Trapez";
  _matrixtype   = "stencil";
  _smoothertype   = "matrix";
  _smoother     = "GS";
  _coarsesolver = "direct";
  _transfertype = "matrix";
  for (std::map <std::string, std::string>::const_iterator p = parameters.begin(); p != parameters.end(); p++)
  {
    if (p->first == "stenciltype")
    {
      _stenciltype = p->second;
    }
    else if (p->first == "matrixtype")
    {
      _matrixtype = p->second;
    }
    else if (parameters.find("smoother") != parameters.end())
    {
      _smoother = p->second;
    }
    else if (parameters.find("smoothertype") != parameters.end())
    {
      _smoothertype = p->second;
    }
    else if (parameters.find("coarsesolver") != parameters.end())
    {
      _coarsesolver = p->second;
    }
    else if (parameters.find("transfertype") != parameters.end())
    {
      _transfertype = p->second;
    }
  }
  if(_matrixtype!="stencil")
  {
    assert(_smoothertype!="stencil");
  }
}
