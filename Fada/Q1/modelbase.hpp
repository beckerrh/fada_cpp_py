//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef modelbase_hpp
#define modelbase_hpp

#include  <map>
#include  "boundary_conditions.hpp"

class UniformGrid;

/*-------------------------------------------------*/
class ModelBase
{
protected:
  std::shared_ptr <UniformGrid const> _ug;
  std::string _stenciltype, _matrixtype, _smoothertype, _smoother, _coarsesolver, _transfertype;
  std::shared_ptr <BoundaryConditions const> _boundaryconditions;

public:
  ~ModelBase()
  {
  }

  ModelBase(const std::map <std::string, std::string>& parameters, std::shared_ptr <BoundaryConditions const> boundaryconditions = nullptr);
  ModelBase(const ModelBase& model) : _ug(model._ug), _stenciltype(model._stenciltype), _boundaryconditions(model._boundaryconditions)
  {
  }
};

#endif
