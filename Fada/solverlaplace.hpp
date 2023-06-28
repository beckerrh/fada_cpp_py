//
//  solverlaplace.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef solverlaplace_hpp
#define solverlaplace_hpp

#include  "modelinterface.hpp"
#include  "multigridinterface.hpp"
#include  "mgsolver.hpp"
#include  "Q1/nodevector.hpp"

/*-------------------------------------------------*/
class SolverLaplace
{
protected:
  std::shared_ptr <ModelInterface> _model;
  std::shared_ptr <MultiGridInterface> _mggrid;
  MgSolver _mgsolver;

public:
  SolverLaplace() : _model(nullptr), _mggrid(nullptr) {}
  // SolverLaplace(std::shared_ptr <MultiGridInterface> mggrid, std::string stenciltype, std::string matrixtype, std::string smoothertype = "", int updatelength = 0);
  // void set_data(std::shared_ptr <MultiGridInterface> mggrid, std::string stenciltype, std::string matrixtype, std::string smoothertype = "", int updatelength = 0);
  SolverLaplace(std::shared_ptr <MultiGridInterface> mggrid, const std::map<std::string,std::string>& parameters);
  void set_data(std::shared_ptr <MultiGridInterface> mggrid, const std::map<std::string,std::string>& parameters);
  std::string toString() const;

  std::shared_ptr <const ModelInterface> getModel() const
  {
    return(_model);
  }

  int testsolve(bool print = true, std::string problem = "DirichletRhsOne");

  const NodeVector& get_solution() const
  {
    std::shared_ptr <const NodeVector> p = std::dynamic_pointer_cast <const NodeVector>(_mgsolver.getU());
    assert(p);
    return(*p);
  }

  NodeVector& get_solution()
  {
    std::shared_ptr <NodeVector> p = std::dynamic_pointer_cast <NodeVector>(_mgsolver.getU());
    assert(p);
    return(*p);
  }
  std::shared_ptr<VectorInterface> get_u()
  {
    return(_mgsolver.getU());
  }

  std::shared_ptr<VectorInterface> get_rhs()
  {
    return(_mgsolver.getF());
  }
};


#endif /* solverlaplace_hpp */
