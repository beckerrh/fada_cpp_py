//
//  solverlaplace.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef solverlaplace_hpp
#define solverlaplace_hpp

#include  "finiteelementinterface.hpp"
#include  "multigridinterface.hpp"
#include  "mgsolver.hpp"
#include  "nodevector.hpp"

/*-------------------------------------------------*/
class SolverLaplace
{
protected:
  std::shared_ptr<FiniteElementInterface> _fem;
  std::shared_ptr<MultiGridInterface> _mggrid;
  MgSolver _mgsolver;
  NodeVector _u, _f;

public:
  SolverLaplace() : _fem(nullptr), _mggrid(nullptr) {}
  SolverLaplace(std::shared_ptr<MultiGridInterface> mggrid, std::string femtype, std::string matrixtype, std::string smoothertype);
  void set_data (std::shared_ptr<MultiGridInterface> mggrid, std::string femtype, std::string matrixtype, std::string smoothertype);
  
  std::string toString() const;

  int testsolve(bool print=true, std::string problem="DirichletRhsOne");
  const NodeVector& get_solution() const {return _u;}
  NodeVector& get_solution() {return _u;}
};


#endif /* solverlaplace_hpp */
