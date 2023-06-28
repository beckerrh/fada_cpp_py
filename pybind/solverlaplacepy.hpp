//
//  solverlaplace.h
//  Fada
//
//  Created by Roland Becker on 29/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef solverlaplacepy_hpp
#define solverlaplacepy_hpp

#include  "../Fada/solverlaplace.hpp"
#include  "uniformmultigridpy.hpp"


/*-------------------------------------------------*/
class SolverLaplacePy : public SolverLaplace
{
protected:
public:
  SolverLaplacePy(const UniformMultiGridPy& umg, const std::map<std::string,std::string>& parameters);

  pybind11::array_t<double> get_solution();
};

#endif /* solverlaplacepy_hpp */
