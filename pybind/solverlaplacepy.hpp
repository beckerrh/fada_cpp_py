//
//  solverlaplace.h
//  Fada
//
//  Created by Roland Becker on 29/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef  solverlaplacepy_hpp
#define  solverlaplacepy_hpp

#include  <carma>
#include  "Fada/solverlaplace.hpp"
#include  "uniformmultigridpy.hpp"


/*-------------------------------------------------*/
class SolverLaplacePy : public SolverLaplace
{
protected:
public:
  SolverLaplacePy(const UniformMultiGridPy& umg, const std::map<std::string,std::string>& parameters)
  {
    auto mggrid = std::make_shared<UniformMultiGrid>(umg);
    set_data(mggrid, parameters);
  }

  pybind11::array_t<double> get_solution()
  {
      return carma::mat_to_arr<double>(SolverLaplace::get_solution());      
  }
};

#endif
