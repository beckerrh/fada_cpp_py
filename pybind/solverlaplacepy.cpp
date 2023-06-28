//
//  solverlaplace.c
//  Fada
//
//  Created by Roland Becker on 29/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "solverlaplacepy.hpp"
#include  <carma>

/*-------------------------------------------------*/
SolverLaplacePy::SolverLaplacePy(const UniformMultiGridPy& umg, const std::map<std::string,std::string>& parameters)
{
  auto mggrid = std::make_shared<UniformMultiGrid>(umg);
  set_data(mggrid, parameters);
}

/*-------------------------------------------------*/
pybind11::array_t<double> SolverLaplacePy::get_solution()
{
  return carma::mat_to_arr<double>(SolverLaplace::get_solution());
}
