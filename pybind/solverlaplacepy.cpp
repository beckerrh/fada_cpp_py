//
//  solverlaplace.c
//  Fada
//
//  Created by Roland Becker on 29/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "solverlaplacepy.hpp"
#include  "carma/carma.h"

/*-------------------------------------------------*/
SolverLaplacePy::SolverLaplacePy(const UniformMultiGridPy& umg, std::string femtype, std::string matrixtype, std::string smoothertype)
{
  auto mggrid = std::make_shared<UniformMultiGrid>(umg);
//  mggrid->set_size(nlevelmax, nlevels, n0);
  set_data(mggrid, femtype, matrixtype, smoothertype);
}

/*-------------------------------------------------*/
pybind11::array_t<double> SolverLaplacePy::get_solution()
{
  return carma::mat_to_arr<double>(_u.data());
}
