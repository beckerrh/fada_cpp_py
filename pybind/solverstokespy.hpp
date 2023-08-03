//
//  solverlaplace.h
//  Fada
//
//  Created by Roland Becker on 29/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef  solverstokespy_hpp
#define  solverstokespy_hpp

#include  <carma>
#include  <pybind11/pybind11.h>
#include  <pybind11/numpy.h>
#include  <pybind11/pytypes.h>
#include  "../Fada/Stokes/solverstokes.hpp"
#include  "uniformmultigridpy.hpp"


/*-------------------------------------------------*/
class SolverStokesPy : public SolverStokes
{
public:
  SolverStokesPy(const std::map<std::string,std::string>& parameters) : SolverStokes(parameters) {}

  py::tuple chorin_sde(bool print=true)
  {
      StokesInfoSde info = SolverStokes::chorin_sde(print);
      return py::make_tuple(
          carma::col_to_arr(info.dt),
          carma::col_to_arr(info.err_p),
          carma::col_to_arr(info.err_v)
            );     
  }

};

#endif
