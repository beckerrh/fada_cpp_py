//
//  operator.cpp
//  Fada
//
//  Created by Roland Becker on 04/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef operatorpy_hpp
#define operatorpy_hpp

#include  <pybind11/pybind11.h>
#include  <pybind11/numpy.h>
#include  "../Fada/operator.hpp"
#include  "uniformmultigridpy.hpp"

/*-------------------------------------------------*/
class MgSolverPy : public MgSolver
{
public:
  MgSolverPy(int nlevelmax, int nlevels, pybind11::array_t<int> n0, std::string femtype="Q1", std::string matrixtype="Full", std::string smoothertype="Jac");
  MgSolverPy(pybind11::array_t<int> n, pybind11::array_t<double> bounds, std::string femtype="Q1", std::string matrixtype="Full", std::string smoothertype="Jac");
  MgSolverPy(UniformMultiGridPy& umg);
  pybind11::array_t<int> get_dimensions() const;
  pybind11::array_t<double> get_solution();
  void set_solution(pybind11::array_t<double>);
  pybind11::array_t<double> get_test();
  void set_test(pybind11::array_t<double>&);
};

#endif /* operatorpy_hpp */
