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
class OperatorPy : public Operator
{
public:
  OperatorPy(int nlevelmax, int nlevels, pybind11::array_t<int> n0, std::string femtype="Q1", std::string matrixtype="Full");
  OperatorPy(pybind11::array_t<int> n, pybind11::array_t<double> bounds, std::string femtype="Q1", std::string matrixtype="Full");
//  OperatorPy(UniformMultiGridPy umg, pybind11::array_t<int> n, std::string femtype="Q1", std::string matrixtype="Full");
  pybind11::array_t<int> get_dimensions() const;
  pybind11::array_t<double> get_solution();
};

#endif /* operatorpy_hpp */
