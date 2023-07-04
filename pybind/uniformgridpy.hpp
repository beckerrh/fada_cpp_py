//
//  uniformmultigrid.hpp
//  Fada
//
//  Created by Roland Becker on 22/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef uniformgridpy_hpp
#define uniformgridpy_hpp

#include  "../Fada/uniformgrid.hpp"
#include  <pybind11/pybind11.h>
#include  <pybind11/numpy.h>

/*-------------------------------------------------*/
class UniformGridPy : public UniformGrid
{
public:
  UniformGridPy() : UniformGrid() {}
  UniformGridPy(pybind11::array_t<int>& n);
  pybind11::array_t<int> get_dimensions() const;
//  pybind11::array_t<int> n();
//  pybind11::array_t<double> dx();
//  pybind11::array_t<double> bounds();
};

#endif
