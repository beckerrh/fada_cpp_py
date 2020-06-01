//
//  uniformmultigrid.hpp
//  Fada
//
//  Created by Roland Becker on 22/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef uniformmultigridpy_hpp
#define uniformmultigridpy_hpp

#include  "../Fada/uniformmultigrid.hpp"
#include  <pybind11/pybind11.h>
#include  <pybind11/numpy.h>

/*-------------------------------------------------*/
class UniformMultiGridPy : public UniformMultiGrid
{
public:
  UniformMultiGridPy(int nlevelmax, int nlevels, pybind11::array_t<int>& n);
  pybind11::array_t<int> get_dimensions() const;
//  pybind11::array_t<int> n();
//  pybind11::array_t<double> dx();
//  pybind11::array_t<double> bounds();
};

#endif /* uniformmultigrid_hpp */
