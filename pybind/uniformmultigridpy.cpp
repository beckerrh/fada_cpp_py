//
//  uniformmultigrid.cpp
//  Fada
//
//  Created by Roland Becker on 22/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "uniformmultigridpy.hpp"
#include  "carma.h"

/*-------------------------------------------------*/
UniformMultiGridPy::UniformMultiGridPy(int nlevelmax, int nlevels, pybind11::array_t<int> n) : UniformMultiGrid()
{
  armaicvec n0(py::len(n));
  auto buf = n.request();
  int* ptr = (int *) buf.ptr;
  for (int i=0;i<py::len(n);i++)
  {
    n0[i] = ptr[i];
  }
  UniformMultiGrid::set_size(nlevelmax, nlevels, n0);
}

/*-------------------------------------------------*/
pybind11::array_t<int> UniformMultiGridPy::get_dimensions() const
{
    arma::Col<int> dims; dims.ones(3);
    for(int i=0;i<UniformMultiGrid::dim();i++) dims[i] =  UniformMultiGrid::nmax(i);
    return carma::col_to_arr<int>(dims);
}

/*-------------------------------------------------*/
pybind11::array_t<int> UniformMultiGridPy::n()
{
//  std::cerr << "UniformMultiGridPy::n() " << UniformMultiGrid::n();
//  arma::Mat<int> nmat = UniformMultiGrid::n();
//  pybind11::array_t<int> pn = carma::mat_to_arr<int>(nmat);
//  std::cerr << "UniformMultiGridPy::n() " << py::len(pn);
  return carma::mat_to_arr<int>(UniformMultiGrid::n(), true);
}

/*-------------------------------------------------*/
pybind11::array_t<double> UniformMultiGridPy::dx()
{
  return carma::mat_to_arr<double>(UniformMultiGrid::dx(), true);
}

/*-------------------------------------------------*/
pybind11::array_t<double> UniformMultiGridPy::bounds() 
{
  return carma::mat_to_arr<double>(&UniformMultiGrid::bounds(), true);
}

