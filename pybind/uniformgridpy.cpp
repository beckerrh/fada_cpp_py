//
//  uniformmultigrid.cpp
//  Fada
//
//  Created by Roland Becker on 22/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "uniformgridpy.hpp"
#include  <carma>

/*-------------------------------------------------*/
UniformGridPy::UniformGridPy(pybind11::array_t<int>& n) : UniformGrid()
{
  armaicvec n0(pybind11::len(n));
  auto buf = n.request();
  int* ptr = (int *) buf.ptr;
  for (int i=0;i<pybind11::len(n);i++)
  {
    n0[i] =  ptr[i];
  }
//  std::cerr << "UniformGridPy::UniformGridPy() " << nlevelmax << " " << nlevels << " " << n0<<"\n";
  UniformGrid::set_size(n0);
//  std::cerr << "UniformGridPy::UniformGridPy() n = " << UniformMultiGrid::n()<<"\n";
}
//
/*-------------------------------------------------*/
pybind11::array_t<int> UniformGridPy::get_dimensions() const
{
    arma::Col<int> dims; dims.ones(3);
    // std::cerr<< "UniformGrid::n()=" << UniformGrid::n() << "\n";
    // std::cerr<< "UniformGrid::dim()=" << UniformGrid::dim() << "\n";
    for(int i=0;i<UniformGrid::dim();i++) dims[i] =  UniformGrid::n(i);
    // std::cerr<< "dims" << dims << "\n";
    return carma::col_to_arr<int>(dims, true);
}
//
///*-------------------------------------------------*/
//pybind11::array_t<int> UniformGridPy::n()
//{
//  return carma::mat_to_arr<int>(UniformMultiGrid::n(), true);
//}
//
///*-------------------------------------------------*/
//pybind11::array_t<double> UniformGridPy::dx()
//{
//  return carma::mat_to_arr<double>(UniformMultiGrid::dx(), true);
//}
//
///*-------------------------------------------------*/
//pybind11::array_t<double> UniformGridPy::bounds()
//{
//  return carma::mat_to_arr<double>(&UniformMultiGrid::bounds(), true);
//}
