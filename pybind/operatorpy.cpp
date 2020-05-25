//
//  operator.cpp
//  Fada
//
//  Created by Roland Becker on 04/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "operatorpy.hpp"
#include  "carma.h"

/*-------------------------------------------------*/
OperatorPy::OperatorPy(int nlevelmax, int nlevels, py::array_t<int> n0, std::string femtype, std::string matrixtype) : Operator(false)
{
  //    std::cerr << "cocou nlevelmax " << nlevelmax << " nlevels " << nlevels << "\n";
//  armaicvec n0av(py::len(n));
//  auto buf = n0.request();
//  int* ptr = (int *) buf.ptr;
//  for (int i=0;i<py::len(n0);i++)
//  {
//    n0av[i] = ptr[i];
//  }
  armaicvec n0av(carma::arr_to_col<int>(n0));
  Operator::set_size(nlevelmax, nlevels, n0av, femtype, matrixtype);
}
/*-------------------------------------------------*/
OperatorPy::OperatorPy(pybind11::array_t<int> n, pybind11::array_t<double> bounds, std::string femtype, std::string matrixtype)
{
  auto p = std::make_shared<armamat>(std::move(carma::arr_to_mat<double>(bounds)));
  UniformMultiGrid umg(carma::arr_to_mat<int>(n), p);
  Operator::set_size(umg, femtype, matrixtype);
}

/*-------------------------------------------------*/
py::array_t<int> OperatorPy::get_dimensions() const
{
  arma::Col<int> dims; dims.ones(3);
  for(int i=0;i<Operator::dim();i++) dims[i] =  Operator::nmax(i);
  return carma::col_to_arr<int>(dims);
}

/*-------------------------------------------------*/
py::array_t<double> OperatorPy::get_solution()
{
  armavec& av = Operator::get_solution().arma();
  const armaicvec& n = Operator::get_solution().n();
  return carma::col_to_arr<double>(av);
}
