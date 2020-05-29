//
//  operator.cpp
//  Fada
//
//  Created by Roland Becker on 04/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "operatorpy.hpp"
#include  "carma/carma.h"
//#include  "arma2py.hpp"

/*-------------------------------------------------*/
MgSolverPy::MgSolverPy(int nlevelmax, int nlevels, pybind11::array_t<int> n0, std::string femtype, std::string matrixtype, std::string smoothertype) : MgSolver(false)
{
  //    std::cerr << "cocou nlevelmax " << nlevelmax << " nlevels " << nlevels << "\n";
//  armaicvec n0av(py::len(n));
//  auto buf = n0.request();
//  int* ptr = (int *) buf.ptr;
//  for (int i=0;i<py::len(n0);i++)
//  {
//    n0av[i] = ptr[i];
//  }
//  armaicvec n0av(carma::arr_to_col<int>(n0));
//  armaicvec n0av(arr2col<int>(n0));
  MgSolver::set_size(nlevelmax, nlevels, carma::arr_to_col<int>(n0, true), femtype, matrixtype, smoothertype);
}
/*-------------------------------------------------*/
MgSolverPy::MgSolverPy(pybind11::array_t<int> n, pybind11::array_t<double> bounds, std::string femtype, std::string matrixtype, std::string smoothertype)
{
//  auto p = std::make_shared<armamat>(std::move(carma::arr_to_mat<double>(bounds)));
//  auto p = std::make_shared<armamat>(std::move(arr2mat<double>(bounds)));
//  UniformMultiGrid umg(carma::arr_to_mat<int>(n), p);
//  UniformMultiGrid umg(arr2mat<int>(n), p);
//  MgSolver::set_size(arr2mat<int>(n), femtype, matrixtype);
  MgSolver::set_size(carma::arr_to_mat<int>(n, true), femtype, matrixtype, smoothertype);
  std::cerr << " MgSolverPy::MgSolverPy() n = "<< _mggrid << "\n";
}
/*-------------------------------------------------*/
MgSolverPy::MgSolverPy(UniformMultiGridPy& umg)
{
  MgSolver::set_size(carma::arr_to_mat<int>(umg.n()));
}

/*-------------------------------------------------*/
pybind11::array_t<int> MgSolverPy::get_dimensions() const
{
  arma::Col<int> dims; dims.ones(3);
  for(int i=0;i<MgSolver::dim();i++) dims[i] =  MgSolver::nmax(i);
  return carma::col_to_arr<int>(dims, true);
//  return col2arr<int>(dims);
}

/*-------------------------------------------------*/
pybind11::array_t<double> MgSolverPy::get_solution()
{
  std::cerr << "MgSolverPy::get_solution() _u " << _u.data().memptr() << "\n";
//  return carma::col_to_arr<double>(_test, false);
  return carma::mat_to_arr<double>(_u.data());
}

/*-------------------------------------------------*/
void MgSolverPy::set_solution(pybind11::array_t<double> v)
{
  _u.data() =  carma::arr_to_mat<double>(v);
  std::cerr << "MgSolverPy::set_solution() _u " << _u.data().memptr() << "\n";
}

/*-------------------------------------------------*/
pybind11::array_t<double> MgSolverPy::get_test()
{
  std::cerr << "MgSolverPy::get_solution() _test " << _test.memptr() << "\n";
//  return carma::col_to_arr<double>(_test, false);
  return carma::mat_to_arr<double>(_test);
  pybind11::array_t<double> pytest = carma::mat_to_arr<double>(_test);
  carma::update_array(_test, pytest);
  return pytest;
}

/*-------------------------------------------------*/
void MgSolverPy::set_test(pybind11::array_t<double>& v)
{
  _test =  carma::arr_to_mat<double>(v);
  carma::update_array(_test, v);
  std::cerr << "MgSolverPy::set_solution() _test " << _test.memptr() << "\n";
}
