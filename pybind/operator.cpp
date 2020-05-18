//
//  operator.c
//  Fada
//
//  Created by Roland Becker on 04/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <iostream>
#include  <pybind11/pybind11.h>
#include  <pybind11/numpy.h>
#include  <armadillo>
#include  "../Fada/operator.hpp"
#include  "carma.h"

namespace py = pybind11;


/*-------------------------------------------------*/
class OperatorPy : public Operator
{
public:
  OperatorPy(int nlevels, py::array_t<int> n, std::string matrixtype="Q1") : Operator()
  {
//    std::cerr << "cocou\n" << py::len(n);
    armaicvec n0(py::len(n));
    auto buf = n.request();
    int* ptr = (int *) buf.ptr;
    for (int i=0;i<py::len(n);i++)
    {
      n0[i] = ptr[i];
    }
    Operator::set_size(nlevels, n0, matrixtype);
//    std::cerr << _n;
  }
  py::array_t<int> get_dimensions() const
  {
    arma::Col<int> dims; dims.ones(3);
//    for(int i=0;i<Operator::dim();i++) dims[i] =  Operator::n()[i];
    for(int i=0;i<Operator::dim();i++) dims[i] =  Operator::nmax(i);
    return carma::col_to_arr<int>(dims);
  }

  py::array_t<double> get_solution()
  {
    armavec& av = Operator::get_solution().arma();
    const armaicvec& n = Operator::get_solution().n();
    return carma::col_to_arr<double>(av);
  }
};

PYBIND11_MODULE(pyfada, m) {
    m.doc() = "fada plugin";
    py::class_<OperatorPy>(m, "Operator")
    .def(py::init<int, py::array_t<int> >())
    .def(py::init<int, py::array_t<int>, std::string >())
    .def("testsolve", &OperatorPy::testsolve, py::arg("print")=true, py::arg("problem")="Random")
    .def("get_solution", &OperatorPy::get_solution)
    .def("get_dimensions", &OperatorPy::get_dimensions)
    .def("nall", &OperatorPy::nall)
    .def("dim", &OperatorPy::dim)
    .def_readwrite("maxiter", &OperatorPy::maxiter)
    .def_readwrite("smoother", &OperatorPy::smoother);
//    .def_readwrite("optmem", &OperatorPy::optmem);
}

