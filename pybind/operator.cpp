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
#include  "../Fada/operateur.hpp"
#include  "carma.h"

namespace py = pybind11;


/*-------------------------------------------------*/
class OperatorPy : public Operateur
{
public:
  OperatorPy(int nlevels, py::array_t<int> n) : Operateur()
  {
    std::cerr << "cocou\n" << py::len(n);
    armaicvec n0(py::len(n));
    auto buf = n.request();
    int* ptr = (int *) buf.ptr;
    for (int i=0;i<py::len(n);i++)
    {
      n0[i] = ptr[i];
    }
    Operateur::set_size(nlevels, n0);
    std::cerr << _n;
  }
  py::array_t<double> get_solution()
  {
    armavec& av = Operateur::get_solution().arma();
    const armaicvec& n = Operateur::get_solution().n();
    return carma::col_to_arr<double>(av);
  }
};

PYBIND11_MODULE(pyfada, m) {
    m.doc() = "fada plugin";
    py::class_<OperatorPy>(m, "Operateur")
    .def(py::init<int, py::array_t<int> >())
    .def("testsolve", &OperatorPy::testsolve)
    .def("get_solution", &OperatorPy::get_solution)
    .def_readwrite("maxiter", &OperatorPy::maxiter)
    .def_readwrite("smoother", &OperatorPy::smoother);
}
