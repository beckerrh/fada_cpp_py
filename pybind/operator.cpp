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

namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}

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
};

PYBIND11_MODULE(pyfada, m) {
    m.doc() = "fada plugin";
    py::class_<armaicvec>(m, "armaivec");
    py::class_<OperatorPy>(m, "Operateur")
    .def(py::init<int, py::array_t<int> >())
    .def("testsolve", &OperatorPy::testsolve)
    .def_readwrite("maxiter", &OperatorPy::maxiter)
    .def_readwrite("smoother", &OperatorPy::smoother);
//        .def("n", &Operateur::n);
    m.def("add", &add, "A function which adds two numbers");
}
