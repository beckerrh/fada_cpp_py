//
//  bindings.cpp
//  Fada
//
//  Created by Roland Becker on 22/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "operatorpy.hpp"
#include  "uniformmultigridpy.hpp"

/*=================================================*/
PYBIND11_MODULE(pyfada, m) {
    m.doc() = "fada plugin";
#
    pybind11::class_<UniformMultiGridPy>(m, "UniformMultiGrid")
      .def(pybind11::init<int, int, pybind11::array_t<int> >())
      .def("__repr__", &UniformMultiGridPy::toString)
      .def("get_dimensions", &UniformMultiGridPy::get_dimensions)
      .def("n", &UniformMultiGridPy::n)
      .def("bounds", &UniformMultiGridPy::bounds)
      .def("dx", &UniformMultiGridPy::dx)
      .def("dim", &UniformMultiGridPy::dim);
#
    pybind11::class_<OperatorPy>(m, "Operator")
//      .def(pybind11::init<UniformMultiGridPy>())
      .def(pybind11::init<int, int, pybind11::array_t<int> >())
      .def(pybind11::init<int, int, pybind11::array_t<int>, std::string, std::string >())
      .def(pybind11::init<pybind11::array_t<int>, pybind11::array_t<double>, std::string, std::string >())
      .def(pybind11::init<pybind11::array_t<int>, pybind11::array_t<double> >())
      .def("testsolve", &OperatorPy::testsolve, pybind11::arg("print")=true, pybind11::arg("problem")="Random")
      .def("get_solution", &OperatorPy::get_solution)
      .def("get_dimensions", &OperatorPy::get_dimensions)
      .def("nall", &OperatorPy::nall)
      .def("dim", &OperatorPy::dim)
      .def_readwrite("maxiter", &OperatorPy::maxiter)
      .def_readwrite("smoother", &OperatorPy::smoother);
}
