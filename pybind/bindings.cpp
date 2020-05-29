//
//  bindings.cpp
//  Fada
//
//  Created by Roland Becker on 22/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "uniformmultigridpy.hpp"
#include  "solverlaplacepy.hpp"

/*=================================================*/
PYBIND11_MODULE(pyfada, m) {
    m.doc() = "fada plugin";
#
    pybind11::class_<UniformMultiGridPy>(m, "UniformMultiGrid")
      .def(pybind11::init<int, int, pybind11::array_t<int>& >())
      .def("get_dimensions", &UniformMultiGridPy::get_dimensions)
//      .def("n", &UniformMultiGridPy::n)
//      .def("bounds", &UniformMultiGridPy::bounds)
//      .def("dx", &UniformMultiGridPy::dx)
      .def("dim", &UniformMultiGridPy::dim)
      .def("nall", &UniformMultiGridPy::nall)
      .def("__repr__", &UniformMultiGridPy::toString);
#
    pybind11::class_<SolverLaplacePy>(m, "SolverLaplace")
      .def(pybind11::init<const UniformMultiGridPy&, std::string, std::string, std::string>())
      .def("testsolve", &SolverLaplacePy::testsolve, pybind11::arg("print")=true, pybind11::arg("problem")="Random")
      .def("get_solution", &SolverLaplacePy::get_solution)
      .def("__repr__", &SolverLaplacePy::toString);
//      .def("set_solution", &MgSolverPy::set_solution)
//      .def("get_test", &MgSolverPy::get_test)
//      .def("set_test", &MgSolverPy::set_test)
//      .def("get_dimensions", &MgSolverPy::get_dimensions)
//      .def("nall", &MgSolverPy::nall)
//      .def("dim", &MgSolverPy::dim)
//      .def_readwrite("maxiter", &MgSolverPy::maxiter);
}