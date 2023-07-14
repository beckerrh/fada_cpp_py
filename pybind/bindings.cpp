//
//  bindings.cpp
//  Fada
//
//  Created by Roland Becker on 22/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <pybind11/stl.h>
#include  "uniformgridpy.hpp"
#include  "uniformmultigridpy.hpp"
#include  "solverlaplacepy.hpp"

// PYBIND11_MAKE_OPAQUE(std::map<std::string, std::string>);

/*=================================================*/
PYBIND11_MODULE(pyfada, m) {
    m.doc() = "fada plugin";
#
    pybind11::class_<UniformGridPy>(m, "UniformGrid")
      .def(pybind11::init<>())
      .def(pybind11::init<pybind11::array_t<int>& >(),pybind11::arg("n0"))
      .def("get_dimensions", &UniformGridPy::get_dimensions)
      .def("dim", &UniformGridPy::dim)
      .def("__repr__", &UniformGridPy::toString)
      .def("savehdf5", &UniformGridPy::savehdf5)
      .def("loadhdf5", &UniformGridPy::loadhdf5);
#
    pybind11::class_<UniformMultiGridPy>(m, "UniformMultiGrid")
      .def(pybind11::init<int, pybind11::array_t<int>& >(),pybind11::arg("nlevels"),pybind11::arg("n0"))
      .def("get_dimensions", &UniformMultiGridPy::get_dimensions)
    //      .def("n", &UniformMultiGridPy::n)
    //      .def("bounds", &UniformMultiGridPy::bounds)
    //      .def("dx", &UniformMultiGridPy::dx)
      .def("dim", &UniformMultiGridPy::dim)
      .def("n_fine", &UniformMultiGridPy::n_fine)
      .def("__repr__", &UniformMultiGridPy::toString);
#
    pybind11::class_<SolverLaplacePy>(m, "SolverLaplace")
      // py::bind_map<std::map<std::string, std::string>>(m, "StringStringMap");
      .def(pybind11::init<const UniformMultiGridPy&, const std::map<std::string, std::string>&>(), pybind11::arg("umg"), pybind11::arg("parameters"))
      .def("testsolve", &SolverLaplacePy::testsolve, pybind11::arg("print")=true, pybind11::arg("problem")="Random")
      .def("get_solution", &SolverLaplacePy::get_solution)
      .def("__repr__", &SolverLaplacePy::toString);
//      .def("set_solution", &MgSolverPy::set_solution)
//      .def("get_test", &MgSolverPy::get_test)
//      .def("set_test", &MgSolverPy::set_test)
//      .def("get_dimensions", &MgSolverPy::get_dimensions)
//      .def("n_fine", &MgSolverPy::n_fine)
//      .def("dim", &MgSolverPy::dim)
//      .def_readwrite("maxiter", &MgSolverPy::maxiter);
}
