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
#include  "solverstokespy.hpp"

// PYBIND11_MAKE_OPAQUE(std::map<std::string, std::string>);

/*=================================================*/
PYBIND11_MODULE(pyfada, m) {
    m.doc() = "fada plugin";
#
    pybind11::class_<StokesInfo>(m, "StokesInfo")
        .def("__repr__", &StokesInfo::toString)
        .def_readwrite("niter_mean_v", &StokesInfo::niter_mean_v)
        .def_readwrite("niter_mean_p", &StokesInfo::niter_mean_p)
        .def_readwrite("err_v", &StokesInfo::err_v)
        .def_readwrite("err_p", &StokesInfo::err_p)
        .def_readwrite("niter", &StokesInfo::niter);

    pybind11::class_<LaplaceInfo>(m, "LaplaceInfo")
        .def("__repr__", &LaplaceInfo::toString)
        .def_readwrite("err", &LaplaceInfo::err)
        .def_readwrite("niter", &LaplaceInfo::niter);

    pybind11::class_<MgSolver::IterationInfo>(m, "MgSolver::IterationInfo")
        .def("__repr__", &MgSolver::IterationInfo::toString)
        .def_readwrite("maxiter", &MgSolver::IterationInfo::maxiter)
        .def_readwrite("tol_rel", &MgSolver::IterationInfo::tol_rel)
        .def_readwrite("tol_abs", &MgSolver::IterationInfo::tol_abs);

    // pybind11::class_<StokesInfoSde>(m, "StokesInfoSde")
    //     .def_readwrite("niter_mean_v", &StokesInfo::niter_mean_v)
    //     .def_readwrite("niter_mean_p", &StokesInfo::niter_mean_p)
    //     .def_readwrite("err_v", &StokesInfoSde::err_v)
    //     .def_readwrite("err_p", &StokesInfoSde::err_p)
    //     .def_readwrite("niter", &StokesInfoSde::niter);

    pybind11::class_<UniformGridPy>(m, "UniformGrid")
      .def(pybind11::init<>())
      .def(pybind11::init<pybind11::array_t<int>& >(),pybind11::arg("n0"))
      .def("get_dimensions", &UniformGridPy::get_dimensions)
      .def("dim", &UniformGridPy::dim)
      .def("__repr__", &UniformGridPy::toString)
      .def("savehdf5", &UniformGridPy::savehdf5)
      .def("n_gridpoints", &UniformGridPy::n_gridpoints)
      .def("n", &UniformGridPy::n)
      .def("dx", &UniformGridPy::dx)
      .def("loadhdf5", &UniformGridPy::loadhdf5,pybind11::arg("filename"));
#
    pybind11::class_<UniformMultiGridPy>(m, "UniformMultiGrid")
      .def(pybind11::init<int, pybind11::array_t<int>& >(),pybind11::arg("nlevels"),pybind11::arg("n0"))
      .def("get_dimensions", &UniformMultiGridPy::get_dimensions)
    //      .def("n", &UniformMultiGridPy::n)
    //      .def("bounds", &UniformMultiGridPy::bounds)
    //      .def("dx", &UniformMultiGridPy::dx)
      .def("dim", &UniformMultiGridPy::dim)
      .def("__repr__", &UniformMultiGridPy::toString);
#
    pybind11::class_<SolverLaplacePy>(m, "SolverLaplace")
      // py::bind_map<std::map<std::string, std::string>>(m, "StringStringMap");
      // .def(pybind11::init<const UniformMultiGridPy&, const std::map<std::string, std::string>&>(), pybind11::arg("umg"), pybind11::arg("parameters"))
      .def(pybind11::init<const std::map<std::string, std::string>&>(), pybind11::arg("parameters"))
      .def("testsolve", &SolverLaplacePy::testsolve, pybind11::arg("print")=true, pybind11::arg("info")=MgSolver::IterationInfo())
      .def("save_for_visu", &SolverLaplacePy::save_for_visu)
      .def("n_gridpoints", &SolverLaplacePy::n_gridpoints)
      .def("get_grid", &SolverLaplacePy::get_grid)
      .def("get_solution_nodes", &SolverLaplacePy::get_solution_nodes)
      .def("__repr__", &SolverLaplacePy::toString);
#
    pybind11::class_<SolverStokesPy>(m, "SolverStokes")
      .def(pybind11::init<const std::map<std::string, std::string>&>(), pybind11::arg("parameters"))
      .def("chorin_stationary", &SolverStokesPy::chorin_stationary, pybind11::arg("print")=true)
      .def("chorin_sde", &SolverStokesPy::chorin_sde, pybind11::arg("print")=true)
      .def("solve_stationary", &SolverStokesPy::solve_stationary, pybind11::arg("print")=true)
      .def("save_for_visu", &SolverStokesPy::save_for_visu)
      .def("n_gridpoints", &SolverStokesPy::n_gridpoints)
      .def("__repr__", &SolverStokesPy::toString);
//      .def("set_solution", &MgSolverPy::set_solution)
//      .def("get_test", &MgSolverPy::get_test)
//      .def("set_test", &MgSolverPy::set_test)
//      .def("get_dimensions", &MgSolverPy::get_dimensions)
//      .def("n_gridpoints", &MgSolverPy::n_gridpoints)
//      .def("dim", &MgSolverPy::dim)
//      .def_readwrite("maxiter", &MgSolverPy::maxiter);
}
