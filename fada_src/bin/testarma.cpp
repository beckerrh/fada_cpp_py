//
//  testarma.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <stdio.h>
#include  <iostream>
#include  <armadillo>
#include  "Fada/gridvector.hpp"
#include  "Fada/umfmatrix.hpp"
#include  <amgcl/solver/lgmres.hpp>
#include  "Fada/uniformmultigrid.hpp"
#include  "Fada/solverlaplace.hpp"
#include  "Fada/analyticalfunctioninterface.hpp"


/*-------------------------------------------------*/
void solve_umf(std::shared_ptr <MatrixInterface const> A, GridVector& x, GridVector& b)
{
  UmfMatrix umf;

  // umf.getSparseMatrix() = A;
  umf.init(A);
  umf.computeLu();
  umf.solve(x, b);
  // arma::sp_mat A2(A.rows(), A.cols(), A.values(), A.nrows(), A.nrows());
  // armavec y(x);
  // y.fill(0);
  // A->dot(y,x);
  // y-=b;
  // std::cerr << "zero ? " << y.is_zero(1e-12) << std::endl;
}

/*-------------------------------------------------*/
int main(int argc, char **argv)
{
  std::map<std::string,std::string> parameters;
  parameters["n0"] = "3";
  parameters["matrixtype"] = "matrix";
  parameters["application"] = "Random";
  auto solver = std::make_shared <SolverLaplace>(parameters);
  auto u = std::make_shared <Vector <GridVector>>();
  auto f = std::make_shared <Vector <GridVector>>();
  auto mggrid = solver->get_mgrid();
  u->set_size(mggrid->get(0)->n());
  u->fill(0);
  f->set_size(mggrid->get(0)->n());
  auto p = std::make_shared<RandomFunction>();
  solver->getModel()->rhs(f, mggrid->get(0), solver->getApplication());

  auto M = solver->getModel()->newMatrix(mggrid->get(0));
  // arma::umat locations;
  // armavec    values;
  // M->get_locations_values(locations, values);
  // typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;
  // auto p = std::make_shared <SparseMatrixDef>(locations, values);
  solve_umf(M, u->get(), f->get());


  // armaicvec n;
  //
  // n = { 4, 2, 3 };
  // // n << 4 << 2 << 3 << arma::endr;
  // GridVector v(n), w(n);
  // v.fill(2);
  // w.data().randu();
  //
  // GridVector u(v);
  //
  // u = 3 * v.data() + w.data();
  // std::cerr << u << std::endl;
  // std::cerr << "u(2,1,2) =" << u.at(2, 1, 2) << std::endl;
  // std::cerr << arma::dot(v.data(), w.data()) << std::endl;
  //
  //
  // int       size = 3;
  // arma::vec b    = arma::randu <arma::vec>(size);
  // arma::vec x(size), y(size);
  //
  // UmfMatrix     umf;
  // SparseMatrix& sp    = umf.getSparseMatrix();
  // int           nelem = 4;
  // arma::umat    locations(2, nelem);
  // armavec       values(nelem);
  // locations(0, 0) = 0; locations(1, 0) = 0; values[0] = 1.1;
  // locations(0, 1) = 1; locations(1, 1) = 1; values[1] = 2.2;
  // locations(0, 2) = 2; locations(1, 2) = 2; values[2] = 3.3;
  // locations(0, 3) = 2; locations(1, 3) = 0; values[3] = 1.3;
  // sp.set_elements(locations, values);
  // sp.save(std::cerr);
  // umf.computeLu();
  // umf.solve(x, b);
  // std::cerr << "b=" << b.t();
  // std::cerr << "x=" << x.t();
  // y.fill(0);
  // sp.dot(y, x);
  // std::cerr << "y=" << y.t();
  //
  // arma::sp_mat A(locations, values);
  // arma::vec    x2 = arma::spsolve(A, b); // solve one system
  //
  // bool status = arma::spsolve(x2, A, b); // use default solver
  // if (status == false)
  // {
  //    std::cout << "no solution" << std::endl;
  // }
  //
  // arma::spsolve(x2, A, b, "lapack");   // use LAPACK  solver
  // std::cerr << "b=" << b.t();
  // std::cerr << "x=" << x.t();
  // y = A * x;
  // std::cerr << "y=" << y.t();
  //
  // arma::spsolve(x2, A, b, "superlu");  // use SuperLU solver
  // std::cerr << "b=" << b.t();
  // std::cerr << "x=" << x.t();
  // y = A * x;
  // std::cerr << "y=" << y.t();
  //
  // arma::superlu_opts opts;
  // opts.allow_ugly  = true;
  // opts.equilibrate = true;
  // arma::spsolve(x2, A, b, "superlu", opts);
  // std::cerr << "b=" << b.t();
  // std::cerr << "x=" << x.t();
  // y = A * x;
  // std::cerr << "y=" << y.t();
}
