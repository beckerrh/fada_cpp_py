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
#include  "Fada/vector.hpp"
#include  "Fada/umfmatrix.hpp"


/*-------------------------------------------------*/
int main(int argc, char** argv)
{
    armaicvec n;
    n << 4 << 2 << 3 << arma::endr;
    Vector v(n), w(n);
    v.fill(2);
    w.randu();
    
    Vector u(v);
    
    u = 3*v.arma() + w.arma();
    std::cerr << u << std::endl;
    std::cerr << "u(2,1,2) =" << u.at(2,1,2) << std::endl;
    std::cerr << arma::dot(v.arma(),w.arma()) << std::endl;

  
  int size = 3;
  arma::vec b = arma::randu<arma::vec>(size);
  arma::vec x(size), y(size);
  
  UmfMatrix umf;
  SparseMatrix& sp = umf.getSparseMatrix();
  int nelem=4;
  arma::umat locations(2, nelem);
  armavec values(nelem);
  locations(0, 0) = 0; locations(1, 0) = 0; values[0] = 1.1;
  locations(0, 1) = 1; locations(1, 1) = 1; values[1] = 2.2;
  locations(0, 2) = 2; locations(1, 2) = 2; values[2] = 3.3;
  locations(0, 3) = 2; locations(1, 3) = 0; values[3] = 1.3;
  sp.set_elements(locations, values);
  sp.save(std::cerr);
  umf.computeLu();
  umf.solve(x, b);
  std::cerr << "b="<<b.t();
  std::cerr << "x="<<x.t();
  y.fill(0);
  sp.dot(y,x);
  std::cerr << "y="<<y.t();

//  umf.computeLu(A);
//  umf.solve(x, b);
//
//    arma::vec x = arma::spsolve(A, b);  // solve one system
//
//    bool status = arma::spsolve(x, A, b);  // use default solver
//    if(status == false)  { std::cout << "no solution" << std::endl; }
//
//    arma::spsolve(x, A, b, "lapack" );  // use LAPACK  solver
//    arma::spsolve(x, A, b, "superlu");  // use SuperLU solver
//
//    arma::superlu_opts opts;
//
//    opts.allow_ugly  = true;
//    opts.equilibrate = true;
//
//    arma::spsolve(x, A, b, "superlu", opts);
}
