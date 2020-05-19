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
    
    arma::sp_mat A = arma::sprandu<arma::sp_mat>(1000, 1000, 0.1);
    
    arma::vec b = arma::randu<arma::vec>(1000);
    arma::mat B = arma::randu<arma::mat>(1000, 5);
    
    arma::vec x = arma::spsolve(A, b);  // solve one system
    arma::mat X = arma::spsolve(A, B);  // solve several systems
    
    bool status = arma::spsolve(x, A, b);  // use default solver
    if(status == false)  { std::cout << "no solution" << std::endl; }
    
    arma::spsolve(x, A, b, "lapack" );  // use LAPACK  solver
    arma::spsolve(x, A, b, "superlu");  // use SuperLU solver
    
    arma::superlu_opts opts;
    
    opts.allow_ugly  = true;
    opts.equilibrate = true;
    
    arma::spsolve(x, A, b, "superlu", opts);
}
