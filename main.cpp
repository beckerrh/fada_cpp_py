//
//  main.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <stdio.h>
#include  <stdlib.h>
#include  <ctime>
#include  "Fada/vecteur.hpp"
#include  "Fada/vecteurmg.hpp"
#include  "Fada/operateur.hpp"

/*-------------------------------------------------*/
inline double seconds(void)
{
    static const double secs_per_tick = 1.0 / CLOCKS_PER_SEC;
    return ( (double) clock() ) * secs_per_tick;
}

/*-------------------------------------------------*/
int main(int argc, char** argv)
{
  armaicvec n0;
  n0 << 3 << 3 << arma::endr;
  int smoother = 1, levels = 12;
  int         maxiter = 100;
  double      tol_rel = 1e-8, tol_abs = 1e-14;

  double t0 = seconds();
  /*
   smoother [0:jacobi 1:gauss_seidel 2:gauss_seidel_symm]
   */
  Operateur   A(smoother, levels, n0);

  Vecteur     u,f;

  u.set_size(A.n());
  f.set_size(u);

  A.right(f);
  A.boundary(u);
  A.boundary(f, u);

  int iter = A.solve(u, f, maxiter, tol_rel, tol_abs);

  std::string filename("solution.hdf");
  u.output(filename);
  
  printf("Vous avez utilise le lisseur ");
  if   (A.smoother()==0  ) printf("[Jac]");
  else if(A.smoother()==1) printf("[GS1]");
  else if(A.smoother()==2) printf("[GS2]");
  printf("\n\nNo. Iterations %3d (N = %6d)\n",iter, A.nall());
  printf("Total time: %6.2f\n", seconds()-t0);
}
