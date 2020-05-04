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
#include  "Fada/vector.hpp"
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
  int levels=12, dim=3;
  if(dim==2)
  {
    n0 << 3 << 3 << arma::endr;
    levels = 6;
  }
  else
  {
    n0 << 3 << 3 << 3 << arma::endr;
    levels = 6;
  }

  double t0 = seconds();
  /*
   smoother [0:jacobi 1:gauss_seidel 2:gauss_seidel_symm]
   */
  Operateur   A(levels, n0);
  A.smoother = "jac";

//  Vecteur     u,f;
  vector u,f;

//  std::cerr << "A.n()="<<A.n();
  
  u.set_size(A.n());
//  std::cerr << "u="<<u;
  f.set_size(u);
//  std::cerr << "f="<<f;
//  exit(1);


  A.right(f);
  A.boundary(u);
  A.boundary(f, u);

//  std::cerr << "f="<<f;

  int iter = A.solve(u, f);

  std::string filename("solution.hdf");
  u.output(filename);
  
  printf("Vous avez utilise le lisseur %s", A.smoother.c_str());
  printf("\n\nNo. Iterations %3d (N = %6d)\n",iter, A.nall());
  printf("Total time: %6.2f\n", seconds()-t0);
}
