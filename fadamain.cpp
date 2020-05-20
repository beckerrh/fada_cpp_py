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
#include  "Fada/operator.hpp"
#include  "Fada/uniformmultigrid.hpp"

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
  int nlevels=12, dim=2;
  if(dim==2)
  {
    n0 << 3 << 3 << arma::endr;
//    nlevels = 3;
  }
  else
  {
    n0 << 3 << 3 << 3 << arma::endr;
    nlevels = 8;
  }

  double t0 = seconds();
//  UniformMultiGrid mggrid;
//  mggrid.set_size(nlevels, n0);

  
//  Operateur   A(nlevels, n0);
  Operator   A;
  A.set_size(nlevels, n0);
//  A.smoother = "gs2";

  int iter = A.testsolve();
  Vector& u = A.get_solution();
  printf("u = %10.4e  %10.4e\n", arma::mean(u.arma()), arma::max(u.arma()));

//  std::string filename("solution.hdf");
//  u.output(filename);
  
  printf("Vous avez utilise le lisseur %s\n", A.smoother.c_str());
  printf("No. Iterations %3d (N = %6d)\n",iter, A.nall());
  printf("Total time: %6.2f\n", seconds()-t0);
}
