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
#include  "Fada/solverlaplace.hpp"
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
  int nlevelmax=12, dim=2;
  if(dim==2)
  {
    n0 << 3 << 3 << arma::endr;
    nlevelmax = 6;
  }
  else
  {
    n0 << 3 << 3 << 3 << arma::endr;
//    nlevels = 8;
//    n0 << 2 << 2 << 2 << arma::endr;
    nlevelmax = 6;
  }
  int nlevels=nlevelmax;
  double t0 = seconds();
  auto mggrid = std::make_shared<UniformMultiGrid>();
  mggrid->set_size(nlevelmax, nlevels, n0);
  std::string smoother = "GS2";
  auto solver = std::make_shared<SolverLaplace>(mggrid, "Q1", "Trapez", smoother);
  int iter = solver->testsolve();
  const NodeVector& u = solver->get_solution();
  printf("u = %10.4e  %10.4e\n", arma::mean(u.data()), arma::max(u.data()));

  std::string filename("solution.hdf");
  u.output(filename);
  
  printf("Vous avez utilise le lisseur %s\n", smoother.c_str());
  printf("No. Iterations %3d (N = %6d dim = %2d)\n",iter, mggrid->get(0)->nall(), (int)mggrid->dim());
  printf("Total time: %6.2f\n", seconds()-t0);
}
