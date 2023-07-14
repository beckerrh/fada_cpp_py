//
//  main.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

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
  std::map<std::string,std::string> parameters;
  armaicvec n0;
  int nlevels=-1, dim=2;
  std::string problem("Random");
  if(argc%2!=1)
  {
    std::cerr << "arguments with '-d' or '-m' or '-st' or '-sm' '-smt' '-p' '-bc -tr'\n";
    std::cerr << "argc="<<argc<<" argv=" << *argv<<"\n";
    exit(1);
  }
  for(int i=1; i<argc;i++)
  {
    std::cerr << "i="<<i<<" argv="<<argv[i]<<"\n";
    if(!strcmp(argv[i], "-d"))
    {
      dim = atoi(argv[i+1]);
    }
    else if(!strcmp(argv[i], "-n"))
    {
      nlevels = atoi(argv[i+1]);
    }
    else if(!strcmp(argv[i], "-st"))
    {
      parameters["stenciltype"] = argv[i+1];
    }
    else if(!strcmp(argv[i], "-m"))
    {
      parameters["matrixtype"] = argv[i+1];
    }
    else if(!strcmp(argv[i], "-smt"))
    {
      parameters["smoothertype"] = argv[i+1];
    }
    else if(!strcmp(argv[i], "-sm"))
    {
      parameters["smoother"] = argv[i+1];
    }
    else if(!strcmp(argv[i], "-p"))
    {
      problem = argv[i+1];
    }
    else if(!strcmp(argv[i], "-bc"))
    {
      parameters["boundary_condition"] = argv[i+1];
    }
    else if(!strcmp(argv[i], "-tr"))
    {
      parameters["transfertype"] = argv[i+1];
    }
    else
    {
      std::cerr << "unknown argument" << argv[i] << "\n";
      exit(1);
    }
    i++;
  }
  if(dim==2)
  {
    n0 = {5,5};
    if(nlevels==-1) nlevels = 9;
  }
  else
  {
    n0 = {5,5,5};
    if(nlevels==-1) nlevels = 6;
  }
  double t0 = seconds();
  auto mggrid = std::make_shared<UniformMultiGrid>();
  mggrid->set_size(nlevels, n0);
  int updatemem = 0;
  auto solver = std::make_shared<SolverLaplace>(mggrid, parameters);
  int iter = solver->testsolve(true, problem);
  // const GridVector& u = solver->get_solution();
  // printf("u = %10.4e  %10.4e\n", arma::mean(u.data()), arma::max(u.data()));
  //
  std::string filename("solution.hdf");
  // u.output(filename);
  arma::hdf5_name spec(filename);
  solver->get_solution().save(spec);
  mggrid->get(0)->savehdf5("grid.hdf");

  printf("Vous avez utilise le lisseur %s\n", parameters["smoother"].c_str());
  printf("No. Iterations %3d (N = %6d dim = %2d)\n",iter, mggrid->get(0)->n_fine(), (int)mggrid->dim());
  printf("Total time: %6.2f\n", seconds()-t0);
}
