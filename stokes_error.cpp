//
//  main.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <ctime>
#include  "Fada/Stokes/solverstokes.hpp"

/*-------------------------------------------------*/
inline double seconds(void)
{
    static const double secs_per_tick = 1.0 / CLOCKS_PER_SEC;
    return ( (double) clock() ) * secs_per_tick;
}

/*-------------------------------------------------*/
int main(int argc, char** argv)
{
    // double dt=0.1;
    std::map<std::string,std::string> parameters;
    parameters["dt"] = "0.01";
    if(argc%2!=1)
    {
        std::cerr << "arguments with '-n' or '-sm' '-p'\n";
        std::cerr << "argc="<<argc<<" argv=" << *argv<<"\n";
        exit(1);
    }
    for(int i=1; i<argc;i++)
    {
        std::cerr << "i="<<i<<" argv="<<argv[i]<<"\n";
        if(!strcmp(argv[i], "-n"))
        {
            parameters["nlevels"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-sm"))
        {
            parameters["smoother"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-n0"))
        {
            parameters["n0"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-p"))
        {
            parameters["application"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-dt"))
        {
            parameters["dt"] = std::atof(argv[i+1]);
        }
        else
        {
            std::cerr << "unknown argument" << argv[i] << "\n";
            exit(1);
        }
        i++;
    }
    // parameters["dt"] = std::to_string(dt);
    double t0 = seconds();

    auto solver = std::make_shared<SolverStokes>(parameters);  
    // auto info = solver->chorin_stationary();
    auto info = solver->solve_stationary();
    
    
    printf("niter: %6d mg(v): %3.1f mg(p): %3.1f errp: %9.3e errv: %9.3e\n", info.niter, info.niter_mean_v, info.niter_mean_p, info.err_p, info.err_v);
    solver->save_for_visu();
  
    // int iter = solver->testsolve(true);
    // armavec err = solver->compute_error();
    // solver->save_for_visu();
    //
    // printf("Number of grid points: %10d niter: %3d\n", solver->get_mgrid()->get(0)->n_gridpoints(), iter);
    // printf("Total time: %6.2f\n", seconds()-t0);
    // std::cerr << "Error: " << err.t()<<"\n";
}
