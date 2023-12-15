//
//  main.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <ctime>
#include  "Fada/solverstokes.hpp"

/*-------------------------------------------------*/
int main(int argc, char** argv)
{
    auto pmap = args2map(argc, argv);
    pmap["model"] = "coupled";
    auto solver = std::make_shared<SolverStokes>(pmap);  
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
