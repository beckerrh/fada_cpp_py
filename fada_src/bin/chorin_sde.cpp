//
//  main.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "Fada/solverstokes.hpp"

/*-------------------------------------------------*/
int main(int argc, char** argv)
{
    auto pmap = args2map(argc, argv);
    pmap["dt"] = "0.0625";
    pmap["nt"] = "8";
    pmap["nt"] = "5";
    pmap["end_time"] = "10";
    pmap["niter_out"] = "5";
    pmap["n0"] = "13";
    pmap["nlevels"] = "4";
    pmap["nsamples"] = "2";
    pmap["q"] = "0.1";
    pmap["model"] = "decoupled";

    auto solver = std::make_shared<SolverStokes>(pmap);  
    StokesInfoSde info = solver->chorin_sde(false);
}
