//
//  main.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <ctime>
#include  "Fada/solverlaplace.hpp"
#include  "Fada/timer.hpp"
#include  "Fada/uniformmultigrid.hpp"
//
// /*-------------------------------------------------*/
// inline double seconds(void)
// {
//     static const double secs_per_tick = 1.0 / CLOCKS_PER_SEC;
//     return ( (double) clock() ) * secs_per_tick;
// }

/*-------------------------------------------------*/
int main(int argc, char** argv)
{
    Timer T;
    T.enrol("all", false);
    T.start("all");
    std::map<std::string,std::string> parameters;
    int nlevels=-1, dim=2, n0=-1;
    parameters["application"] = "Random_dir";
    if(argc%2!=1)
    {
        std::cerr << "arguments with '-d' or '-m' or '-st' or '-sm' '-smt' '-p' '-bc -tr'\n";
        std::cerr << "argc="<<argc<<" argv=" << *argv<<"\n";
        exit(1);
    }
    for(int i=1; i<argc;i++)
    {
        std::cerr << "i="<<i<<" argv="<<argv[i]<<" "<<argv[i+1]<<"\n";
        if(!strcmp(argv[i], "-d"))
        {
            dim = atoi(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-n"))
        {
            nlevels = atoi(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-n0"))
        {
            n0 = atoi(argv[i+1]);
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
        else if(!strcmp(argv[i], "-app"))
        {
            parameters["application"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-tr"))
        {
            parameters["transfertype"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-method"))
        {
            parameters["method"] = argv[i+1];
        }
        else
        {
            std::cerr << "unknown argument" << argv[i] << "\n";
            exit(1);
        }
        i++;
    }
    if(n0==-1) n0 = 5;
    if(dim==2)
    {
        if(nlevels==-1) nlevels = 9;
    }
    else
    {
        if(nlevels==-1) nlevels = 6;
    }
    parameters["nlevels"] = std::to_string(nlevels);
    parameters["n0"] = std::to_string(n0);
    // std::cerr << "n0="<<n0<<"\n";
    // auto mggrid = std::make_shared<UniformMultiGrid>();
    // mggrid->set_size(nlevels, n0);
    // int updatemem = 0;
    auto solver = std::make_shared<SolverLaplace>(parameters);
    LaplaceInfo info = solver->testsolve(true, MgSolver::IterationInfo(30,1e-10));
    solver->save_for_visu();

    T.stop("all");
    printf("No. Iterations %3d (N = %6d dim = %2d)\n",info.niter,solver->n_gridpoints(), (int)solver->get_mgrid()->dim());
    printf("Total time: %6.2f\n", T.get("all"));
    if(info.has_error) printf("Error: %12.3e\n", info.err);
}
