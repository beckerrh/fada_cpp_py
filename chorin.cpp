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
    // params = {'dt':1/16, 'nt':8, 'end_time':1, 'niter_out':2, 'n0':13, 'nlevels':4, 'mgtimer':'false', 'nsamples':10, 'q':0.1}

    std::map<std::string,std::string> parameters;
    // default values
    parameters["dt"] = "0.0625"; // 1/16
    parameters["nt"] = "8";
    parameters["end_time"] = "1";
    parameters["niter_out"] = "2";
    parameters["n0"] = "13";
    parameters["nlevels"] = "4";
    parameters["nsamples"] = "1000";
    parameters["q"] = "0.1";
    parameters["mgtimer"] = "false";
    if(argc%2!=1)
    {
        std::cerr << "arguments with '-n' or '-sm' '-p'\n";
        std::cerr << "argc="<<argc<<" argv=" << *argv<<"\n";
        exit(1);
    }
    for(int i=1; i<argc;i++)
    {
        std::cerr << "i="<<i<<" argv="<<argv[i]<<"\n";
        if(!strcmp(argv[i], "-nt"))
        {
            parameters["nt"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-dt"))
        {
            parameters["dt"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-nt"))
        {
            parameters["nt"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-end_time"))
        {
            parameters["end_time"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-p"))
        {
            parameters["application"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-n0"))
        {
            parameters["n0"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-nlevels"))
        {
            parameters["nlevels"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-nsamples"))
        {
            parameters["nsamples"] = argv[i+1];
        }
        else if(!strcmp(argv[i], "-q"))
        {
            parameters["q"] = argv[i+1];
        }
        else
        {
            std::cerr << "unknown argument" << argv[i] << "\n";
            exit(1);
        }
        i++;
    }
    double t0 = seconds();

    auto solver = std::make_shared<SolverStokes>(parameters);  
    StokesInfoSde info = solver->chorin_sde(false);
}
