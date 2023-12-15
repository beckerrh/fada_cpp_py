//
//  solverlaplace.cpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

// #include  <omp.h>
#include  <sstream>
#include  <iomanip>
#include  "solver.hpp"
#include  "sparsematrix.hpp"
#include  "boundary_conditions.hpp"
#include  "uniformgrid.hpp"
#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;

/*-------------------------------------------------*/
void Solver::save_for_visu() const
{
    std::string datadir = _parameters.get<std::string>("datadir");
    const fs::path p{datadir.c_str()};
    if(not fs::exists(p))
    {
        std::string command = "mkdir -p " + datadir;
        system( command.c_str() );
    }    
    std::string filename(datadir + "/solution.hdf");
    PointDataMap u_intp = _model->to_point_data(_mgsolver->getU(), _mggrid->get(0));
    armaicvec n(_parameters.get<int>("dim")+1);
    int count=0;        
    for(auto p:u_intp)
    {
        n[count] = p.second->n_elem;
        count++;        
    }
    arma::hdf5_name spec(filename, "n");
    n.save(spec);
    for(auto p:u_intp)
    {
        arma::hdf5_name spec(filename, p.first, arma::hdf5_opts::append);
        p.second->save(spec);        
    }
    _mggrid->get(0)->savehdf5(datadir + "/grid.hdf");
}

/*-------------------------------------------------*/
std::string Solver::toString() const
{
    std::stringstream ss;
    ss << "application=" << _application->toString();
    ss << "model=" << _model->toString();
    ss << "mggrid=" << _mggrid->toString();
    return(ss.str());
}

/*-------------------------------------------------*/
void Solver::define_parameters()
{
    _parameters.add("application", "Sinus_per");
    _parameters.add("datadir", "none");
    _parameters.add("mgtimer", 0);
    _parameters.add("mgdebug", 0);
    _parameters.add("n0", 4);
    _parameters.add("dim", 2);
    _parameters.add("nlevels", 2);
    _parameters.add("dt", 0.0);
    _parameters.add("nt", 1000);
    _parameters.add("niter_out", 100);
    _parameters.add("nsamples", 100);
    _parameters.add("end_time", 1.0);
    _parameters.add("q", 0.0);
}

/*-------------------------------------------------*/
void Solver::set_data()
{
    if(_parameters.get<std::string>("datadir")=="none")
    {
        _parameters.set("datadir", "datadir_" +get_name()+"_"+_parameters.get<std::string>("application"));
    }
    int n00 = _parameters.get<int>("n0");
    armaicvec n0 = {n00, n00};
    int ref_factor = 2;
    _mggrid = std::make_shared<UniformMultiGrid>(_parameters.get<int>("nlevels"), n0, nullptr, ref_factor);
}
