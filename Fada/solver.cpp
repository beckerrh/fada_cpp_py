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
    const fs::path p{_datadir.c_str()};
    if(not fs::exists(p))
    {
        std::string command = "mkdir -p " + _datadir;
        system( command.c_str() );
    }    
    std::string filename(_datadir + "/solution.hdf");
    PointDataMap u_intp = _model->to_point_data(_mgsolver->getU(), _mggrid->get(0));
    armaicvec n(_dim+1);
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
    _mggrid->get(0)->savehdf5(_datadir + "/grid.hdf");
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
void Solver::set_data(const std::map<std::string,std::string>& parameters)
{
    int nlevels(8);
    _applicationname = "Sinus_per";
    _datadir = "";
    _dim=2;
    _dt=0.0;
    _mgtimer=true; _mgdebug=false;
    int n00 = 4;
    for(std::map<std::string,std::string>::const_iterator p = parameters.begin(); p != parameters.end(); p++)
    {
        if(p->first=="nlevels")
        {
            nlevels = std::stoi(p->second);
        }
        else if(p->first=="dim")
        {
            _dim = std::stoi(p->second);
        }
        else if(p->first=="application")
        {
            _applicationname = p->second.c_str();
        }
        else if(p->first=="datadir")
        {
            _datadir = p->second.c_str();
        }
        else if(p->first=="n0")
        {
            n00 = std::stoi(p->second);
        }
        else if(p->first=="dt")
        {
            _dt = std::stod(p->second);
        }
        else if(p->first=="dt")
        {
            _dt = std::stod(p->second);
        }
        else if(p->first=="mgdebug")
        {
            _mgdebug = p->second=="true";
        }
        else if(p->first=="mgtimer")
        {
            _mgtimer = p->second=="true";
        }
    }
    _parameters["dt"] = std::to_string(_dt);
    if(_datadir=="")
    {
        _datadir = "datadir_" +_applicationname;
    }
    try{
        _end_time = std::stod(_parameters.at("end_time"));
        _nt = std::stoi(_parameters.at("nt"));
        _niter_out = std::stoi(_parameters.at("niter_out"));
        _nsamples = std::stoi(_parameters.at("nsamples"));
        _q = std::stod(_parameters.at("q"));        
    }
    catch(...){
        // std::cerr << ("no sde!");
        _end_time = 1.0; _nt=0; _niter_out=0; _nsamples=0; _q=0.0;
    }         
    armaicvec n0 = {n00, n00};
    int ref_factor = 2;
    _mggrid = std::make_shared<UniformMultiGrid>(nlevels, n0, nullptr, ref_factor);
}
