//
//  solverlaplace.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef  solver_hpp
#define  solver_hpp

#include  <sstream>
#include  "gridvector.hpp"
#include  "mgsolver.hpp"
#include  "modelinterface.hpp"
#include  "sparsematrix.hpp"
#include  "uniformmultigrid.hpp"

/*-------------------------------------------------*/
class Solver
{
protected:
    bool _mgtimer, _mgdebug;
    double _end_time, _q, _dt;
    int _dim, _niter_out, _nt, _nsamples;
    std::string _applicationname, _datadir;
    std::map<std::string,std::string> _parameters;
    std::shared_ptr<UniformMultiGrid> _mggrid;
    std::shared_ptr<ApplicationInterface> _application;    
    std::shared_ptr<ModelInterface> _model;
    std::shared_ptr<MgSolver> _mgsolver;

public:
    Solver() : _mggrid(nullptr), _application(nullptr), _model(nullptr), _mgsolver(nullptr), _parameters() {}
    Solver(const std::map<std::string,std::string>& parameters) : _mggrid(nullptr), _application(nullptr), _model(nullptr), _mgsolver(nullptr)
    {
        _parameters = parameters;
        set_data(parameters);
    }      
    void set_data(const std::map<std::string,std::string>& parameters);
    std::string toString() const;

    int n_gridpoints() const
    {
        return _mggrid->get(0)->n_gridpoints();        
    }
    void save_for_visu() const;
    std::shared_ptr<UniformMultiGrid const> get_mgrid() const {return _mggrid;}
    std::shared_ptr<ApplicationInterface const> getApplication() const {return _application;}
    std::shared_ptr<ModelInterface const> getModel() const {return _model;}
};


#endif
