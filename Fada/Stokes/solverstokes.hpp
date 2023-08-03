//
//  solverlaplace.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef  solverstokes_hpp
#define  solverstokes_hpp

#include  <sstream>
#include  "../gridvector.hpp"
#include  "../mgsolver.hpp"
#include  "../modelinterface.hpp"
#include  "../sparsematrix.hpp"
#include  "../uniformmultigrid.hpp"
#include  "model_p.hpp"
#include  "model_v.hpp"
#include  "model_stokes.hpp"
#include  "stokesapplication.hpp"
#include  "stokesvector.hpp"

struct StokesInfo {
    int niter;
    double niter_mean_v;
    double niter_mean_p;
    double err_v;
    double err_p;
    std::string toString(){std::stringstream ss; ss << "errors: "<<err_p<< " "<<err_v<<"\n"; ss<<"niter="<<niter<<"\n"; return ss.str();}
};

struct StokesInfoSde {
    int niter;
    double niter_mean_v;
    double niter_mean_p;
    armavec dt;
    armavec err_v;
    armavec err_p;
    StokesInfoSde(int n) : dt(n), err_v(n), err_p(n){}    
};

/*-------------------------------------------------*/
class SolverStokes
{
protected:
    double _end_time, _q, _dt;
    int _dim, _niter_out, _nt, _nsamples;
    std::string _applicationname, _datadir;
    std::map<std::string,std::string> _parameters;
    std::shared_ptr<UniformMultiGrid> _mggrid;
    std::shared_ptr<StokesApplication> _application;
    
    std::shared_ptr<ModelInterface> _model_all;
    std::shared_ptr<ModelInterface> _model_p;
    std::vector<std::shared_ptr<ModelInterface>> _model_v;
    
    // MgSolver _mgsolver_p;
    std::shared_ptr<MgSolver> _mgsolver_all;
    std::shared_ptr<MgSolver> _mgsolver_p;
    std::vector<std::shared_ptr<MgSolver>> _mgsolver_v;

    std::vector<std::shared_ptr<SparseMatrix const>> _B;

    void make_chorin_solver();
    void make_coupled_solver();
    
    const ModelP2d& getModelP() const
    {
        auto p = std::dynamic_pointer_cast<Model<ModelP2d, Vector<GridVector>> const>(_model_p);
        assert(p);
        return(p->get());
    }
    const ModelV2d& getModelV(int i) const
    {
        auto p = std::dynamic_pointer_cast<Model<ModelV2d, Vector<GridVector>> const>(_model_v[i]);
        assert(p);
        return(p->get());
    }
    // bool _directoryExists( std::string pzPath ) const;

    typedef std::vector<std::shared_ptr<GridVector>> VelocityVector;
    typedef std::shared_ptr<GridVector> PressureVector;
    VelocityVector newVelocityVector();
    PressureVector newPressureVector();

public:
    SolverStokes() : _model_p(nullptr), _model_v(), _mggrid(nullptr), _mgsolver_all(nullptr), _mgsolver_p(nullptr), _mgsolver_v(), _parameters() {}
    SolverStokes(const std::map<std::string,std::string>& parameters) : _model_p(nullptr), _model_v(), _mggrid(nullptr), _mgsolver_all(nullptr), _mgsolver_p(nullptr), _mgsolver_v()
    {
        _parameters = parameters;
        set_data(parameters);
    }      
    void set_data(const std::map<std::string,std::string>& parameters);
    std::string toString() const;
    std::shared_ptr<UniformMultiGrid const> get_mgrid() const {return _mggrid;} 

    int n_gridpoints() const
    {
        return _mggrid->get(0)->n_gridpoints();        
    }

    const GridVector& get_solution_p() const
    {
        auto p = std::dynamic_pointer_cast <const GridVector>(_mgsolver_p->getU());
        assert(p);
        return(*p);
    }
    GridVector& get_solution_p()
    {
        auto p = std::dynamic_pointer_cast <GridVector>(_mgsolver_p->getU());
        assert(p);
        return(*p);
    }
    const GridVector& get_rhs_p() const
    {
        auto p = std::dynamic_pointer_cast <const GridVector>(_mgsolver_p->getF());
        assert(p);
        return(*p);
    }
    GridVector& get_rhs_p()
    {
        auto p = std::dynamic_pointer_cast <GridVector>(_mgsolver_p->getF());
        assert(p);
        return(*p);
    }
    const GridVector& get_solution_v(int i) const
    {
        auto p = std::dynamic_pointer_cast <const GridVector>(_mgsolver_v[i]->getU());
        assert(p);
        return(*p);
    }
    GridVector& get_solution_v(int i)
    {
        auto p = std::dynamic_pointer_cast <GridVector>(_mgsolver_v[i]->getU());
        assert(p);
        return(*p);
    }
    const GridVector& get_rhs_v(int i) const
    {
        auto p = std::dynamic_pointer_cast <const GridVector>(_mgsolver_v[i]->getF());
        assert(p);
        return(*p);
    }
    GridVector& get_rhs_v(int i)
    {
        auto p = std::dynamic_pointer_cast <GridVector>(_mgsolver_v[i]->getF());
        assert(p);
        return(*p);
    }
    void save(SolverStokes::PressureVector p, SolverStokes::VelocityVector v, const armaicvec& iters) const;
    void load(std::vector<SolverStokes::VelocityVector>& vs, int& nend) const;
    void save_for_visu() const;
    std::map<std::string,double> compute_error() const;
    StokesInfo solve_stationary(bool print = true);
    StokesInfoSde chorin_sde(bool print=true);
    StokesInfo chorin_stationary(bool print=true);
};


#endif
