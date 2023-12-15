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
#include  "gridvector.hpp"
#include  "mgsolver.hpp"
#include  "modelinterface.hpp"
#include  "sparsematrix.hpp"
#include  "uniformmultigrid.hpp"
#include  "q0.hpp"
#include  "q1_shifted.hpp"
#include  "model_stokes.hpp"
#include  "solver.hpp"
#include  "stokesapplication.hpp"
#include  "stokesvector.hpp"


/*-------------------------------------------------*/
struct StokesInfo {
    int niter;
    double niter_mean_v, niter_mean_p;
    double err_v, err_p;
    std::string toString(){std::stringstream ss; ss << "errors p, v: "<<err_p<< " "<<err_v<<"\n"; ss<<"niter="<<niter<<"\n"; return ss.str();}
};

struct StokesInfoDynamic {
    int niter;
    double niter_mean_v, niter_mean_p;
    armavec dt;
    armavec err_v, err_p, diff_p, diff_v;
    StokesInfoSde(int n) : dt(n), err_v(n), err_p(n), diff_p(n), diff_v(n){}    
};

/*-------------------------------------------------*/
class SolverStokes : public Solver
{
protected:
    std::shared_ptr<ModelInterface> _model_p;
    std::vector<std::shared_ptr<ModelInterface>> _model_v;
    
    std::shared_ptr<MgSolver> _mgsolver_p;
    std::vector<std::shared_ptr<MgSolver>> _mgsolver_v;

    std::vector<std::shared_ptr<SparseMatrix const>> _B;
    
    const Q02d& getQ0() const
    {
        auto p = std::dynamic_pointer_cast<Model<Q02d, Vector<GridVector>> const>(_model_p);
        assert(p);
        return(p->get());
    }
    const Q1shifted2d& getModelV(int i) const
    {
        auto p = std::dynamic_pointer_cast<Model<Q1shifted2d, Vector<GridVector>> const>(_model_v[i]);
        assert(p);
        return(p->get());
    }

    typedef std::vector<std::shared_ptr<GridVector>> VelocityVector;
    typedef std::shared_ptr<GridVector> PressureVector;
    VelocityVector newVelocityVector();
    PressureVector newPressureVector();

public:
    SolverStokes() : Solver(), _model_p(nullptr), _model_v(), _mgsolver_p(nullptr), _mgsolver_v() {}
    // SolverStokes(const std::map<std::string,std::string>& parameters) : Solver(parameters) {}
    SolverStokes(const std::map<std::string,std::string>& parameters) : Solver() 
    {
       define_parameters();
       _parameters.fill_from_map(parameters);
       set_data();        
    }
    // SolverStokes(int argc, char** argv) : Solver()
    // {
    //    define_parameters();
    //    _parameters.fill_from_args(argc, argv);
    //    set_data();
    // }
    std::string  get_name() const {return "SolverStokes";}
    void define_parameters();
    void set_data();
    // SolverStokes(const std::map<std::string,std::string>& parameters);
    // void set_data(const std::map<std::string,std::string>& parameters);
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
    // void save_for_visu() const;
    std::map<std::string,double> compute_error() const;
    StokesInfo solve_stationary(bool print = true);
    StokesInfoSde chorin(bool print=true);
    StokesInfoSde chorin_sde(bool print=true);
    StokesInfo chorin_stationary(bool print=true);
};


#endif
