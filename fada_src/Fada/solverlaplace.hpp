//
//  solverlaplace.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef solverlaplace_hpp
#define solverlaplace_hpp

#include  "applicationinterface.hpp"
#include  "modelinterface.hpp"
#include  "multigridinterface.hpp"
#include  "mgsolver.hpp"
#include  "gridvector.hpp"
#include  "solver.hpp"


/*-------------------------------------------------*/
class Linear2D: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Sinus2D: public AnalyticalFunctionInterface
{
    double A,B;
public:
    Sinus2D(double A_=2, double B_=-1) : AnalyticalFunctionInterface(), A(A_), B(B_) {}
    double operator()(double x, double y) const;
};
class Sinus2D_rhs: public AnalyticalFunctionInterface
{
    double A,B;
public:
    Sinus2D_rhs(double A_=2, double B_=-1) : AnalyticalFunctionInterface(), A(A_), B(B_) {}
    double operator()(double x, double y) const;
};
class Sinus_per2D: public AnalyticalFunctionInterface
{
public:
    double A,B;
    Sinus_per2D(double A_, double B_) : AnalyticalFunctionInterface(), A(A_), B(B_) {}
    double operator()(double x, double y) const;
};
class Sinus_per2D_rhs: public AnalyticalFunctionInterface
{
public:
    double A,B;
    Sinus_per2D_rhs(double A_, double B_) : AnalyticalFunctionInterface(), A(A_), B(B_) {}
    double operator()(double x, double y) const;
};

class LaplaceApplication: public ApplicationInterface
{
protected:
    void set_application(int dim, std::string name);
public:
    LaplaceApplication(int dim, std::string name) : ApplicationInterface("Laplace", dim)
    {
        set_application(dim, name);
    }
};

/*-------------------------------------------------*/
struct LaplaceInfo {
    bool has_error;
    int niter;
    double err;
    std::string toString(){std::stringstream ss; if(has_error) ss << "error: "<<err<< "\n"; ss<<"niter="<<niter<<"\n"; return ss.str();}
    LaplaceInfo() : has_error(false) {}
};


/*-------------------------------------------------*/
class SolverLaplace : public Solver
{
protected:

public:
    SolverLaplace() : Solver() {}
    SolverLaplace(const SolverLaplace& solver) : Solver(solver) {}
    SolverLaplace(const std::map<std::string,std::string>& parameters) : Solver() 
    {
       define_parameters();
       _parameters.fill_from_map(parameters);
       set_data();        
    }
    std::string  get_name() const {return "SolverLaplace";}
    // SolverLaplace(int argc, char** argv) : Solver()
    // {
    //    define_parameters();
    //    _parameters.fill_from_args(argc, argv);
    //    set_data();
    // }
    // Solver(const std::map<std::string,std::string>& parameters) : _mggrid(nullptr), _application(nullptr), _model(nullptr), _mgsolver(nullptr), _parameters()
    // {
    //    // define_parameters();
    //    _parameters.fill_from_map(parameters);
    //    set_data();
    // }
    // SolverLaplace(const std::map<std::string,std::string>& parameters);
    void define_parameters();
    void set_data();

    LaplaceInfo testsolve(bool print = true, MgSolver::IterationInfo info=MgSolver::IterationInfo());
    std::shared_ptr<armavec const> get_solution_nodes() const 
    {
        PointDataMap u_intp = _model->to_point_data(_mgsolver->getU(), _mggrid->get(0));
        assert(u_intp["u"]);
        return std::dynamic_pointer_cast<armavec const>(u_intp["u"]);        
    }
};


#endif
