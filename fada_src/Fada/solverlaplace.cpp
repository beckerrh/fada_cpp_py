//
//  solverlaplace.cpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "solverlaplace.hpp"
#include  "q0.hpp"
#include  "q1.hpp"
#include  "q1_shifted.hpp"
#include  "tokenize.hpp"


/*-------------------------------------------------*/
double Linear2D::operator()(double x, double y) const {return 3.3+2.2*x-1.1*y;}

double Sinus2D::operator()(double x, double y) const {return sin(A*x*x+B*y*y);}
double Sinus2D_rhs::operator()(double x, double y) const {return (4*A*A*x*x+4*B*B*y*y)*sin(A*x*x+B*y*y)-2*(A+B)*cos(A*x*x+B*y*y);}

// double Sinus2D::operator()(double x, double y) const {return sin(1.1*x-2.2*y);}
// double Sinus2D_rhs::operator()(double x, double y) const {return (2.2*2.2+1.1*1.1)*sin(1.1*x-2.2*y);}

double Sinus_per2D::operator()(double x, double y) const {return sin(2*A*M_PI*x)*sin(2*B*M_PI*y);}
double Sinus_per2D_rhs::operator()(double x, double y) const {return 4*(A*A+B*B)*M_PI*M_PI*sin(2*A*M_PI*x)*sin(2*B*M_PI*y);}

/*-------------------------------------------------*/
void LaplaceApplication::set_application(int dim, std::string name)
{
    std::vector<std::string> tokens = tokenize(name,"_");
    assert(tokens.size()==2);
    for(auto& b: *_boundaryconditions)
    {
        b[0] = tokens[1];
        b[1] = tokens[1];
    }        
    if(tokens[0]=="Linear")
    {
        _rhs["u"] = std::make_shared<ConstantFunction>(0.0);
        _sol["u"] = std::make_shared<Linear2D>();
        for(int i=0;i<dim;i++)
            for(int j=0;j<2;j++)
                _boundaryconditions->get_bf(i,j)["u"] = _sol["u"];
        // std::cerr << "_boundaryconditions = " << *_boundaryconditions <<"\n";
    }
    else if(tokens[0]=="Sinus")
    {
        if(tokens[1]=="dir")
        {
            _rhs["u"] = std::make_shared<Sinus2D_rhs>();
            _sol["u"] = std::make_shared<Sinus2D>();
            
        }
        else
        {
            _rhs["u"] = std::make_shared<Sinus_per2D_rhs>(1.,1.);
            _sol["u"] = std::make_shared<Sinus_per2D>(1.,1.);
            
        }
        for(int i=0;i<dim;i++)
            for(int j=0;j<2;j++)
                _boundaryconditions->get_bf(i,j)["u"] = _sol["u"];
        // std::cerr << "_boundaryconditions = " << *_boundaryconditions <<"\n";
    }
    else if(tokens[0]=="RhsOne")
    {
        _rhs["u"] = std::make_shared<ConstantFunction>();
    }
    else if(tokens[0]=="Random")
    {
        _rhs["u"] = std::make_shared<RandomFunction>();
    }
    else
    {
        _not_written_(" *** unknown name " + name);
    }
}
/*-------------------------------------------------*/
void SolverLaplace::define_parameters()
{
    Solver::define_parameters();
    _parameters.add("method", "Q1");
}

/*-------------------------------------------------*/
void SolverLaplace::set_data()
{
    std::cerr << "parameters:\n" << _parameters << "\n";
    Solver::set_data();
    int dim = _parameters.get<int>("dim");
    _application = std::make_shared<LaplaceApplication>(dim, _parameters.get<std::string>("application"));
    if (dim == 2)
    {
        if(_parameters.get<std::string>("method")=="Q1")
        {
            _model = std::make_shared<Model<Q12d, Vector<GridVector>>>("u", _parameters, _application);
        }
        else if(_parameters.get<std::string>("method")=="Q0")
        {
            _model = std::make_shared<Model<Q02d, Vector<GridVector>>>("u", _parameters, _application);
        }
        else if(_parameters.get<std::string>("method")=="Q1_0")
        {
            ParameterMap parameters2(_parameters);
            parameters2.set("direction", 0); 
            parameters2.set("dt", 0.0); 
            // std::map<std::string,std::string> parameters2(parameters.begin(), parameters.end());
            // parameters2["dt"] = "0.0";
            // parameters2["direction"] = "0";
            _model = std::make_shared<Model<Q1shifted2d, Vector<GridVector>>>("u", parameters2, _application);
        }
        else if(_parameters.get<std::string>("method")=="Q1_1")
        {
            ParameterMap parameters2(_parameters);
            parameters2.set("direction", 1); 
            parameters2.set("dt", 0.0); 
            _model = std::make_shared<Model<Q1shifted2d, Vector<GridVector>>>("u", parameters2, _application);
        }
        else
        {
            _not_written_("unknown method " +_parameters.get<std::string>("method"));
        }
    }
    else if (dim == 3)
    {
        _model = std::make_shared<Model<Q13d, Vector<GridVector>>>("u", _parameters, _application);
    }
    _mgsolver = std::make_shared<MgSolver>(_parameters.get<int>("mgtimer"), _parameters.get<int>("mgdebug"));
    _mgsolver->set_sizes(_mggrid, _model);    
}

// /*-------------------------------------------------*/
// SolverLaplace::SolverLaplace(const std::map<std::string,std::string>& parameters) : Solver(parameters)
// {
//     // set_data(mggrid, parameters);
//     _application = std::make_shared<LaplaceApplication>(_dim, parameters.at("application"));
//     if(parameters.find("method")==parameters.end())
//     {
//         _not_written_("needs paramater 'method'");
//     }
//     if (_dim == 2)
//     {
//         if(parameters.at("method")=="Q1")
//         {
//             _model = std::make_shared<Model<Q12d, Vector<GridVector>>>("u", parameters, _application);
//         }
//         else if(parameters.at("method")=="Q0")
//         {
//             _model = std::make_shared<Model<Q02d, Vector<GridVector>>>("u", parameters, _application);
//         }
//         else if(parameters.at("method")=="Q1_0")
//         {
//             std::map<std::string,std::string> parameters2(parameters.begin(), parameters.end());
//             parameters2["dt"] = "0.0";
//             parameters2["direction"] = "0";
//             _model = std::make_shared<Model<Q1shifted2d, Vector<GridVector>>>("u", parameters2, _application);
//         }
//         else if(parameters.at("method")=="Q1_1")
//         {
//             std::map<std::string,std::string> parameters2(parameters.begin(), parameters.end());
//             parameters2["dt"] = "0.0";
//             parameters2["direction"] = "1";
//             _model = std::make_shared<Model<Q1shifted2d, Vector<GridVector>>>("u", parameters2, _application);
//         }
//         else
//         {
//             _not_written_("unknown method " +parameters.at("method"));
//         }
//     }
//     else if (_dim == 3)
//     {
//         _model = std::make_shared<Model<Q13d, Vector<GridVector>>>("u", parameters, _application);
//     }
//     _mgsolver = std::make_shared<MgSolver>(_mgtimer, _mgdebug);
//     _mgsolver->set_sizes(_mggrid, _model);
//     // _mgsolver.set_sizes(_mggrid, _model, updatelength);
// }

/*-------------------------------------------------*/
LaplaceInfo SolverLaplace::testsolve(bool print, MgSolver::IterationInfo info)
{
    LaplaceInfo result;
    _model->rhs(_mgsolver->getF(), _mggrid->get(0), _application);
    _model->boundary(_mgsolver->getF(), _mggrid->get(0), _application->get_bc());
    std::shared_ptr <const GridVector> pf = std::dynamic_pointer_cast <const GridVector>(_mgsolver->getF());
    std::shared_ptr <const GridVector> pu = std::dynamic_pointer_cast <const GridVector>(_mgsolver->getU());
    // std::cerr << "f " << pf->min() << " " << pf->max() << "\n";
    std::cerr << "f " << pf->min() << " " << pf->max() << "\n";
    // std::cerr << "f " << *pf << "\n";
    result.niter = _mgsolver->solve(print, info);
    std::cerr << "u " << pu->min() << " " << pu->max() << "\n";
    // std::cerr << "f " << *pu << "\n";
    if(_application->has_solution())
    {
        result.has_error = true;
        result.err = _model->compute_error(_mgsolver->getU(), _mggrid->get(0), _application).at("u");
    }
    return result;
}
