//
//  solverlaplace.cpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "solverlaplace.hpp"
#include  "q1.hpp"
#include  "../tokenize.hpp"
// #include  "../analyticalfunctioninterface.hpp"


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
std::string SolverLaplace::toString() const
{
    std::stringstream ss;
    ss << "fem=" << _model->toString();
    ss << "mggrid=" << _mggrid->toString();
    return(ss.str());
}

/*-------------------------------------------------*/
SolverLaplace::SolverLaplace(std::shared_ptr <MultiGridInterface> mggrid, const std::map<std::string,std::string>& parameters)
{
    set_data(mggrid, parameters);
}

/*-------------------------------------------------*/
void SolverLaplace::set_data(std::shared_ptr <MultiGridInterface> mggrid, const std::map<std::string,std::string>& parameters)
{
    int updatelength(0);
    for(std::map<std::string,std::string>::const_iterator p = parameters.begin(); p != parameters.end(); p++)
    {
        if(p->first=="updatelength")
        {
            updatelength = std::atoi(p->second.c_str());
        }
    }
    _mggrid = mggrid;
    size_t dim = mggrid->dim();
    _application = std::make_shared<LaplaceApplication>(dim, parameters.at("application"));
    if (dim == 2)
    {
        _model = std::make_shared<Model<Q12d, Vector<GridVector>>>("u", parameters, _application);
    }
    else if (dim == 3)
    {
        _model = std::make_shared<Model<Q13d, Vector<GridVector>>>("u", parameters, _application);
    }
    // _model->set_grid(mggrid->get(0));
    _mgsolver.set_sizes(_mggrid, _model, updatelength);
}

/*-------------------------------------------------*/
LaplaceInfo SolverLaplace::testsolve(bool print, std::string application)
{
    LaplaceInfo info;
    _model->rhs(_mgsolver.getF(), _mggrid->get(0), _application);
    _model->boundary_zero(_mgsolver.getF(), _mggrid->get(0));
    // // _u.set_size(_mggrid->get(0)->n());
    // // _u.fill(0);
    // // _f.set_size(_u);
    // if (application == "DirichletRhsOne")
    // {
    //     auto p = std::make_shared<ConstantFunction>();
    //     _model->rhs(get_rhs(), _mggrid->get(0), p);
    //     _model->boundary_zero(get_rhs(), _mggrid->get(0));
    //     _model->boundary_zero(get_u(), _mggrid->get(0));
    // }
    // else if (application == "Random")
    // {
    //     auto p = std::make_shared<RandomFunction>();
    //     _model->rhs(get_rhs(), _mggrid->get(0), p);
    //     _model->boundary_zero(get_rhs(), _mggrid->get(0));
    //     _model->boundary_zero(get_u(), _mggrid->get(0));
    // }
    // else if (application == "Linear")
    // {
    //     get_rhs()->fill(0);
    //     _model->boundary_linear(get_u(), _mggrid->get(0));
    //     _model->boundary_linear(get_rhs(), _mggrid->get(0));
    // }
    // else
    // {
    //     std::cerr << "unknwon application " << application << "\n";
    //     assert(0);
    //     exit(1);
    // }
    //
    std::shared_ptr <const GridVector> pf = std::dynamic_pointer_cast <const GridVector>(_mgsolver.getF());
    std::cerr << "f " << pf->min() << " " << pf->max() << "\n";
    info.niter = _mgsolver.solve(print);
    std::cerr << "u " << get_solution().min() << " " << get_solution().max() << "\n";
    if(_application->has_solution())
    {
        // std::shared_ptr <const GridVector> pu = std::dynamic_pointer_cast <const GridVector>(_mgsolver.getU());
        info.err = _model->compute_error(_mgsolver.getU(), _mggrid->get(0), _application).at("u");
    }
    return info;
}
