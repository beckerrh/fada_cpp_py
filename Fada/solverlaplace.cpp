//
//  solverlaplace.cpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "solverlaplace.hpp"
#include  "Q1/q1.hpp"
#include  "boundary_conditions.hpp"

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
  std::string bc("dirichlet");
  int updatelength(0);
  for(std::map<std::string,std::string>::const_iterator p = parameters.begin(); p != parameters.end(); p++)
  {
    if(p->first=="updatelength")
    {
      updatelength = std::atoi(p->second.c_str());
    }
    else if(p->first=="boundary_condition")
    {
      bc = p->second.c_str();
    }
  }
  _mggrid = mggrid;
  size_t dim = mggrid->dim();
  auto boundaryconditions = std::make_shared<BoundaryConditions>(dim);
  if(bc=="periodic")
  {
    for(auto& b: *boundaryconditions)
    {
      b[0] = "periodic";
      b[1] = "periodic";
    }
  }
  else if(bc=="dirichlet")
  {
    for(auto& b: *boundaryconditions)
    {
      b[0] = "dir";
      b[1] = "dir";
    }
  }
  else
  {
    std::cerr << bc << " *** unknown bc " << bc << "\n";
    exit(1);
  }
  if (dim == 2)
  {
    _model = std::make_shared<Model<Q12d, Vector<GridVector>>>(parameters, boundaryconditions);
  }
  else if (dim == 3)
  {
    _model = std::make_shared<Model<Q13d, Vector<GridVector>>>(parameters, boundaryconditions);
  }
  _model->set_grid(mggrid->get(0));
  _mgsolver.set_sizes(_mggrid, _model, updatelength);
}

/*-------------------------------------------------*/
int SolverLaplace::testsolve(bool print, std::string problem)
{
  // _u.set_size(_mggrid->get(0)->n());
  // _u.fill(0);
  // _f.set_size(_u);
  if (problem == "DirichletRhsOne")
  {
    _model->rhs_one(get_rhs());
    _model->boundary_zero(get_rhs());
    _model->boundary_zero(get_u());
  }
  else if (problem == "Random")
  {
    _model->rhs_random(get_rhs());
    _model->boundary_zero(get_rhs());
    _model->boundary_zero(get_u());
  }
  else if (problem == "Linear")
  {
    get_rhs()->fill(0);
    _model->boundary_linear(get_u());
    _model->boundary_linear(get_rhs());
  }
  else
  {
    std::cerr << "unknwon problem " << problem << "\n";
    assert(0);
    exit(1);
  }
//
  std::shared_ptr <const GridVector> p = std::dynamic_pointer_cast <const GridVector>(get_rhs());
  std::cerr << "f " << p->min() << " " << p->max() << "\n";
  int iter = _mgsolver.solve(print);
  std::cerr << "u " << get_solution().min() << " " << get_solution().max() << "\n";
  return iter;
//  return 0;
}
