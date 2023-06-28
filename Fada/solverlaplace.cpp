//
//  solverlaplace.cpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "solverlaplace.hpp"
#include  "Q1/q1.hpp"

/*-------------------------------------------------*/
std::string SolverLaplace::toString() const
{
  std::stringstream ss;
  ss << "fem=" << _model->toString();
  ss << "mggrid=" << _mggrid->toString();
  return(ss.str());
}

// /*-------------------------------------------------*/
// SolverLaplace::SolverLaplace(std::shared_ptr <MultiGridInterface> mggrid, std::string stenciltype, std::string matrixtype, std::string smoothertype, int updatelength)
// {
//   set_data(mggrid, stenciltype, matrixtype, smoothertype, updatelength);
// }

// /*-------------------------------------------------*/
// void SolverLaplace::set_data(std::shared_ptr <MultiGridInterface> mggrid, std::string stenciltype, std::string matrixtype, std::string smoothertype, int updatelength)
// {
//   _mggrid = mggrid;
//   size_t dim = mggrid->dim();
//   if (dim == 2)
//   {
//     _model = std::make_shared<Model<Q12d, Vector<NodeVector>>>(stenciltype, matrixtype);
//     // _model = std::shared_ptr <ModelInterface>(new Model <Q12d, Vector<NodeVector>>(matrixtype));
//   }
//   else if (dim == 3)
//   {
//     _model = std::make_shared<Model<Q13d, Vector<NodeVector>>>(stenciltype, matrixtype);
//     // _model = std::shared_ptr <ModelInterface>(new Model <Q13d, Vector<NodeVector>>(matrixtype));
//   }
//   _model->set_grid(mggrid->get(0));
//   _mgsolver.set_sizes(_mggrid, _model, smoothertype, updatelength);
// }

/*-------------------------------------------------*/
SolverLaplace::SolverLaplace(std::shared_ptr <MultiGridInterface> mggrid, const std::map<std::string,std::string>& parameters)
{
  set_data(mggrid, parameters);
}

/*-------------------------------------------------*/
void SolverLaplace::set_data(std::shared_ptr <MultiGridInterface> mggrid, const std::map<std::string,std::string>& parameters)
{
  std::string stenciltype("Trapez"), matrixtype("stencil"), smoothertype("GS");
  int updatelength(0);
  for(std::map<std::string,std::string>::const_iterator p = parameters.begin(); p != parameters.end(); p++)
  {
    if(p->first=="stenciltype")
    {
      stenciltype = p->second;
    }
    else if(p->first=="matrixtype")
    {
      matrixtype = p->second;
    }
    else if(p->first=="smoothertype")
    {
      smoothertype = p->second;
    }
    else if(p->first=="updatelength")
    {
      updatelength = std::atoi(p->second.c_str());
    }
    else
    {
      std::cerr<< "unknwon parameter " << p->first << " value=" << p->second << "\n";
      exit(1);
    }
  }
  _mggrid = mggrid;
  size_t dim = mggrid->dim();
  if (dim == 2)
  {
    _model = std::make_shared<Model<Q12d, Vector<NodeVector>>>(stenciltype, matrixtype);
  }
  else if (dim == 3)
  {
    _model = std::make_shared<Model<Q13d, Vector<NodeVector>>>(stenciltype, matrixtype);
  }
  _model->set_grid(mggrid->get(0));
  _mgsolver.set_sizes(_mggrid, _model, smoothertype, updatelength);
}

/*-------------------------------------------------*/
int SolverLaplace::testsolve(bool print, std::string problem)
{
  // _u.set_size(_mggrid->get(0)->n());
  // _u.fill(0);
  // _f.set_size(_u);
  // std::cerr << "u " << _u.data().n_elem << " " << _mggrid->get(0)->nall() << "\n";
  if (problem == "DirichletRhsOne")
  {
    _model->rhs_one(get_rhs());
    _model->boundary(get_rhs());
    _model->boundary(get_u());
  }
  else if (problem == "Random")
  {
    _model->rhs_random(get_rhs());
    _model->boundary(get_rhs());
    _model->boundary(get_u());
  }
  else if (problem == "Linear")
  {
    get_rhs()->fill(0);
    _model->boundary(get_u());
    _model->boundary(get_rhs());
  }
//
//   _model->vector2vectormg(_mggrid, 0, fmg(0), _f);
//   _model->vector2vectormg(_mggrid, 0, umg(0), _u);
  int iter = _mgsolver.solve(print);
  std::cerr << "u " << get_solution().min() << " " << get_solution().max() << "\n";
//   _model->vectormg2vector(_mggrid, 0, _u, umg(0));
  return iter;
//  return 0;
}
