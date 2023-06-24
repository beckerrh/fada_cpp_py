//
//  solverlaplace.cpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <sstream>
#include  "solverlaplace.hpp"
#include  "q1.hpp"

/*-------------------------------------------------*/
SolverLaplace::SolverLaplace(std::shared_ptr<MultiGridInterface> mggrid, std::string femtype, std::string matrixtype, std::string smoothertype, int updatelength)
{
  set_data(mggrid, femtype, matrixtype, smoothertype, updatelength);
}

/*-------------------------------------------------*/
std::string SolverLaplace::toString() const
{
  std::stringstream ss;
  ss << "fem="<<_fem->toString();
  ss << "mggrid="<<_mggrid->toString();
  return ss.str();
}

/*-------------------------------------------------*/
void SolverLaplace::set_data (std::shared_ptr<MultiGridInterface> mggrid, std::string femtype, std::string matrixtype, std::string smoothertype, int updatelength)
{
  _mggrid = mggrid;
  size_t dim = mggrid->dim();
  if(dim==2)
  {
    if(femtype=="Q1")
    {
//      _fem = std::unique_ptr<FiniteElementInterface>(new Q12d(matrixtype));
      _fem = std::unique_ptr<FiniteElementInterface>(new FiniteElement<Q12d,NodeVector>(matrixtype));
//      matrixtype
    }
    else
    {
      std::cerr << "unknown fem '" << femtype<<"'\n";
    }
  }
  else if(dim==3)
  {
    if(femtype=="Q1")
    {
//      _fem = std::unique_ptr<FiniteElementInterface>(new Q13d(matrixtype));
      _fem = std::unique_ptr<FiniteElementInterface>(new FiniteElement<Q13d,NodeVector>(matrixtype));
    }
    else
    {
      std::cerr << "unknown fem '" << femtype<<"'\n";
    }
  }
  _fem->set_grid(mggrid->get(0));
  _mgsolver.set_sizes(_mggrid, _fem, smoothertype, updatelength);
}
/*-------------------------------------------------*/
int SolverLaplace::testsolve(bool print, std::string problem)
{
  _u.set_size(_mggrid->get(0)->n());
   _u.fill(0);
   _f.set_size(_u);
   std::cerr << " u " << _u.data().n_elem << " " << _mggrid->get(0)->nall() << "\n";
   if(problem=="DirichletRhsOne")
   {
     _fem->rhs_one(_f);
   }
   else if(problem=="Random")
   {
     _fem->rhs_random(_f);
   }
   else if(problem=="Linear")
   {
     _f.fill(0);
     _fem->boundary(_u);
     _fem->boundary(_f);
   }
//
//   _fem->vector2vectormg(_mggrid, 0, fmg(0), _f);
//   _fem->vector2vectormg(_mggrid, 0, umg(0), _u);
  int iter = _mgsolver.solve(_u, _f, print);
//   _fem->vectormg2vector(_mggrid, 0, _u, umg(0));
   return iter;
//  return 0;
}
