//
//  uniformmultigrid.cpp
//  Fada
//
//  Created by Roland Becker on 15/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "uniformgrid.hpp"
#include  "uniformmultigrid.hpp"
#include  <sstream>

/*-------------------------------------------------*/
UniformMultiGrid& UniformMultiGrid::operator=(const UniformMultiGrid& umg)
{
  _n = umg._n;
  _dx = umg._dx;
  _bounds = umg._bounds;
  _nall = umg._nall;
  return *this;
}

/*-------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const UniformMultiGrid& mg)
{
  os << "nlevels = " << mg.nlevels() << "\n";
  for(int l=0;l<mg.nlevels();l++)
  {
    os << mg.n(l).t();
  }
  return os;
}
/*-------------------------------------------------*/
std::string UniformMultiGrid::toString() const
{
  std::stringstream ss;
  ss << *this;
  return ss.str();
}
/*-------------------------------------------------*/
void UniformMultiGrid::_set_size(std::shared_ptr<armamat> bounds)
{
  int dim = (int) _n.n_rows;
  _nall = arma::prod(_n, 0);
  if(bounds==nullptr)
  {
    _bounds.set_size(2, dim);
    for(int i=0;i<dim;i++)
    {
      _bounds(0, i) = 0.0;
      _bounds(1, i) = 1.0;
    }
  }
  else _bounds = *bounds;
  _dx.set_size(dim, _n.n_cols);
  for(int l=0;l<_n.n_cols;l++)
  {
    for(int i=0;i<dim;i++) _dx(i,l) = (_bounds(1,i)-_bounds(0,i))/(_n(i,l)-1);
  }
}
/*-------------------------------------------------*/
void UniformMultiGrid::set_size(const armaimat&  n, std::shared_ptr<armamat> bounds)
{
  _n = n;
  _set_size(bounds);
}

/*-------------------------------------------------*/
void UniformMultiGrid::set_size(int nlevelmax, int nlevels, const armaicvec& n0, std::shared_ptr<armamat> bounds)
{
  if(nlevels>nlevelmax)
  {
    std::cerr << nlevels << " = nlevels > nlevelmax = " << nlevelmax << "\n";
    exit(1);
  }
  int dim = (int) n0.n_elem;
  _n.set_size(dim, nlevels);
  for(int i=0;i<dim;i++)
  {
    for(int l=0;l<nlevels;l++)
    {
//      _n(i,l) = int(pow(2,l))*(n0[i]-1)+1;
      _n(i,l) = int(pow(2,nlevelmax-1-l))*(n0[i]-1)+1;
    }
  }
  _set_size(bounds);
}

/*-------------------------------------------------*/
std::shared_ptr<UniformGrid> UniformMultiGrid::grid(int l) const
{
  std::shared_ptr<armamat> p = std::make_shared<armamat>(std::move(_bounds));
  return std::shared_ptr<UniformGrid>(new UniformGrid(n(l), p));
  // ne marche pas à cause du premier argument ('expects an lvalue')
//  return std::make_shared<UniformGrid>(n(l), p);
}
