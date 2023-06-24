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
  assert(0);
  return *this;
}

/*-------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const UniformMultiGrid& mg)
{
  os << mg.toString();
  return os;
}
/*-------------------------------------------------*/
std::string UniformMultiGrid::toString() const
{
  std::stringstream ss;
  ss << "nlevels = " << nlevels() << "\n";
  for(int l=0;l<nlevels();l++)
  {
    ss << get(l)->toString();
  }
  return ss.str();
}
/*-------------------------------------------------*/
void UniformMultiGrid::set_size(int nlevelmax, int nlevels, const armaicvec& n0, std::shared_ptr<armamat> bp)
{
  if(nlevels>nlevelmax)
  {
    std::cerr << nlevels << " = nlevels > nlevelmax = " << nlevelmax << "\n";
    exit(1);
  }
  if(nlevelmax==0)
  {
    std::cerr<< " nlevelmax = " << nlevelmax << "\n";
    exit(1);
  }
  arma::uword dim = n0.n_elem;
  armaicvec n(dim);
  std::shared_ptr<armamat> bounds;
  if(bp==nullptr)
  {
    bounds = std::unique_ptr<armamat>(new armamat(2, dim));
    for(int i=0;i<dim;i++)
    {
      (*bounds)(0, i) = 0.0;
      (*bounds)(1, i) = 1.0;
    }
  }
  else
  {
    bounds = bp;
  }
  _grids.resize(nlevels);
  for(int l=0;l<nlevels;l++)
  {
    for(int i=0;i<dim;i++)
    {
      n[i] = int(pow(2,nlevelmax-1-l))*(n0[i]-1)+1;
    }
    _grids[l] = std::unique_ptr<GridInterface>(new UniformGrid(n, bounds));
  }
}
