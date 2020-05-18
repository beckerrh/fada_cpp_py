//
//  uniformmultigrid.cpp
//  Fada
//
//  Created by Roland Becker on 15/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "uniformgrid.hpp"
#include  "uniformmultigrid.hpp"

/*-------------------------------------------------*/
void UniformMultiGrid::set_size(int nlevels, const armaicvec& n0, const armamat* bounds)
{
  int dim = n0.n_elem;
  _n.set_size(dim, nlevels);
  
  for(int i=0;i<dim;i++)
  {
    for(int l=0;l<nlevels;l++)
    {
      _n(i,l) = int(pow(2,l))*(n0[i]-1)+1;
    }
  }
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
  _dx.set_size(dim, nlevels);
  for(int l=0;l<nlevels;l++)
  {
    for(int i=0;i<dim;i++) _dx(i,l) = (_bounds(1,i)-_bounds(0,i))/(_n(i,l)-1);
  }
//  std::cerr << "_dx = " << _dx << "\n";
}

/*-------------------------------------------------*/
UniformGrid UniformMultiGrid::grid(int l) const
{
  if(l==-1) l = maxlevel();
  return UniformGrid(n(l), &_bounds);
}
