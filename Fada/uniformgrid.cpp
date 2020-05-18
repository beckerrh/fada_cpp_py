//
//  uniformgrid.cpp
//  Fada
//
//  Created by Roland Becker on 15/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include "uniformgrid.hpp"

/*-------------------------------------------------*/
void UniformGrid::set_size(const armaicvec& n, const armamat* bounds)
{
  _n = n;
  int dim = _n.n_elem;
  _dx.set_size(dim);
  if(bounds==nullptr)
  {
    _bounds = nullptr;
    for(int i=0;i<dim;i++) _dx[i] = 1.0/(_n[i]-1);
  }
  else
  {
    _bounds = *bounds;
    for(int i=0;i<dim;i++) _dx[i] = (_bounds(1,i)-_bounds(0,i))/(_n[i]-1);
  }
}
