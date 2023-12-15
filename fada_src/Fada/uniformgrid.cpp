//
//  uniformgrid.cpp
//  Fada
//
//  Created by Roland Becker on 15/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <sstream>
#include "uniformgrid.hpp"

/*-------------------------------------------------*/
void UniformGrid::savehdf5(const std::string& filename) const
{
  _n.save(arma::hdf5_name(filename,"n"));
  _dx.save(arma::hdf5_name(filename,"dx",arma::hdf5_opts::append));
  _bounds.save(arma::hdf5_name(filename,"bounds",arma::hdf5_opts::append));
}
/*-------------------------------------------------*/
void UniformGrid::loadhdf5(const std::string& filename)
{
  // std::cerr <<"@@@@@ loadhdf5\n";
  _n.load(arma::hdf5_name(filename,"n"));
  _dx.load(arma::hdf5_name(filename,"dx"));
  _bounds.load(arma::hdf5_name(filename,"bounds"));
  // std::cerr <<"_n\n" << _n << "\n";
  // std::cerr <<"loadhdf5\n" << this->toString() << "\n";
}

/*-------------------------------------------------*/
std::string UniformGrid::toString() const
{
  std::stringstream ss;
  ss << "n="<<_n.t();
  ss << "dx="<<_dx.t();
  ss << "bounds="<<_bounds.t();
  return ss.str();
}

/*-------------------------------------------------*/
void UniformGrid::set_size(const armaicvec& n, std::shared_ptr<armamat> bounds)
{
  _n = n;
  arma::uword dim = _n.n_elem;
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
