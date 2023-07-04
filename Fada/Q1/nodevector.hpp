//
//  NodeVector.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef NodeVector_h
#define NodeVector_h

#include  <armadillo>
#include  "../typedefs.hpp"
#include  "../vectorinterface.hpp"

/*-------------------------------------------------*/
class NodeVector : public armavec
{
protected:
  armaicvec _n, _ofs;
  int _ofsp;
  void set_ofs()
  {
    _ofs.set_size(_n.n_elem);
    _ofs.fill(1);
    for(int i=0;i<_n.n_elem;i++)
    {
      for(int j=i+1;j<_n.n_elem;j++)
      {
        _ofs[i] *= _n[j];
      }
    }
    _ofsp = 1;
    for(int i=0;i<_n.n_elem-1;i++)
    {
      _ofsp += _ofs[i];
    }
  }

public:
  NodeVector() : armavec(), _n(), _ofs() {}
  NodeVector(const armaicvec& n) : armavec(arma::prod(n)), _n(n)
  {
    set_ofs();
  }
  NodeVector& operator=(const NodeVector& v)
  {
    armavec::operator=(v);
    assert(arma::all(_n==v.n()));
    return *this;
  }
  NodeVector& operator=(const armavec& v)
  {
    armavec::operator=(v);
    return *this;
  }
  void set_size(const armaicvec& n)
  {
    //    std::cerr << "NodeVector::set_size() n = " << n.t();
    _n = n;
    set_ofs();
    armavec::set_size(arma::prod(n));
  }
  void set_size(const NodeVector& u)
  {
    set_size(u.n());
  }
  void fill_bdry(double d=0);
  int dim() const {return (int) _n.n_elem;}
  const armaicvec& n() const {return _n;}
  const armaicvec& ofs() const {return _ofs;}
  double& at(int ix, int iy)
  {
    return (*this)[_ofs[0]*ix+iy];
  }
  const double& at(int ix, int iy) const
  {
    return (*this)[_ofs[0]*ix+iy];
  }
  double& at(int ix, int iy, int iz)
  {
    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz];
  }
  const double& at(int ix, int iy, int iz) const
  {
    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz];
  }
  void scale(double d)
  {
    //    this->arma() *= d;
    // _data *= d;
    *this *= d;
  }
  double norm(double p=2) const
  {
    const armavec& tarma = static_cast<const armavec&>(*this);
    return arma::norm(tarma,p);
  }
  void add(double d, const NodeVector& v)
  {
    armavec& tarma = static_cast<armavec&>(*this);
    const armavec& varma = static_cast<const armavec&>(v);
    tarma += d*varma;
  }
  void equal(const NodeVector& v)
  {
    armavec::operator=(v);
  }
  double dot(const NodeVector& v) const
  {
    const armavec& tarma = static_cast<const armavec&>(*this);
    const armavec& varma = static_cast<const armavec&>(v);
    return arma::dot(tarma, varma);
  }
  void savehdf5(const std::string& filename) const
  {
    save(arma::hdf5_name(filename));
  }
};
std::ostream& operator<<(std::ostream& os, const NodeVector& v);

#endif /* NodeVector_h */
