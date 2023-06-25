//
//  SeamVector.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef SeamVector_h
#define SeamVector_h

#include  <armadillo>
#include  "../typedefs.hpp"

class NodeVector;

/*-------------------------------------------------*/
class SeamVector : public armavec
{
protected:
  // armavec _data;
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
  // SeamVector() : _data(), _n(), _ofs() {}
  // SeamVector(const armaicvec& n) : _data(arma::prod(n)), _n(n)
  SeamVector() : armavec(), _n(), _ofs() {}
  SeamVector(const armaicvec& n) : armavec(arma::prod(n)), _n(n)
  {
    set_ofs();
  }
  SeamVector& operator=(const SeamVector& v)
  {
    assert(arma::all(_n==v.n()));
    return *this;
  }
  void set_size(const armaicvec& n)
  {
    //    std::cerr << "SeamVector::set_size() n = " << n.t();
    _n = n;
    set_ofs();
    armavec::set_size(arma::prod(n));
  }
  void set_size(const SeamVector& u)
  {
    set_size(u.n());
  }
  int dim() const {return (int) _n.n_elem;}
  const armaicvec& n() const {return _n;}
  const armaicvec& ofs() const {return _ofs;}
  double& atp(int ix, int iy)
  {
    // _ofs[0]*(ix+1)+iy+1 = _ofs[0]*ix+iy + _ofs[0]+1
    return (*this)[_ofs[0]*ix+iy+_ofsp];
  }
  const double& atp(int ix, int iy) const
  {
    return (*this)[_ofs[0]*ix+iy+_ofsp];
  }
  double& atp(int ix, int iy, int iz)
  {
    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz+_ofsp];
  }
  const double& atp(int ix, int iy, int iz) const
  {
    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz+_ofsp];
  }
  void tovector(NodeVector& u) const;
  void fromvector(const NodeVector& u);
};
std::ostream& operator<<(std::ostream& os, const SeamVector& v);

#endif /* SeamVector_h */