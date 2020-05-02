//
//  vector.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef vector_h
#define vector_h

#include  "typedefs.hpp"

/*-------------------------------------------------*/
class vector : public armavec
{
private:
  armaicvec _n, _ofs;

public:
  vector(int dim) : armavec(), _n(dim), _ofs(dim) {}
  vector(const armaicvec& n) : armavec(arma::prod(n)), _n(n), _ofs(arma::cumprod(arma::reverse(n))) {}
  vector& operator=(const vector& v)
  {
    armavec::operator=(v);
    _n = v._n;
    _ofs = arma::cumprod(arma::reverse(_n));
    return *this;
  }
  vector& operator=(const armavec& v)
  {
    armavec::operator=(v);
    return *this;
  }
  const armaicvec& n() const {return _n;}
  const armaicvec& ofs() const {return _ofs;}
  armavec& arma()
  {
    armavec& t = static_cast<armavec&>(*this);
    return t;
  }
  const armavec& arma() const
  {
    const armavec& t = static_cast<const armavec&>(*this);
    return t;
  }
  double& operator()(int ix, int iy)
  {
    return (*this)[_ofs[0]*ix+iy];
  }
  const double& operator()(int ix, int iy) const
  {
    return (*this)[_ofs[0]*ix+iy];
  }
  double& operator()(int ix, int iy, int iz)
  {
    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz];
  }
  const double& operator()(int ix, int iy, int iz) const
  {
    return (*this)[_ofs[0]*ix+_ofs[2]*iy+iz];
  }

};
std::ostream& operator<<(std::ostream& os, const vector& v)
{
  const armavec& tarma =static_cast<const armavec&>(v);
  os << tarma.t()<< "n=" << v.n().t()<< "ofs=" << v.ofs().t();
  return os;
}

#endif /* vector_h */
