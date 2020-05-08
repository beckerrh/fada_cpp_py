//
//  Vector.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef Vector_h
#define Vector_h

#include  "typedefs.hpp"

/*-------------------------------------------------*/
class Vector : public armavec
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
    for(int i=1;i<_n.n_elem;i++)
    {
      _ofsp += _ofs[i];
    }
    //    std::cerr << "_n = " << _n.t();
    //    std::cerr << "_ofs = " << _ofs.t();
  }

public:
  Vector() : armavec(), _n(), _ofs() {}
  Vector(const armaicvec& n) : armavec(arma::prod(n)), _n(n)
  {
    set_ofs();
    assert(0);
  }
  Vector& operator=(const Vector& v)
  {
    armavec::operator=(v);
    assert(arma::all(_n==v.n()));
//    _n = v._n;
//    set_ofs();
    return *this;
  }
  Vector& operator=(const armavec& v)
  {
    armavec::operator=(v);
    return *this;
  }
  void set_size(const armaicvec& n)
  {
//    std::cerr << "Vector::set_size() n = " << n.t();
    _n = n;
    set_ofs();
    armavec::set_size(arma::prod(n));
  }
  void set_size(const Vector& u)
  {
    set_size(u.n());
  }
  int dim() const {return _n.n_elem;}
//  int n(int i) const {return _n[i];}
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
//  double& operator()(int ix, int iy)
//  {
//    return (*this)[_ofs[0]*ix+iy];
//  }
//  const double& operator()(int ix, int iy) const
//  {
//    return (*this)[_ofs[0]*ix+iy];
//  }
//  double& operator()(int ix, int iy, int iz)
//  {
//    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz];
//  }
//  const double& operator()(int ix, int iy, int iz) const
//  {
//    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz];
//  }
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
  double& atp(int ix, int iy)
  {
    return (*this)[_ofs[0]*ix+iy+_ofsp];
  }
  const double& atp(int ix, int iy) const
  {
//    return (*this)[_ofs[0]*ix+iy+_ofs[0]+1];
    return (*this)[_ofs[0]*ix+iy+_ofsp];
  }
  double& atp(int ix, int iy, int iz)
  {
    // apparemment c'est plus rapide comme ça (??)
//    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz+1+_ofs[0]+_ofs[1]];
    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz+_ofsp];
  }
  const double& atp(int ix, int iy, int iz) const
  {
    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz+_ofsp];
//    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz+1+_ofs[0]+_ofs[1]];
//    return (*this)[_ofs[0]*(ix+1)+_ofs[1]*(iy+1)+iz+1];
  }
  void add(double d, const Vector& v)
  {
    this->arma() += d*v.arma();
  }
  double dot(const Vector& v) const
  {
    return arma::dot(this->arma(),v.arma());
  }
  void scale(double d)
  {
    this->arma() *= d;
  }
  double norm(double p=2) const
  {
    return arma::norm(this->arma(),p);
  }
  void output(const std::string& filename) const
  {
    arma::hdf5_name spec(filename);
    this->arma().save(spec);
  }
};
std::ostream& operator<<(std::ostream& os, const Vector& v);

#endif /* Vector_h */
