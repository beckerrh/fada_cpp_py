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
protected:
  armaicvec _n, _ofs;
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
//    std::cerr << "_n = " << _n.t();
//    std::cerr << "_ofs = " << _ofs.t();
  }
public:
  vector() : armavec(), _n(), _ofs() {}
  vector(const armaicvec& n) : armavec(arma::prod(n)), _n(n)
  {
    set_ofs();
    assert(0);
  }
  vector& operator=(const vector& v)
  {
    armavec::operator=(v);
    _n = v._n;
    set_ofs();
    return *this;
  }
  vector& operator=(const armavec& v)
  {
    armavec::operator=(v);
    return *this;
  }
  int dim() const {return _n.n_elem;}
  int n(int i) const {return _n[i];}
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
    return (*this)[_ofs[0]*ix+_ofs[1]*iy+iz];
  }
  void set_size(const armaicvec& n)
  {
    _n = n;
    set_ofs();
    armavec::set_size(arma::prod(n));
  }
  void set_size(const vector& u)
  {
    set_size(u.n());
  }
  void add(double d, const vector& v)
  {
    this->arma() += d*v.arma();
  }
  double dot(const vector& v) const
  {
    return arma::dot(this->arma(),v.arma());
  }
  void output(const std::string& filename) const
  {
    arma::hdf5_name spec(filename);
    this->arma().save(spec);
  }
};
std::ostream& operator<<(std::ostream& os, const vector& v);

#endif /* vector_h */
