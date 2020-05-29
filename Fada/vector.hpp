//
//  Vector.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef Vector_h
#define Vector_h

#include  <armadillo>

typedef arma::Col<double> armavec;
typedef arma::Col<int> armaicvec;
//#include  "typedefs.hpp"

#include  "vectorinterface.hpp"

/*-------------------------------------------------*/
class Vector : public VectorInterface
{
protected:
  armavec _data;
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
//        std::cerr << "_n = " << _n.t();
//        std::cerr << "_ofs = " << _ofs.t();
//    std::cerr << "_ofsp = " << _ofsp << "\n";
  }

public:
  Vector() : _data(), _n(), _ofs() {}
  Vector(const armaicvec& n) : _data(arma::prod(n)), _n(n)
  {
    set_ofs();
//    assert(0);
  }
  Vector& operator=(const Vector& v)
  {
    _data = v.data();
//    armavec::operator=(v);
    assert(arma::all(_n==v.n()));
//    _n = v._n;
//    set_ofs();
    return *this;
  }
  Vector& operator=(const armavec& v)
  {
    _data = v;
//    armavec::operator=(v);
    return *this;
  }
    Vector& operator*=(double d)
    {
      _data *= d;
  //    armavec::operator=(v);
      return *this;
    }
  void set_size(const armaicvec& n)
  {
//    std::cerr << "Vector::set_size() n = " << n.t();
    _n = n;
    set_ofs();
//    armavec::set_size(arma::prod(n));
    _data.set_size(arma::prod(n));
  }
  void set_size(const Vector& u)
  {
    set_size(u.n());
  }
  void fill_bdry(double d=0);
  void fill_bdry2(double d=0);
  int dim() const {return (int) _n.n_elem;}
  const armaicvec& n() const {return _n;}
  const armaicvec& ofs() const {return _ofs;}
  armavec& data()
  {
    return _data;
//    armavec& t = static_cast<armavec&>(*this);
//    return t;
  }
  const armavec& data() const
  {
    return _data;
//    const armavec& t = static_cast<const armavec&>(*this);
//    return t;
  }
  double& at(int ix, int iy)
  {
    return _data[_ofs[0]*ix+iy];
  }
  const double& at(int ix, int iy) const
  {
    return _data[_ofs[0]*ix+iy];
  }
  double& atp(int ix, int iy)
  {
    // _ofs[0]*(ix+1)+iy+1 = _ofs[0]*ix+iy + _ofs[0]+1
    return _data[_ofs[0]*ix+iy+_ofsp];
  }
  const double& atp(int ix, int iy) const
  {
    return _data[_ofs[0]*ix+iy+_ofsp];
  }
  double& at(int ix, int iy, int iz)
  {
    return _data[_ofs[0]*ix+_ofs[1]*iy+iz];
  }
  const double& at(int ix, int iy, int iz) const
  {
    return _data[_ofs[0]*ix+_ofs[1]*iy+iz];
  }
  double& atp(int ix, int iy, int iz)
  {
    return _data[_ofs[0]*ix+_ofs[1]*iy+iz+_ofsp];
  }
  const double& atp(int ix, int iy, int iz) const
  {
    return _data[_ofs[0]*ix+_ofs[1]*iy+iz+_ofsp];
  }
  void fill(double d)
  {
    _data.fill(d);
  }
  void add(double d, const Vector& v)
  {
    _data += d*v.data();
//    this->arma() += d*v.arma();
  }
  double dot(const Vector& v) const
  {
//    return arma::dot(this->arma(),v.arma());
    return arma::dot(_data,v.data());
  }
  void scale(double d)
  {
//    this->arma() *= d;
    _data *= d;
  }
  double norm(double p=2) const
  {
    return arma::norm(_data,p);
//    return arma::norm(this->arma(),p);
  }
  void output(const std::string& filename) const
  {
    arma::hdf5_name spec(filename);
//    this->arma().save(spec);
    _data.save(spec);
  }
};
std::ostream& operator<<(std::ostream& os, const Vector& v);

#endif /* Vector_h */
