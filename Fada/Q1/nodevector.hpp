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
class NodeVector : public VectorInterface
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
  NodeVector() : _data(), _n(), _ofs() {}
  NodeVector(const armaicvec& n) : _data(arma::prod(n)), _n(n)
  {
    set_ofs();
    //    assert(0);
  }
  NodeVector& operator=(const NodeVector& v)
  {
    _data = v.data();
    assert(arma::all(_n==v.n()));
    return *this;
  }
  NodeVector& operator=(const armavec& v)
  {
    _data = v;
    return *this;
  }
  NodeVector& operator*=(double d)
  {
    _data *= d;
    return *this;
  }
  void set_size(const armaicvec& n)
  {
    //    std::cerr << "NodeVector::set_size() n = " << n.t();
    _n = n;
    set_ofs();
    _data.set_size(arma::prod(n));
  }
  void set_size(const NodeVector& u)
  {
    set_size(u.n());
  }
  void fill_bdry(double d=0);
  int dim() const {return (int) _n.n_elem;}
  const armaicvec& n() const {return _n;}
  const armaicvec& ofs() const {return _ofs;}
  armavec& data()
  {
    return _data;
  }
  const armavec& data() const
  {
    return _data;
  }
  double& at(int ix, int iy)
  {
    return _data[_ofs[0]*ix+iy];
  }
  const double& at(int ix, int iy) const
  {
    return _data[_ofs[0]*ix+iy];
  }
  double& at(int ix, int iy, int iz)
  {
    return _data[_ofs[0]*ix+_ofs[1]*iy+iz];
  }
  const double& at(int ix, int iy, int iz) const
  {
    return _data[_ofs[0]*ix+_ofs[1]*iy+iz];
  }
  void fill(double d)
  {
    _data.fill(d);
  }
  void add(double d, const VectorInterface& v)
  {
//    std::cerr << "??? " << data().n_elem<<"\n";
//    std::cerr << "??? " << v.data().n_elem<<"\n";
    _data += d*v.data();
    //    this->arma() += d*v.arma();
  }
  void equal(const VectorInterface& v)
  {
    _data = v.data();
  }  double dot(const VectorInterface& v) const
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
std::ostream& operator<<(std::ostream& os, const NodeVector& v);

#endif /* NodeVector_h */
