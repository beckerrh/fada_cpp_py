#ifndef __vecteur_h
#define __vecteur_h

#include  <stdio.h>
#include  <armadillo>
#include  "typedefs.hpp"

typedef arma::Col<double> armavec;

/*-------------------------------------------------*/

class Vecteur
{
private:
  armavec  _val;
  armaicvec _n;

public:
  const int& nx() const {return _n[0];}
  const int& ny() const {return _n[1];}
  const armaicvec& n() const {return _n;}
  const armavec& val() const {return _val;}

  Vecteur() : _val(), _n() {}
  Vecteur(const Vecteur& v) : _val(v.val()), _n(v._n) {assert(0);}
  Vecteur& operator=(const Vecteur& a)
  {
    _val = a.val();
    return *this;
  }
  Vecteur& operator=(double a)
  {
    _val.fill(a);
    return *this;
  }
  void set_size(const armaicvec& n)
  {
    _n = n;
    _val.set_size(arma::prod(n));
  }
  void set_size(const Vecteur& u)
  {
    set_size(u.n());
  }
  double& operator()(int i, int j);
  const double& operator()(int i, int j) const;  
  double operator*(const Vecteur& v) const;
  double dot(const Vecteur& v) const {return (*this)*v;}
  void add(double d, const Vecteur& v);
  
  void output(const std::string& filename) const
  {
    arma::hdf5_name spec(filename);
    _val.save(spec);
  }  
};

/*-------------------------------------------------*/
inline double& Vecteur::operator()(int i, int j)
{
  return _val[_n[1]*i+j];
}
inline const double& Vecteur::operator()(int i, int j) const
{
  return _val[_n[1]*i+j];
}
/*-------------------------------------------------*/
inline double Vecteur::operator* (const Vecteur& v) const
{
  return arma::dot(_val,v.val());
}
inline void Vecteur::add(double d, const Vecteur& v)
{
  _val += d * v.val();
}

/*-------------------------------------------------*/

#endif
