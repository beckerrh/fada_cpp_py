#ifndef __vecteur_h
#define __vecteur_h

#include  "array.h"
#include  <stdio.h>
#include  <armadillo>

typedef arma::Col<double> armavec;

/**************************************************/

class Vecteur
{
private:
//  Array<double>  _val;
  armavec  _val;
  int      _nx, _ny;
  
//  Array<double>& val() {return _val;}
  armavec& val() {return _val;}

public:
  const int& nx() const {return _nx;}
  const int& ny() const {return _ny;}
//  const Array<double>& val() const {return _val;}
  const armavec& val() const {return _val;}

  Vecteur() : _nx(0), _ny(0), _val() {}
  Vecteur(int nx, int ny) : _nx(nx), _ny(ny), _val(nx*ny)
  {
    _val.fill(0);
//    (*this).reinit(nx, ny);
  }
  
  Vecteur& operator=(const Vecteur& a)
  {
    _nx = a.nx(); _ny = a.ny();
    _val = a.val();
    return *this;
  }
  
  Vecteur& operator=(double a)
  {
    _val.fill(a);
//    val() = a;
    return *this;
  }
  
  void reinit(int nx, int ny)
  {
    _nx = nx; _ny = ny;
    _val.set_size(nx*ny);
    _val.fill(0);
//    val().reinit(nx*ny);
  }
  void reinit(const Vecteur& u)
  {
    reinit(u.nx(),u.ny());
  }
  
  double& operator()(int i, int j)
  {
    return _val[_ny*i+j];
    return _val(_ny*i+j);
  }
  const double& operator()(int i, int j) const
  {
    return _val[_ny*i+j];
    return _val(_ny*i+j);
  }
  const double& operator()(int i) const
  {
    return _val[i];
    return _val(i);
  }
  
  double operator* (const Vecteur& v) const;
//  {
////    return arma::dot(_val, v.val());
//    double d = 0.;
//    for(int i=0;i<_nx*_ny;i++)
//    {
//      d +=  val()(i)*v.val()(i);
//    }
//    return d;
//  }
  
  void add(double d, const Vecteur& v)
  {
    _val += d * v.val();
//    for(int i=0;i<_nx*_ny;i++)
//    {
//      val()(i) += d * v.val()(i);
//    }
  }
  
  void boundary(const Vecteur& v, double d=1.0)
  {
    for(int i=0;i<_nx;i++)
    {
      (*this)(i,0)     = d*v(i,0);
      (*this)(i,_ny-1) = d*v(i,_ny-1);
    }
    for(int j=0;j<_ny;j++)
    {
      (*this)(0,j)     = d*v(0,j);
      (*this)(_nx-1,j) = d*v(_nx-1,j);
    }
  }
  
  void boundary()
  {
    for(int i=0;i<_nx;i++)
    {
      (*this)(i,0)     = 0;
      (*this)(i,_ny-1) = 0;
    }
    for(int j=0;j<_ny;j++)
    {
      (*this)(0,j)     = 0;
      (*this)(_nx-1,j) = 0;
    }
  }
  
  void right()
  {
    for(int i=0;i<_nx;i++)
    {
      for(int j=0;j<_ny;j++)
      {
        (*this)(i,j) = 1.;
      }
    }
  }
  
};

#endif
