#ifndef __vecteur_h
#define __vecteur_h

#include "array.h"
#include <stdio.h>

/**************************************************/

class Vecteur
{
private:
  Array<double>  _val;
  int            _nx, _ny;
  
  Array<double>& val() {return _val;}
  
public:
  const int& nx() const {return _nx;}
  const int& ny() const {return _ny;}
  const Array<double>& val() const {return _val;}

  Vecteur() : _nx(0), _ny(0) {}
  Vecteur(int nx, int ny)
  {
    (*this).reinit(nx, ny);
  }
  
  Vecteur(const Vecteur& a)
  {
    (*this).reinit(a.nx(),a.ny());
    val() = a.val();
  }
  Vecteur& operator=(const Vecteur& a)
  {
    val() = a.val();
    return *this;
  }
  
  Vecteur& operator=(double a)
  {
    val() = a;
    return *this;
  }
  
  void reinit(int nx, int ny)
  {
    _nx = nx; _ny = ny;
    val().reinit(nx*ny);
  }
  void reinit(const Vecteur& u)
  {
    reinit(u.nx(),u.ny());
  }
  
  double& operator()(int i, int j)
  {
    return _val(_ny*i+j);
  }
  const double& operator()(int i, int j) const
  {
    return _val(_ny*i+j);
  }
  const double& operator()(int i) const { return val()(i); }
  
  double operator* (const Vecteur& v) const
  {
    double d = 0.;
    for(int i=0;i<_nx*_ny;i++)
    {
      d +=  val()(i)*v.val()(i);
    }
    return d;
  }
  
  void equ(double d, const Vecteur& v)
  {
    for(int i=0;i<_nx*_ny;i++)
    {
      val()(i) = d * v.val()(i);
    }
  }
  
  void add(double d, const Vecteur& v)
  {
    for(int i=0;i<_nx*_ny;i++)
    {
      val()(i) += d * v.val()(i);
    }
  }
  
  void sadd(double e, double d, const Vecteur& v)
  {
    for(int i=0;i<_nx*_ny;i++)
    {
      val()(i) = e * val()(i)  +  d * v.val()(i);
    }
  }
  
  void boundary(const Vecteur& v)
  {
    for(int i=0;i<_nx;i++)
    {
      (*this)(i,0)     = v(i,0);
      (*this)(i,_ny-1) = v(i,_ny-1);
    }
    for(int j=0;j<_ny;j++)
    {
      (*this)(0,j)     = v(0,j);
      (*this)(_nx-1,j) = v(_nx-1,j);
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
