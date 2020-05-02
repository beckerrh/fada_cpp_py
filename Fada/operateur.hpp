//
//  operateur.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef __operateur_h
#define __operateur_h

#include  <armadillo>
#include  "typedefs.hpp"
#include  "vecteur.hpp"
#include  "vecteurmg.hpp"

/*-------------------------------------------------*/
class Operateur
{
private:
  int       _smoother;
  armairvec _nall;
  armaimat  _n;
  Array<VecteurMG> omgmem;
  
  void restrict(int l, VecteurMG& out, const VecteurMG& in) const;
  void prolongate(int l, VecteurMG& out, const VecteurMG& in) const;
  void residual(int l, VecteurMG& r, const VecteurMG& u, const VecteurMG& f) const;

  void jacobi       (Vecteur& out, const Vecteur& in) const;
  void gauss_seidel1(Vecteur& out, const Vecteur& in) const;
  void gauss_seidel2(Vecteur& out, const Vecteur& in) const;
  void smooth       (Vecteur& out, const Vecteur& in) const;
  void mgstep(int l, VecteurMG& u, VecteurMG& f, VecteurMG& d, VecteurMG& w, VecteurMG& Aw, double tol);
  void update(int l, VecteurMG& u, VecteurMG& f, VecteurMG& d, VecteurMG& w, VecteurMG& Aw) const;

public:
  Operateur(int smoother, int nlevels, const armaicvec& n0);

  void set_size();
  int minlevel() const   { return 0;}
  int maxlevel() const   { return nlevels()-1;}
  int nlevels()  const   { return _n.n_cols;}
  const armaicvec&  n(int l) const { return _n.col(l);}
  const armaicvec&  n() const { return _n.col(nlevels()-1);}
  int nall()     const   { return _nall(nlevels()-1);}

  int smoother() const   { return _smoother;}

  void set_size(VecteurMG& v) const;

  void dot(Vecteur& out, const Vecteur& in, double d) const;
  int solve(Vecteur& out, const Vecteur& in, int maxiter, double tol_rel, double tol_abs);

  void boundary(Vecteur& v, const Vecteur& w, double d=1.0) const
  {
    int nx = v.nx(), ny = v.ny();
    for(int i=0;i<nx;i++)
    {
      v(i,0)     = d*w(i,0);
      v(i,ny-1) = d*w(i,ny-1);
    }
    for(int j=0;j<ny;j++)
    {
      v(0,j)     = d*w(0,j);
      v(nx-1,j) = d*w(nx-1,j);
    }
  }
  
  void boundary(Vecteur& v) const
  {
    int nx = v.nx(), ny = v.ny();
    for(int i=0;i<nx;i++)
    {
      v(i,0)     = 0;
      v(i,ny-1) = 0;
    }
    for(int j=0;j<ny;j++)
    {
      v(0,j)     = 0;
      v(nx-1,j) = 0;
    }
  }
  
  void right(Vecteur& v) const
  {
    int nx = v.nx(), ny = v.ny();
    double d = 1.0/nx/ny;
    for(int i=0;i<nx;i++)
    {
      for(int j=0;j<ny;j++)
      {
        v(i,j) = d;
      }
    }
  }
};

#endif
