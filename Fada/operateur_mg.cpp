//
//  operateur_mg.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include "operateur.hpp"

/*-------------------------------------------------*/
int Operateur::solve(Vecteur& out, const Vecteur& in, int maxiter, double tol_rel, double tol_abs)
{
  VecteurMG& u = omgmem(0);
  VecteurMG& f = omgmem(1);
  VecteurMG& d = omgmem(2);
  VecteurMG& w = omgmem(3);
  VecteurMG& Aw = omgmem(4);
  f(nlevels()-1) = in;
  
  int maxlevel = nlevels()-1;
  double res, tol;
  for(int iter=0; iter<maxiter+1; iter++)
  {
    residual(maxlevel, d, u, f);
    res = sqrt(d(maxlevel)*d(maxlevel));
    if(iter==0)
    {
      tol = fmax(tol_abs, tol_rel*res);
    }
    printf("%3d %10.3e %10.3e\n", iter, res, tol);
    if(res <= tol)
    {
      out = u(nlevels()-1);
      return iter;
    }
    mgstep(maxlevel, u, f, d, w, Aw, tol);
  }
  return -1;
}
/*-------------------------------------------------*/
void Operateur::mgstep(int l, VecteurMG& u, VecteurMG& f, VecteurMG& d, VecteurMG& w, VecteurMG& Aw, double tol)
{
  if(l==minlevel())
  {
    residual(l, d, u, f);
    smooth(w(l), d(l));
    update(l, u, f, d, w, Aw);
  }
  else
  {
    smooth(w(l), d(l));
    update(l, u, f, d, w, Aw);
    
    restrict(l-1, d, d);
    
    f(l-1) = d(l-1);
    u(l-1) = 0.0;
    mgstep(l-1, u, f, d, w, Aw, tol);
    
    prolongate(l, w, u);
    residual(l, d, u, f);
    update(l, u, f, d, w, Aw);
    smooth(w(l), d(l));
    update(l, u, f, d, w, Aw);
  }
}
/*-------------------------------------------------*/
void Operateur::update(int l, VecteurMG& u, VecteurMG& f, VecteurMG& d, VecteurMG& w, VecteurMG& Aw) const
{
  /*
   input:
      d contains the residual (d = f - Au)
   output:
      d contgains the new residual
   */
  Aw(l) = 0.0;
  dot(Aw(l), w(l), 1.0);
  double d1 = w(l).dot(d(l));
  double d2 = w(l).dot(Aw(l));
  double omega = d1/d2;
  omega = fmax(fmin(omega, 2.0),0.5);
//  printf("l=%2d '%4.2f'\n",l, omega);
//  omega = 1.0;
  u(l).add(omega, w(l));
  d(l).add(-omega, Aw(l));
}
