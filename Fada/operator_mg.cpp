//
//  operateur_mg.cpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include "operator.hpp"

/*-------------------------------------------------*/
int Operator::solve(bool print)
{
  VectorMG& u = _mgmem(0);
  VectorMG& f = _mgmem(1);
  VectorMG& d = _mgmem(2);
  VectorMG& w = _mgmem(3);
  VectorMG& Aw = _mgmem(4);

  int maxlevel = nlevels()-1;
  double res, tol;
//  std::cerr << "smoother " << smoother << std::endl;
//  std::cerr << "maxlevel " << maxlevel << std::endl;
  for(int iter=0; iter<this->maxiter+1; iter++)
  {
//    std::cerr <<arma::norm(d(maxlevel).arma())<<" "<<arma::norm(f(maxlevel).arma())<<" "<<arma::norm(u(maxlevel).arma()) << "\n";
    residual(maxlevel, d, u, f);
//    std::cerr <<arma::norm(d(maxlevel).arma())<<" "<<arma::norm(f(maxlevel).arma())<<" "<<arma::norm(u(maxlevel).arma()) << "\n";

    res = arma::norm(d(maxlevel).arma());
    if(iter==0)
    {
      tol = fmax(this->tol_abs, this->tol_rel*res);
      if(print) printf("--- %10.3e\n", tol);
    }
    if(print) printf("%3d %10.3e\n", iter, res);
    if(res <= tol)
    {
      return iter;
    }
    mgstep(maxlevel, u, f, d, w, Aw, tol);
  }
  return -1;
}

/*-------------------------------------------------*/
void Operator::mgstep(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, VectorMG& Aw, double tol)
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
    u(l-1).fill(0.0);
    mgstep(l-1, u, f, d, w, Aw, tol);
    
    prolongate(l, w, u);
    residual(l, d, u, f);
    update(l, u, f, d, w, Aw);
    smooth(w(l), d(l));
    update(l, u, f, d, w, Aw);
  }
}
/*-------------------------------------------------*/
void Operator::update(int l, VectorMG& u, VectorMG& f, VectorMG& r, VectorMG& w, VectorMG& Aw, bool print) const
{
  /*
   input:
      d contains the residual (d = f - Au)
   output:
      d contgains the new residual
   */
  Aw(l).fill(0.0);
  dot(Aw(l), w(l), 1.0);
  double omega = 1.0;
  if(optmem>=0)
  {
    double d1 = w(l).dot(r(l));
    double d2 = w(l).dot(Aw(l));
//    double d1 = Aw(l).dot(r(l));
//    double d2 = Aw(l).dot(Aw(l));
    omega = d1/d2;
    omega = fmax(fmin(omega, 10.0),0.1);
//    printf("l=%2d '%4.2f' %10.3e\n",l, omega, arma::norm(w(l).arma()));
  }
  if(print) printf("l=%2d '%4.2f' %10.3e\n",l, omega, arma::norm(w(l).arma()));
  u(l).add(omega, w(l));
  r(l).add(-omega, Aw(l));
}
