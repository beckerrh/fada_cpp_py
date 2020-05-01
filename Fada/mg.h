#ifndef __mg_h
#define __mg_h

#include  <math.h>
//#include  "info.h"

/**************************************************************************/
template<class MATRIX, class VECTOR>
inline int mg_solve(MATRIX& A, VECTOR& u, VECTOR& f, VECTOR& d,
                    int maxiter, double tol_rel, double tol_abs)
{
  int it, reached=0;
  int maxlevel = A.levels()-1;
  double res = sqrt(f(maxlevel)*f(maxlevel));
  double tol = fmax(tol_abs, tol_rel*res);
  printf("%3d %8.2f\n", 0, res);
//  reached = info.check_residual(0,res);
  for (it=1;!maxiter+1;it++)
  {
    res = mgstep(maxlevel, A, u, f, d, tol);
    reached = (res <= tol);
    printf("%3d %10.3e\n", it, res);
    if(reached) return it;
  }
  return -1;
}


/**************************************************************************/
template<class MATRIX, class VECTOR>
inline void smooth(int l, MATRIX& A, VECTOR& u, VECTOR& f, VECTOR& d, bool compute_residual)
{
  if(compute_residual) A.residual(l, d, u, f);
  A.smooth_pre(l, d);
  u(l).add(1., d(l));
}

/**************************************************************************/
template<class MATRIX, class VECTOR>
inline double mgstep(int l, MATRIX& A,VECTOR& u,VECTOR& f,VECTOR& d, double tol)
{
  double res=-1;
  if(l==A.minlevel())
  {
    smooth(l, A, u, f, d, true);
    if(l==A.levels()-1)
    {
      res = sqrt(d(l)*d(l));
      if(res<=tol) return res;
    }
  }
  else
  {
    if(l==A.levels()-1) smooth(l, A, u, f, d, true);
    else smooth(l, A, u, f, d, false);
    A.residual(l, d, u, f);
    if(l==A.levels()-1)
    {
      res = sqrt(d(l)*d(l));
      if(res<=tol) return res;
    }
    d(l-1) = 0.0;
    A.restrict(l-1, d);
    
    f(l-1) = d(l-1);
    u(l-1) = 0.0;
    mgstep(l-1, A, u, f, d, tol);
    
    A.prolongate(l,u);
    
    smooth(l, A, u, f, d, true);
    if(l==A.levels()-1)
    {
      if(res<=tol) return res;
    }
  }
  return res;
}

#endif
