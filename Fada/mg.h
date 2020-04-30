#ifndef __mg_h
#define __mg_h

#include  <math.h>
#include  "info.h"

/**************************************************************************/
template<class MATRIX, class VECTOR>
inline int mg_solve(MATRIX& A, VECTOR& u, VECTOR& f, VECTOR& d, Info& info)
{
  int it, reached=0;
  double res;
  int maxlevel = A.levels()-1;
  res = sqrt(f(maxlevel)*f(maxlevel));
  reached = info.check_residual(0,res);
  for (it=1;!reached;it++)
  {
    res = mgstep2(maxlevel, A, u, f, d);
    reached = info.check_residual(it,res);
    if (reached) continue;
  }
  if (reached<0) return 1;
  return 0;
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
inline double mgstep2(int l, MATRIX& A,VECTOR& u,VECTOR& f,VECTOR& d)
{
  double res=1.0;
  if(l==A.minlevel())
  {
    smooth(l, A, u, f, d, true);
    if(l==A.levels()-1)
    {
      res = sqrt(d(l)*d(l));
    }
  }
  else
  {
    if(l==A.levels()-1) smooth(l, A, u, f, d, true);
    else smooth(l, A, u, f, d, false);
    A.residual(l  , d, u, f);
    if(l==A.levels()-1)
    {
      res = sqrt(d(l)*d(l));
    }
    d(l-1) = 0.0;
    A.restrict(l-1, d);
    
    f(l-1) = d(l-1);
    u(l-1) = 0.0;
    mgstep2(l-1, A, u, f, d);
    
    A.prolongate(l,u);
    
    smooth(l, A, u, f, d, true);
    if(l==A.levels()-1)
    {
      return res;
    }
  }
  return res;
}

#endif
