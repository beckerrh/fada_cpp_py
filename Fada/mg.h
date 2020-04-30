#ifndef __mg_h
#define __mg_h

#include  <math.h>
#include  "info.h"

/**************************************************************************/
template<class MATRIX, class VECTOR, class MEM>
inline int mg_solve(MATRIX& A, VECTOR& in, MEM& mem, Info& info)
{
  VECTOR& d = mem(0);
  VECTOR& u = mem(1);
  
  int it, reached=0;
  double res;
  int maxlevel = A.levels()-1;
  res = sqrt(in(maxlevel)*in(maxlevel));
  reached = info.check_residual(0,res);
  for (it=1;!reached;it++)
  {
    res = mgstep2(maxlevel, A, u, in, d);
    reached = info.check_residual(it,res);
    if (reached) continue;
  }
  in(maxlevel) = u(maxlevel);
  if (reached<0) return 1;
  return 0;
}


/**************************************************************************/
template<class MATRIX, class VECTOR>
inline void smooth(int l, MATRIX& A,VECTOR& u,VECTOR& f,VECTOR& d)
{
  A.residual(l, d, u, f);
  A.smooth_pre(l, d);
  u(l).add(1., d(l));
}

/**************************************************************************/
template<class MATRIX, class VECTOR>
inline double mgstep2(int l, MATRIX& A,VECTOR& u,VECTOR& f,VECTOR& d)
{
  smooth(l, A, u, f, d);
  double res;
  if(l>A.minlevel())
    {
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
      
      smooth(l, A, u, f, d);
      if(l==A.levels()-1)
      {
        return res;
      }
    }
  return res;
}

#endif
