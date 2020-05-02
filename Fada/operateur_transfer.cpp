#include <math.h>
#include "operateur.hpp"

/*-------------------------------------------------*/
void Operateur::restrict(int l, VecteurMG& out, const VecteurMG& in) const
//   V_{l+1} --> V_l
{
  int n = out(l).nx(), m = out(l).ny();
  for(int i=1;i<n-1;i++)
  {
    for(int j=1;j<m-1;j++)
    {
      out(l)(i,j) = in(l+1)(2*i  ,2*j  )
      +  0.5 *( in(l+1)(2*i-1,2*j  ) + in(l+1)(2*i+1,2*j  )
               + in(l+1)(2*i  ,2*j-1) + in(l+1)(2*i  ,2*j+1) )
      +  0.25*( in(l+1)(2*i-1,2*j-1) + in(l+1)(2*i-1,2*j+1)
               + in(l+1)(2*i+1,2*j-1) + in(l+1)(2*i+1,2*j+1) );
    }
  }
//  boundary(out(l));
}

/*-------------------------------------------------*/
void Operateur::prolongate(int l, VecteurMG& out, const VecteurMG& in) const
//   V_{l-1} --> V_l
{
  int n = in(l-1).nx(), m = in(l-1).ny();
  out(l) = 0.0;
  for(int i=1;i<n-1;i++)
  {
    for(int j=1;j<m-1;j++)
    {
      out(l)(2*i  ,2*j  ) +=        in(l-1)(i,j);
      out(l)(2*i-1,2*j  ) += 0.5  * in(l-1)(i,j);
      out(l)(2*i+1,2*j  ) += 0.5  * in(l-1)(i,j);
      out(l)(2*i  ,2*j-1) += 0.5  * in(l-1)(i,j);
      out(l)(2*i  ,2*j+1) += 0.5  * in(l-1)(i,j);
      out(l)(2*i-1,2*j-1) += 0.25 * in(l-1)(i,j);
      out(l)(2*i-1,2*j+1) += 0.25 * in(l-1)(i,j);
      out(l)(2*i+1,2*j-1) += 0.25 * in(l-1)(i,j);
      out(l)(2*i+1,2*j+1) += 0.25 * in(l-1)(i,j);
    }
  }
}
