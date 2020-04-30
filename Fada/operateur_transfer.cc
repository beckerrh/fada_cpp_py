#include <math.h>
#include "operateur.h"

/**************************************************/
void Operateur::restrict(int l, VecteurMG& v)    //   V_{l+1} --> V_l
{
  int n = v(l).n(), m = v(l).m();
  for(int i=1;i<n-1;i++)
  {
    for(int j=1;j<m-1;j++)
    {
      v(l)(i,j) = v(l+1)(2*i  ,2*j  )
      +  0.5 *( v(l+1)(2*i-1,2*j  ) + v(l+1)(2*i+1,2*j  )
               + v(l+1)(2*i  ,2*j-1) + v(l+1)(2*i  ,2*j+1) )
      +  0.25*( v(l+1)(2*i-1,2*j-1) + v(l+1)(2*i-1,2*j+1)
               + v(l+1)(2*i+1,2*j-1) + v(l+1)(2*i+1,2*j+1) );
    }
  }
  v(l).boundary();
}

/**************************************************/
void Operateur::prolongate(int l, VecteurMG& v)   //   V_{l-1} --> V_l
{
  int n = v(l-1).n(), m = v(l-1).m();
  for(int i=1;i<n-1;i++)
  {
    for(int j=1;j<m-1;j++)
    {
      v(l)(2*i  ,2*j  ) +=        v(l-1)(i,j);
      v(l)(2*i-1,2*j  ) += 0.5  * v(l-1)(i,j);
      v(l)(2*i+1,2*j  ) += 0.5  * v(l-1)(i,j);
      v(l)(2*i  ,2*j-1) += 0.5  * v(l-1)(i,j);
      v(l)(2*i  ,2*j+1) += 0.5  * v(l-1)(i,j);
      v(l)(2*i-1,2*j-1) += 0.25 * v(l-1)(i,j);
      v(l)(2*i-1,2*j+1) += 0.25 * v(l-1)(i,j);
      v(l)(2*i+1,2*j-1) += 0.25 * v(l-1)(i,j);
      v(l)(2*i+1,2*j+1) += 0.25 * v(l-1)(i,j);
    }
  }
}
