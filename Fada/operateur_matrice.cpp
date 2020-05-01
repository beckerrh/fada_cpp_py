#include <math.h>
#include "operateur.h"

/**************************************************/

void Operateur::vmult(Vecteur& out, const Vecteur& in, double d) const
{
  // Laplacien   elements finis q1  (9-point-stencil)
  int n = out.nx(), m = out.ny();
    
  for(int i=1;i<n-1;i++)
  {
    for(int j=1;j<m-1;j++)
    {
      out(i,j) += d*(8. * in(i,j) - in(i-1,j )  - in(i+1,j)
      - in(i  ,j-1) - in(i  ,j+1)
      - in(i-1,j-1) - in(i-1,j+1)
      - in(i+1,j-1) - in(i+1,j+1));
    } 
  }
  
  // Conditions aux limites ( Dirichlet
  
  out.boundary(in, d);
}


/**************************************************/

void Operateur::jacobi(int l, VecteurMG& out, double omega)
{
  int n = out(l).nx(), m = out(l).ny();
  for(int i=1;i<n-1;i++)
  {
    for(int j=1;j<m-1;j++)
    {
      out(l)(i,j) = omega * out(l)(i,j);
    }
  }
}

/**************************************************/

void Operateur::gauss_seidel_pre(int l, VecteurMG& out)
{
  int n = out(l).nx(), m = out(l).ny();
  for(int i=1;i<n-1;i++)
  {
    for(int j=1;j<m-1;j++)
    {
      out(l)(i,j) = 0.125 * ( out(l)(i  ,j  ) 
                             + out(l)(i-1,j-1) + out(l)(i-1,j  )
                             + out(l)(i-1,j+1) + out(l)(i  ,j-1) );
    }
  }
}

/**************************************************/

void Operateur::gauss_seidel_post(int l, VecteurMG& out)
{
  int n = out(l).nx(), m = out(l).ny();
  for(int i=n-2;i>=1;i--)
  {
    for(int j=m-2;j>=1;j--)
    {
      out(l)(i,j) = 0.125 * ( out(l)(i  ,j  ) 
                             + out(l)(i+1,j-1) + out(l)(i+1,j  )
                             + out(l)(i+1,j+1) + out(l)(i  ,j+1) );
    }
  }
}
