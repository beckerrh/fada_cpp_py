#include <math.h>
#include "operateur.hpp"

/*-------------------------------------------------*/
void Operateur::dot(Vecteur& out, const Vecteur& in, double d) const
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
  
  // Conditions aux limites ( Dirichlet)
  
  //  out.boundary(in, d);
}


/*-------------------------------------------------*/
void Operateur::jacobi(Vecteur& out, const Vecteur& in) const
{
  double omega = 0.1;
  int n = out.nx(), m = out.ny();
  for(int i=1;i<n-1;i++)
  {
    for(int j=1;j<m-1;j++)
    {
      out(i,j) = omega * in(i,j);
    }
  }
}

/*-------------------------------------------------*/
void Operateur::gauss_seidel1(Vecteur& out, const Vecteur& in) const
{
  int n = out.nx(), m = out.ny();
//  out = in;
  for(int i=1;i<n-1;i++)
  {
    for(int j=1;j<m-1;j++)
    {
      out(i,j) = 0.125 * ( in(i  ,j  )+ out(i-1,j-1) + out(i-1,j  )+ out(i-1,j+1) + out(i  ,j-1) );
    }
  }
}

/*-------------------------------------------------*/
void Operateur::gauss_seidel2(Vecteur& out, const Vecteur& in) const
{
  int n = out.nx(), m = out.ny();
  for(int i=n-2;i>=1;i--)
  {
    for(int j=m-2;j>=1;j--)
    {
      out(i,j) = 0.125 * ( in(i  ,j  )
                             + out(i+1,j-1) + out(i+1,j  )
                             + out(i+1,j+1) + out(i  ,j+1) );
    }
  }
}
