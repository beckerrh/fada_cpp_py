#include <math.h>
#include <stdio.h>
#include "operateur.hpp"


/*-------------------------------------------------*/
Operateur::Operateur(int smoother, int nlevels, const armaicvec& n0) : _smoother(smoother)
{
  omgmem.set_size(5);
  int dim = 2;
  _n.set_size(dim, nlevels);
  
  for(int i=0;i<n0.n_elem;i++)
  {
    for(int l=0;l<nlevels;l++)
    {
      _n(i,l) = int(pow(2,l))*(n0[i]-1)+1;
    }
  }
  _nall = arma::prod(_n, 0);
  for(int i=0;i<omgmem.n();i++)
  {
    set_size(omgmem(i));
  }
}

/*-------------------------------------------------*/
void Operateur::set_size(VecteurMG& v) const
{
  int lev = nlevels();
  v.val().set_size(lev);
  for(int l=0;l<lev;l++)
  {
//    v(l).set_size(nx(l),ny(l));
    v(l).set_size(n(l));
  }
}

/*-------------------------------------------------*/
void Operateur::smooth(Vecteur& out, const Vecteur& in) const
{
  if(smoother()==0)
  {
    jacobi(out, in);
  }
  else if(smoother()==1)
  {
    gauss_seidel1(out, in);
  }
  else if(smoother()==2)
  {
    gauss_seidel2(out, in);
  }
}

/*-------------------------------------------------*/
void Operateur::residual(int l, VecteurMG& r, const VecteurMG& u, const VecteurMG& f) const
{
  r(l) =  f(l);
  dot(r(l),u(l), -1.0);
}
