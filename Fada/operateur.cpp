#include <math.h>
#include <stdio.h>
#include "operateur.hpp"

/*-------------------------------------------------*/
void Operateur::set_parameters()
{
  maxiter = 100;
  smoother = "jac";
  tol_rel = 1e-10;
  tol_abs = 1e-13;
}

/*-------------------------------------------------*/
Operateur::Operateur()
{
  set_parameters();
}

/*-------------------------------------------------*/
Operateur::Operateur(int nlevels, const armaicvec& n0)
{
  set_parameters();
  set_size(nlevels, n0);
}

/*-------------------------------------------------*/
void Operateur::set_size(int nlevels, const armaicvec& n0)
{
  _mgmem.set_size(5);
  int dim = n0.n_elem;
  _n.set_size(dim, nlevels);
  
  for(int i=0;i<n0.n_elem;i++)
  {
    for(int l=0;l<nlevels;l++)
    {
      _n(i,l) = int(pow(2,l))*(n0[i]-1)+1;
    }
  }
  _nall = arma::prod(_n, 0);
  for(int i=0;i<_mgmem.n();i++)
  {
    set_size(_mgmem(i));
  }
}

/*-------------------------------------------------*/
void Operateur::set_size(VecteurMG& v) const
{
  int lev = nlevels();
  v.val().set_size(lev);
  for(int l=0;l<lev;l++)
  {
    v(l).set_size(n(l));
  }
}

/*-------------------------------------------------*/
void Operateur::smooth(vector& out, const vector& in) const
{
  if(this->smoother=="jac")
  {
    jacobi(out, in);
  }
  else if(this->smoother=="gs1")
  {
    gauss_seidel1(out, in);
  }
  else if(this->smoother=="gs2")
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


/*-------------------------------------------------*/
int Operateur::solve(vector& out, const vector& in)
{
  VecteurMG& u = _mgmem(0);
  VecteurMG& f = _mgmem(1);

  f(nlevels()-1) = in;
  int iter = solve();
  out = u(nlevels()-1);
  return iter;
}

/*-------------------------------------------------*/
int Operateur::testsolve()
{
  VecteurMG& u = _mgmem(0);
  VecteurMG& f = _mgmem(1);

  right(f(nlevels()-1));
  u(nlevels()-1).fill(0.0);
  boundary(u(nlevels()-1));
  boundary(f(nlevels()-1), u(nlevels()-1));
  return solve();
}
