#include <math.h>
#include <stdio.h>
#include "operateur.h"
#include "mg.h"


/**************************************************/
Operateur::Operateur(int s, int l, int n, int m)
{
  omgmem.reinit(2);
  ocgmem.reinit(2);
  o_smoother = s;
  o_n0 = n;
  o_m0 = m;
  o_levels = l;
  o_n.reinit(l);
  o_m.reinit(l);
  
  for(int i=0;i<l;i++)
  {
    o_n(i) = int(pow(2,i))*(n-1)+1;
    o_m(i) = int(pow(2,i))*(m-1)+1;
  }
}


/**************************************************/
void Operateur::reinit()
{
  reinit(umg);
  
  for(int i=0;i<ocgmem.n();i++)
  {
    ocgmem(i).reinit(umg(levels()-1));
  }
  for(int i=0;i<omgmem.n();i++)
  {
    omgmem(i).reinit(umg);
  }
}

/**************************************************/

void Operateur::reinit(VecteurMG& v) const
{
  int lev = levels();
  v.val().reinit(lev);
  for(int l=0;l<lev;l++)
  {
    v(l).reinit(n(l),m(l));
  }
}

/************************************************/
void Operateur::solve(Vecteur& out, const Vecteur& in, Info& info)
{
  umg(levels()-1).equ(1., in);
  mg_solve(*this, umg, omgmem, info);
  out.equ(1.,umg(levels()-1));
}

/**************************************************/
void Operateur::smooth(int l, VecteurMG& out)
{
  if(smoother()==0)
  {
    double omega = 0.1;
    jacobi(l,out,omega);
  }
  else if(smoother()==1)
  {
    gauss_seidel_pre(l,out);
  }
  else if(smoother()==2)
  {
    gauss_seidel_pre (l,out);
    gauss_seidel_post(l,out);
  }
}

/**************************************************/
void Operateur::smooth_pre(int l, VecteurMG& out)
{
  if(smoother()==0)
  {
    double omega = 0.1;
    jacobi(l,out,omega);
  }
  else if(smoother()==1)
  {
    gauss_seidel_pre(l,out);
  }
  else if(smoother()==2)
  {
    gauss_seidel_pre(l,out);
  }
}

/**************************************************/
void Operateur::smooth_post(int l, VecteurMG& out)
{
  if(smoother()==0)
  {
    double omega = 0.05;
    jacobi(l,out,omega);
  }
  else if(smoother()==1)
  {
    gauss_seidel_pre(l,out);
  }
  else if(smoother()==2)
  {
    gauss_seidel_post(l,out);
  }
}

/**************************************************/
void Operateur::residual(int l, VecteurMG& r, VecteurMG& u, VecteurMG& f)
{
  vmult(r(l),u(l));
  r(l).sadd(-1.,1.,f(l));
}
