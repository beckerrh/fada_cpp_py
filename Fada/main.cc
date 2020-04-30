#include  <stdio.h>
#include  <stdlib.h>
#include  <ctime>
#include  "vecteur.h"
#include  "vecteurmg.h"
#include  "operateur.h"
#include  "info.h"

inline double seconds(void)
{
    static const double secs_per_tick = 1.0 / CLOCKS_PER_SEC;
    return ( (double) clock() ) * secs_per_tick;
}

/**************************************************/

int main(int argc, char** argv)
{
  // preblabla(argc,argv);
  int n0 = 2,  m0 = 2;
  int smoother = 1, levels = 12;

  double t0 = seconds();
  /*
   smoother [0:jacobi 1:gauss_seidel 2:gauss_seidel_symm]
   */
  Operateur   A(smoother, levels, n0, m0);

  Vecteur     u,f;

  int         iterations = 100, print = 1;
  double      tol_rel = 1e-10, tol_abs = 1e-15;;
  Info        info(tol_rel,tol_abs,print,iterations,"Fada");

  int n = A.n(), m = A.m();
  u.reinit(n,m);
  f.reinit(u);
  A.reinit();

  u.boundary();
  f.right();
  f.boundary(u);

  A.solve(u, f, info);

//  u.output_plotmtv();
  
  printf("Vous avez utilise le lisseur ");
  if   (A.smoother()==0  ) printf("[Jac]");
  else if(A.smoother()==1) printf("[GS ]");
  else if(A.smoother()==2) printf("[GSS]");
  printf("\n\nNo. Iterations %3d (N = %6d)  %7.2e\n\n",
   info.iteration(),A.n()*A.m(),info.residual());
  printf("Total time: %6.2f\n", seconds()-t0);
}
