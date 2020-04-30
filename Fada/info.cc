#include <math.h>
#include <stdio.h>
#include "info.h"

const char* Info::std_txt = "It:";

#define MAX(a,b) ( ((a)>(b)) ? (a) : (b) )

Info::Info(double f, double t, int p, int m, const char* txt)
  : text(txt)
{
  reduction() = f;
  globaltol() = t;
  printstep() = p;
  maxiter  () = m;
  residual () = 0.;
  iteration() = 0;
  reduction_rate() = 1.;
  last_reduction_rate() = 1.;
  lastit = -1;
}

void Info::reset()
{
  iteration() = 0;
  putada()    = 0;
}

double Info::compute_reduction_rate() const
{ 
  return pow(residual()/first_res,1./MAX(1,iteration()));
}

int Info::check_residual(int i, double res)
{
  if (!i)
    {
      first_res = res;
      last_res = res;
      tol() = res*reduction();
    }
  residual () = res;
  iteration() = i;
  if(i)
    {
      last_reduction_rate() = res/last_res;
      last_res = res;
      reduction_rate() = compute_reduction_rate();
    }

  ierg = 0;
  if ( (res<tol()) || (res<globaltol()) )    
    {
      ierg =  1;
    }
  else if ((i>=maxiter()) || (res>1.e+8)) 
    {
      putada() = 1;
      ierg = -1;
    }
  
  lastit = i;
  if (printstep() && !(i%printstep()) )
    printf("%s %4d  Res: %8.2e\n",text,i,res);
  return ierg;
}
