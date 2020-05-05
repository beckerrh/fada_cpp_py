#include <math.h>
#include <stdio.h>
#include <armadillo>
#include "operator.hpp"

/*-------------------------------------------------*/
void Operator::set_parameters()
{
  maxiter = 100;
  smoother = "jac";
  tol_rel = 1e-10;
  tol_abs = 1e-13;
  optmem = 0;
}

/*-------------------------------------------------*/
Operator::Operator()
{
  set_parameters();
}

/*-------------------------------------------------*/
Operator::Operator(int nlevels, const armaicvec& n0)
{
  set_parameters();
  set_size(nlevels, n0);
}

/*-------------------------------------------------*/
void Operator::set_size(int nlevels, const armaicvec& n0)
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
void Operator::set_size(VectorMG& v) const
{
  int lev = nlevels();
  v.val().set_size(lev);
  for(int l=0;l<lev;l++)
  {
    v(l).set_size(n(l));
  }
}

/*-------------------------------------------------*/
void Operator::smooth(vector& out, const vector& in) const
{
  out.fill(0.0);
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
void Operator::residual(int l, VectorMG& r, const VectorMG& u, const VectorMG& f) const
{
  r(l) =  f(l);
  dot(r(l),u(l), -1.0);
}

/*-------------------------------------------------*/
int Operator::testsolve(bool print)
{
  _u.set_size(n());
  _u.fill(0);
  _f.set_size(_u);
  right(_f);
  boundary(_u);
  boundary(_f, _u);

  VectorMG& umg = _mgmem(0);
  VectorMG& fmg = _mgmem(1);

  fmg(nlevels()-1) = _f;
  umg(nlevels()-1) = _u;
  int iter = solve(print);
  _u = umg(nlevels()-1);
  return iter;
}
/*-------------------------------------------------*/
void Operator::boundary(vector& v, const vector& w, double d) const
 {
   int dim = v.dim();
   if(dim==2)
   {
     int nx = v.n(0), ny = v.n(1);
     for(int ix=0;ix<nx;ix++)
     {
       v(ix,0)    = d*w(ix,0);
       v(ix,ny-1) = d*w(ix,ny-1);
     }
     for(int iy=0;iy<ny;iy++)
     {
       v(0,iy)    = d*w(0,iy);
       v(nx-1,iy) = d*w(nx-1,iy);
     }
   }
   else if(dim==3)
   {
     int nx = v.n(0), ny = v.n(1), nz = v.n(2);
     for(int ix=0;ix<nx;ix++)
     {
       for(int iy=0;iy<ny;iy++)
       {
         v(ix,iy,0)    = d*w(ix,iy,0);
         v(ix,iy,nz-1) = d*w(ix,iy,nz-1);
       }
     }
     for(int ix=0;ix<nx;ix++)
     {
       for(int iz=0;iz<nz;iz++)
       {
         v(ix,0,   iz) = d*w(ix,0,   iz);
         v(ix,ny-1,iz) = d*w(ix,ny-1,iz);
       }
     }
     for(int iy=0;iy<ny;iy++)
     {
       for(int iz=0;iz<nz;iz++)
       {
         v(0,iy,iz)    = d*w(0,iy,iz);
         v(nx-1,iy,iz) = d*w(nx-1,iy,iz);
       }
     }
   }
 }
 
 /*-------------------------------------------------*/
 void Operator::boundary(vector& v) const
 {
   int dim = v.dim();
   if(dim==2)
   {
     int nx = v.n(0), ny = v.n(1);
     for(int ix=0;ix<nx;ix++)
     {
       v(ix,0)    = 0.0;
       v(ix,ny-1) = 0.0;
     }
     for(int iy=0;iy<ny;iy++)
     {
       v(0,iy)    = 0.0;
       v(nx-1,iy) = 0.0;
     }
   }
   else if(dim==3)
   {
     int nx = v.n(0), ny = v.n(1), nz = v.n(2);
     for(int ix=0;ix<nx;ix++)
     {
       for(int iy=0;iy<ny;iy++)
       {
         v(ix,iy,0)    = 0.0;
         v(ix,iy,nz-1) = 0.0;
       }
     }
     for(int ix=0;ix<nx;ix++)
     {
       for(int iz=0;iz<nz;iz++)
       {
         v(ix,0,   iz) = 0.0;
         v(ix,ny-1,iz) = 0.0;
       }
     }
     for(int iy=0;iy<ny;iy++)
     {
       for(int iz=0;iz<nz;iz++)
       {
         v(0,   iy,iz) = 0.0;
         v(nx-1,iy,iz) = 0.0;
       }
     }
   }
 }

/*-------------------------------------------------*/
void Operator::right(vector& v) const
{
  int dim = v.dim();
  double d = (4.0*dim+pow(2.0,dim))/arma::prod(v.n());
  if(dim==2)
  {
    int nx = v.n(0), ny = v.n(1);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        v(ix,iy) = d;
      }
    }
  }
  else if(dim==3)
  {
    int nx = v.n(0), ny = v.n(1), nz = v.n(2);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<nz;iz++)
        {
          v(ix,iy, iz) = d;
        }
      }
    }
  }
  else
  {
    std::cerr << "wrong dimension " << dim << std::endl;
    exit(1);
  }
}
