#include  <math.h>
#include  <stdio.h>
#include  <armadillo>
#include  <algorithm>
#include  "operator.hpp"
#include  "updater.hpp"

/*-------------------------------------------------*/
void Operator::set_parameters()
{
  maxiter = 100;
  smoother = "jac";
  tol_rel = 1e-8;
  tol_abs = 1e-13;
  optmem = 0;
}

Operator::~Operator()
{
  for(int l=0;l<nlevels();l++)
  {
    delete _mgupdate(l);
  }
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
  _mgupdate.set_size(nlevels);
  std::string type="cyc";
//  type="coef";
//  type="ortho";
//  type="restart";
  for(int l=0;l<nlevels;l++)
  {
    if(optmem==0) _mgupdate(l) = new UpdaterSimple;
    else _mgupdate(l) = new Updater;
  }
  for(int l=0;l<nlevels;l++)
  {
    _mgupdate(l)->setParameters(l, this, std::max(0,optmem), type);
    _mgupdate(l)->set_size(n(l)+2);
  }
}

/*-------------------------------------------------*/
void Operator::set_size(VectorMG& v) const
{
  v.set_size(nlevels());
  for(int l=0;l<nlevels();l++)
  {
//    v(l).set_size(n(l));
    v(l).set_size(n(l)+2);
  }
}

/*-------------------------------------------------*/
void Operator::smooth(int l, Vector& out, const Vector& in) const
{
  out.fill(0.0);
  if(this->smoother=="jac")
  {
    jacobi(l, out, in);
  }
  else if(this->smoother=="gs1")
  {
    gauss_seidel1(l, out, in);
  }
  else if(this->smoother=="gs2")
  {
    gauss_seidel2(l, out, in);
  }
}

/*-------------------------------------------------*/
void Operator::residual(int l, VectorMG& r, const VectorMG& u, const VectorMG& f) const
{
  r(l) =  f(l);
  dot(l, r(l),u(l), -1.0);
}

/*-------------------------------------------------*/
void Operator::vectormg2vector(int l, Vector& u, const VectorMG& umg) const
{
//  u = umg(l);
    if(dim()==2)
    {
      int nx = _n(0,l), ny = _n(1,l);
      for(int ix=0;ix<nx;ix++)
      {
        for(int iy=0;iy<ny;iy++)
        {
          u.at(ix,iy) = umg(l).atp(ix,iy);
        }
      }
    }
  else if(dim()==3)
  {
    int nx = _n(0,l), ny = _n(1,l);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<ny;iz++)
        {
          u.at(ix,iy,iz) = umg(l).atp(ix,iy,iz);
        }
      }
    }
  }
}
/*-------------------------------------------------*/
void Operator::vector2vectormg(int l, VectorMG& umg, const Vector& u) const
{
//  umg(l) = u;
    if(dim()==2)
    {
      int nx = _n(0,l), ny = _n(1,l);
      for(int ix=0;ix<nx;ix++)
      {
        for(int iy=0;iy<ny;iy++)
        {
          umg(l).atp(ix,iy) = u.at(ix,iy);
        }
      }
    }
  else if(dim()==3)
  {
    int nx = _n(0,l), ny = _n(1,l);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<ny;iz++)
        {
          umg(l).atp(ix,iy,iz) = u.at(ix,iy,iz);
        }
      }
    }
  }
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

  vector2vectormg(nlevels()-1, fmg, _f);
  vector2vectormg(nlevels()-1, umg, _u);
//  fmg(nlevels()-1) = _f;
//  umg(nlevels()-1) = _u;
  int iter = solve(print);
  vectormg2vector(nlevels()-1, _u, umg);
//  _u = umg(nlevels()-1);
  return iter;
}
/*-------------------------------------------------*/
void Operator::boundary(Vector& v, const Vector& w, double d) const
 {
   if(dim()==2)
   {
     int nx = _n(0,_n.n_cols-1), ny = _n(1,_n.n_cols-1);
     for(int ix=0;ix<nx;ix++)
     {
       v.at(ix,0)    = d*w.at(ix,0);
       v.at(ix,ny-1) = d*w.at(ix,ny-1);
     }
     for(int iy=0;iy<ny;iy++)
     {
       v.at(0,iy)    = d*w.at(0,iy);
       v.at(nx-1,iy) = d*w.at(nx-1,iy);
     }
   }
   else if(dim()==3)
   {
     int nx = _n(0,_n.n_cols-1), ny = _n(1,_n.n_cols-1), nz = _n(2,_n.n_cols-1);
     for(int ix=0;ix<nx;ix++)
     {
       for(int iy=0;iy<ny;iy++)
       {
         v.at(ix,iy,0)    = d*w.at(ix,iy,0);
         v.at(ix,iy,nz-1) = d*w.at(ix,iy,nz-1);
       }
     }
     for(int ix=0;ix<nx;ix++)
     {
       for(int iz=0;iz<nz;iz++)
       {
         v.at(ix,0,   iz) = d*w.at(ix,0,   iz);
         v.at(ix,ny-1,iz) = d*w.at(ix,ny-1,iz);
       }
     }
     for(int iy=0;iy<ny;iy++)
     {
       for(int iz=0;iz<nz;iz++)
       {
         v.at(0,iy,iz)    = d*w.at(0,iy,iz);
         v.at(nx-1,iy,iz) = d*w.at(nx-1,iy,iz);
       }
     }
   }
 }
 
 /*-------------------------------------------------*/
 void Operator::boundary(Vector& v) const
 {
   if(dim()==2)
   {
     int nx = _n(0,_n.n_cols-1), ny = _n(1,_n.n_cols-1);
     for(int ix=0;ix<nx;ix++)
     {
       v.at(ix,0)    = 0.0;
       v.at(ix,ny-1) = 0.0;
     }
     for(int iy=0;iy<ny;iy++)
     {
       v.at(0,iy)    = 0.0;
       v.at(nx-1,iy) = 0.0;
     }
   }
   else if(dim()==3)
   {
     int nx = _n(0,_n.n_cols-1), ny = _n(1,_n.n_cols-1), nz = _n(2,_n.n_cols-1);
     for(int ix=0;ix<nx;ix++)
     {
       for(int iy=0;iy<ny;iy++)
       {
         v.at(ix,iy,0)    = 0.0;
         v.at(ix,iy,nz-1) = 0.0;
       }
     }
     for(int ix=0;ix<nx;ix++)
     {
       for(int iz=0;iz<nz;iz++)
       {
         v.at(ix,0,   iz) = 0.0;
         v.at(ix,ny-1,iz) = 0.0;
       }
     }
     for(int iy=0;iy<ny;iy++)
     {
       for(int iz=0;iz<nz;iz++)
       {
         v.at(0,   iy,iz) = 0.0;
         v.at(nx-1,iy,iz) = 0.0;
       }
     }
   }
 }

/*-------------------------------------------------*/
void Operator::right(Vector& v) const
{
  double d = (4.0*dim()+pow(2.0,dim()))/arma::prod(n());
  if(dim()==2)
  {
    int nx = _n(0,_n.n_cols-1), ny = _n(1,_n.n_cols-1);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        v.at(ix,iy) = d;
      }
    }
  }
  else if(dim()==3)
  {
    int nx = _n(0,_n.n_cols-1), ny = _n(1,_n.n_cols-1), nz = _n(2,_n.n_cols-1);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<nz;iz++)
        {
          v.at(ix,iy, iz) = d;
        }
      }
    }
  }
  else
  {
    std::cerr << "wrong dimension " << dim() << std::endl;
    exit(1);
  }
}
