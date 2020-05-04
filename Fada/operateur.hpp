//
//  operateur.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef __operateur_h
#define __operateur_h

#include  "typedefs.hpp"
#include  "vector.hpp"
#include  "vecteurmg.hpp"

/*-------------------------------------------------*/
class Operateur
{
protected:
  armairvec _nall;
  armaimat  _n;
  Array<VecteurMG> omgmem;
  
  void restrict(int l, VecteurMG& out, const VecteurMG& in) const;
  void prolongate(int l, VecteurMG& out, const VecteurMG& in) const;
  void residual(int l, VecteurMG& r, const VecteurMG& u, const VecteurMG& f) const;

  void jacobi       (vector& out, const vector& in) const;
  void gauss_seidel1(vector& out, const vector& in) const;
  void gauss_seidel2(vector& out, const vector& in) const;
  void smooth       (vector& out, const vector& in) const;
  void mgstep(int l, VecteurMG& u, VecteurMG& f, VecteurMG& d, VecteurMG& w, VecteurMG& Aw, double tol);
  void update(int l, VecteurMG& u, VecteurMG& f, VecteurMG& d, VecteurMG& w, VecteurMG& Aw, bool print=false) const;
  int solve();

public:
  std::string smoother;
  int maxiter;
  double tol_rel, tol_abs;
  
  Operateur();
  Operateur(int nlevels, const armaicvec& n0);
  void set_size(int nlevels, const armaicvec& n0);

//  void set_size();
  int minlevel() const   { return 0;}
  int maxlevel() const   { return nlevels()-1;}
  int nlevels()  const   { return _n.n_cols;}
  arma::subview_col<int>  n(int l) const { return _n.col(l);}
  arma::subview_col<int>  n() const { return _n.col(nlevels()-1);}
  int nall()     const   { return _nall(nlevels()-1);}

  void set_size(VecteurMG& v) const;

  void dot(vector& out, const vector& in, double d) const;
  int testsolve();
  int solve(vector& out, const vector& in);
  void boundary(vector& v, const vector& w, double d=1.0) const
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
   
   void boundary(vector& v) const
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

  void right(vector& v) const
  {
    int dim = v.dim();
    double d = 1.0/arma::prod(v.n());
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
};

#endif
