//
//  trapezmatrix.cpp
//  Fada
//
//  Created by Roland Becker on 15/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "trapezmatrix.hpp"
#include  <math.h>
#include  "sparsematrix.hpp"

TrapezMatrix2d::~TrapezMatrix2d() {}
TrapezMatrix3d::~TrapezMatrix3d() {}

/*-------------------------------------------------*/
void TrapezMatrix2d::set_grid(const armaicvec& n, const armavec& dx)
{
  assert(n.n_elem==2);
  _nx = n[0];
  _ny = n[1];
  _vol = arma::prod(dx);
  _dx = dx[0];
  _dy = dx[1];
  assert(_dx==_dy);
}

/*-------------------------------------------------*/
//void TrapezMatrix3d::set_grid(std::shared_ptr<GridInterface> grid)
void TrapezMatrix3d::set_grid(const armaicvec& n, const armavec& dx)
{
  assert(n.n_elem==3);
  _nx = n[0];
  _ny = n[1];
  _nz = n[2];
  _vol = arma::prod(dx);
  _dx = dx[0];
  _dy = dx[1];
  _dz = dx[2];
  assert(_dx==_dy);
  assert(_dx==_dz);
  assert(_dy==_dz);
}
/*-------------------------------------------------*/
void TrapezMatrix2d::get_sparse_matrix(SparseMatrix& sp) const
{
  std::cerr << "Not written!\n"; exit(1);
  assert(0);
}
/*-------------------------------------------------*/
void TrapezMatrix3d::get_sparse_matrix(SparseMatrix& sp) const
{
  std::cerr << "Not written!\n"; exit(1);
  assert(0);
}

/*-------------------------------------------------*/
void TrapezMatrix2d::_boundary(NodeVector& out) const
{
  for(int ix=0;ix<_nx;ix++)
  {
    out.atp(ix,0)    = 0.0;
    out.atp(ix,_ny-1) = 0.0;
  }
  for(int iy=0;iy<_ny;iy++)
  {
    out.atp(0,iy)    = 0.0;
    out.atp(_nx-1,iy) = 0.0;
  }
}

/*-------------------------------------------------*/
void TrapezMatrix2d::dot(NodeVector& out, const NodeVector& in, double d) const
{
  // Laplacien   elements finis q1  (5-point-stencil)
  double d0 = 4.0 * d;
  double d1 = -1.0 * d;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      out.atp(ix,iy) += d0 * in.atp(ix,iy)
      + d1* in.atp(ix-1,iy  ) + d1* in.atp(ix+1,iy)
      + d1* in.atp(ix  ,iy-1) + d1* in.atp(ix  ,iy+1);
    }
  }
  // Conditions aux limites ( Dirichlet)
  _boundary(out);
}


/*-------------------------------------------------*/
void TrapezMatrix2d::jacobi(NodeVector& out, const NodeVector& in) const
{
  double d0inv = 0.25;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      out.atp(ix,iy) = d0inv * in.atp(ix,iy);
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void TrapezMatrix2d::gauss_seidel1(NodeVector& out, const NodeVector& in) const
{
  /*
   (ix+p)*ny + iy+q < ix*ny + iy
   p*ny +q < 0
   p=-1 q=-1,0,1
   p= 0 q=-1
   */
  double d0inv = 0.25;
  double d1 = -1.0;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      out.atp(ix,iy) = d0inv * (
                                in.atp(ix,iy)
                                -d1* out.atp(ix-1,iy  )
                                -d1* out.atp(ix  ,iy-1)
                                );
    }
  }
}
/*-------------------------------------------------*/
void TrapezMatrix2d::gauss_seidel2(NodeVector& out, const NodeVector& in) const
{
  double d0inv = 0.25;
  double d1 = -1.0;
  for(int ix=_nx-1;ix>=0;ix--)
  {
    for(int iy=_ny-1;iy>=0;iy--)
    {
      out.atp(ix,iy) = d0inv * (
                                in.atp(ix,iy)
                                -d1* out.atp(ix+1,iy  )
                                -d1* out.atp(ix  ,iy+1)
                                );
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void TrapezMatrix3d::_boundary(NodeVector& out) const
{
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      out.atp(ix,iy,0)    = 0.0;
      out.atp(ix,iy,_nz-1) = 0.0;
    }
  }
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iz=0;iz<_nz;iz++)
    {
      out.atp(ix,0,   iz) = 0.0;
      out.atp(ix,_ny-1,iz) = 0.0;
    }
  }
  for(int iy=0;iy<_ny;iy++)
  {
    for(int iz=0;iz<_nz;iz++)
    {
      out.atp(0,   iy,iz) = 0.0;
      out.atp(_nx-1,iy,iz) = 0.0;
    }
  }
}

/*-------------------------------------------------*/
void TrapezMatrix3d::dot(NodeVector& out, const NodeVector& in, double d) const
{
  double e = d/arma::mean(out.n());
  double d0 = 6.0 * e;
  double d1 = -1.0 * e;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        out.atp(ix,iy,iz) +=
        d0* in.atp(ix,iy,iz)

        +d1* in.atp(ix-1,iy  ,iz  )
        +d1* in.atp(ix+1,iy  ,iz  )
        +d1* in.atp(ix  ,iy-1,iz  )
        +d1* in.atp(ix  ,iy+1,iz  )
        +d1* in.atp(ix  ,iy  ,iz-1)
        +d1* in.atp(ix+1,iy  ,iz+1)
        ;
      }
    }
  }
  // Conditions aux limites ( Dirichlet)
  _boundary(out);
}


/*-------------------------------------------------*/
void TrapezMatrix3d::jacobi(NodeVector& out, const NodeVector& in) const
{
  double e = 1.0/arma::mean(out.n());
  double d0 = 7.0 * e;
  double d0inv = 1.0/d0;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        out.atp(ix,iy,iz) = d0inv * in.atp(ix,iy,iz);
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void TrapezMatrix3d::gauss_seidel1(NodeVector& out, const NodeVector& in) const
{

  /*
   (ix+p)*ny*nz + (iy+q)*nz + iz+r < ix*ny*nz + iy*nz + iz
   p*ny*nz +q*nz +r < 0
   p=-1 q=-1,0,1  r=-1,0,1
   p= 0 q=-1 r=-1,0,1 q=0 r=-1
   */
  double e = 1.0/arma::mean(out.n());
  double d0 = 7.0 * e;
  double d1 = -1.0 * e;
  double d0inv = 1.0/d0;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        out.atp(ix,iy,iz) = d0inv*(
                                   in.atp(ix,iy,iz)
                                   -d1* out.atp(ix-1,iy  ,iz  )
                                   -d1* out.atp(ix  ,iy-1,iz  )
                                   -d1* out.atp(ix  ,iy  ,iz-1)
                                   );
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void TrapezMatrix3d::gauss_seidel2(NodeVector& out, const NodeVector& in) const
{
  double e = 1.0/arma::mean(out.n());
  double d0 = 7.0 * e;
  double d1 = -1.0 * e;
  double d0inv = 1.0/d0;
  for(int ix=_nx-1;ix>=0;ix--)
  {
    for(int iy=_ny-1;iy>=0;iy--)
    {
      for(int iz=_nz-1;iz>=0;iz--)
      {
        out.atp(ix,iy,iz) = d0inv*(
                                   in.atp(ix,iy,iz)+
                                   -d1* out.atp(ix+1,iy  ,iz  )
                                   -d1* out.atp(ix  ,iy+1,iz  )
                                   -d1* out.atp(ix  ,iy  ,iz+1)
                                   );
      }
    }
  }
  _boundary(out);
}
