//
//  fullmatrix.cpp
//  Fada
//
//  Created by Roland Becker on 15/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <math.h>
#include  "fullmatrix.hpp"
#include  "uniformgrid.hpp"
#include  "typedefs.hpp"
#include  "sparsematrix.hpp"

FullMatrix2d::~FullMatrix2d() {}
FullMatrix3d::~FullMatrix3d() {}

/*-------------------------------------------------*/
void FullMatrix2d::set_grid(std::shared_ptr<GridInterface> grid)
{
  std::shared_ptr<UniformGrid> ug = std::dynamic_pointer_cast<UniformGrid>(grid);
  assert(ug);
  assert(ug->dim()==2);
  _nx = ug->nx();
  _ny = ug->ny();
  _vol = arma::prod(ug->dx());
  _dx =ug->dx(0);
  _dy =ug->dx(1);
  assert(_dx==_dy);
}

/*-------------------------------------------------*/
void FullMatrix3d::set_grid(std::shared_ptr<GridInterface> grid)
{
  std::shared_ptr<UniformGrid> ug = std::dynamic_pointer_cast<UniformGrid>(grid);
  assert(ug);
  assert(ug->dim()==3);
  _nx = ug->nx();
  _ny = ug->ny();
  _nz = ug->nz();
  _vol = arma::prod(ug->dx());
  _dx =ug->dx(0);
  _dy =ug->dx(1);
  _dz =ug->dx(2);
  assert(_dx==_dy);
  assert(_dx==_dz);
}
/*-------------------------------------------------*/
void FullMatrix2d::get_sparse_matrix(SparseMatrix& sp) const
{
  double d0 = 8.0/3.0;
  double d1 = -1.0/3.0;
  // (nx+2)*(ny+2) - nx*ny + nx*ny - (nx-2)*(ny-2) = 4*(nx+ny)
  int size = 9*(_nx-2)*(_ny-2) + 4*(_nx+_ny);
  arma::umat locations(2, size);
  armavec values(size);
  int ofsx = _ny + 2;
  int ofsp = ofsx+1;
//  std::cerr << "_nx _ny " << _nx << " " << _ny << " ofsx " << ofsx << " ofsp " << ofsp << "\n";
  int count=0, i, j;
  for(int ix=1;ix<_nx-1;ix++)
  {
    for(int iy=1;iy<_ny-1;iy++)
    {
      i = ofsx*ix + iy + ofsp;
      j = ofsx*ix + iy + ofsp;
      locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d0;
      count++;
      j = ofsx*(ix-1) + iy + ofsp;
      locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
      count++;
      j = ofsx*(ix+1) + iy + ofsp;
      locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
      count++;
      j = ofsx*ix + iy-1 + ofsp;
      locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
      count++;
      j = ofsx*ix + iy+1 + ofsp;
      locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
      count++;
      j = ofsx*(ix-1) + iy-1 + ofsp;
      locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
      count++;
      j = ofsx*(ix-1) + iy+1 + ofsp;
      locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
      count++;
      j = ofsx*(ix+1) + iy-1 + ofsp;
      locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
      count++;
      j = ofsx*(ix+1) + iy+1 + ofsp;
      locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
      count++;
    }
  }
  //bdry
  for(int ix=0;ix<_nx;ix++)
  {
    i = ofsx*ix + 0 + ofsp;
    locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
    count++;
    i = ofsx*ix + _ny-1 + ofsp;
    locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
    count++;
  }
  for(int iy=1;iy<_ny-1;iy++)
  {
    i = ofsx*0 + iy + ofsp;
    locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
    count++;
    i = ofsx*(_nx-1) + iy + ofsp;
    locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
    count++;
  }
  //aux
  for(int ix=0;ix<_nx+2;ix++)
  {
    i = ofsx*ix + 0;
    locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
    count++;
    i = ofsx*ix + _ny+1;
    locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
    count++;
  }
  for(int iy=1;iy<_ny+1;iy++)
  {
    i = ofsx*0 + iy;
    locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
    count++;
    i = ofsx*(_nx+1) + iy;
    locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
    count++;
  }
//  std::cerr << "locations i " << locations.row(0);
//  std::cerr << "locations j " << locations.row(1);
//  std::cerr << "values " << values.t();
  sp.set_elements(locations, values);
}
/*-------------------------------------------------*/
void FullMatrix3d::get_sparse_matrix(SparseMatrix& sp) const
{
  double e = _dx;
  double d0 = 8.0/3.0 * e;
  double d1 = -0.0/3.0 * e;
  double d2 = -1.0/6.0 * e;
  double d3 = -1.0/12.0 * e;
  // (nx+2)*(ny+2)*(nz+2) - nx*ny*nz + nx*ny*z - (nx-2)*(ny-2)*(nz-2)
  //= 8 + 4*(nx+ny+nz)+2*(nx*ny+nx*nz+ny*nz) -( -8 + 4*(nx+ny+nz) - 2*(nx*ny+nx*nz+ny*nz)  )
  //= 16 + 4*(nx*ny+nx*nz+ny*nz)
  int size = 27*(_nx-2)*(_ny-2)*(_nz-2) + 16 + 4*(_nx*_ny+_nx*_nz+_ny*_nz);
  arma::umat locations(2, size);
  armavec values(size);
  int ofsy = _nz + 2;
  int ofsx = ofsy*(_ny + 2);
  int ofsp = ofsx+ofsy+1;
//  std::cerr << "ofsx " << ofsx << " ofsy " << ofsy << " ofsp " << ofsp<< " size " << size << "\n";
  int count=0, i, j;
  for(int ix=1;ix<_nx-1;ix++)
  {
    for(int iy=1;iy<_ny-1;iy++)
    {
      for(int iz=1;iz<_nz-1;iz++)
      {
        i = ofsx*ix + ofsy*iy + iz + ofsp;
        j = ofsx*ix + ofsy*iy + iz + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d0;
        count++;
        j = ofsx*(ix-1) + ofsy*iy + iz + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
        count++;
        j = ofsx*(ix+1) + ofsy*iy + iz + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
        count++;
        j = ofsx*ix + ofsy*(iy-1) + iz + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
        count++;
        j = ofsx*ix + ofsy*(iy+1) + iz + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
        count++;
        j = ofsx*ix + ofsy*iy + iz-1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
        count++;
        j = ofsx*ix + ofsy*iy + iz+1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d1;
        count++;

        j = ofsx*(ix-1) + ofsy*(iy-1) + iz + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*(ix-1) + ofsy*(iy+1) + iz + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*(ix+1) + ofsy*(iy-1) + iz + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*(ix+1) + ofsy*(iy+1) + iz + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*(ix-1) + ofsy*iy + iz-1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*(ix-1) + ofsy*iy + iz+1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*(ix+1) + ofsy*iy + iz-1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*(ix+1) + ofsy*iy + iz+1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*ix + ofsy*(iy-1) + iz-1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*ix + ofsy*(iy-1) + iz+1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*ix + ofsy*(iy+1) + iz-1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;
        j = ofsx*ix + ofsy*(iy+1) + iz+1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d2;
        count++;

        j = ofsx*(ix-1) + ofsy*(iy-1) + iz-1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d3;
        count++;
        j = ofsx*(ix-1) + ofsy*(iy-1) + iz+1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d3;
        count++;
        j = ofsx*(ix-1) + ofsy*(iy+1) + iz-1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d3;
        count++;
        j = ofsx*(ix-1) + ofsy*(iy+1) + iz+1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d3;
        count++;
        j = ofsx*(ix+1) + ofsy*(iy-1) + iz-1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d3;
        count++;
        j = ofsx*(ix+1) + ofsy*(iy-1) + iz+1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d3;
        count++;
        j = ofsx*(ix+1) + ofsy*(iy+1) + iz-1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d3;
        count++;
        j = ofsx*(ix+1) + ofsy*(iy+1) + iz+1 + ofsp;
        locations(0, count ) = i;      locations(1, count ) = j;      values[count] = d3;
        count++;
      }
    }
  }
  //bdry
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      i = ofsx*ix + ofsy*iy + 0 + ofsp;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
      i = ofsx*ix + ofsy*iy + _nz-1 + ofsp;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
    }
  }
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iz=1;iz<_nz-1;iz++)
    {
      i = ofsx*ix + ofsy*0 + iz + ofsp;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
      i = ofsx*ix + ofsy*(_ny-1) + iz + ofsp;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
    }
  }
  for(int iy=1;iy<_ny-1;iy++)
  {
    for(int iz=1;iz<_nz-1;iz++)
    {
      i = ofsx*0 + ofsy*iy + iz + ofsp;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
      i = ofsx*(_nx-1) + ofsy*iy + iz + ofsp;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
    }
  }
  //aux
  for(int ix=0;ix<_nx+2;ix++)
  {
    for(int iy=0;iy<_ny+2;iy++)
    {
      i = ofsx*ix + ofsy*iy + 0;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
      i = ofsx*ix + ofsy*iy + _nz+1;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
    }
  }
  for(int ix=0;ix<_nx+2;ix++)
  {
    for(int iz=1;iz<_nz+1;iz++)
    {
      i = ofsx*ix + ofsy*0 + iz;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
      i = ofsx*ix + ofsy*(_ny+1) + iz;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
    }
  }
  for(int iy=1;iy<_ny+1;iy++)
  {
    for(int iz=1;iz<_nz+1;iz++)
    {
      i = ofsx*0 + ofsy*iy + iz;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
      i = ofsx*(_nx+1) + ofsy*iy + iz;
      locations(0, count ) = i;      locations(1, count ) = i;      values[count] = 1;
      count++;
    }
  }
//  for(int i=0;i<size;i++)
//  {
//    std::cerr <<locations(0, i) << " : " << locations(1, i) << " -> " << values[i]<<"\n";
//  }
  sp.set_elements(locations, values);
}

/*-------------------------------------------------*/
void FullMatrix2d::_boundary(Vector& out) const
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
void FullMatrix2d::dot(Vector& out, const Vector& in, double d) const
{
  // Laplacien   elements finis q1  (9-point-stencil)
  double d0 = 8.0/3.0 * d;
  double d1 = -1.0/3.0 * d;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      out.atp(ix,iy) += d0 * in.atp(ix,iy)
      + d1* in.atp(ix-1,iy  ) + d1* in.atp(ix+1,iy)
      + d1* in.atp(ix  ,iy-1) + d1* in.atp(ix  ,iy+1)
      + d1* in.atp(ix-1,iy-1) + d1* in.atp(ix-1,iy+1)
      + d1* in.atp(ix+1,iy-1) + d1* in.atp(ix+1,iy+1);
    }
  }
  // Conditions aux limites ( Dirichlet)
  _boundary(out);
}


/*-------------------------------------------------*/
void FullMatrix2d::jacobi(Vector& out, const Vector& in) const
{
  double omega = 0.8;
  double d0inv = 3.0/8.0 * omega;
  d0inv = 0.1;
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
void FullMatrix2d::gauss_seidel1(Vector& out, const Vector& in) const
{
  /*
   (ix+p)*ny + iy+q < ix*ny + iy
   p*ny +q < 0
   p=-1 q=-1,0,1
   p= 0 q=-1
   */
//  double omega = 0.8;
  double d0inv = 3.0/8.0;
  double d1 = -1.0/3.0;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      out.atp(ix,iy) = d0inv * (
                                in.atp(ix,iy)
                                -d1* out.atp(ix-1,iy-1)
                                -d1* out.atp(ix-1,iy  )
                                -d1* out.atp(ix-1,iy+1)
                                -d1* out.atp(ix  ,iy-1)
                                );
    }
  }
  _boundary(out);
}
/*-------------------------------------------------*/
void FullMatrix2d::gauss_seidel2(Vector& out, const Vector& in) const
{
//  double omega = 0.8;
  double d0inv = 3.0/8.0;
  double d1 = -1.0/3.0;
  for(int ix=_nx-1;ix>=0;ix--)
  {
    for(int iy=_ny-1;iy>=0;iy--)
    {
      out.atp(ix,iy) = d0inv * (
                                in.atp(ix,iy)
                                -d1* out.atp(ix+1,iy-1)
                                -d1* out.atp(ix+1,iy  )
                                -d1* out.atp(ix+1,iy+1)
                                -d1* out.atp(ix  ,iy+1)
                                );
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void FullMatrix3d::_boundary(Vector& out) const
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
void FullMatrix3d::dot(Vector& out, const Vector& in, double d) const
{
//  double e = d/arma::mean(out.n());
  double e = d*_dx;
  double d0 = 8.0/3.0 * e;
  double d1 = -0.0/3.0 * e;
  double d2 = -1.0/6.0 * e;
  double d3 = -1.0/12.0 * e;
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
        +d1* in.atp(ix  ,iy  ,iz+1)
        
        +d2* in.atp(ix-1,iy-1,iz  )
        +d2* in.atp(ix-1,iy+1,iz  )
        +d2* in.atp(ix+1,iy-1,iz  )
        +d2* in.atp(ix+1,iy+1,iz  )
        +d2* in.atp(ix-1,iy  ,iz-1)
        +d2* in.atp(ix-1,iy  ,iz+1)
        +d2* in.atp(ix+1,iy  ,iz-1)
        +d2* in.atp(ix+1,iy  ,iz+1)
        +d2* in.atp(ix  ,iy-1,iz-1)
        +d2* in.atp(ix  ,iy-1,iz+1)
        +d2* in.atp(ix  ,iy+1,iz-1)
        +d2* in.atp(ix  ,iy+1,iz+1)
        
        +d3* in.atp(ix-1,iy-1,iz-1)
        +d3* in.atp(ix-1,iy-1,iz+1)
        +d3* in.atp(ix-1,iy+1,iz-1)
        +d3* in.atp(ix-1,iy+1,iz+1)
        +d3* in.atp(ix+1,iy-1,iz-1)
        +d3* in.atp(ix+1,iy-1,iz+1)
        +d3* in.atp(ix+1,iy+1,iz-1)
        +d3* in.atp(ix+1,iy+1,iz+1)
        ;
      }
    }
  }
  // Conditions aux limites ( Dirichlet)
  _boundary(out);
}


/*-------------------------------------------------*/
void FullMatrix3d::jacobi(Vector& out, const Vector& in) const
{
  double omega = 0.8;
  double d0inv = 3.0/8.0 * arma::mean(out.n()) * omega;
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
void FullMatrix3d::gauss_seidel1(Vector& out, const Vector& in) const
{
  
  /*
   (ix+p)*ny*nz + (iy+q)*nz + iz+r < ix*ny*nz + iy*nz + iz
   p*ny*nz +q*nz +r < 0
   p=-1 q=-1,0,1  r=-1,0,1
   p= 0 q=-1 r=-1,0,1 q=0 r=-1
   */
  double e = 1.0/arma::mean(out.n());
  double d0 = 8.0/3.0 * e;
  double d1 = -0.0/3.0 * e;
  double d2 = -1.0/6.0 * e;
  double d3 = -1.0/12.0 * e;
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
                                   
                                   -d2* out.atp(ix-1,iy-1,iz  )
                                   -d2* out.atp(ix-1,iy+1,iz  )
                                   -d2* out.atp(ix-1,iy  ,iz-1)
                                   -d2* out.atp(ix-1,iy  ,iz+1)
                                   -d2* out.atp(ix  ,iy-1,iz-1)
                                   -d2* out.atp(ix  ,iy-1,iz+1)
                                   
                                   -d3* out.atp(ix-1,iy-1,iz-1)
                                   -d3* out.atp(ix-1,iy-1,iz+1)
                                   -d3* out.atp(ix-1,iy+1,iz-1)
                                   -d3* out.atp(ix-1,iy+1,iz+1)
                                   );
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void FullMatrix3d::gauss_seidel2(Vector& out, const Vector& in) const
{
  double e = 1.0/arma::mean(out.n());
  double d0 = 8.0/3.0 * e;
  double d1 = -0.0/3.0 * e;
  double d2 = -1.0/6.0 * e;
  double d3 = -1.0/12.0 * e;
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
                                   
                                   -d2* out.atp(ix+1,iy  ,iz+1)
                                   -d2* out.atp(ix+1,iy-1,iz  )
                                   -d2* out.atp(ix+1,iy+1,iz  )
                                   -d2* out.atp(ix+1,iy  ,iz-1)
                                   -d2* out.atp(ix  ,iy+1,iz-1)
                                   -d2* out.atp(ix  ,iy+1,iz+1)
                                   
                                   -d3* out.atp(ix+1,iy-1,iz-1)
                                   -d3* out.atp(ix+1,iy-1,iz+1)
                                   -d3* out.atp(ix+1,iy+1,iz-1)
                                   -d3* out.atp(ix+1,iy+1,iz+1)
                                   );
      }
    }
  }
  _boundary(out);
}
