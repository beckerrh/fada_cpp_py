//
//  stencil3d.cpp
//  Fada
//
//  Created by Roland Becker on 02/06/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  "stencil3d.hpp"
#include  "sparsematrix.hpp"

/*-------------------------------------------------*/
void Stencil3d::_boundary(NodeVector& out) const
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
void Stencil3d27::set_grid(const armaicvec& n, const armavec& coef)
{
  assert(n.n_elem==3);
  _nx = n[0];
  _ny = n[1];
  _nz = n[2];
  _coef = coef;
//  std::cerr << "Stencil3d27() _coef="<<arma::sum(_coef)<<"\n";
}
/*-------------------------------------------------*/
void Stencil3d27::dot(NodeVector& out, const NodeVector& in, double d) const
{
  arma::vec::fixed<27> coef = d*_coef;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        out.atp(ix,iy,iz) +=
           coef[ 0]* in.atp(ix-1,iy-1,iz-1)
          +coef[ 1]* in.atp(ix-1,iy-1,iz  )
          +coef[ 2]* in.atp(ix-1,iy-1,iz+1)
          +coef[ 3]* in.atp(ix-1,iy  ,iz-1)
          +coef[ 4]* in.atp(ix-1,iy  ,iz  )
          +coef[ 5]* in.atp(ix-1,iy  ,iz+1)
          +coef[ 6]* in.atp(ix-1,iy+1,iz-1)
          +coef[ 7]* in.atp(ix-1,iy+1,iz  )
          +coef[ 8]* in.atp(ix-1,iy+1,iz+1)
          +coef[ 9]* in.atp(ix  ,iy-1,iz-1)
          +coef[10]* in.atp(ix  ,iy-1,iz  )
          +coef[11]* in.atp(ix  ,iy-1,iz+1)
          +coef[12]* in.atp(ix  ,iy  ,iz-1)
          +coef[13]* in.atp(ix  ,iy  ,iz  )
          +coef[14]* in.atp(ix  ,iy  ,iz+1)
          +coef[15]* in.atp(ix  ,iy+1,iz-1)
          +coef[16]* in.atp(ix  ,iy+1,iz  )
          +coef[17]* in.atp(ix  ,iy+1,iz+1)
          +coef[18]* in.atp(ix+1,iy-1,iz-1)
          +coef[19]* in.atp(ix+1,iy-1,iz  )
          +coef[20]* in.atp(ix+1,iy-1,iz+1)
          +coef[21]* in.atp(ix+1,iy  ,iz-1)
          +coef[22]* in.atp(ix+1,iy  ,iz  )
          +coef[23]* in.atp(ix+1,iy  ,iz+1)
          +coef[24]* in.atp(ix+1,iy+1,iz-1)
          +coef[25]* in.atp(ix+1,iy+1,iz  )
          +coef[26]* in.atp(ix+1,iy+1,iz+1)
        ;
      }
    }
  }
  // Conditions aux limites ( Dirichlet)
  _boundary(out);
}


/*-------------------------------------------------*/
void Stencil3d27::jacobi(NodeVector& out, const NodeVector& in) const
{
  double d0inv = 1.0/_coef[13];
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
void Stencil3d27::gauss_seidel1(NodeVector& out, const NodeVector& in) const
{
  /*
   (ix+p)*ny*nz + (iy+q)*nz + iz+r < ix*ny*nz + iy*nz + iz
   p*ny*nz +q*nz +r < 0
   p=-1 q=-1,0,1  r=-1,0,1
   p= 0 q=-1 r=-1,0,1 q=0 r=-1
   */
  double d0inv = 1.0/_coef[13];
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        out.atp(ix,iy,iz) = d0inv*(
                                   in.atp(ix,iy,iz)
                                   -_coef[ 0]* out.atp(ix-1,iy-1,iz-1)
                                   -_coef[ 1]* out.atp(ix-1,iy-1,iz  )
                                   -_coef[ 2]* out.atp(ix-1,iy-1,iz+1)
                                   -_coef[ 3]* out.atp(ix-1,iy  ,iz-1)
                                   -_coef[ 4]* out.atp(ix-1,iy  ,iz  )
                                   -_coef[ 5]* out.atp(ix-1,iy  ,iz+1)
                                   -_coef[ 6]* out.atp(ix-1,iy+1,iz-1)
                                   -_coef[ 7]* out.atp(ix-1,iy+1,iz  )
                                   -_coef[ 8]* out.atp(ix-1,iy+1,iz+1)
                                   -_coef[ 9]* out.atp(ix  ,iy-1,iz-1)
                                   -_coef[10]* out.atp(ix  ,iy-1,iz  )
                                   -_coef[11]* out.atp(ix  ,iy-1,iz+1)
                                   -_coef[12]* out.atp(ix  ,iy  ,iz-1)
                                   );
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil3d27::gauss_seidel2(NodeVector& out, const NodeVector& in) const
{
  double d0inv = 1.0/_coef[13];
  for(int ix=_nx-1;ix>=0;ix--)
  {
    for(int iy=_ny-1;iy>=0;iy--)
    {
      for(int iz=_nz-1;iz>=0;iz--)
      {
        out.atp(ix,iy,iz) = d0inv*(
                                   in.atp(ix,iy,iz)
                                   -_coef[14]* out.atp(ix  ,iy  ,iz+1)
                                   -_coef[15]* out.atp(ix  ,iy+1,iz-1)
                                   -_coef[16]* out.atp(ix  ,iy+1,iz  )
                                   -_coef[17]* out.atp(ix  ,iy+1,iz+1)
                                   -_coef[18]* out.atp(ix+1,iy-1,iz-1)
                                   -_coef[19]* out.atp(ix+1,iy-1,iz  )
                                   -_coef[20]* out.atp(ix+1,iy-1,iz+1)
                                   -_coef[21]* out.atp(ix+1,iy  ,iz-1)
                                   -_coef[22]* out.atp(ix+1,iy  ,iz  )
                                   -_coef[23]* out.atp(ix+1,iy  ,iz+1)
                                   -_coef[24]* out.atp(ix+1,iy+1,iz-1)
                                   -_coef[25]* out.atp(ix+1,iy+1,iz  )
                                   -_coef[26]* out.atp(ix+1,iy+1,iz+1)
                                   );
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil3d27::get_sparse_matrix(SparseMatrix& sp) const
{
  // (nx+2)*(ny+2)*(nz+2) - nx*ny*nz + nx*ny*z - (nx-2)*(ny-2)*(nz-2)
  //= 8 + 4*(nx+ny+nz)+2*(nx*ny+nx*nz+ny*nz) -( -8 + 4*(nx+ny+nz) - 2*(nx*ny+nx*nz+ny*nz)  )
  //= 16 + 4*(nx*ny+nx*nz+ny*nz)
  int size = 27*(_nx-2)*(_ny-2)*(_nz-2) + 16 + 4*(_nx*_ny+_nx*_nz+_ny*_nz);
  int ofsy = _nz + 2;
  int ofsx = ofsy*(_ny + 2);
  int ofsp = ofsx+ofsy+1;
//  std::cerr << "ofsx " << ofsx << " ofsy " << ofsy << " ofsp " << ofsp<< " size " << size << "\n";
  int i, j;
  Construct_Elements ce(size);
  for(int ix=1;ix<_nx-1;ix++)
  {
    for(int iy=1;iy<_ny-1;iy++)
    {
      for(int iz=1;iz<_nz-1;iz++)
      {
        i = ofsx*ix + ofsy*iy + iz + ofsp;

        j = ofsx*(ix-1) + ofsy*(iy-1) + iz-1 + ofsp;
        ce.add(i, j, _coef[0]);
        j = ofsx*(ix-1) + ofsy*(iy-1) + iz   + ofsp;
        ce.add(i, j, _coef[1]);
        j = ofsx*(ix-1) + ofsy*(iy-1) + iz+1 + ofsp;
        ce.add(i, j, _coef[2]);
        j = ofsx*(ix-1) + ofsy*(iy  ) + iz-1 + ofsp;
        ce.add(i, j, _coef[3]);
        j = ofsx*(ix-1) + ofsy*(iy  ) + iz   + ofsp;
        ce.add(i, j, _coef[4]);
        j = ofsx*(ix-1) + ofsy*(iy  ) + iz+1 + ofsp;
        ce.add(i, j, _coef[5]);
        j = ofsx*(ix-1) + ofsy*(iy+1) + iz-1 + ofsp;
        ce.add(i, j, _coef[6]);
        j = ofsx*(ix-1) + ofsy*(iy+1) + iz   + ofsp;
        ce.add(i, j, _coef[7]);
        j = ofsx*(ix-1) + ofsy*(iy+1) + iz+1 + ofsp;
        ce.add(i, j, _coef[8]);

        j = ofsx*(ix  ) + ofsy*(iy-1) + iz-1 + ofsp;
        ce.add(i, j, _coef[9]);
        j = ofsx*(ix  ) + ofsy*(iy-1) + iz   + ofsp;
        ce.add(i, j, _coef[10]);
        j = ofsx*(ix  ) + ofsy*(iy-1) + iz+1 + ofsp;
        ce.add(i, j, _coef[11]);
        j = ofsx*(ix  ) + ofsy*(iy  ) + iz-1 + ofsp;
        ce.add(i, j, _coef[12]);
        j = ofsx*(ix  ) + ofsy*(iy  ) + iz   + ofsp;
        ce.add(i, j, _coef[13]);
        j = ofsx*(ix  ) + ofsy*(iy  ) + iz+1 + ofsp;
        ce.add(i, j, _coef[14]);
        j = ofsx*(ix  ) + ofsy*(iy+1) + iz-1 + ofsp;
        ce.add(i, j, _coef[15]);
        j = ofsx*(ix  ) + ofsy*(iy+1) + iz   + ofsp;
        ce.add(i, j, _coef[16]);
        j = ofsx*(ix  ) + ofsy*(iy+1) + iz+1 + ofsp;
        ce.add(i, j, _coef[17]);

        j = ofsx*(ix+1) + ofsy*(iy-1) + iz-1 + ofsp;
        ce.add(i, j, _coef[18]);
        j = ofsx*(ix+1) + ofsy*(iy-1) + iz   + ofsp;
        ce.add(i, j, _coef[19]);
        j = ofsx*(ix+1) + ofsy*(iy-1) + iz+1 + ofsp;
        ce.add(i, j, _coef[20]);
        j = ofsx*(ix+1) + ofsy*(iy  ) + iz-1 + ofsp;
        ce.add(i, j, _coef[21]);
        j = ofsx*(ix+1) + ofsy*(iy  ) + iz   + ofsp;
        ce.add(i, j, _coef[22]);
        j = ofsx*(ix+1) + ofsy*(iy  ) + iz+1 + ofsp;
        ce.add(i, j, _coef[23]);
        j = ofsx*(ix+1) + ofsy*(iy+1) + iz-1 + ofsp;
        ce.add(i, j, _coef[24]);
        j = ofsx*(ix+1) + ofsy*(iy+1) + iz   + ofsp;
        ce.add(i, j, _coef[25]);
        j = ofsx*(ix+1) + ofsy*(iy+1) + iz+1 + ofsp;
        ce.add(i, j, _coef[26]);
      }
    }
  }
  //bdry
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      i = ofsx*ix + ofsy*iy + 0 + ofsp;
      ce.add(i, i, 1);
      i = ofsx*ix + ofsy*iy + _nz-1 + ofsp;
      ce.add(i, i, 1);
    }
  }
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iz=1;iz<_nz-1;iz++)
    {
      i = ofsx*ix + ofsy*0 + iz + ofsp;
      ce.add(i, i, 1);
      i = ofsx*ix + ofsy*(_ny-1) + iz + ofsp;
      ce.add(i, i, 1);
    }
  }
  for(int iy=1;iy<_ny-1;iy++)
  {
    for(int iz=1;iz<_nz-1;iz++)
    {
      i = ofsx*0 + ofsy*iy + iz + ofsp;
      ce.add(i, i, 1);
      i = ofsx*(_nx-1) + ofsy*iy + iz + ofsp;
      ce.add(i, i, 1);
    }
  }
  //aux
  for(int ix=0;ix<_nx+2;ix++)
  {
    for(int iy=0;iy<_ny+2;iy++)
    {
      i = ofsx*ix + ofsy*iy + 0;
      ce.add(i, i, 1);
      i = ofsx*ix + ofsy*iy + _nz+1;
      ce.add(i, i, 1);
    }
  }
  for(int ix=0;ix<_nx+2;ix++)
  {
    for(int iz=1;iz<_nz+1;iz++)
    {
      i = ofsx*ix + ofsy*0 + iz;
      ce.add(i, i, 1);
      i = ofsx*ix + ofsy*(_ny+1) + iz;
      ce.add(i, i, 1);
    }
  }
  for(int iy=1;iy<_ny+1;iy++)
  {
    for(int iz=1;iz<_nz+1;iz++)
    {
      i = ofsx*0 + ofsy*iy + iz;
      ce.add(i, i, 1);
      i = ofsx*(_nx+1) + ofsy*iy + iz;
      ce.add(i, i, 1);
    }
  }
  sp.set_elements(ce.locations(), ce.values());
}

/*-------------------------------------------------*/
void Stencil3d7::set_grid(const armaicvec& n, const armavec& coef)
{
  assert(n.n_elem==3);
  _nx = n[0];
  _ny = n[1];
  _nz = n[2];
  _coef = coef;
//  std::cerr << "Stencil3d27() _coef="<<arma::sum(_coef)<<"\n";
}
/*-------------------------------------------------*/
void Stencil3d7::dot(NodeVector& out, const NodeVector& in, double d) const
{
  arma::vec::fixed<7> coef = d*_coef;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        out.atp(ix,iy,iz) +=
           coef[ 0]* in.atp(ix-1,iy  ,iz  )
          +coef[ 1]* in.atp(ix  ,iy-1,iz  )
          +coef[ 2]* in.atp(ix  ,iy  ,iz-1)
          +coef[ 3]* in.atp(ix  ,iy  ,iz  )
          +coef[ 4]* in.atp(ix  ,iy  ,iz+1)
          +coef[ 5]* in.atp(ix  ,iy+1,iz  )
          +coef[ 6]* in.atp(ix+1,iy  ,iz  )
        ;
      }
    }
  }
  // Conditions aux limites ( Dirichlet)
  _boundary(out);
}


/*-------------------------------------------------*/
void Stencil3d7::jacobi(NodeVector& out, const NodeVector& in) const
{
  double d0inv = 1.0/_coef[13];
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
void Stencil3d7::gauss_seidel1(NodeVector& out, const NodeVector& in) const
{
  /*
   (ix+p)*ny*nz + (iy+q)*nz + iz+r < ix*ny*nz + iy*nz + iz
   p*ny*nz +q*nz +r < 0
   p=-1 q=-1,0,1  r=-1,0,1
   p= 0 q=-1 r=-1,0,1 q=0 r=-1
   */
  double d0inv = 1.0/_coef[3];
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        out.atp(ix,iy,iz) = d0inv*(
                                   in.atp(ix,iy,iz)
                                   -_coef[0]* out.atp(ix-1,iy  ,iz  )
                                   -_coef[1]* out.atp(ix  ,iy-1,iz  )
                                   -_coef[2]* out.atp(ix  ,iy  ,iz-1)
                                   );
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil3d7::gauss_seidel2(NodeVector& out, const NodeVector& in) const
{
  double d0inv = 1.0/_coef[3];
  for(int ix=_nx-1;ix>=0;ix--)
  {
    for(int iy=_ny-1;iy>=0;iy--)
    {
      for(int iz=_nz-1;iz>=0;iz--)
      {
        out.atp(ix,iy,iz) = d0inv*(
                                   in.atp(ix,iy,iz)
                                   -_coef[4]* out.atp(ix  ,iy  ,iz+1)
                                   -_coef[5]* out.atp(ix  ,iy+1,iz  )
                                   -_coef[6]* out.atp(ix+1,iy  ,iz  )
                                   );
      }
    }
  }
  _boundary(out);
}

/*-------------------------------------------------*/
void Stencil3d7::get_sparse_matrix(SparseMatrix& sp) const
{
  // (nx+2)*(ny+2)*(nz+2) - nx*ny*nz + nx*ny*z - (nx-2)*(ny-2)*(nz-2)
  //= 8 + 4*(nx+ny+nz)+2*(nx*ny+nx*nz+ny*nz) -( -8 + 4*(nx+ny+nz) - 2*(nx*ny+nx*nz+ny*nz)  )
  //= 16 + 4*(nx*ny+nx*nz+ny*nz)
  int size = 7*(_nx-2)*(_ny-2)*(_nz-2) + 16 + 4*(_nx*_ny+_nx*_nz+_ny*_nz);
  int ofsy = _nz + 2;
  int ofsx = ofsy*(_ny + 2);
  int ofsp = ofsx+ofsy+1;
//  std::cerr << "ofsx " << ofsx << " ofsy " << ofsy << " ofsp " << ofsp<< " size " << size << "\n";
  int i, j;
  Construct_Elements ce(size);
  for(int ix=1;ix<_nx-1;ix++)
  {
    for(int iy=1;iy<_ny-1;iy++)
    {
      for(int iz=1;iz<_nz-1;iz++)
      {
        i = ofsx*ix + ofsy*iy + iz + ofsp;

        j = ofsx*(ix-1) + ofsy*(iy  ) + iz   + ofsp;
        ce.add(i, j, _coef[0]);

        j = ofsx*(ix  ) + ofsy*(iy-1) + iz   + ofsp;
        ce.add(i, j, _coef[1]);
        j = ofsx*(ix  ) + ofsy*(iy  ) + iz-1 + ofsp;
        ce.add(i, j, _coef[2]);
        j = ofsx*(ix  ) + ofsy*(iy  ) + iz   + ofsp;
        ce.add(i, j, _coef[3]);
        j = ofsx*(ix  ) + ofsy*(iy) + iz+1   + ofsp;
        ce.add(i, j, _coef[4]);
        j = ofsx*(ix  ) + ofsy*(iy+1) + iz   + ofsp;
        ce.add(i, j, _coef[5]);
        j = ofsx*(ix+1) + ofsy*(iy  ) + iz   + ofsp;
        ce.add(i, j, _coef[6]);
      }
    }
  }
  //bdry
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      i = ofsx*ix + ofsy*iy + 0 + ofsp;
      ce.add(i, i, 1);
      i = ofsx*ix + ofsy*iy + _nz-1 + ofsp;
      ce.add(i, i, 1);
    }
  }
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iz=1;iz<_nz-1;iz++)
    {
      i = ofsx*ix + ofsy*0 + iz + ofsp;
      ce.add(i, i, 1);
      i = ofsx*ix + ofsy*(_ny-1) + iz + ofsp;
      ce.add(i, i, 1);
    }
  }
  for(int iy=1;iy<_ny-1;iy++)
  {
    for(int iz=1;iz<_nz-1;iz++)
    {
      i = ofsx*0 + ofsy*iy + iz + ofsp;
      ce.add(i, i, 1);
      i = ofsx*(_nx-1) + ofsy*iy + iz + ofsp;
      ce.add(i, i, 1);
    }
  }
  //aux
  for(int ix=0;ix<_nx+2;ix++)
  {
    for(int iy=0;iy<_ny+2;iy++)
    {
      i = ofsx*ix + ofsy*iy + 0;
      ce.add(i, i, 1);
      i = ofsx*ix + ofsy*iy + _nz+1;
      ce.add(i, i, 1);
    }
  }
  for(int ix=0;ix<_nx+2;ix++)
  {
    for(int iz=1;iz<_nz+1;iz++)
    {
      i = ofsx*ix + ofsy*0 + iz;
      ce.add(i, i, 1);
      i = ofsx*ix + ofsy*(_ny+1) + iz;
      ce.add(i, i, 1);
    }
  }
  for(int iy=1;iy<_ny+1;iy++)
  {
    for(int iz=1;iz<_nz+1;iz++)
    {
      i = ofsx*0 + ofsy*iy + iz;
      ce.add(i, i, 1);
      i = ofsx*(_nx+1) + ofsy*iy + iz;
      ce.add(i, i, 1);
    }
  }
  sp.set_elements(ce.locations(), ce.values());
}
