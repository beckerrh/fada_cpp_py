//
//  vector.c
//  Fada
//
//  Created by Roland Becker on 03/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <stdio.h>
#include  "seamvector.hpp"
#include  "uniformgrid.hpp"
#include  "nodevector.hpp"

/*-------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const SeamVector& v)
{
//  const armavec& tarma =static_cast<const armavec&>(v);
//  os << tarma.t()<< "n=" << v.n().t()<< "ofs=" << v.ofs().t();
  os << v.t() << "n=" << v.n().t() << "ofs=" << v.ofs().t();
  return(os);
}

/*-------------------------------------------------*/
void SeamVector::tovector(NodeVector& u) const
{
  armaicvec n = u.n();
  if (dim() == 2)
  {
    for (int ix = 0; ix < n[0]; ix++)
    {
      for (int iy = 0; iy < n[1]; iy++)
      {
        u.at(ix, iy) = this->atp(ix, iy);
      }
    }
  }
  else if (dim() == 3)
  {
    for (int ix = 0; ix < n[0]; ix++)
    {
      for (int iy = 0; iy < n[1]; iy++)
      {
        for (int iz = 0; iz < n[2]; iz++)
        {
          u.at(ix, iy, iz) = this->atp(ix, iy, iz);
        }
      }
    }
  }
}

/*-------------------------------------------------*/
void SeamVector::fromvector(const NodeVector& u)
{
  // std::cerr<<"_n="<<_n.t()<<"\n" << "u.n()"<<u.n()<<"\n";
  // std::cerr << "this " << data().t() << "\n";
  // std::cerr << "u " << u.data().t() << "\n";
  // data().fill(0);
  armaicvec n = u.n();
  if (dim() == 2)
  {
    for (int ix = 0; ix < n[0]; ix++)
    {
      for (int iy = 0; iy < n[1]; iy++)
      {
        this->atp(ix, iy) = u.at(ix, iy);
      }
    }
  }
  else if (dim() == 3)
  {
    for (int ix = 0; ix < n[0]; ix++)
    {
      for (int iy = 0; iy < n[1]; iy++)
      {
        for (int iz = 0; iz < n[2]; iz++)
        {
          this->atp(ix, iy, iz) = u.at(ix, iy, iz);
        }
      }
    }
  }
  // std::cerr << "this " << data().t() << "\n";
}

//
//
// /*-------------------------------------------------*/
// void SeamVector::fill_bdry(double d)
// {
//   if(dim()==2)
//   {
//     int nx = _n[0], ny = _n[1];
//     for(int ix=0;ix<nx;ix++)
//     {
//       this->at(ix,0)    = d;
//       this->at(ix,ny-1) = d;
//     }
//     for(int iy=0;iy<ny;iy++)
//     {
//       this->at(0,iy)    = d;
//       this->at(nx-1,iy) = d;
//     }
//   }
//   else if(dim()==3)
//   {
//     int nx = _n[0], ny = _n[1], nz = _n[2];
//     for(int ix=0;ix<nx;ix++)
//     {
//       for(int iy=0;iy<ny;iy++)
//       {
//         this->at(ix,iy,0)    = d;
//         this->at(ix,iy,nz-1) = d;
//       }
//     }
//     for(int ix=0;ix<nx;ix++)
//     {
//       for(int iz=0;iz<nz;iz++)
//       {
//         this->at(ix,0,   iz) = d;
//         this->at(ix,ny-1,iz) = d;
//       }
//     }
//     for(int iy=0;iy<ny;iy++)
//     {
//       for(int iz=0;iz<nz;iz++)
//       {
//         this->at(0,   iy,iz) = d;
//         this->at(nx-1,iy,iz) = d;
//       }
//     }
//   }
// }
