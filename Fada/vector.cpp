//
//  vector.c
//  Fada
//
//  Created by Roland Becker on 03/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <stdio.h>
#include  "vector.hpp"

/*-------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const Vector& v)
{
  const armavec& tarma =static_cast<const armavec&>(v);
  os << tarma.t()<< "n=" << v.n().t()<< "ofs=" << v.ofs().t();
  return os;
}

/*-------------------------------------------------*/
void Vector::fill_bdry(double d)
{
  if(dim()==2)
  {
    int nx = _n[0], ny = _n[1];
    for(int ix=0;ix<nx;ix++)
    {
      this->at(ix,0)    = d;
      this->at(ix,ny-1) = d;
    }
    for(int iy=0;iy<ny;iy++)
    {
      this->at(0,iy)    = d;
      this->at(nx-1,iy) = d;
    }
  }
  else if(dim()==3)
  {
    int nx = _n[0], ny = _n[1], nz = _n[2];
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        this->at(ix,iy,0)    = d;
        this->at(ix,iy,nz-1) = d;
      }
    }
    for(int ix=0;ix<nx;ix++)
    {
      for(int iz=0;iz<nz;iz++)
      {
        this->at(ix,0,   iz) = d;
        this->at(ix,ny-1,iz) = d;
      }
    }
    for(int iy=0;iy<ny;iy++)
    {
      for(int iz=0;iz<nz;iz++)
      {
        this->at(0,   iy,iz) = d;
        this->at(nx-1,iy,iz) = d;
      }
    }
  }
}

/*-------------------------------------------------*/
void Vector::fill_bdry2(double d)
{
  if(dim()==2)
  {
    int nx = _n[0], ny = _n[1];
    std::cerr << nx << " " << ny << "\n";
    for(int ix=1;ix<nx-1;ix++)
    {
      this->at(ix,1)    = d;
      this->at(ix,ny-2) = d;
    }
    for(int iy=1;iy<ny-1;iy++)
    {
      this->at(1,iy)    = d;
      this->at(nx-2,iy) = d;
    }
  }
  else if(dim()==3)
  {
    int nx = _n[0], ny = _n[1], nz = _n[2];
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iy=1;iy<ny-1;iy++)
      {
        this->at(ix,iy,1)    = d;
        this->at(ix,iy,nz-2) = d;
      }
    }
    for(int ix=1;ix<nx-1;ix++)
    {
      for(int iz=1;iz<nz-1;iz++)
      {
        this->at(ix,1,   iz) = d;
        this->at(ix,ny-2,iz) = d;
      }
    }
    for(int iy=1;iy<ny-1;iy++)
    {
      for(int iz=1;iz<nz-1;iz++)
      {
        this->at(1,   iy,iz) = d;
        this->at(nx-2,iy,iz) = d;
      }
    }
  }
}
