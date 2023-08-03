//
//  gridvector.cpp
//  Fada
//
//  Created by Roland Becker on 03/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <stdio.h>
#include  "gridvector.hpp"
#include  "../uniformgrid.hpp"

/*-------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const GridVector& v)
{
    const armavec& tarma =static_cast<const armavec&>(v);
    os << tarma.t()<< "n=" << v.n().t()<< "ofs=" << v.ofs().t();
    // os << v.data().t()<< "n=" << v.n().t()<< "ofs=" << v.ofs().t();
    return os;
}

/*-------------------------------------------------*/
void GridVector::after_prolongate()
{
//     std::cerr << "** after_prolongate\n";
//     if(_mean)
//     {
//         project_mean();
//         return;
//     }
}

/*-------------------------------------------------*/
void GridVector::after_restrict()
{
    // std::cerr << "** after_restrict " << "\n";
    if(_mean)
    {
        project_mean();
        return;
    }
    boundary_zero();
}

/*-------------------------------------------------*/
void GridVector::after_smooth()
{
    // std::cerr << "** after_smooth\n";
    if(_mean)
    {
        project_mean();
    }    
}

/*-------------------------------------------------*/
void GridVector::project_mean()
{
    armavec& tarma = static_cast<armavec&>(*this);
    double mean = arma::accu(tarma)/this->n_elem;
    // std::cerr << "mean="<<mean<<"\n";
    *this -= mean;
    // mean = arma::accu(tarma)/this->n_elem;
    // std::cerr << "mean="<<mean<<"\n";
}

/*-------------------------------------------------*/
void GridVector::boundary_zero()
{
    if(dim()==2)
    {
        int nx = _n[0], ny = _n[1];
        if(_bc[0][0])
        {
            for(int iy=0;iy<ny;iy++)
            {
                this->at(0,iy)    = 0;
            }
          
        }
        if(_bc[0][1])
        {
            for(int iy=0;iy<ny;iy++)
            {
                this->at(nx-1,iy) = 0;
            }          
        }
        if(_bc[1][0])
        {
            for(int ix=0;ix<nx;ix++)
            {
                this->at(ix,0)    = 0;
            }
          
        }
        if(_bc[1][1])
        {
            for(int ix=0;ix<nx;ix++)
            {
                this->at(ix,ny-1) = 0;
            }
          
        }
    }
    else if(dim()==3)
    {
        int nx = _n[0], ny = _n[1], nz = _n[2];
        if(_bc[0][0])
        {
            for(int iy=0;iy<ny;iy++)
            {
                for(int iz=0;iz<nz;iz++)
                {
                    this->at(0,   iy,iz) = 0;
                }
            }
        }
        if(_bc[0][1])
        {
            for(int iy=0;iy<ny;iy++)
            {
                for(int iz=0;iz<nz;iz++)
                {
                    this->at(nx-1,iy,iz) = 0;
                }
            }
        }
        if(_bc[1][0])
        {
            for(int ix=0;ix<nx;ix++)
            {
                for(int iz=0;iz<nz;iz++)
                {
                    this->at(ix,0,   iz) = 0;
                }
            }
        }
        if(_bc[1][1])
        {
            for(int ix=0;ix<nx;ix++)
            {
                for(int iz=0;iz<nz;iz++)
                {
                    this->at(ix,ny-1,iz) = 0;
                }
            }
        }
        if(_bc[2][0])
        {
            for(int ix=0;ix<nx;ix++)
            {
                for(int iy=0;iy<ny;iy++)
                {
                    this->at(ix,iy,0)    = 0;
                }
            }
        }
        if(_bc[2][1])
        {
            for(int ix=0;ix<nx;ix++)
            {
                for(int iy=0;iy<ny;iy++)
                {
                    this->at(ix,iy,nz-1) = 0;
                }
            }
        }
    }
}
