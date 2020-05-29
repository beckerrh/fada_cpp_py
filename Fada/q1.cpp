//
//  q1.cpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <cassert>
#include  <armadillo>
#include  "q1.hpp"
#include  "uniformgrid.hpp"
#include  "vector.hpp"
#include  "fullmatrix.hpp"
#include  "trapezmatrix.hpp"
#include  "transferq1.hpp"
#include  "smoothersimple.hpp"
#include  "smootherumf.hpp"

double lin2d(double x, double y) {return 3.0*x+2.0*y;}
double lin3d(double x, double y, double z) {return 3.0*x+2.0*y+z;}

Q12d::~Q12d(){}
Q13d::~Q13d(){}

/*-------------------------------------------------*/
void Q12d::set_grid(std::shared_ptr<GridInterface> grid)
{
  std::shared_ptr<UniformGrid> ug = std::dynamic_pointer_cast<UniformGrid>(grid);
//  const UniformGrid* ug = dynamic_cast<const UniformGrid*>(grid);
  assert(ug);
  assert(ug->dim()==2);
  _nx = ug->nx();
  _ny = ug->ny();
  _vol = arma::prod(ug->dx());
  _ug = ug;
//  std::cerr << "Q12d::set_grid() " <<_nx << " " << _ny << "\n";
//  std::cerr << "_ug->n()="<<_ug->n().t();
}

/*-------------------------------------------------*/
void Q13d::set_grid(std::shared_ptr<GridInterface> grid)
{
  std::shared_ptr<UniformGrid> ug = std::dynamic_pointer_cast<UniformGrid>(grid);
  assert(ug);
  assert(ug->dim()==3);
  _nx = ug->nx();
  _ny = ug->ny();
  _nz = ug->nz();
  _vol = arma::prod(ug->dx());
  _ug = ug;
}
/*-------------------------------------------------*/
void Q12d::set_size_mgvector(const GridInterface& grid, Vector& u) const
{
  u.set_size(grid.n()+2);
  u.fill_bdry(0);
  u.fill_bdry2(0);
}
/*-------------------------------------------------*/
void Q13d::set_size_mgvector(const GridInterface& grid, Vector& u) const
{
  u.set_size(grid.n()+2);
  u.fill_bdry(0);
  u.fill_bdry2(0);
}

/*-------------------------------------------------*/
std::unique_ptr<MatrixInterface> Q12d::newMatrix() const
{
  if(_matrixtype=="Full")
  {
    return std::unique_ptr<MatrixInterface>(new FullMatrix2d);
  }
  else if(_matrixtype=="Trapez")
  {
    return std::unique_ptr<MatrixInterface>(new TrapezMatrix2d);
  }
  else
  {
    std::cerr << "unknown matrix '" << _matrixtype<<"'\n";
    std::exit(1);
  }
}
/*-------------------------------------------------*/
std::unique_ptr<SmootherInterface> Q12d::newSmoother(std::string type) const
{
  return std::unique_ptr<SmootherInterface>(new SmootherSimple(type));
}
/*-------------------------------------------------*/
std::unique_ptr<SmootherInterface> Q12d::newCoarseSolver(std::string type) const
{
  return std::unique_ptr<SmootherInterface>(new SmootherUmf);
}

/*-------------------------------------------------*/
std::unique_ptr<TransferInterface> Q12d::newTransfer() const
{
  return std::unique_ptr<TransferInterface>(new TransferQ12d);
}

/*-------------------------------------------------*/
std::unique_ptr<MatrixInterface> Q13d::newMatrix() const
{
  if(_matrixtype=="Full")
  {
    return std::unique_ptr<MatrixInterface>(new FullMatrix3d);
  }
  else if(_matrixtype=="Trapez")
  {
    return std::unique_ptr<MatrixInterface>(new TrapezMatrix3d);
  }
  else
  {
    std::cerr << "unknown matrix '" << _matrixtype<<"'\n";
    std::exit(1);
  }
}
/*-------------------------------------------------*/
std::unique_ptr<SmootherInterface> Q13d::newSmoother(std::string type) const
{
  return std::unique_ptr<SmootherInterface>(new SmootherSimple(type));
}
/*-------------------------------------------------*/
std::unique_ptr<SmootherInterface> Q13d::newCoarseSolver(std::string type) const
{
  return std::unique_ptr<SmootherInterface>(new SmootherUmf);
}

/*-------------------------------------------------*/
std::unique_ptr<TransferInterface> Q13d::newTransfer() const
{
  return std::unique_ptr<TransferInterface>(new TransferQ13d);
}

/*-------------------------------------------------*/
//void  Q12d::vectormg2vector(const UniformMultiGrid& mggrid, int l, Vector& u, const Vector& umg) const
void  Q12d::vectormg2vector(Vector& u, const Vector& umg) const
{
//  int nx = mggrid.nx(l), ny = mggrid.ny(l);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      u.at(ix,iy) = umg.atp(ix,iy);
    }
  }
}
/*-------------------------------------------------*/
//void  Q12d::vector2vectormg(const UniformMultiGrid& mggrid, int l, Vector& umg, const Vector& u) const
void  Q12d::vector2vectormg(Vector& umg, const Vector& u) const
{
//  int nx = mggrid.nx(l), ny = mggrid.ny(l);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      umg.atp(ix,iy) = u.at(ix,iy);
    }
  }
}

/*-------------------------------------------------*/
//void  Q13d::vectormg2vector(const UniformMultiGrid& mggrid, int l, Vector& u, const Vector& umg) const
void  Q13d::vectormg2vector(Vector& u, const Vector& umg) const
{
//  int nx = mggrid.nx(l), ny = mggrid.ny(l), nz = mggrid.nz(l);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        u.at(ix,iy,iz) = umg.atp(ix,iy,iz);
      }
    }
  }
}
/*-------------------------------------------------*/
//void  Q13d::vector2vectormg(const UniformMultiGrid& mggrid, int l, Vector& umg, const Vector& u) const
void  Q13d::vector2vectormg(Vector& umg, const Vector& u) const
{
//  int nx = mggrid.nx(l), ny = mggrid.ny(l), nz = mggrid.nz(l);
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        umg.atp(ix,iy,iz) = u.at(ix,iy,iz);
      }
    }
  }
}

/*-------------------------------------------------*/
void Q12d::rhs_one(Vector& v) const
{
//  std::cerr << "v.n()="<<v.n().t();
//  assert(_ug);
//  std::cerr << "_ug->n()="<<_ug->nx();
//  assert(arma::all(v.n()==_ug->n()));
  double d = _vol*12.0;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      v.at(ix,iy) = d;
    }
  }
}
/*-------------------------------------------------*/
void Q12d::rhs_random(Vector& v) const
{
  arma::arma_rng::set_seed_random();
  v.data().randu();
  v *= 100*_vol;
}
/*-------------------------------------------------*/
void Q12d::boundary(Vector& v) const
{
  for(int ix=0;ix<_nx;ix++)
  {
    v.at(ix,0)    = lin2d(_ug->x(ix,0), _ug->y(ix,0));
    v.at(ix,_ny-1) = lin2d(_ug->x(ix,_ny-1), _ug->y(ix,_ny-1));
  }
  for(int iy=0;iy<_ny;iy++)
  {
    v.at(0,iy)    = lin2d(_ug->x(0,iy), _ug->y(0,iy));
    v.at(_nx-1,iy) = lin2d(_ug->x(_nx-1,iy), _ug->y(_nx-1,iy));
  }
}

/*-------------------------------------------------*/
void Q13d::rhs_one(Vector& v) const
{
  double d = _vol*20.0;
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      for(int iz=0;iz<_nz;iz++)
      {
        v.at(ix,iy, iz) = d;
      }
    }
  }
}
/*-------------------------------------------------*/
void Q13d::rhs_random(Vector& v) const
{
  arma::arma_rng::set_seed_random();
  v.data().randu();
  v *= 100*_vol;
}
/*-------------------------------------------------*/
void Q13d::boundary(Vector& v) const
{
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      v.at(ix,iy,0)    = lin3d(_ug->x(ix,iy,0), _ug->y(ix,iy,0), _ug->z(ix,iy,0));
      v.at(ix,iy,_nz-1) = lin3d(_ug->x(ix,iy,_nz-1), _ug->y(ix,iy,_nz-1), _ug->z(ix,iy,_nz-1));
    }
  }
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iz=0;iz<_nz;iz++)
    {
      v.at(ix,0,   iz) = lin3d(_ug->x(ix,0,iz), _ug->y(ix,0,iz), _ug->z(ix,0,iz));
      v.at(ix,_ny-1,iz) = lin3d(_ug->x(ix,_ny-1,iz), _ug->y(ix,_ny-1,iz), _ug->z(ix,_ny-1,iz));
    }
  }
  for(int iy=0;iy<_ny;iy++)
  {
    for(int iz=0;iz<_nz;iz++)
    {
      v.at(0,   iy,iz) = lin3d(_ug->x(0,iy,iz), _ug->y(0,iy,iz), _ug->z(0,iy,iz));
      v.at(_nx-1,iy,iz) = lin3d(_ug->x(_nx-1,iy,iz), _ug->y(_nx-1,iy,iz), _ug->z(_nx-1,iy,iz));
    }
  }
}
