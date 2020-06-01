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
#include  "nodevector.hpp"
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
std::unique_ptr<VectorInterface> Q12d::newMgvector(const GridInterface& grid) const
{
  std::unique_ptr<VectorInterface> p = std::make_unique<NodeVector>();
  p->set_size(grid.n()+2);
  p->fill_bdry(0);
  p->fill_bdry2(0);
  return p;
}
/*-------------------------------------------------*/
std::unique_ptr<VectorInterface> Q13d::newMgvector(const GridInterface& grid) const
{
  std::unique_ptr<VectorInterface> p = std::make_unique<NodeVector>();
  p->set_size(grid.n()+2);
  p->fill_bdry(0);
  p->fill_bdry2(0);
  return p;
}

/*-------------------------------------------------*/
std::unique_ptr<MatrixInterface> Q12d::newMatrix(const GridInterface& grid) const
{
  const UniformGrid* ug = dynamic_cast<const UniformGrid*>(&grid);
  assert(ug);
  const armaicvec& n = ug->n();
  const armavec& dx = ug->dx();
  if(_type=="Full")
  {
    return std::unique_ptr<MatrixInterface>(new Matrix<FullMatrix2d,NodeVector>(n, dx));
  }
  else if(_type=="Trapez")
  {
//    return std::unique_ptr<MatrixInterface>(new TrapezMatrix2d(n, dx));
    return std::unique_ptr<MatrixInterface>(new Matrix<TrapezMatrix2d,NodeVector>(n, dx));
  }
  else
  {
    std::cerr << "unknown matrix '" << _type<<"'\n";
    std::exit(1);
  }
}
/*-------------------------------------------------*/
std::unique_ptr<SmootherInterface> Q12d::newSmoother(std::string type, const GridInterface& grid) const
{
  return std::unique_ptr<SmootherInterface>(new SmootherSimple(type));
}
/*-------------------------------------------------*/
std::unique_ptr<SmootherInterface> Q12d::newCoarseSolver(std::string type, const GridInterface& grid) const
{
  return std::unique_ptr<SmootherInterface>(new SmootherUmf);
}

/*-------------------------------------------------*/
std::unique_ptr<TransferInterface> Q12d::newTransfer(const GridInterface& grid) const
{
  const UniformGrid* ug = dynamic_cast<const UniformGrid*>(&grid);
  assert(ug);
  const armaicvec& n = ug->n();
  const armavec& dx = ug->dx();
//  return std::unique_ptr<TransferInterface>(new TransferQ12d(n, dx));
  return std::unique_ptr<TransferInterface>(new Transfer<TransferQ12d,NodeVector>(n, dx));
}

/*-------------------------------------------------*/
std::unique_ptr<MatrixInterface> Q13d::newMatrix(const GridInterface& grid) const
{
  const UniformGrid* ug = dynamic_cast<const UniformGrid*>(&grid);
  assert(ug);
  const armaicvec& n = ug->n();
  const armavec& dx = ug->dx();
  if(_type=="Full")
  {
//    return std::unique_ptr<MatrixInterface>(new FullMatrix3d(n,dx));
    return std::unique_ptr<MatrixInterface>(new Matrix<FullMatrix3d,NodeVector>(n, dx));
  }
  else if(_type=="Trapez")
  {
//    return std::unique_ptr<MatrixInterface>(new TrapezMatrix3d(n,dx));
    return std::unique_ptr<MatrixInterface>(new Matrix<TrapezMatrix3d,NodeVector>(n, dx));
  }
  else
  {
    std::cerr << "unknown matrix '" << _type<<"'\n";
    std::exit(1);
  }
}
/*-------------------------------------------------*/
std::unique_ptr<SmootherInterface> Q13d::newSmoother(std::string type, const GridInterface& grid) const
{
  return std::unique_ptr<SmootherInterface>(new SmootherSimple(type));
}
/*-------------------------------------------------*/
std::unique_ptr<SmootherInterface> Q13d::newCoarseSolver(std::string type, const GridInterface& grid) const
{
  return std::unique_ptr<SmootherInterface>(new SmootherUmf);
}

/*-------------------------------------------------*/
std::unique_ptr<TransferInterface> Q13d::newTransfer(const GridInterface& grid) const
{
  const UniformGrid* ug = dynamic_cast<const UniformGrid*>(&grid);
  assert(ug);
  const armaicvec& n = ug->n();
  const armavec& dx = ug->dx();
//  return std::unique_ptr<TransferInterface>(new TransferQ13d(n, dx));
  return std::unique_ptr<TransferInterface>(new Transfer<TransferQ13d,NodeVector>(n, dx));
}

/*-------------------------------------------------*/
void  Q12d::vectormg2vector(NodeVector& u, const NodeVector& umg) const
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
void  Q12d::vector2vectormg(NodeVector& umg, const NodeVector& u) const
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
void  Q13d::vectormg2vector(NodeVector& u, const NodeVector& umg) const
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
void  Q13d::vector2vectormg(NodeVector& umg, const NodeVector& u) const
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
void Q12d::rhs_one(NodeVector& v) const
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
void Q12d::rhs_random(NodeVector& v) const
{
  arma::arma_rng::set_seed_random();
  v.data().randu();
  v *= 100*_vol;
}
/*-------------------------------------------------*/
void Q12d::boundary(NodeVector& v) const
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
void Q13d::rhs_one(NodeVector& v) const
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
void Q13d::rhs_random(NodeVector& v) const
{
  arma::arma_rng::set_seed_random();
  v.data().randu();
  v *= 100*_vol;
}
/*-------------------------------------------------*/
void Q13d::boundary(NodeVector& v) const
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
