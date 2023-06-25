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
#include  "../uniformgrid.hpp"
#include  "nodevector.hpp"
#include  "../matrixinterface.hpp"
#include  "transferq1.hpp"
#include  "../smoothersimple.hpp"
#include  "../smootherumf.hpp"
#include  "stencil2d.hpp"
#include  "stencil3d.hpp"

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
  // p->set_size(grid.n()+2);
  p->set_size(grid.n());
  p->fill_bdry(0);
  return p;
}
/*-------------------------------------------------*/
std::unique_ptr<VectorInterface> Q13d::newMgvector(const GridInterface& grid) const
{
  std::unique_ptr<VectorInterface> p = std::make_unique<NodeVector>();
  // p->set_size(grid.n()+2);
  p->set_size(grid.n());
  p->fill_bdry(0);
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
    arma::vec::fixed<9> coef;
    coef[0] = -1.0/3.0;
    coef[1] = -1.0/3.0;
    coef[2] = -1.0/3.0;
    coef[3] = -1.0/3.0;
    coef[4] =  8.0/3.0;
    coef[5] = -1.0/3.0;
    coef[6] = -1.0/3.0;
    coef[7] = -1.0/3.0;
    coef[8] = -1.0/3.0;
    return std::unique_ptr<MatrixInterface>(new Matrix<Stencil2d9,NodeVector>(n, coef));
  }
  else if(_type=="Trapez")
  {
    arma::vec::fixed<5> coef;
    coef[0] = -1.0;
    coef[1] = -1.0;
    coef[2] =  4.0;
    coef[3] = -1.0;
    coef[4] = -1.0;
    return std::unique_ptr<MatrixInterface>(new Matrix<Stencil2d5,NodeVector>(n, coef));
    // return std::unique_ptr<MatrixInterface>(new Matrix<TrapezMatrix2d,NodeVector>(n, dx));
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
  // return std::unique_ptr<SmootherInterface>(new SmootherSimple("Jac"));
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
    assert(dx[0]==dx[1]);
    double e = dx[0];
    double d0 = 8.0/3.0 * e;
    double d1 = -0.0/3.0 * e;
    double d2 = -1.0/6.0 * e;
    double d3 = -1.0/12.0 * e;
    arma::vec::fixed<27> coef;

    coef[13] = d0;

    coef[ 4] = coef[10] = coef[12] = coef[14] = coef[16] = coef[22] = d1;

    coef[ 1] = coef[ 3] = coef[ 5] = coef[ 7] = coef[ 9] = coef[11] = d2;
    coef[15] = coef[17] = coef[19] = coef[21] = coef[23] = coef[25] = d2;

    coef[ 0] = coef[ 2] = coef[ 6] = coef[ 8] = coef[18] = coef[20] = coef[24] = coef[26] = d3;
    return std::unique_ptr<MatrixInterface>(new Matrix<Stencil3d27,NodeVector>(n, coef));
//    return std::unique_ptr<MatrixInterface>(new Matrix<FullMatrix3d,NodeVector>(n, dx));
  }
  else if(_type=="Trapez")
  {
    assert(dx[0]==dx[1]);
    assert(dx[0]==dx[2]);
    arma::vec::fixed<7> coef;
    coef[0] = -1.0;
    coef[1] = -1.0;
    coef[2] = -1.0;
    coef[3] = 6.0;
    coef[4] = -1.0;
    coef[5] = -1.0;
    coef[6] = -1.0;
    coef * dx[0];
    return std::unique_ptr<MatrixInterface>(new Matrix<Stencil3d7,NodeVector>(n, coef));
//    return std::unique_ptr<MatrixInterface>(new TrapezMatrix3d(n,dx));
    // return std::unique_ptr<MatrixInterface>(new Matrix<TrapezMatrix3d,NodeVector>(n, dx));
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
      // u.at(ix,iy) = umg.atp(ix,iy);
      u.at(ix,iy) = umg.at(ix,iy);
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
      // umg.atp(ix,iy) = u.at(ix,iy);
      umg.at(ix,iy) = u.at(ix,iy);
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
        // u.at(ix,iy,iz) = umg.atp(ix,iy,iz);
        u.at(ix,iy,iz) = umg.at(ix,iy,iz);
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
        // umg.atp(ix,iy,iz) = u.at(ix,iy,iz);
        umg.at(ix,iy,iz) = u.at(ix,iy,iz);
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
