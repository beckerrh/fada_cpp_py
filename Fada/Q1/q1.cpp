//
//  q1.cpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <cassert>
#include  <armadillo>
#include  "../matrixinterface.hpp"
#include  "../uniformgrid.hpp"
#include  "../smoothersimple.hpp"
#include  "../smootherumf.hpp"
#include  "q1.hpp"
#include  "nodevector.hpp"
#include  "transferq1.hpp"
#include  "stencil2d.hpp"
#include  "stencil3d.hpp"

double lin2d(double x, double y) {return 3.0*x+2.0*y;}
double lin3d(double x, double y, double z) {return 3.0*x+2.0*y+z;}


/*-------------------------------------------------*/
void Q1::set_grid(std::shared_ptr<GridInterface const> grid)
{
  std::shared_ptr<UniformGrid const> ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
  assert(ug);
  assert(ug->dx().n_elem == ug->dim());
  assert(ug->n().n_elem == ug->dim());
  _nx = ug->nx();
  _ny = ug->ny();
  _nz = -1;
  if(ug->dim()==3) _nz = ug->nz();
  _vol = arma::prod(ug->dx());
  _ug = ug;
}
/*-------------------------------------------------*/
std::shared_ptr<VectorInterface> Q1::newVector(std::shared_ptr<GridInterface const> grid) const
{
  // std::shared_ptr<VectorInterface> p = std::make_unique<NodeVector>();
  std::shared_ptr<VectorInterface> p = std::make_unique<Vector<NodeVector>>();
  // p->set_size(grid.n()+2);
  p->set_size(grid->n());
  p->fill_bdry(0);
  return p;
}
/*-------------------------------------------------*/
std::shared_ptr<SmootherInterface> Q1::newSmoother(std::string type, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const
{
  return std::shared_ptr<SmootherInterface>(new SmootherSimple(type, matrix));
}
/*-------------------------------------------------*/
std::shared_ptr<SmootherInterface> Q1::newCoarseSolver(std::string type, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const
{
  return std::shared_ptr<SmootherInterface>(new SmootherSimple("GS", matrix));
  std::shared_ptr<MatrixInterface const> matrixforumf;
  std::shared_ptr<FemAndMatrixInterface const> stencil = std::dynamic_pointer_cast<FemAndMatrixInterface const>(matrix);
  if(stencil)
  {
    std::shared_ptr<FemAndMatrixInterface const> stencil = std::dynamic_pointer_cast<FemAndMatrixInterface const>(matrix);
    arma::umat locations;
    armavec values;
    stencil->get_locations_values(locations, values);
    typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;
    matrixforumf = std::make_shared<SparseMatrixDef>(locations, values);
  }
  else
  {
    matrixforumf = matrix;
  }
  // return std::shared_ptr<SmootherInterface>(new Smoother<arma::sp_mat, Vector<armavec>>(matrix));
  return std::shared_ptr<SmootherInterface>(new Smoother<SmootherUmf, Vector<armavec>>(matrixforumf));
}

/*-------------------------------------------------*/
std::shared_ptr<MatrixInterface> Q1::newMatrix(std::shared_ptr<GridInterface const> grid) const
{
  std::shared_ptr<FemAndMatrixInterface> p = newStencil(grid);
  if(_matrixtype=="stencil")
  {
    return p;
  }
  else
  {
    arma::umat locations;
    armavec values;
    p->get_locations_values(locations, values);
    // return std::shared_ptr<MatrixInterface>(new SparseMatrix(locations, values));
    typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;
    return std::shared_ptr<MatrixInterface>(new SparseMatrixDef(locations, values));
  }
}
/*-------------------------------------------------*/
void Q1::rhs_one(NodeVector& v) const
{
  v.fill(_vol*12.0);
}
/*-------------------------------------------------*/
void Q1::rhs_random(NodeVector& v) const
{
  arma::arma_rng::set_seed_random();
  // v.data().randu();
  v.randu();
  v *= 100*_vol;
}
/*-------------------------------------------------*/
std::shared_ptr<FemAndMatrixInterface> Q12d::newStencil(std::shared_ptr<GridInterface const> grid) const
{
  std::shared_ptr<UniformGrid const> ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
  assert(ug);
  const armaicvec& n = ug->n();
  const armavec& dx = ug->dx();
  if(_stenciltype=="Full")
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
    // return std::shared_ptr<MatrixInterface>(new Matrix<Stencil2d9,Vector<NodeVector>>(n, coef));
    return std::shared_ptr<FemAndMatrixInterface>(new FemAndMatrix<Stencil2d9,Vector<NodeVector>>(n, coef));
  }
  else if(_stenciltype=="Trapez")
  {
    arma::vec::fixed<5> coef;
    coef[0] = -1.0;
    coef[1] = -1.0;
    coef[2] =  4.0;
    coef[3] = -1.0;
    coef[4] = -1.0;
    // return std::shared_ptr<MatrixInterface>(new Matrix<Stencil2d5,Vector<NodeVector>>(n, coef));
    return std::shared_ptr<FemAndMatrixInterface>(new FemAndMatrix<Stencil2d5,Vector<NodeVector>>(n, coef));
  }
  else
  {
    std::cerr << "unknown matrix '" << _stenciltype<<"'\n";
    std::exit(1);
  }
}

/*-------------------------------------------------*/
std::shared_ptr<TransferInterface> Q12d::newTransfer(std::shared_ptr<GridInterface const> grid) const
{
  std::shared_ptr<UniformGrid const> ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
  // const UniformGrid* ug = dynamic_cast<const UniformGrid*>(grid);
  assert(ug);
  const armaicvec& n = ug->n();
  const armavec& dx = ug->dx();
// return std::shared_ptr<TransferInterface>(new Transfer<TransferQ12d,NodeVector>(n, dx));
return std::shared_ptr<TransferInterface>(new Transfer<TransferQ12d,Vector<NodeVector>>(n, dx));
}

/*-------------------------------------------------*/
std::shared_ptr<FemAndMatrixInterface> Q13d::newStencil(std::shared_ptr<GridInterface const> grid) const
{
  std::shared_ptr<UniformGrid const> ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
  // const UniformGrid* ug = dynamic_cast<const UniformGrid*>(grid);
  assert(ug);
  const armaicvec& n = ug->n();
  const armavec& dx = ug->dx();
  if(_stenciltype=="Full")
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
    // return std::shared_ptr<MatrixInterface>(new Matrix<Stencil3d27,Vector<NodeVector>>(n, coef));
    return std::shared_ptr<FemAndMatrixInterface>(new FemAndMatrix<Stencil3d27,Vector<NodeVector>>(n, coef));
  }
  else if(_stenciltype=="Trapez")
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
    // return std::shared_ptr<MatrixInterface>(new Matrix<Stencil3d7,Vector<NodeVector>>(n, coef));
    return std::shared_ptr<FemAndMatrixInterface>(new FemAndMatrix<Stencil3d7,Vector<NodeVector>>(n, coef));
  }
  else
  {
    std::cerr << "unknown matrix '" << _stenciltype<<"'\n";
    std::exit(1);
  }
}

/*-------------------------------------------------*/
std::shared_ptr<TransferInterface> Q13d::newTransfer(std::shared_ptr<GridInterface const> grid) const
{
  std::shared_ptr<UniformGrid const> ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
  // const UniformGrid* ug = dynamic_cast<const UniformGrid*>(grid);
  assert(ug);
  const armaicvec& n = ug->n();
  const armavec& dx = ug->dx();
  // return std::shared_ptr<TransferInterface>(new Transfer<TransferQ13d,NodeVector>(n, dx));
  return std::shared_ptr<TransferInterface>(new Transfer<TransferQ13d,Vector<NodeVector>>(n, dx));
}


/*-------------------------------------------------*/
void Q12d::boundary(NodeVector& u) const
{
  for(int ix=0;ix<_nx;ix++)
  {
    u.at(ix,0)    = 0;
    u.at(ix,_ny-1) = 0;
  }
  for(int iy=0;iy<_ny;iy++)
  {
    u.at(0,iy)    = 0;
    u.at(_nx-1,iy) = 0;
  }
}


/*-------------------------------------------------*/
void Q13d::boundary(NodeVector& u) const
{
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iy=0;iy<_ny;iy++)
    {
      u.at(ix,iy,0)     = 0;
      u.at(ix,iy,_nz-1) = 0;
    }
  }
  for(int ix=0;ix<_nx;ix++)
  {
    for(int iz=0;iz<_nz;iz++)
    {
      u.at(ix,0,   iz)   = 0;
      u.at(ix,_ny-1,iz)  = 0;
    }
  }
  for(int iy=0;iy<_ny;iy++)
  {
    for(int iz=0;iz<_nz;iz++)
    {
      u.at(0,   iy,iz)  = 0;
      u.at(_nx-1,iy,iz) = 0;
    }
  }
}
