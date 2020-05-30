#include  <math.h>
#include  <stdio.h>
#include  <armadillo>
#include  <algorithm>
#include  "mgsolver.hpp"
#include  "updater.hpp"
#include  "uniformgrid.hpp"

/*-------------------------------------------------*/
std::string MgSolver::toString() const
{
  std::stringstream ss;
  ss << "nlevels = " << _nlevels;
  return ss.str();
}
/*-------------------------------------------------*/
void MgSolver::set_parameters()
{
  maxiter = 100;
  tol_rel = 1e-8;
  tol_abs = 1e-12;
}
/*-------------------------------------------------*/
void MgSolver::set_sizes(std::shared_ptr<MultiGridInterface> mgrid, std::shared_ptr<FiniteElementInterface> fem, std::string smoothertype, int updatemem)
{
  _timer.enrol("smooth");
  _timer.enrol("transfer");
  _timer.enrol("residual");
  _timer.enrol("update");
  _timer.enrol("solvecoarse");
  //  std::cerr << "mggrid" << *mggrid << "\n";
  _nlevels = mgrid->nlevels();
  _fem = fem;
  _mgmem.resize(4);
  for(int i=0;i<_mgmem.size();i++)
  {
    _set_size_vectormg(mgrid, _mgmem[i]);
  }
  _mgmatrix.resize(_nlevels);
  for(int l=0;l<_nlevels;l++)
  {
    _mgmatrix[l] = _fem->newMatrix(*mgrid->get(l));
//    _mgmatrix[l]->set_grid(mgrid->get(l));
  }
  _mgsmoother.resize(_nlevels);
  for(int l=0;l<_nlevels;l++)
  {
    if(l<_nlevels-1)
    {
      _mgsmoother[l] = _fem->newSmoother(smoothertype, *mgrid->get(l));
    }
    else
    {
      _mgsmoother[l] = _fem->newCoarseSolver("direct", *mgrid->get(l));
    }
    _mgsmoother[l]->set_matrix(_mgmatrix[l]);
  }
  _mgtransfer.resize(_nlevels-1);
  for(int l=0;l<_nlevels-1;l++)
  {
    _mgtransfer[l] = _fem->newTransfer(*mgrid->get(l+1));
//    _mgtransfer[l]->set_grid(mgrid->get(l+1));
  }
  _mgupdate.resize(_nlevels);
  _mgupdatesmooth.resize(_nlevels);
  std::string type="cyc";
  //  type="coef";
  //  type="ortho";
  //  type="restart";
  for(int l=0;l<_nlevels;l++)
  {
    if(updatemem==0)
    {
      _mgupdate[l] = std::unique_ptr<UpdaterInterface>(new UpdaterSimple);
    }
    else
    {
      _mgupdate[l] = std::unique_ptr<UpdaterInterface>(new Updater);
    }
    _mgupdatesmooth[l] = std::unique_ptr<UpdaterInterface>(new UpdaterSimple);
  }
  for(int l=0;l<_nlevels;l++)
  {
    _mgupdate[l]->setParameters(*_fem, *mgrid->get(l), _mgmatrix[l], updatemem, type);
    //    _mgupdate(l)->setParameters(l, this, updatemem, type);
//    _mgupdate[l]->set_size(mgrid->get(l)->n()+2);
    //    _mgupdatesmooth(l)->setParameters(l, this, 0, type);
    _mgupdatesmooth[l]->setParameters(*_fem, *mgrid->get(l), _mgmatrix[l], 0, type);
//    _mgupdatesmooth[l]->set_size(mgrid->get(l)->n()+2);
  }
}
/*-------------------------------------------------*/
void MgSolver::_set_size_vectormg(std::shared_ptr<MultiGridInterface> mgrid, VectorMG& v) const
{
  v.resize(_nlevels);
  for(int l=0;l<_nlevels;l++)
  {
    v[l] = _fem->newMgvector(*mgrid->get(l));
  }
}

/*-------------------------------------------------*/
void MgSolver::residual(int l, VectorInterface& r, const VectorInterface& u, const VectorInterface& f) const
{
//  r =  f;
  r.equal(f);
  _mgmatrix[l]->dot(r, u, -1.0);
}
/*-------------------------------------------------*/
int MgSolver::solve(VectorInterface& u, const VectorInterface& f, bool print)
{
  VectorMG& umg = _mgmem[0];
  VectorMG& fmg = _mgmem[1];
  VectorMG& d   = _mgmem[2];
  VectorMG& w   = _mgmem[3];
  
  _fem->vector2vectormg(*fmg[0], f);
  _fem->vector2vectormg(*umg[0], u);
  
  int maxlevel = 0;
  double res, tol=0;
  for(int iter=0; iter<this->maxiter+1; iter++)
  {
    _timer.start("residual");
    residual(maxlevel, *d[maxlevel], *umg[maxlevel], *fmg[maxlevel]);
    d[maxlevel]->fill_bdry(0);
    d[maxlevel]->fill_bdry2(0);
    _timer.stop("residual");
    res = d[maxlevel]->norm();
    if(iter==0)
    {
      tol = fmax(this->tol_abs, this->tol_rel*res);
      if(print) printf("-mg- --- %10.3e ---\n", tol);
    }
    if(print) printf("-mg- %3d %10.3e\n", iter, res);
    if(res <= tol)
    {
      _fem->vectormg2vector(u, *umg[0]);
      return iter;
    }
    mgstep(maxlevel, umg, fmg, d, w, tol);
  }
  return -1;
}
/*-------------------------------------------------*/
void MgSolver::mgstep(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, double tol)
{
  if(l==_nlevels-1)
  {
    _timer.start("solvecoarse");
    _mgsmoother[l]->solve(*u[l], *f[l]);
    _timer.stop("solvecoarse");
  }
  else
  {
    _timer.start("smooth");
    _mgsmoother[l]->pre(*w[l], *d[l]);
    _timer.stop("smooth");
    _timer.start("update");
    //    _mgupdatesmooth[l]->addUpdate(w[l], u[l], d[l]);
    _mgupdate[l]->addUpdate(*w[l], *u[l], *d[l]);
    _timer.stop("update");
    
    _timer.start("transfer");
    _mgtransfer[l]->restrict(*d[l+1], *d[l]);
    _timer.stop("transfer");
    
//    f[l+1] = d[l+1];
    f[l+1]->equal(*d[l+1]);
    u[l+1]->fill(0.0);
    mgstep(l+1, u, f, d, w, tol);
    
    _timer.start("transfer");
    _mgtransfer[l]->prolongate(*w[l], *u[l+1]);
    _timer.stop("transfer");
    _timer.start("residual");
    residual(l, *d[l], *u[l], *f[l]);
    _timer.stop("residual");
    _timer.start("update");
    _mgupdate[l]->addUpdate(*w[l], *u[l], *d[l]);
    _timer.stop("update");
    _timer.start("smooth");
    _mgsmoother[l]->post(*w[l], *d[l]);
    _timer.stop("smooth");
    _timer.start("update");
    _mgupdate[l]->addUpdate(*w[l], *u[l], *d[l]);
    //    _mgupdatesmooth[l]->addUpdate(w[l], u[l], d[l]);
    _timer.stop("update");
  }
}
