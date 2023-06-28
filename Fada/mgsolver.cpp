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
void MgSolver::set_parameters(int maxiter, double tol_rel, double tol_abs)
{
  _maxiter = maxiter;
  _tol_rel = tol_rel;
  _tol_abs = tol_abs;
}
/*-------------------------------------------------*/
void MgSolver::set_sizes(std::shared_ptr<MultiGridInterface> mgrid, std::shared_ptr<ModelInterface> fem, std::string smoothertype, int updatemem)
{
  _timer.enrol("smooth");
  _timer.enrol("transfer");
  _timer.enrol("residual");
  _timer.enrol("update");
  _timer.enrol("solvecoarse");
  //  std::cerr << "mggrid" << *mggrid << "\n";
  _nlevels = mgrid->nlevels();
  _model = fem;
  _mgmem.resize(4);
  for(int i=0;i<_mgmem.size();i++)
  {
    _set_size_vectormg(mgrid, _mgmem[i]);
  }
  _mgmatrix.resize(_nlevels);
  for(int l=0;l<_nlevels;l++)
  {
    _mgmatrix[l] = _model->newMatrix(mgrid->get(l));
  }
  _mgsmoother.resize(_nlevels);
  for(int l=0;l<_nlevels;l++)
  {
    if(l<_nlevels-1)
    {
      _mgsmoother[l] = _model->newSmoother(smoothertype, mgrid->get(l), _mgmatrix[l]);
    }
    else
    {
      _mgsmoother[l] = _model->newCoarseSolver("direct", mgrid->get(l), _mgmatrix[l]);
    }
    // _mgsmoother[l]->set_matrix(_mgmatrix[l]);
  }
  _mgtransfer.resize(_nlevels-1);
  for(int l=0;l<_nlevels-1;l++)
  {
    _mgtransfer[l] = _model->newTransfer(mgrid->get(l+1));
  }
  _mgupdate.resize(_nlevels);
  _mgupdatesmooth.resize(_nlevels);
  for(int l=0;l<_nlevels;l++)
  {
    if(updatemem==0)
    {
      // _mgupdate[l] = std::shared_ptr<UpdaterInterface>(new UpdaterConstant);
      _mgupdate[l] = std::shared_ptr<UpdaterInterface>(new UpdaterSimple);
    }
    else
    {
      _mgupdate[l] = std::shared_ptr<UpdaterInterface>(new Updater);
    }
    // _mgupdatesmooth[l] = std::shared_ptr<UpdaterInterface>(new UpdaterSimple);
    _mgupdatesmooth[l] = std::shared_ptr<UpdaterInterface>(new UpdaterConstant(0.8));
  }
  for(int l=0;l<_nlevels;l++)
  {
    _mgupdate[l]->setParameters(*_model, mgrid->get(l), _mgmatrix[l], updatemem);
    _mgupdatesmooth[l]->setParameters(*_model, mgrid->get(l), _mgmatrix[l], 0);
  }
}
/*-------------------------------------------------*/
void MgSolver::_set_size_vectormg(std::shared_ptr<MultiGridInterface> mgrid, VectorMG& v) const
{
  v.resize(_nlevels);
  for(int l=0;l<_nlevels;l++)
  {
    v[l] = _model->newVector(mgrid->get(l));
  }
}

/*-------------------------------------------------*/
void MgSolver::residual(int l, std::shared_ptr<VectorInterface> r, std::shared_ptr<VectorInterface const> u, std::shared_ptr<VectorInterface const>  f) const
{
  // std::cerr << "\nu\n";
  // u->save(std::cerr);

  r->equal(*f);
  _mgmatrix[l]->dot(r, u, -1.0);

  // std::cerr << "\nr\n";
  // r->save(std::cerr);
}
/*-------------------------------------------------*/
// int MgSolver::solve(VectorInterface& u, const VectorInterface& f, bool print)
int MgSolver::solve(bool print)
{
  VectorMG& umg = _mgmem[0];
  VectorMG& fmg = _mgmem[1];
  VectorMG& d   = _mgmem[2];
  VectorMG& w   = _mgmem[3];

  // _model->vector2vectormg(*fmg[0], f);
  // _model->vector2vectormg(*umg[0], u);
  // fmg[0]->equal(f);
  // umg[0]->equal(u);

  int maxlevel = 0;
  double res, tol=0;
  for(int iter=0; iter<this->_maxiter+1; iter++)
  {
    _timer.start("residual");
    // residual(maxlevel, *d[maxlevel], *umg[maxlevel], *fmg[maxlevel]);
    residual(maxlevel, d[maxlevel], umg[maxlevel], fmg[maxlevel]);
    d[maxlevel]->fill_bdry(0);
    _timer.stop("residual");
    res = d[maxlevel]->norm();
    if(iter==0)
    {
      tol = fmax(this->_tol_abs, this->_tol_rel*res);
      if(print) printf("-mg- ---tol= %10.3e ---\n", tol);
    }
    if(print) printf("-mg- %3d %10.3e\n", iter, res);
    if(res <= tol)
    {
      // _model->vectormg2vector(u, *umg[0]);
      // u.equal(*umg[0]);
      return iter;
    }
    mgstep(maxlevel, umg, fmg, d, w, tol);
  }
  return -1;
}
/*-------------------------------------------------*/
void MgSolver::mgstep(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, double tol)
{
  // std::cerr << "MgSolver::mgstep() l=" << l << "_nlevels="<<_nlevels<<"\n";
  if(l==_nlevels-1)
  {
    _timer.start("solvecoarse");
    // std::cerr << "\n BEFORE u\n";
    // u[l]->save(std::cerr);
    // std::cerr << "\n BEFORE f\n";
    // f[l]->save(std::cerr);
    _mgsmoother[l]->solve(u[l], f[l]);
    // std::cerr << "\n AFTER u\n";
    // u[l]->save(std::cerr);
    _timer.stop("solvecoarse");
  }
  else
  {
    _timer.start("smooth");
    _mgsmoother[l]->pre(w[l], d[l]);
    _timer.stop("smooth");
    _timer.start("update");
    //    _mgupdatesmooth[l]->addUpdate(w[l], u[l], d[l]);
    _mgupdate[l]->addUpdate(w[l], u[l], d[l]);
    _timer.stop("update");

    _timer.start("transfer");
    _mgtransfer[l]->restrict(d[l+1], d[l]);
    _timer.stop("transfer");

    f[l+1]->equal(*d[l+1]);
    u[l+1]->fill(0.0);
    mgstep(l+1, u, f, d, w, tol);

    _timer.start("transfer");
    _mgtransfer[l]->prolongate(w[l], u[l+1]);
    _timer.stop("transfer");
    // _timer.start("residual");
    // residual(l, *d[l], *u[l], *f[l]);
    // _timer.stop("residual");
    _timer.start("update");
    _mgupdate[l]->addUpdate(w[l], u[l], d[l]);
    _timer.stop("update");
    _timer.start("smooth");
    _mgsmoother[l]->post(w[l], d[l]);
    _timer.stop("smooth");
    _timer.start("update");
    _mgupdate[l]->addUpdate(w[l], u[l], d[l]);
    //    _mgupdatesmooth[l]->addUpdate(w[l], u[l], d[l]);
    _timer.stop("update");
  }
}
