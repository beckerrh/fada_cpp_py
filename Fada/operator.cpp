#include  <math.h>
#include  <stdio.h>
#include  <armadillo>
#include  <algorithm>
#include  "q1.hpp"
#include  "operator.hpp"
#include  "updater.hpp"
#include  "uniformgrid.hpp"

/*-------------------------------------------------*/
void Operator::set_parameters()
{
  maxiter = 100;
  smoother = "jac";
  tol_rel = 1e-8;
  tol_abs = 1e-12;
  _optmem = 0;
}
/*-------------------------------------------------*/
Operator::~Operator() {}
/*-------------------------------------------------*/
Operator::Operator(bool printtimer) : _fem(nullptr), _timer(printtimer)
{
  set_parameters();
}
/*-------------------------------------------------*/
void Operator::_set_size(std::string femtype, std::string matrixtype)
{
  _timer.enrol("smooth");
  _timer.enrol("transfer");
  _timer.enrol("residual");
  _timer.enrol("update");
  _timer.enrol("solvecoarse");
  //  std::cerr << "_mggrid" << _mggrid << "\n";
  int nlevels = _mggrid.nlevels();
  if(dim()==2)
  {
    if(femtype=="Q1")
    {
      _fem = std::unique_ptr<FiniteElementInterface>(new Q12d);
    }
    else
    {
      std::cerr << "unknown fem '" << femtype<<"'\n";
    }
  }
  else if(dim()==3)
  {
    if(femtype=="Q1")
    {
      _fem = std::unique_ptr<FiniteElementInterface>(new Q13d);
    }
    else
    {
      std::cerr << "unknown fem '" << femtype<<"'\n";
    }
  }
  _fem->set_grid(_mggrid.grid());

  _mgmem.set_size(4);
  for(int i=0;i<_mgmem.n();i++)
  {
    set_size(_mgmem(i));
  }
  _mgmatrix.set_size(nlevels);
  for(int l=0;l<nlevels;l++)
  {
    _mgmatrix(l) = _fem->newMatrix(matrixtype);
    _mgmatrix(l)->set_grid(_mggrid.grid(l));
  }
//  _mgmatrix(0)->get_sparse_matrix(_spmat.getSparseMatrix());
  _mgmatrix(nlevels-1)->get_sparse_matrix(_spmat.getSparseMatrix());
  _spmat.computeLu();
  _mgtransfer.set_size(nlevels-1);
  for(int l=0;l<nlevels-1;l++)
  {
    _mgtransfer(l) = _fem->newTransfer(matrixtype);
//    _mgtransfer(l)->set_grid(_mggrid.grid(l));
    _mgtransfer(l)->set_grid(_mggrid.grid(l+1));
  }
  _mgupdate.set_size(nlevels);
  _mgupdatesmooth.set_size(nlevels);
  std::string type="cyc";
  //  type="coef";
  //  type="ortho";
  //  type="restart";
  for(int l=0;l<nlevels;l++)
  {
    if(_optmem==0)
    {
      _mgupdate(l) = std::unique_ptr<UpdaterInterface>(new UpdaterSimple);
    }
    else
    {
      _mgupdate(l) = std::unique_ptr<UpdaterInterface>(new Updater);
    }
    _mgupdatesmooth(l) = std::unique_ptr<UpdaterInterface>(new UpdaterSimple);
  }
  for(int l=0;l<nlevels;l++)
  {
    _mgupdate(l)->setParameters(l, this, std::max(0,_optmem), type);
    _mgupdate(l)->set_size(_mggrid.n(l)+2);
    _mgupdatesmooth(l)->setParameters(l, this, std::max(0,_optmem), type);
    _mgupdatesmooth(l)->set_size(_mggrid.n(l)+2);
  }
}
/*-------------------------------------------------*/
void Operator::set_size(const UniformMultiGrid& umg, std::string femtype, std::string matrixtype)
{
  _mggrid = umg;
  _set_size(femtype, matrixtype);
}

/*-------------------------------------------------*/
void Operator::set_size(int nlevelmax, int nlevels, const armaicvec& n0, std::string femtype, std::string matrixtype)
{
  _mggrid.set_size(nlevelmax, nlevels, n0);
  _set_size(femtype, matrixtype);
}

/*-------------------------------------------------*/
void Operator::set_size(VectorMG& v) const
{
  v.set_size(_mggrid.nlevels());
  for(int l=0;l<_mggrid.nlevels();l++)
  {
    v(l).set_size(_mggrid.n(l)+2);
    v(l).fill_bdry(0);
    v(l).fill_bdry2(0);
    //    std::cerr << "_mggrid.n(l) " << _mggrid.n(l).t() << "\n";
    //    v(l).fill(7);
    //    v(l).fill_bdry(1);
    //    v(l).fill_bdry2(2);
    //    std::cerr << v(l) << "\n";
    //    _boundary(l, v(l));
    //    std::cerr << v(l) << "\n";
    //    assert(0);
  }
}

/*-------------------------------------------------*/
void Operator::smoothpre(int l, Vector& out, const Vector& in) const
  {
  out.fill(0.0);
  if(this->smoother=="jac")
  {
    _mgmatrix(l)->jacobi(out, in);
  }
  else if(this->smoother=="gs")
  {
    _mgmatrix(l)->gauss_seidel1(out, in);
  }
  else if(this->smoother=="gs1")
  {
    _mgmatrix(l)->gauss_seidel1(out, in);
  }
  else if(this->smoother=="gs2")
  {
    _mgmatrix(l)->gauss_seidel2(out, in);
  }
}

/*-------------------------------------------------*/
void Operator::smoothpost(int l, Vector& out, const Vector& in) const
{
  out.fill(0.0);
  if(this->smoother=="jac")
  {
    _mgmatrix(l)->jacobi(out, in);
  }
  else if(this->smoother=="gs")
  {
    _mgmatrix(l)->gauss_seidel2(out, in);
  }
  else if(this->smoother=="gs1")
  {
    _mgmatrix(l)->gauss_seidel1(out, in);
  }
  else if(this->smoother=="gs2")
  {
    _mgmatrix(l)->gauss_seidel2(out, in);
  }
}

/*-------------------------------------------------*/
void Operator::residual(int l, Vector& r, const Vector& u, const Vector& f) const
{
  r =  f;
  _mgmatrix(l)->dot(r, u, -1.0);
}

/*-------------------------------------------------*/
void Operator::vectormg2vector(int l, Vector& u, const VectorMG& umg) const
{
  if(_mggrid.dim()==2)
  {
    int nx = _mggrid.nx(l), ny = _mggrid.ny(l);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        u.at(ix,iy) = umg(l).atp(ix,iy);
      }
    }
  }
  else if(_mggrid.dim()==3)
  {
    int nx = _mggrid.nx(l), ny = _mggrid.ny(l), nz = _mggrid.nz(l);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<nz;iz++)
        {
          u.at(ix,iy,iz) = umg(l).atp(ix,iy,iz);
        }
      }
    }
  }
}
/*-------------------------------------------------*/
void Operator::vector2vectormg(int l, VectorMG& umg, const Vector& u) const
{
  if(_mggrid.dim()==2)
  {
    int nx = _mggrid.nx(l), ny = _mggrid.ny(l);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        umg(l).atp(ix,iy) = u.at(ix,iy);
      }
    }
  }
  else if(_mggrid.dim()==3)
  {
    int nx = _mggrid.nx(l), ny = _mggrid.ny(l), nz = _mggrid.nz(l);
    for(int ix=0;ix<nx;ix++)
    {
      for(int iy=0;iy<ny;iy++)
      {
        for(int iz=0;iz<nz;iz++)
        {
          umg(l).atp(ix,iy,iz) = u.at(ix,iy,iz);
        }
      }
    }
  }
}

/*-------------------------------------------------*/
int Operator::testsolve(bool print, std::string problem)
{
  _u.set_size(_mggrid.n(0));
  _u.fill(0);
  _f.set_size(_u);
  if(problem=="DirichletRhsOne")
  {
    _fem->rhs_one(_f);
  }
  else if(problem=="Random")
  {
    _fem->rhs_random(_f);
  }
  else if(problem=="Linear")
  {
    _f.fill(0);
    _fem->boundary(_u);
    _fem->boundary(_f);
  }
  VectorMG& umg = _mgmem(0);
  VectorMG& fmg = _mgmem(1);

//  vector2vectormg(_mggrid.nlevels()-1, fmg, _f);
//  vector2vectormg(_mggrid.nlevels()-1, umg, _u);
  vector2vectormg(0, fmg, _f);
  vector2vectormg(0, umg, _u);
  int iter = solve(print);
//  vectormg2vector(_mggrid.nlevels()-1, _u, umg);
  vectormg2vector(0, _u, umg);
  return iter;
}
