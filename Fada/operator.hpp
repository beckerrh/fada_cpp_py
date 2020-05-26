//
//  operateur.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef __operateur_h
#define __operateur_h

#include  "typedefs.hpp"
#include  "finiteelementinterface.hpp"
#include  "updaterinterface.hpp"
#include  "matrixinterface.hpp"
#include  "transferinterface.hpp"
#include  "uniformmultigrid.hpp"
#include  "timer.hpp"
#include  "umfmatrix.hpp"

/*-------------------------------------------------*/
class Operator
{
protected:
  Timer _timer;
  int _optmem;
  UniformMultiGrid _mggrid;
  std::shared_ptr<FiniteElementInterface> _fem;
  UmfMatrix _spmat;
  Array<std::shared_ptr<MatrixInterface> > _mgmatrix;
  Array<std::shared_ptr<TransferInterface> > _mgtransfer;
  Array<VectorMG> _mgmem;
  Vector _u, _f;
  mutable Array<std::shared_ptr<UpdaterInterface> > _mgupdate, _mgupdatesmooth;

  void _set_size(std::string femtype, std::string matrixtype);

  void vectormg2vector(int l, Vector& u, const VectorMG& umg) const;
  void vector2vectormg(int l, VectorMG& umg, const Vector& u) const;
  void residual(int l, Vector& r, const Vector& u, const Vector& f) const;

  void smoothpre(int l, Vector& out, const Vector& in) const;
  void smoothpost(int l, Vector& out, const Vector& in) const;
  void mgstep(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, double tol);
  int solve(bool print=true);
  void solve_coarse(int l, Vector& u, const Vector& f, Vector& d, Vector& w) const;

public:
  std::string smoother;
  int maxiter;
  double tol_rel, tol_abs;

  void set_parameters();
  ~Operator();
  Operator(bool printtimer=true);
  void set_size(int nlevelmax, int nlevels, const armaicvec& n0, std::string femtype="Q1", std::string matrixtype="Full");
  void set_size(const UniformMultiGrid& umg, std::string femtype="Q1", std::string matrixtype="Full");

  int nall() const { return _mggrid.nall();}
  int dim() const { return _mggrid.dim();}
  int nmax(int i) const {return _mggrid.nmax(i);}

  const Vector& get_solution() const {return _u;}
  Vector& get_solution() {return _u;}

  void set_size(VectorMG& v) const;

  void dot(int l, Vector& out, const Vector& in, double d) const {_mgmatrix(l)->dot(out, in, d);}
  int testsolve(bool print=true, std::string problem="DirichletRhsOne");
};

#endif
