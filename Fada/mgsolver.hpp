//
//  operateur.hpp
//  Fada
//
//  Created by Roland Becker on 02/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef __operateur_h
#define __operateur_h

#include  <vector>

#include  "typedefs.hpp"

#include  "coarsesolverinterface.hpp"
#include  "modelinterface.hpp"
#include  "updaterinterface.hpp"
#include  "matrixinterface.hpp"
#include  "multigridinterface.hpp"
#include  "smootherinterface.hpp"
#include  "transferinterface.hpp"
#include  "vectorinterface.hpp"
#include  "timer.hpp"

typedef std::vector<std::shared_ptr<VectorInterface>> VectorMG;

/*-------------------------------------------------*/
class MgSolver
{
protected:
  Timer _timer;
  size_t _nlevels;
  std::shared_ptr<ModelInterface const> _model;
  std::vector<std::shared_ptr<MatrixInterface const> > _mgmatrix;
  std::vector<std::shared_ptr<SmootherInterface const> > _mgsmoother;
  std::vector<std::shared_ptr<TransferInterface const> > _mgtransfer;
  std::shared_ptr<CoarseSolverInterface const> _mgcoarsesolver;
  mutable std::vector<std::shared_ptr<UpdaterInterface> > _mgupdate, _mgupdatesmooth;
  int _maxiter;
  double _tol_rel, _tol_abs;
  std::vector<VectorMG> _mgmem;
  void _set_size_vectormg(std::shared_ptr<MultiGridInterface> mgrid, VectorMG& v) const;
  void residual(int l, std::shared_ptr<VectorInterface> r, std::shared_ptr<VectorInterface const> u, std::shared_ptr<VectorInterface const>  f) const;
  void mgstep(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, double tol);

public:

  void set_parameters(int maxiter=30, double tol_rel = 1e-8, double tol_abs = 1e-12);
  ~MgSolver() {}
  MgSolver(bool printtimer=true, bool debug=false): _model(nullptr), _timer(printtimer, debug)
  {
    set_parameters();
  }
  void set_sizes(std::shared_ptr<MultiGridInterface> mgrid, std::shared_ptr<ModelInterface> model, int updatemem=0);
  std::string toString() const;
  int solve(bool print=true);
  std::shared_ptr<VectorInterface const> getU() const {return _mgmem[0][0];}
  std::shared_ptr<VectorInterface> getU() {return _mgmem[0][0];}
  std::shared_ptr<VectorInterface> getF() {return _mgmem[1][0];}
};

#endif
