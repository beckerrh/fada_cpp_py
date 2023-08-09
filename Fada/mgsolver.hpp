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
  mutable Timer _timer;
  size_t _nlevels;
  std::shared_ptr<MultiGridInterface> _mgrid;
  std::shared_ptr<ModelInterface> _model;
  std::vector<std::shared_ptr<MatrixInterface> > _mgmatrix;
  std::vector<std::shared_ptr<SmootherInterface> > _mgsmoother;
  std::vector<std::shared_ptr<TransferInterface> > _mgtransfer;
  std::shared_ptr<CoarseSolverInterface> _mgcoarsesolver;
  mutable std::vector<std::shared_ptr<UpdaterInterface> > _mgupdate, _mgupdatesmooth;
  int _niter_post, _niter_pre;
  // int _maxiter, _niter_post, _niter_pre;
  // double _tol_rel, _tol_abs;
  std::vector<VectorMG> _mgmem;
  void _set_size_vectormg(std::shared_ptr<MultiGridInterface> mgrid, VectorMG& v) const;
  // void residual(int l, std::shared_ptr<VectorInterface> r, std::shared_ptr<VectorInterface const> u, std::shared_ptr<VectorInterface const>  f) const;
  void mgstep(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, double tol);

public:
    /*-------------------------------------------------*/
    struct IterationInfo {
        int maxiter;
        double tol_rel, tol_abs;
        std::string toString(){std::stringstream ss; ss << "maxiter: "<<maxiter<< "\n"; ss<<"tol_rel="<<tol_rel<<"\ttol_abs="<<tol_abs<<"\n"; return ss.str();}
        IterationInfo(int maxiter_=30, double tol_rel_ = 1e-8, double tol_abs_ = 1e-12) : maxiter(maxiter_), tol_rel(tol_rel_), tol_abs(tol_abs_) {}
    };
  // void set_parameters(int maxiter=30, double tol_rel = 1e-8, double tol_abs = 1e-12);
  ~MgSolver() {}
  MgSolver(bool printtimer=true, bool debug=false): _model(nullptr), _timer(printtimer, debug)
  {
    // set_parameters();
  }
  void set_sizes(std::shared_ptr<MultiGridInterface> mgrid, std::shared_ptr<ModelInterface> model);
  std::string toString() const;
  int solve(bool print=true, IterationInfo info=IterationInfo());
  std::shared_ptr<MatrixInterface const> getMatrix(int level=0) const {return _mgmatrix[level];}
  std::shared_ptr<VectorInterface const> getU() const {return _mgmem[0][0];}
  std::shared_ptr<VectorInterface> getU() {return _mgmem[0][0];}
  std::shared_ptr<VectorInterface const> getF() const {return _mgmem[1][0];}
  std::shared_ptr<VectorInterface> getF() {return _mgmem[1][0];}
  void update_coefficients(double dt);
};

#endif
