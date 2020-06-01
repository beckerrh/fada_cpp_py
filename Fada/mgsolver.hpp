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

//#include  "array.hpp"
#include  "typedefs.hpp"

#include  "finiteelementinterface.hpp"
#include  "updaterinterface.hpp"
#include  "matrixinterface.hpp"
#include  "multigridinterface.hpp"
#include  "smootherinterface.hpp"
#include  "transferinterface.hpp"
#include  "vectorinterface.hpp"

//#include  "uniformmultigrid.hpp"
#include  "timer.hpp"
//#include  "vector.hpp"

typedef std::vector<std::shared_ptr<VectorInterface>> VectorMG;
//typedef std::vector<Vector> VectorMG;

/*-------------------------------------------------*/
class MgSolver
{
protected:
  Timer _timer;
  size_t _nlevels;
  std::shared_ptr<FiniteElementInterface> _fem;
  std::vector<std::shared_ptr<MatrixInterface> > _mgmatrix;
  std::vector<std::shared_ptr<SmootherInterface> > _mgsmoother;
  std::vector<std::shared_ptr<TransferInterface> > _mgtransfer;
  mutable std::vector<std::shared_ptr<UpdaterInterface> > _mgupdate, _mgupdatesmooth;

  std::vector<VectorMG> _mgmem;

  void _set_size_vectormg(std::shared_ptr<MultiGridInterface> mgrid, VectorMG& v) const;
  void residual(int l, VectorInterface& r, const VectorInterface& u, const VectorInterface& f) const;
  void mgstep(int l, VectorMG& u, VectorMG& f, VectorMG& d, VectorMG& w, double tol);

public:
  int maxiter;
  double tol_rel, tol_abs;

  void set_parameters();
  ~MgSolver() {}
  MgSolver(bool printtimer=true): _fem(nullptr), _timer(printtimer)
  {
    set_parameters();
  }
  void set_sizes(std::shared_ptr<MultiGridInterface> mgrid, std::shared_ptr<FiniteElementInterface> fem, std::string smoothertype, int updatemem=0);
  std::string toString() const;
  int solve(VectorInterface& u, const VectorInterface& f, bool print=true);
};

#endif
