//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef q1_hpp
#define q1_hpp

#include  "../feminterface.hpp"
#include  "../modelinterface.hpp"
#include  "nodevector.hpp"

class UniformGrid;

/*-------------------------------------------------*/
class Q1
{
protected:
  std::shared_ptr<UniformGrid const> _ug;
  std::string _stenciltype, _matrixtype;
  size_t _nx, _ny, _nz;
  double _vol;

public:
  ~Q1() {}
  // Q1() : _ug(nullptr), _stenciltype(), _matrixtype() {}
  Q1(std::string stenciltype, std::string matrixtype="stencil") : _ug(nullptr), _stenciltype(stenciltype), _matrixtype(matrixtype) {}
  Q1(const Q1& model) : _ug(model._ug), _stenciltype(model._stenciltype) {}

  void set_grid(std::shared_ptr<GridInterface const> grid);
  std::shared_ptr<VectorInterface> newVector(std::shared_ptr<GridInterface const>grid) const;
  std::shared_ptr<SmootherInterface> newSmoother(std::string type, std::shared_ptr<GridInterface const>grid, std::shared_ptr<MatrixInterface const> matrix) const;
  std::shared_ptr<SmootherInterface> newCoarseSolver(std::string type,std::shared_ptr<GridInterface const>grid, std::shared_ptr<MatrixInterface const> matrix) const;
  std::shared_ptr<MatrixInterface> newMatrix(std::shared_ptr<GridInterface const>grid) const;
  // virtual std::shared_ptr<MatrixInterface> newStencil(std::shared_ptr<GridInterface const>grid) const=0;
  virtual std::shared_ptr<FemAndMatrixInterface> newStencil(std::shared_ptr<GridInterface const>grid) const=0;
  void rhs_one(NodeVector& v) const;
  void rhs_random(NodeVector& v) const;
};

/*-------------------------------------------------*/
class Q12d : public Q1
{
public:
  ~Q12d() {}
  // Q12d() : Q1() {}
  Q12d(std::string stenciltype, std::string matrixtype="stencil") : Q1(stenciltype, matrixtype) {}
  Q12d(const Q12d& model) : Q1(model) {}

  std::string toString() const {return "Q12d";}
  void boundary(NodeVector& v) const;
  // std::shared_ptr<MatrixInterface> newStencil(std::shared_ptr<GridInterface const>grid) const;
  std::shared_ptr<FemAndMatrixInterface> newStencil(std::shared_ptr<GridInterface const>grid) const;
  std::shared_ptr<TransferInterface> newTransfer(std::shared_ptr<GridInterface const>grid) const;
};

/*-------------------------------------------------*/
class Q13d : public Q1
{
public:
  ~Q13d() {}
  // Q13d() : Q1() {}
  Q13d(std::string stenciltype, std::string matrixtype="stencil") : Q1(stenciltype, matrixtype) {}
  Q13d(const Q13d& model) : Q1(model) {}

  std::string toString() const {return "Q13d";}
  void boundary(NodeVector& v) const;
  // std::shared_ptr<MatrixInterface> newStencil(std::shared_ptr<GridInterface const> grid) const;
  std::shared_ptr<FemAndMatrixInterface> newStencil(std::shared_ptr<GridInterface const>grid) const;
  std::shared_ptr<TransferInterface> newTransfer(std::shared_ptr<GridInterface const>grid) const;
};

#endif /* q1_hpp */
