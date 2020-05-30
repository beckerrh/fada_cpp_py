//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef q1_hpp
#define q1_hpp

#include  "finiteelementinterface.hpp"
#include  "nodevector.hpp"

class UniformGrid;
/*-------------------------------------------------*/
class Q12d
{
protected:
  std::string _matrixtype;
  size_t _nx, _ny;
  double _vol;
  std::shared_ptr<UniformGrid> _ug;

public:
  ~Q12d();
  Q12d() : _ug(nullptr), _matrixtype() {}
  Q12d(std::string matrixtype) : _ug(nullptr), _matrixtype(matrixtype) {}
  Q12d(const Q12d& fem) : _ug(fem._ug), _matrixtype(fem._matrixtype) {}

  std::string toString() const {return "Q12d";}
  void set_grid(std::shared_ptr<GridInterface> grid);
  void rhs_one(NodeVector& v) const;
  void rhs_random(NodeVector& v) const;
  void boundary(NodeVector& v) const;
  std::unique_ptr<MatrixInterface> newMatrix(const GridInterface& grid) const;
  std::unique_ptr<SmootherInterface> newSmoother(std::string type, const GridInterface& grid) const;
  std::unique_ptr<SmootherInterface> newCoarseSolver(std::string type, const GridInterface& grid) const;
  std::unique_ptr<TransferInterface> newTransfer(const GridInterface& grid) const;
  void vectormg2vector(NodeVector& u, const NodeVector& umg) const;
  void vector2vectormg(NodeVector& umg, const NodeVector& u) const;
  std::unique_ptr<VectorInterface> newMgvector(const GridInterface& grid) const;
};

/*-------------------------------------------------*/
class Q13d
{
protected:
  std::string _matrixtype;
  size_t _nx, _ny, _nz;
  double _vol;
  std::shared_ptr<UniformGrid> _ug;

public:
  ~Q13d();
  Q13d() : _ug(nullptr), _matrixtype() {}
  Q13d(std::string matrixtype) : _ug(nullptr), _matrixtype(matrixtype) {}
  Q13d(const Q13d& fem) : _ug(fem._ug), _matrixtype(fem._matrixtype) {}

  std::string toString() const {return "Q13d";}
  void set_grid(std::shared_ptr<GridInterface> grid);
  void rhs_one(NodeVector& v) const;
  void rhs_random(NodeVector& v) const;
  void boundary(NodeVector& v) const;
  std::unique_ptr<MatrixInterface> newMatrix(const GridInterface& grid) const;
  std::unique_ptr<SmootherInterface> newSmoother(std::string type, const GridInterface& grid) const;
  std::unique_ptr<SmootherInterface> newCoarseSolver(std::string type, const GridInterface& grid) const;
  std::unique_ptr<TransferInterface> newTransfer(const GridInterface& grid) const;
  void vectormg2vector(NodeVector& u, const NodeVector& umg) const;
  void vector2vectormg(NodeVector& umg, const NodeVector& u) const;
  std::unique_ptr<VectorInterface> newMgvector(const GridInterface& grid) const;
};

#endif /* q1_hpp */
