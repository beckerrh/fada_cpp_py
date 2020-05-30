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

class UniformGrid;
/*-------------------------------------------------*/
class Q12d : public FiniteElementInterface
{
protected:
  std::string _matrixtype;
  size_t _nx, _ny;
  double _vol;
  std::shared_ptr<UniformGrid> _ug;

public:
  ~Q12d();
  Q12d(std::string matrixtype) : FiniteElementInterface(), _ug(nullptr), _matrixtype(matrixtype) {}
  Q12d(const Q12d& fem) : FiniteElementInterface(fem), _ug(fem._ug), _matrixtype(fem._matrixtype) {}

  std::string toString() const {return "Q12d";}
  void set_grid(std::shared_ptr<GridInterface> grid);
  void rhs_one(Vector& v) const;
  void rhs_random(Vector& v) const;
  void boundary(Vector& v) const;
  std::unique_ptr<MatrixInterface> newMatrix(const GridInterface& grid) const;
  std::unique_ptr<SmootherInterface> newSmoother(std::string type, const GridInterface& grid) const;
  std::unique_ptr<SmootherInterface> newCoarseSolver(std::string type, const GridInterface& grid) const;
  std::unique_ptr<TransferInterface> newTransfer(const GridInterface& grid) const;
  void vectormg2vector(Vector& u, const Vector& umg) const;
  void vector2vectormg(Vector& umg, const Vector& u) const;
  void set_size_mgvector(const GridInterface& grid, Vector& u) const;
};

/*-------------------------------------------------*/
class Q13d : public FiniteElementInterface
{
protected:
  std::string _matrixtype;
  size_t _nx, _ny, _nz;
  double _vol;
  std::shared_ptr<UniformGrid> _ug;

public:
  ~Q13d();
  Q13d(std::string matrixtype) : FiniteElementInterface(), _ug(nullptr), _matrixtype(matrixtype) {}
  Q13d(const Q13d& fem) : FiniteElementInterface(fem), _ug(fem._ug), _matrixtype(fem._matrixtype) {}

  std::string toString() const {return "Q13d";}
  void set_grid(std::shared_ptr<GridInterface> grid);
  void rhs_one(Vector& v) const;
  void rhs_random(Vector& v) const;
  void boundary(Vector& v) const;
  std::unique_ptr<MatrixInterface> newMatrix(const GridInterface& grid) const;
  std::unique_ptr<SmootherInterface> newSmoother(std::string type, const GridInterface& grid) const;
  std::unique_ptr<SmootherInterface> newCoarseSolver(std::string type, const GridInterface& grid) const;
  std::unique_ptr<TransferInterface> newTransfer(const GridInterface& grid) const;
  void vectormg2vector(Vector& u, const Vector& umg) const;
  void vector2vectormg(Vector& umg, const Vector& u) const;
  void set_size_mgvector(const GridInterface& grid, Vector& u) const;
};

#endif /* q1_hpp */
