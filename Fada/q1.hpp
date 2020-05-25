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
  int _nx, _ny;
  double _vol;
  std::shared_ptr<UniformGrid> _ug;

public:
  ~Q12d();
  Q12d() : FiniteElementInterface(), _ug(nullptr) {}
  Q12d(const Q12d& fem) : FiniteElementInterface(fem), _ug(fem._ug) {}

  void set_grid(std::shared_ptr<GridInterface> grid);
  void rhs_one(Vector& v) const;
  void rhs_random(Vector& v) const;
  void boundary(Vector& v) const;
  std::unique_ptr<MatrixInterface> newMatrix(std::string matrixtype) const;
  std::unique_ptr<TransferInterface> newTransfer(std::string matrixtype) const;
};

/*-------------------------------------------------*/
class Q13d : public FiniteElementInterface
{
protected:
  int _nx, _ny, _nz;
  double _vol;
  std::shared_ptr<UniformGrid> _ug;

public:
  ~Q13d();
  Q13d() : FiniteElementInterface(), _ug(nullptr) {}
  Q13d(const Q13d& fem) : FiniteElementInterface(fem), _ug(fem._ug) {}

  void set_grid(std::shared_ptr<GridInterface> grid);
  void rhs_one(Vector& v) const;
  void rhs_random(Vector& v) const;
  void boundary(Vector& v) const;
  std::unique_ptr<MatrixInterface> newMatrix(std::string matrixtype) const;
  std::unique_ptr<TransferInterface> newTransfer(std::string matrixtype) const;
};

#endif /* q1_hpp */
