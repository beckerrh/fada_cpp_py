
//
//  finiteelementinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef finiteelementinterface_h
#define finiteelementinterface_h

#include  <memory>

class GridInterface;
class MatrixInterface;
class SmootherInterface;
class TransferInterface;
class UniformMultiGrid;
class Vector;
/*-------------------------------------------------*/
class FiniteElementInterface
{
public:
  virtual ~FiniteElementInterface() {}
  FiniteElementInterface() {}
  FiniteElementInterface(const FiniteElementInterface& fem) {}

  virtual std::string toString() const=0;
  virtual void set_grid(std::shared_ptr<GridInterface> grid)=0;
  virtual void rhs_one(Vector& v) const=0;
  virtual void rhs_random(Vector& v) const=0;
  virtual void boundary(Vector& v) const=0;
  virtual std::unique_ptr<MatrixInterface> newMatrix(const GridInterface& grid) const=0;
  virtual std::unique_ptr<SmootherInterface> newSmoother(std::string type,const GridInterface& grid) const=0;
  virtual std::unique_ptr<SmootherInterface> newCoarseSolver(std::string type,const GridInterface& grid) const=0;
  virtual std::unique_ptr<TransferInterface> newTransfer(const GridInterface& grid) const=0;
  virtual void vectormg2vector(Vector& u, const Vector& umg) const=0;
  virtual void vector2vectormg(Vector& umg, const Vector& u) const=0;
  virtual void set_size_mgvector(const GridInterface& grid, Vector& u) const=0;
};


#endif /* finiteelementinterface_h */
