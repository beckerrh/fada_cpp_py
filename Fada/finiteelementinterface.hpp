
//
//  finiteelementinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef finiteelementinterface_h
#define finiteelementinterface_h

class GridInterface;
class MatrixInterface;
class TransferInterface;
class Vector;
/*-------------------------------------------------*/
class FiniteElementInterface
{
public:
  virtual ~FiniteElementInterface() {}
  FiniteElementInterface() {}
  FiniteElementInterface(const FiniteElementInterface& fem) {}

  virtual void set_grid(const GridInterface& grid)=0;
  virtual void rhs_one(Vector& v) const=0;
  virtual void rhs_random(Vector& v) const=0;
  virtual void boundary(Vector& v) const=0;
  virtual std::unique_ptr<MatrixInterface> newMatrix(std::string matrixtype) const=0;
  virtual std::unique_ptr<TransferInterface> newTransfer(std::string matrixtype) const=0;
};


#endif /* finiteelementinterface_h */
