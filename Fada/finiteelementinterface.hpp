
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
#include  <string>

class GridInterface;
class MatrixInterface;
class SmootherInterface;
class TransferInterface;
class VectorInterface;
class UniformMultiGrid;
/*-------------------------------------------------*/
class FiniteElementInterface
{
public:
  virtual ~FiniteElementInterface() {}
  FiniteElementInterface() {}
  FiniteElementInterface(const FiniteElementInterface& fem) {}
  
  virtual std::string toString() const=0;
  virtual void set_grid(std::shared_ptr<GridInterface> grid)=0;
  virtual void rhs_one(VectorInterface& v) const=0;
  virtual void rhs_random(VectorInterface& v) const=0;
  virtual void boundary(VectorInterface& v) const=0;
  virtual std::unique_ptr<MatrixInterface> newMatrix(const GridInterface& grid) const=0;
  virtual std::unique_ptr<SmootherInterface> newSmoother(std::string type,const GridInterface& grid) const=0;
  virtual std::unique_ptr<SmootherInterface> newCoarseSolver(std::string type,const GridInterface& grid) const=0;
  virtual std::unique_ptr<TransferInterface> newTransfer(const GridInterface& grid) const=0;
  virtual void vectormg2vector(VectorInterface& u, const VectorInterface& umg) const=0;
  virtual void vector2vectormg(VectorInterface& umg, const VectorInterface& u) const=0;
  virtual std::unique_ptr<VectorInterface> newMgvector(const GridInterface& grid) const=0;
};

/*-------------------------------------------------*/
template<typename FEM, class VECTOR>
class FiniteElement : public FEM, public FiniteElementInterface
{
protected:
  FEM& get() { return static_cast<FEM&>(*this); }
  FEM const& get() const { return static_cast<FEM const&>(*this); }
  const VECTOR& getVector(const VectorInterface& u) const {return static_cast<const VECTOR&>(u);}
  VECTOR& getVector(VectorInterface& u) const{return static_cast<VECTOR&>(u);}
public:
  FiniteElement<FEM,VECTOR>() : FEM(), FiniteElementInterface() {}
  FiniteElement<FEM,VECTOR>(std::string matrixtype) : FEM(matrixtype), FiniteElementInterface() {}
  std::string toString() const {return get().toString();}
  void set_grid(std::shared_ptr<GridInterface> grid){get().set_grid(grid);}
  void rhs_one(VectorInterface& v) const{get().rhs_one(getVector(v));}
  void rhs_random(VectorInterface& v) const{get().rhs_random(getVector(v));}
  void boundary(VectorInterface& v) const{get().boundary(getVector(v));}
  std::unique_ptr<MatrixInterface> newMatrix(const GridInterface& grid) const{return get().newMatrix(grid);};
  std::unique_ptr<SmootherInterface> newSmoother(std::string type,const GridInterface& grid) const{return get().newSmoother(type,grid);}
  std::unique_ptr<SmootherInterface> newCoarseSolver(std::string type,const GridInterface& grid) const{return get().newCoarseSolver(type,grid);}
  std::unique_ptr<TransferInterface> newTransfer(const GridInterface& grid) const {return get().newTransfer(grid);}
  void vectormg2vector(VectorInterface& u, const VectorInterface& umg) const
  {
    get().vectormg2vector(getVector(u), getVector(umg));
  }
  void vector2vectormg(VectorInterface& umg, const VectorInterface& u) const
  {
    get().vector2vectormg(getVector(umg), getVector(u));
  }
  std::unique_ptr<VectorInterface> newMgvector(const GridInterface& grid) const {return get().newMgvector(grid);}
};

#endif /* finiteelementinterface_h */
