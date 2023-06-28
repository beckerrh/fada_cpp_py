
//
//  modelinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef modelinterface_h
#define modelinterface_h

#include  <memory>
#include  <string>

class GridInterface;
class MatrixInterface;
class SmootherInterface;
class TransferInterface;
class VectorInterface;
class UniformMultiGrid;
/*-------------------------------------------------*/
class ModelInterface
{
public:
  virtual ~ModelInterface() {}
  ModelInterface() {}
  ModelInterface(const ModelInterface& model) {}

  virtual std::string toString() const=0;
  virtual void set_grid(std::shared_ptr<GridInterface const> grid)=0;
  virtual void rhs_one(std::shared_ptr<VectorInterface> v) const=0;
  virtual void rhs_random(std::shared_ptr<VectorInterface> v) const=0;
  virtual void boundary(std::shared_ptr<VectorInterface> v) const=0;
  virtual std::shared_ptr<MatrixInterface> newMatrix(std::shared_ptr<GridInterface const> grid) const=0;
  virtual std::shared_ptr<SmootherInterface> newSmoother(std::string type,std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const=0;
  virtual std::shared_ptr<SmootherInterface> newCoarseSolver(std::string type,std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const=0;
  virtual std::shared_ptr<TransferInterface> newTransfer(std::shared_ptr<GridInterface const> grid) const=0;
  virtual std::shared_ptr<VectorInterface> newVector(std::shared_ptr<GridInterface const> grid) const=0;
};

/*-------------------------------------------------*/
template<typename MODEL, class VECTOR>
class Model : public MODEL, public ModelInterface
{
protected:
  MODEL& get() { return static_cast<MODEL&>(*this); }
  MODEL const& get() const { return static_cast<MODEL const&>(*this); }
  const VECTOR& getVector(std::shared_ptr<VectorInterface const> u) const {return static_cast<const VECTOR&>(*u);}
  VECTOR& getVector(std::shared_ptr<VectorInterface> u) const{return static_cast<VECTOR&>(*u);}
public:
  Model<MODEL,VECTOR>() : MODEL(), ModelInterface() {}
  Model<MODEL,VECTOR>(std::string stenciltype, std::string matrixtype) : MODEL(stenciltype, matrixtype), ModelInterface() {}
  std::string toString() const {return get().toString();}
  void set_grid(std::shared_ptr<GridInterface const> grid){get().set_grid(grid);}
  void rhs_one(std::shared_ptr<VectorInterface> v) const{get().rhs_one(getVector(v));}
  void rhs_random(std::shared_ptr<VectorInterface> v) const{get().rhs_random(getVector(v));}
  void boundary(std::shared_ptr<VectorInterface> v) const{get().boundary(getVector(v));}
  std::shared_ptr<MatrixInterface> newMatrix(std::shared_ptr<GridInterface const> grid) const{return get().newMatrix(grid);};
  std::shared_ptr<SmootherInterface> newSmoother(std::string type,std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const{return get().newSmoother(type,grid, matrix);}
  std::shared_ptr<SmootherInterface> newCoarseSolver(std::string type,std::shared_ptr<GridInterface const>grid, std::shared_ptr<MatrixInterface const> matrix) const{return get().newCoarseSolver(type,grid, matrix);}
  std::shared_ptr<TransferInterface> newTransfer(std::shared_ptr<GridInterface const> grid) const {return get().newTransfer(grid);}
  std::shared_ptr<VectorInterface> newVector(std::shared_ptr<GridInterface const> grid) const {return get().newVector(grid);}
};

#endif /* modelinterface_h */
