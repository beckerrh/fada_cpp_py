
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
#include  <map>
#include  "typedefs.hpp"
#include  "applicationinterface.hpp"

class AnalyticalFunctionInterface;
class BoundaryConditions;
// class ApplicationInterface;
class GridInterface;
class MatrixInterface;
class SmootherInterface;
class CoarseSolverInterface;
class TransferInterface;
class VectorInterface;
class UniformMultiGrid;
/*-------------------------------------------------*/
class ModelInterface
{
protected:
    std::string _varname;
public:
  virtual ~ModelInterface() {}
  ModelInterface(std::string varname) : _varname(varname) {}
  ModelInterface(const ModelInterface& model) : _varname(model._varname) {}

  virtual std::string toString() const=0;
  
  virtual std::shared_ptr<MatrixInterface> newMatrix(std::shared_ptr<GridInterface const> grid) const=0;
  virtual std::shared_ptr<SmootherInterface> newSmoother(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix) const=0;
  virtual std::shared_ptr<CoarseSolverInterface> newCoarseSolver(std::string type,std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const=0;
  virtual std::shared_ptr<TransferInterface> newTransfer(std::shared_ptr<GridInterface const> grid, int ref_factor) const=0;
  virtual std::shared_ptr<VectorInterface> newVector(std::shared_ptr<GridInterface const> grid) const=0;

  virtual PointDataMap to_point_data(std::shared_ptr<VectorInterface const> v, std::shared_ptr<GridInterface const> grid) const=0;
  virtual void rhs(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<ApplicationInterface const> app) const=0;
  virtual void boundary_zero(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid) const=0;
  virtual void boundary(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const=0;
  virtual void update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt)=0;
  virtual std::map<std::string,double> compute_error(std::shared_ptr<VectorInterface const> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<ApplicationInterface const> app) const=0;
};

/*-------------------------------------------------*/
template<typename MODEL, class VECTOR>
class Model : public MODEL, public ModelInterface
{
protected:
  MODEL& get() { return static_cast<MODEL&>(*this); }
  const VECTOR& getVector(std::shared_ptr<VectorInterface const> u) const {return static_cast<const VECTOR&>(*u);}
  VECTOR& getVector(std::shared_ptr<VectorInterface> u) const{return static_cast<VECTOR&>(*u);}
  std::shared_ptr<VECTOR const> getVectorPointer(std::shared_ptr<VectorInterface const> u) const {return std::static_pointer_cast<VECTOR const>(u);}
  std::shared_ptr<VECTOR> getVectorPointer(std::shared_ptr<VectorInterface> u) const{return std::static_pointer_cast<VECTOR>(u);}
  
public:
  Model<MODEL,VECTOR>(std::string varname) : MODEL(), ModelInterface(varname) {}
  Model<MODEL,VECTOR>(std::string varname, const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : MODEL(varname, parameters, app), ModelInterface(varname) {}
  
  MODEL const& get() const { return static_cast<MODEL const&>(*this); }
  std::string toString() const {return get().toString();}
  void rhs(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<ApplicationInterface const> app) const{get().rhs(getVector(v), grid, app->rhs(_varname));}
  void boundary_zero(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid) const{get().boundary_zero(getVector(v), grid);}
  void boundary(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const{get().boundary(getVector(v), grid, bc);}
  std::shared_ptr<MatrixInterface> newMatrix(std::shared_ptr<GridInterface const> grid) const{return get().newMatrix(grid);};
  std::shared_ptr<SmootherInterface> newSmoother(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix) const{return get().newSmoother(grid, matrix);}
  std::shared_ptr<CoarseSolverInterface> newCoarseSolver(std::string type,std::shared_ptr<GridInterface const>grid, std::shared_ptr<MatrixInterface const> matrix) const{return get().newCoarseSolver(type,grid, matrix);}
  std::shared_ptr<TransferInterface> newTransfer(std::shared_ptr<GridInterface const> grid, int ref_factor) const {return get().newTransfer(grid, ref_factor);}
  std::shared_ptr<VectorInterface> newVector(std::shared_ptr<GridInterface const> grid) const {return get().newVector(grid);}
  // PointDataMap to_point_data(std::shared_ptr<VectorInterface const> v, std::shared_ptr<GridInterface const> grid) const {return get().to_point_data(getVector(v), grid);}
  PointDataMap to_point_data(std::shared_ptr<VectorInterface const> v, std::shared_ptr<GridInterface const> grid) const {return get().to_point_data(getVectorPointer(v), grid);}
  void update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt){get().update_coefficients(grid, matrix, dt);}
  std::map<std::string,double> compute_error(std::shared_ptr<VectorInterface const> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<ApplicationInterface const> app) const {return get().compute_error(getVector(v), grid, app->solution(_varname));}
};

#endif /* modelinterface_h */
