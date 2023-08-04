//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef model_p_hpp
#define model_p_hpp


#include  "modelbase.hpp"

class BoundaryConditions;
class GridVector;
class GridInterface;
class FemAndMatrixAndSmootherInterface;
class UniformGrid;

/*-------------------------------------------------*/
class ModelP : public ModelBase
{
protected:
  virtual void get_locations_values_transfer_2(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
  virtual void get_locations_values_transfer_3(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
  virtual void get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;

public:
  ~ModelP()
  {
  }

  ModelP(const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app=nullptr) : ModelBase(parameters, app)
  {
  }

  ModelP(const ModelP& model) : ModelBase(model)
  {
  }

  std::shared_ptr <VectorInterface>   newVector(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <SmootherInterface> newSmoother(std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface> matrix) const;
  std::shared_ptr <CoarseSolverInterface> newCoarseSolver(std::string type, std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface const> matrix) const;
  std::shared_ptr <MatrixInterface>   newMatrix(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <TransferInterface> newTransfer(std::shared_ptr <GridInterface const>grid, int ref_factor) const;

  void boundary_zero(GridVector& v, std::shared_ptr<GridInterface const> grid) const{_not_written_();}
  void boundary(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const{_not_written_();}
  void update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt) {}
};

/*-------------------------------------------------*/
class ModelP2d : public ModelP
{
protected:
    void get_locations_values_transfer_2(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_transfer_3(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    
public:
  ~ModelP2d()
  {
  }

  ModelP2d(const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : ModelP(parameters, app)
  {
  }

  ModelP2d(const ModelP2d& model) : ModelP(model)
  {
  }

  std::string toString() const
  {
    return("ModelP2d");
  }
  PointDataMap to_point_data(const GridVector& v, std::shared_ptr<GridInterface const> grid) const;
  void rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const;
  std::map<std::string,double> compute_error(const GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> sol) const;
};

/*-------------------------------------------------*/
class ModelP3d : public ModelP
{
protected:
    void get_locations_values_transfer_2(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_transfer_3(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    
public:
  ~ModelP3d()
  {
  }
  ModelP3d(const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : ModelP(parameters, app)
  {
  }

  ModelP3d(const ModelP3d& model) : ModelP(model)
  {
  }

  std::string toString() const
  {
    return("ModelP3d");
  }
};

#endif  
