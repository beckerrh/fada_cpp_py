//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef q1_hpp
#define q1_hpp

#include  "modelbase.hpp"

class BoundaryConditions;
class GridVector;
class GridInterface;
class FemAndMatrixAndSmootherInterface;
class UniformGrid;

/*-------------------------------------------------*/
class Q1 : public ModelBase
{
protected:
  virtual void get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
  virtual std::shared_ptr <FemAndMatrixAndSmootherInterface> newStencil(std::shared_ptr <GridInterface const>grid) const = 0;

public:
  Q1(std::string varname, const ParameterMap& parameters, std::shared_ptr<ApplicationInterface const> app=nullptr);

  void set_grid(std::shared_ptr <GridInterface const> grid);
  std::shared_ptr <VectorInterface>   newVector(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <SmootherInterface> newSmoother(std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface> matrix) const;
  std::shared_ptr <CoarseSolverInterface> newCoarseSolver(std::string type, std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface const> matrix) const;
  std::shared_ptr <MatrixInterface>   newMatrix(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <TransferInterface> newTransfer(std::shared_ptr <GridInterface const>grid, int ref_factor) const;

  PointDataMap to_point_data(std::shared_ptr<GridVector const> v, std::shared_ptr<GridInterface const> grid) const;
  void update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt) {}
};

/*-------------------------------------------------*/
class Q12d : public Q1
{
protected:
    void get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    std::shared_ptr <FemAndMatrixAndSmootherInterface> newStencil(std::shared_ptr <GridInterface const>grid) const;
    
public:
  ~Q12d()
  {
  }

  Q12d(std::string varname, const ParameterMap& parameters, std::shared_ptr<ApplicationInterface const> app) : Q1(varname, parameters, app)
  {
  }

  std::string toString() const
  {
    return("Q12d");
  }

  void boundary_zero(GridVector& v, std::shared_ptr<GridInterface const> grid) const;
  void boundary(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const;
  std::map<std::string,double> compute_error(const GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> sol) const;
  void rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const;
};

/*-------------------------------------------------*/
class Q13d : public Q1
{
protected:
    void get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    std::shared_ptr <FemAndMatrixAndSmootherInterface> newStencil(std::shared_ptr <GridInterface const>grid) const;
    
public:
  ~Q13d()
  {
  }
  Q13d(std::string varname, const ParameterMap& parameters, std::shared_ptr<ApplicationInterface const> app) : Q1(varname, parameters, app)
  {
  }

  std::string toString() const
  {
    return("Q13d");
  }

  void boundary_zero(GridVector& v, std::shared_ptr<GridInterface const> grid) const;
  void boundary(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const;
  std::map<std::string,double> compute_error(const GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> sol) const;

};

#endif /* q1_hpp */
