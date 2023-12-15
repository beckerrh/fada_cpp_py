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
class Q0 : public ModelBase
{
protected:
  virtual void get_locations_values_transfer_2(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
  virtual void get_locations_values_transfer_3(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
  virtual void get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;

public:
  ~Q0()
  {
  }

  Q0(std::string varname, const ParameterMap& parameters, std::shared_ptr<ApplicationInterface const> app=nullptr) : ModelBase(varname, parameters, app)
  {
  }

  std::shared_ptr <VectorInterface>   newVector(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <SmootherInterface> newSmoother(std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface> matrix) const;
  std::shared_ptr <CoarseSolverInterface> newCoarseSolver(std::string type, std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface const> matrix) const;
  std::shared_ptr <MatrixInterface>   newMatrix(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <TransferInterface> newTransfer(std::shared_ptr <GridInterface const>grid, int ref_factor) const;

  void boundary_zero(GridVector& v, std::shared_ptr<GridInterface const> grid) const{_not_written_();}
  void boundary(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const{}
  void update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt) {}
};

/*-------------------------------------------------*/
class Q02d : public Q0
{
protected:
    void get_locations_values_transfer_2(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_transfer_3(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    
public:
  ~Q02d()
  {
  }

  Q02d(std::string varname, const ParameterMap& parameters, std::shared_ptr<ApplicationInterface const> app) : Q0(varname, parameters, app)
  {
  }

  std::string toString() const
  {
    return("Q02d");
  }
  PointDataMap to_point_data(std::shared_ptr<GridVector const> v, std::shared_ptr<GridInterface const> grid) const;
  void rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const;
  std::map<std::string,double> compute_error(const GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> sol) const;
};

/*-------------------------------------------------*/
class Q03d : public Q0
{
protected:
    void get_locations_values_transfer_2(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_transfer_3(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    
public:
  ~Q03d()
  {
  }
  Q03d(std::string varname, const ParameterMap& parameters, std::shared_ptr<ApplicationInterface const> app) : Q0(varname, parameters, app)
  {
  }

  std::string toString() const
  {
    return("Q03d");
  }
};

#endif  
