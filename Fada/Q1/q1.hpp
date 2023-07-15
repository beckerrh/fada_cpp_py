//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef q1_hpp
#define q1_hpp

#include  "../modelbase.hpp"

class BoundaryConditions;
class GridVector;
class GridInterface;
class FemAndMatrixAndSmootherInterface;
class UniformGrid;

/*-------------------------------------------------*/
class Q1 : public ModelBase
{
protected:
  size_t _nx, _ny, _nz;
  double _vol;
  virtual void get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;

public:
  ~Q1()
  {
  }

  Q1(const std::map <std::string, std::string>& parameters, std::shared_ptr<BoundaryConditions const> boundaryconditions=nullptr) : ModelBase(parameters, boundaryconditions)
  {
  }

  Q1(const Q1& model) : ModelBase(model)
  {
  }

  void set_grid(std::shared_ptr <GridInterface const> grid);
  std::shared_ptr <VectorInterface>   newVector(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <SmootherInterface const> newSmoother(std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface const> matrix) const;
  std::shared_ptr <CoarseSolverInterface const> newCoarseSolver(std::string type, std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface const> matrix) const;
  std::shared_ptr <MatrixInterface const>   newMatrix(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <TransferInterface const> newTransfer(std::shared_ptr <GridInterface const>grid) const;

  virtual std::shared_ptr <FemAndMatrixAndSmootherInterface> newStencil(std::shared_ptr <GridInterface const>grid) const = 0;
  void rhs_one(GridVector& v) const;
  void rhs_random(GridVector& v) const;
};

/*-------------------------------------------------*/
class Q12d : public Q1
{
protected:
    void get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    
public:
  ~Q12d()
  {
  }

  Q12d(const std::map <std::string, std::string>& parameters, std::shared_ptr<BoundaryConditions const> bc) : Q1(parameters, bc)
  {
  }

  Q12d(const Q12d& model) : Q1(model)
  {
  }

  std::string toString() const
  {
    return("Q12d");
  }

  void boundary_zero(GridVector& v) const;
  void boundary_linear(GridVector& v) const;

  std::shared_ptr <FemAndMatrixAndSmootherInterface> newStencil(std::shared_ptr <GridInterface const>grid) const;
};

/*-------------------------------------------------*/
class Q13d : public Q1
{
protected:
    void get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    
public:
  ~Q13d()
  {
  }
  Q13d(const std::map <std::string, std::string>& parameters, std::shared_ptr<BoundaryConditions const> bc) : Q1(parameters, bc)
  {
  }

  Q13d(const Q13d& model) : Q1(model)
  {
  }

  std::string toString() const
  {
    return("Q13d");
  }

  void boundary_zero(GridVector& v) const;
  void boundary_linear(GridVector& v) const;

  std::shared_ptr <FemAndMatrixAndSmootherInterface> newStencil(std::shared_ptr <GridInterface const>grid) const;
};

#endif /* q1_hpp */
