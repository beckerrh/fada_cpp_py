//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef model_stokes_hpp
#define model_stokes_hpp


#include  "gridvector.hpp"
#include  "modelinterface.hpp"
#include  "model_p.hpp"
#include  "model_v.hpp"

class BoundaryConditions;
class StokesVector;
class GridInterface;
class FemAndMatrixAndSmootherInterface;
class UniformGrid;

/*-------------------------------------------------*/
class ModelStokes : public ModelInterface
{
protected:
    std::shared_ptr<ModelP2d> _model_p;
    std::vector<std::shared_ptr<ModelV2d>> _model_v;

public:
  ~ModelStokes()
  {
  }
  ModelStokes(const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app=nullptr);
  ModelStokes(const ModelStokes& model) : ModelInterface("v-p"), _model_p(model._model_p), _model_v(model._model_v)
  {
  }
  int dim() const {return _model_v.size();}
  std::shared_ptr <VectorInterface>   newVector(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <SmootherInterface> newSmoother(std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface> matrix) const;
  std::shared_ptr <CoarseSolverInterface> newCoarseSolver(std::string type, std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface const> matrix) const;
  std::shared_ptr <MatrixInterface>   newMatrix(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <TransferInterface> newTransfer(std::shared_ptr <GridInterface const>grid, int ref_factor) const;

  // void boundary_zero(StokesVector& v, std::shared_ptr<GridInterface const> grid) const{_not_written_();}
  // void boundary_linear(StokesVector& v, std::shared_ptr<GridInterface const> grid) const{_not_written_();}
  // void update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt) {}
  // void rhs(StokesVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const{_not_written_();}
  // PointDataMap to_point_data(const StokesVector& v, std::shared_ptr<GridInterface const> grid) const{_not_written_();}

  std::string toString() const {return "ModelStokes";}
  PointDataMap to_point_data(std::shared_ptr<VectorInterface const> v, std::shared_ptr<GridInterface const> grid) const;
  void rhs(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<ApplicationInterface const> app) const;
  void boundary_zero(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid) const;
  void boundary(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const;
  void update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt);
  std::map<std::string,double> compute_error(std::shared_ptr<VectorInterface const> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<ApplicationInterface const> app) const;
};

/*-------------------------------------------------*/
class ModelStokes2d : public ModelStokes
{
protected:
    
public:
  ~ModelStokes2d()
  {
  }

  ModelStokes2d(const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : ModelStokes(parameters, app)
  {
  }

  ModelStokes2d(const ModelStokes2d& model) : ModelStokes(model)
  {
  }

  std::string toString() const
  {
    return("ModelStokes2d");
  }
};

/*-------------------------------------------------*/
class ModelStokes3d : public ModelStokes
{
protected:
    
public:
  ~ModelStokes3d()
  {
  }
  ModelStokes3d(const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : ModelStokes(parameters, app)
  {
  }

  ModelStokes3d(const ModelStokes3d& model) : ModelStokes(model)
  {
  }

  std::string toString() const
  {
    return("ModelStokes3d");
  }
};

#endif  
