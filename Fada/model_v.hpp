//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef model_v_hpp
#define model_v_hpp


#include  "modelbase.hpp"

class BoundaryConditions;
class GridVector;
class GridInterface;
class FemAndMatrixAndSmootherInterface;
class UniformGrid;

/*-------------------------------------------------*/
class ModelV : public ModelBase
{
protected:
    double _dt;
    int _idir;
    virtual void get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
    virtual void get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
    virtual void get_locations_values_matrix_divergence(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
    typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;

public:
    ~ModelV()
    {
    }

    ModelV(const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app=nullptr) : ModelBase(parameters, app)
    {
        _idir = std::stoi(parameters.at("direction"));
        _dt = std::stod(parameters.at("dt"));
        // std::cerr << "_idir = "  << _idir << " _dt = "  << _dt << "\n";
    }

    ModelV(const ModelV& model) : ModelBase(model), _idir(model._idir), _dt(model._dt)
    {
    }

    std::shared_ptr <VectorInterface>   newVector(std::shared_ptr <GridInterface const>grid) const;
    std::shared_ptr <SmootherInterface> newSmoother(std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface> matrix) const;
    std::shared_ptr <CoarseSolverInterface> newCoarseSolver(std::string type, std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface const> matrix) const;
    std::shared_ptr <MatrixInterface>   newMatrix(std::shared_ptr <GridInterface const>grid) const;
    std::shared_ptr <TransferInterface> newTransfer(std::shared_ptr <GridInterface const>grid, int ref_factor) const;

    std::shared_ptr <MatrixInterface>   newMatrixDivergence(std::shared_ptr <GridInterface const>grid) const;

    void boundary_zero(GridVector& v, std::shared_ptr<GridInterface const> grid) const{_not_written_();}
    void boundary(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const{_not_written_();}
    void update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt);
};

/*-------------------------------------------------*/
class ModelV2d : public ModelV
{
protected:
    void get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix_divergence(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    
public:
    ~ModelV2d()
    {
    }

    ModelV2d(const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : ModelV(parameters, app)
    {
    }

    ModelV2d(const ModelV2d& model) : ModelV(model)
    {
    }

    std::string toString() const
    {
        return("ModelV2d");
    }
    PointDataMap to_point_data(const GridVector& v, std::shared_ptr<GridInterface const> grid) const;
    void rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const;
    std::map<std::string,double> compute_error(const GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> sol) const;
};

/*-------------------------------------------------*/
class ModelV3d : public ModelV
{
protected:
    void get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix_divergence(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    
public:
    ~ModelV3d()
    {
    }
    ModelV3d(const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : ModelV(parameters, app)
    {
    }

    ModelV3d(const ModelV3d& model) : ModelV(model)
    {
    }

    std::string toString() const
    {
        return("ModelV3d");
    }
};

#endif  
