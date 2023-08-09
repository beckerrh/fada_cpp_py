//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef q1shifted_hpp
#define q1shifted_hpp


#include  "modelbase.hpp"

class BoundaryConditions;
class GridVector;
class GridInterface;
class FemAndMatrixAndSmootherInterface;
class UniformGrid;
class SparseMatrix;

/*-------------------------------------------------*/
class Q1shifted : public ModelBase
{
protected:
    double _dt;
    int _idir;
    bool _periodic;
    virtual void get_locations_values_transfer_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
    virtual void get_locations_values_matrix_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
    virtual void get_locations_values_matrix_divergence_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
    virtual void get_locations_values_transfer_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
    virtual void get_locations_values_matrix_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
    virtual void get_locations_values_matrix_divergence_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const=0;
    virtual PointDataMap to_point_data_per(const GridVector& v, std::shared_ptr<GridInterface const> grid) const=0;
    virtual PointDataMap to_point_data_dir(const GridVector& v, std::shared_ptr<GridInterface const> grid) const=0;

    typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;

public:
    ~Q1shifted()
    {
    }

    Q1shifted(std::string varname, const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app=nullptr);

    std::shared_ptr <VectorInterface>   newVector(std::shared_ptr <GridInterface const>grid) const;
    std::shared_ptr <SmootherInterface> newSmoother(std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface> matrix) const;
    std::shared_ptr <CoarseSolverInterface> newCoarseSolver(std::string type, std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface const> matrix) const;
    std::shared_ptr <MatrixInterface>   newMatrix(std::shared_ptr <GridInterface const>grid) const;
    std::shared_ptr <TransferInterface> newTransfer(std::shared_ptr <GridInterface const>grid, int ref_factor) const;

    std::shared_ptr <MatrixInterface>   newMatrixDivergence(std::shared_ptr <GridInterface const>grid) const;

    PointDataMap to_point_data(std::shared_ptr<GridVector const> v, std::shared_ptr<GridInterface const> grid) const;
    void boundary_zero(GridVector& v, std::shared_ptr<GridInterface const> grid) const{_not_written_();}
    void update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt);
};

/*-------------------------------------------------*/
class Q1shifted2d : public Q1shifted
{
protected:
    void get_locations_values_transfer_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix_divergence_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_transfer_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix_divergence_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    
public:
    ~Q1shifted2d()
    {
    }

    Q1shifted2d(std::string varname, const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : Q1shifted(varname, parameters, app)
    {
    }

    std::string toString() const
    {
        return("Q1shifted2d");
    }
    PointDataMap to_point_data_per(const GridVector& v, std::shared_ptr<GridInterface const> grid) const;
    PointDataMap to_point_data_dir(const GridVector& v, std::shared_ptr<GridInterface const> grid) const;
    void boundary(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const;
    void rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const;
    std::map<std::string,double> compute_error(const GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> sol) const;
};

/*-------------------------------------------------*/
class Q1shifted3d : public Q1shifted
{
protected:
    void get_locations_values_transfer_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix_divergence_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_transfer_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    void get_locations_values_matrix_divergence_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const;
    
public:
    ~Q1shifted3d()
    {
    }
    Q1shifted3d(std::string varname, const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : Q1shifted(varname, parameters, app){}
    std::string toString() const
    {
        return("Q1shifted3d");
    }
    void boundary(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const{_not_written_();}
    PointDataMap to_point_data_per(const GridVector& v, std::shared_ptr<GridInterface const> grid) const{_not_written_(); return PointDataMap();}
    PointDataMap to_point_data_dir(const GridVector& v, std::shared_ptr<GridInterface const> grid) const{_not_written_(); return PointDataMap();}
};

#endif  
