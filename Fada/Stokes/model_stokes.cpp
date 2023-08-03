//
//  p0.cpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <cassert>
#include  <cmath>
#include  <armadillo>
#include  <sstream>
#include  "../analyticalfunctioninterface.hpp"
#include  "../coarsesolverinterface.hpp"
#include  "../coarsesolver_arma.hpp"
#include  "../construct_elements_matrix.hpp"
#include  "../feminterface.hpp"
#include  "../gridindex.hpp"
#include  "../gridvector.hpp"
#include  "../matrixinterface.hpp"
#include  "../uniformgrid.hpp"
#include  "../smoothersimple.hpp"
#include  "../solverumf.hpp"
#include  "../sparsematrix.hpp"
#include  "../sparsematrix_arma.hpp"
#include  "../transferbymatrix.hpp"
#include  "../transferinterface.hpp"
#include  "model_stokes.hpp"
#include  "model_p.hpp"
#include  "model_v.hpp"
#include  "stokesvector.hpp"
#include  "smoother_chorin.hpp"


/*-------------------------------------------------*/
ModelStokes::ModelStokes(const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : ModelInterface("v-p")
{
    assert(app);
    int dim = app->get_dim();
    _model_v.resize(dim);
    for(int i=0; i< dim;i++)
    {
        std::map<std::string,std::string> parameters2(parameters.begin(), parameters.end());
        parameters2["direction"] =  std::to_string(i);
        _model_v[i] = std::make_shared<ModelV2d>(parameters2, app);
    }
    _model_p = std::make_shared<ModelP2d>(parameters, app);    
}

/*-------------------------------------------------*/
std::shared_ptr<VectorInterface> ModelStokes::newVector(std::shared_ptr<GridInterface const> grid) const
{
    int dim = grid->dim();
    auto p = std::make_shared<StokesVector>(dim);
    p->get(dim) = _model_p->newVector(grid);
    for(int i=0;i<dim;i++)
    {
        p->get(i) = _model_v[i]->newVector(grid);
    }
    // for(int i=0;i<p->n_vectors();i++)
    // {
    //     std::cerr << "### " << i << "\n";
    //     p->get(i)->save(std::cerr, arma::arma_ascii);
    // }
    // std::cerr << "### p->get_size()=" << p->get_size()<<"\n";
    return p;
}

/*-------------------------------------------------*/
std::shared_ptr<MatrixInterface> ModelStokes::newMatrix(std::shared_ptr<GridInterface const> grid) const
{
    std::map<std::string,std::shared_ptr<MatrixInterface>> matrices;
    std::map<std::pair<int,int>,std::pair<std::string, std::string>> patterns;
    int dim = grid->dim();
    for(int i=0;i<dim;i++)
    {
        std::stringstream ss;
        ss << "A" << i; 
        matrices[ss.str()] = _model_v[i]->newMatrix(grid);
        patterns[std::pair<int,int>(i,i)] = std::pair<std::string, std::string>(ss.str(),"");
    }
    for(int i=0;i<dim;i++)
    {
        std::stringstream ss;
        ss << "B" << i; 
        matrices[ss.str()] = _model_v[i]->newMatrixDivergence(grid);
        patterns[std::pair<int,int>(dim,i)] = std::pair<std::string, std::string>(ss.str(),"");
        patterns[std::pair<int,int>(i,dim)] = std::pair<std::string, std::string>(ss.str(),"MT");
    }
    // UMF needs diagonals OR peridoic
    int np = matrices.at("B0")->nrows();
    arma::umat locations(2,np);
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    double h = ug->dx(0);
    armavec values(np,arma::fill::value(0.1*h*h));
    locations.row(0) = arma::linspace<arma::uvec>(0,np-1,np).t(); 
    locations.row(1) = arma::linspace<arma::uvec>(0,np-1,np).t();     
    matrices["C"] = std::make_shared<Matrix<SparseMatrix,Vector<armavec>>>(locations, values);
    patterns[std::pair<int,int>(dim,dim)] = std::pair<std::string, std::string>("C","");
    return std::make_shared<SystemMatrix>(matrices, patterns);
}

/*-------------------------------------------------*/
std::shared_ptr<SmootherInterface> ModelStokes::newSmoother(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix) const
{
    return std::make_shared<Smoother_Chorin>(grid, matrix);
}
/*-------------------------------------------------*/
std::shared_ptr<CoarseSolverInterface> ModelStokes::newCoarseSolver(std::string type, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const
{
    return std::make_shared<SolverUmfSystem>(matrix);
}

/*-------------------------------------------------*/
std::shared_ptr<TransferInterface> ModelStokes::newTransfer(std::shared_ptr<GridInterface const> grid, int ref_factor) const
{
    int dim = grid->dim();
    auto p = std::make_shared<SystemTransfer>(dim+1);
    p->get(dim) = _model_p->newTransfer(grid, ref_factor);
    for(int i=0;i<dim;i++)
    {
        p->get(i) = _model_v[i]->newTransfer(grid, ref_factor);
    }
    return p;
    
}
/*-------------------------------------------------*/
void ModelStokes::update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt)
{
    int dim = grid->dim();
    auto pm = std::dynamic_pointer_cast<SystemMatrix>(matrix);  
    const std::map<std::string,std::shared_ptr<MatrixInterface>>&  matrices  = pm->get_matrices();
    for(int i=0;i<dim;i++)
    {
        std::stringstream ss;
        ss << "A" << i; 
        _model_v[i]->update_coefficients(grid, matrices.at(ss.str()), dt);
    }    
}

/*-------------------------------------------------*/
PointDataMap ModelStokes::to_point_data(std::shared_ptr<VectorInterface const> v, std::shared_ptr<GridInterface const> grid) const
{
    int dim = grid->dim();
    auto pv = std::dynamic_pointer_cast<StokesVector const >(v);  
    PointDataMap pdmap = _model_p->to_point_data(*(pv->get_p()), grid);
    for(int i=0;i<dim;i++)
    {
      PointDataMap pdv = _model_v[i]->to_point_data(*(pv->get_v(i)), grid);
      pdmap.insert(pdv.begin(), pdv.end());  
    }
    return pdmap;
}
/*-------------------------------------------------*/
std::map<std::string,double> ModelStokes::compute_error(std::shared_ptr<VectorInterface const> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<ApplicationInterface const> app) const
{
    std::map<std::string,double> err;
    int dim = grid->dim();
    // auto pf = std::dynamic_pointer_cast<SystemAnalyticalFunction const>(sol);
    // assert(pf);
    auto pv = std::dynamic_pointer_cast<StokesVector const>(v);  
    assert(pv);
    // armavec err(dim+1);
    err["p"] = _model_p->compute_error(*pv->get_p(), grid, app->solution("p")).at("p");
    for(int i=0;i<dim;i++)
    {
        std::stringstream ss;
        ss << "v"<<i;
        err[ss.str()] = _model_v[i]->compute_error(*pv->get_v(i), grid, app->solution(ss.str())).at(ss.str());
    }  
    return err;  
}

/*-------------------------------------------------*/
void ModelStokes::rhs(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<ApplicationInterface const> app) const
{
    int dim = grid->dim();
    // auto pf = std::dynamic_pointer_cast<SystemAnalyticalFunction const>(fct);
    // assert(pf);
    auto pv = std::dynamic_pointer_cast<StokesVector>(v);  
    assert(pv);
    // std::cerr << pv->n_vectors()<<"\n";
    // std::cerr << pv->get_size()<<"\n";
    assert(pv->get_p());
    // std::cerr << pv->get_p()->get_arma().n_elem<<"\n";
    _model_p->rhs(*pv->get_p(), grid, app->rhs("p"));
    for(int i=0;i<dim;i++)
    {
        std::stringstream ss;
        ss << "v"<<i;
      _model_v[i]->rhs(*pv->get_v(i), grid, app->rhs(ss.str()));
    }    
}
/*-------------------------------------------------*/
void ModelStokes::boundary_zero(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid) const
{
    int dim = grid->dim();
    auto pv = std::dynamic_pointer_cast<StokesVector>(v);  
    _model_p->boundary_zero(*pv->get_p(), grid);
    for(int i=0;i<dim;i++)
    {
      _model_v[i]->boundary_zero(*pv->get_v(i), grid);
    }    
}
/*-------------------------------------------------*/
void ModelStokes::boundary_linear(std::shared_ptr<VectorInterface> v, std::shared_ptr<GridInterface const> grid) const
{
    int dim = grid->dim();
    auto pv = std::dynamic_pointer_cast<StokesVector>(v);  
    _model_p->boundary_linear(*pv->get_p(), grid);
    for(int i=0;i<dim;i++)
    {
      _model_v[i]->boundary_linear(*pv->get_v(i), grid);
    }        
}
