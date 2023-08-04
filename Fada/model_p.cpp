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
#include  "analyticalfunctioninterface.hpp"
#include  "coarsesolverinterface.hpp"
#include  "coarsesolver_arma.hpp"
#include  "construct_elements_matrix.hpp"
#include  "feminterface.hpp"
#include  "gridindex.hpp"
#include  "gridvector.hpp"
#include  "matrixinterface.hpp"
#include  "uniformgrid.hpp"
#include  "smoothersimple.hpp"
#include  "solverumf.hpp"
#include  "sparsematrix.hpp"
#include  "sparsematrix_arma.hpp"
#include  "transferbymatrix.hpp"
#include  "transferinterface.hpp"
#include  "model_p.hpp"


/*-------------------------------------------------*/
std::shared_ptr<VectorInterface> ModelP::newVector(std::shared_ptr<GridInterface const> grid) const
{
    // return std::make_shared<Vector<GridVector>>(grid->n()-1, _boundaryconditions, false);
    return std::make_shared<Vector<GridVector>>(grid->n()-1, _app->get_bc(), true);
}
/*-------------------------------------------------*/
std::shared_ptr<SmootherInterface> ModelP::newSmoother(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix) const
{
    return std::make_shared<Smoother<SmootherSimple<SparseMatrix>,Vector<armavec>>>(_smoother, matrix);
}
/*-------------------------------------------------*/
std::shared_ptr<CoarseSolverInterface> ModelP::newCoarseSolver(std::string type, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const
{
    return std::make_shared<CoarseSolver<SolverUmf, Vector<armavec>>>(matrix,_coarsesolver);
    // return std::shared_ptr<CoarseSolverInterface>(new CoarseSolver<SolverUmf, Vector<armavec>>(matrix,_coarsesolver));
}

/*-------------------------------------------------*/
std::shared_ptr<MatrixInterface> ModelP::newMatrix(std::shared_ptr<GridInterface const> grid) const
{
    arma::umat locations;
    armavec values;
    get_locations_values_matrix(locations, values, grid);
    // for(int i=0;i<values.n_elem;i++)
    // {
    //     std::cerr << locations.at(0,i) << " " << locations.at(1,i) << " " << values.at(i) << "\n";
    // }
    // std::cerr << "locations=" << locations.t() <<"\n";
    // std::cerr << "values=" << values.t() <<"\n";
    typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;
    auto p = std::make_shared<SparseMatrixDef>(locations, values);
    // armavec diag = p->SparseMatrix::getDiag();
    // std::cerr << "diag="<<diag.n_elem<<"\n";
    // std::cerr << "diag="<<diag.t()<<"\n";
    // std::cerr << "mat="; p->save(std::cerr); std::cerr<<"\n";
    // armavec u(p->nrows(), arma::fill::ones), v(p->nrows(), arma::fill::zeros);
    // p->get().dot(v,u);
    // std::cerr << "v=" << v.t()<<"\n";
    return p;
}

/*-------------------------------------------------*/
std::shared_ptr<TransferInterface> ModelP::newTransfer(std::shared_ptr<GridInterface const> grid, int ref_factor) const
{
    arma::umat locations;
    armavec values;
    if(ref_factor==3)
    {
        get_locations_values_transfer_3(locations, values, grid);        
    }
    else
    {
        get_locations_values_transfer_2(locations, values, grid);        
    }
    // std::cerr << "locations=" << locations.t() <<"\n";
    // std::cerr << "values=" << values.t() <<"\n";
    auto p = std::make_shared<Transfer<TransferByMatrix,Vector<GridVector>>>(locations, values);
    // std::cerr << "trans="; p->save(std::cerr); std::cerr<<"\n";
    //
    // std::cerr << "ncols=" << p->getMatrix().ncols()<<"\n";
    // std::cerr << "nrows=" << p->getMatrix().nrows()<<"\n";
    // armavec u(p->getMatrix().nrows(), arma::fill::ones), v(p->getMatrix().ncols(), arma::fill::zeros);
    // p->getMatrix().Tdot(v,u);
    // std::cerr << "v=" << v.t()<<"\n";
    
    return p;    
}
/*-------------------------------------------------*/
PointDataMap ModelP2d::to_point_data(const GridVector& v, std::shared_ptr<GridInterface const> grid) const
{
    const armaicvec& n = grid->n();
    int nx(n[0]), ny(n[1]);
    int nxc(n[0]-1), nyc(n[1]-1);
    GridIndex IC(n-1), IN(n);
    PointDataMap pdmap;
    pdmap["p"] = std::make_shared<armavec>(arma::prod(n), arma::fill::zeros);
    // armavec p_intp(arma::prod(n), arma::fill::zeros);
    for (int ix = 0; ix < nx; ix++)
    {
        int ixc1 = ix%nxc;
        int ixc2 = ((ix-1)%nxc+nxc)%nxc;
        for (int iy = 0; iy < ny; iy++)
        {
            int in = IN(ix,iy);
            int iyc1 = iy%nyc;
            int iyc2 = ((iy-1)%nyc+nyc)%nyc;
            
            int i1 = IC(ixc1, iyc1);
            int i2 = IC(ixc1, iyc2);
            int i3 = IC(ixc2, iyc1);
            int i4 = IC(ixc2, iyc2);
            
            // p_intp.at(in) = 0.25 * (v[i1]+v[i2]+v[i3]+v[i4]);
            pdmap["p"]->at(in) = 0.25 * (v[i1]+v[i2]+v[i3]+v[i4]);
        }
    }
    return pdmap;
    // return p_intp;
}


/*-------------------------------------------------*/
void ModelP2d::get_locations_values_transfer_3(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    // std::cerr << "grid " << grid->toString() << "\n";
    const armaicvec& n = grid->n();
    int nx(n[0]-1), ny(n[1]-1);
    // std::cerr << "nx " << nx << " ny"<< ny << "\n";
    // int size = 9*(nx-2)*(ny-2) + 6*2*( nx-2 + ny-2) + 4*4;
    int size = 9*nx*ny;
    GridIndex IC(n-1), IF(3*(n-1));
    
    int i, j;
    Construct_Elements_Matrix ce(locations, values, size);
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            i = IC(ix,iy);
            ce.add(i, IF(3*ix+1, 3*iy+1), 1);
            
            ce.add(i, IF(3*ix+1, 3*iy  ), 0.5);
            ce.add(i, IF(3*ix+1, 3*iy+2), 0.5);
            ce.add(i, IF(3*ix  , 3*iy+1), 0.5);
            ce.add(i, IF(3*ix+2, 3*iy+1), 0.5);
            
            ce.add(i, IF(3*ix  , 3*iy  ), 0.25);
            ce.add(i, IF(3*ix+2, 3*iy  ), 0.25);
            ce.add(i, IF(3*ix  , 3*iy+2), 0.25);
            ce.add(i, IF(3*ix+2, 3*iy+2), 0.25);
        }
    }
}
/*-------------------------------------------------*/
void ModelP2d::get_locations_values_transfer_2(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    // std::cerr << "grid " << grid->toString() << "\n";
    const armaicvec& n = grid->n();
    int nx(n[0]-1), ny(n[1]-1);
    int nxf(2*n[0]-2), nyf(2*n[1]-2);
    // std::cerr << "nx " << nx << " ny"<< ny << "\n";
    int size = 16*nx*ny;    
    GridIndex IC(n-1), IF(2*n-2);
    
    int i, j;
    Construct_Elements_Matrix ce(locations, values, size);
    double d1(9./16.), d2(3./16.), d3(1./16.);
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            i = IC(ix,iy);
            ce.add(i, IF(2*ix  , 2*iy  ), d1);
            ce.add(i, IF(2*ix+1, 2*iy  ), d1);
            ce.add(i, IF(2*ix  , 2*iy+1), d1);
            ce.add(i, IF(2*ix+1, 2*iy+1), d1);
            
            ce.add(i, IF(2*ix  , (nyf+(2*iy-1)%nyf)%nyf), d2);
            ce.add(i, IF(2*ix+1, (nyf+(2*iy-1)%nyf)%nyf), d2);
            ce.add(i, IF(2*ix  , (nyf+(2*iy+2)%nyf)%nyf), d2);
            ce.add(i, IF(2*ix+1, (nyf+(2*iy+2)%nyf)%nyf), d2);

            ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, 2*iy  ), d2);
            ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, 2*iy+1), d2);
            ce.add(i, IF((nxf+(2*ix+2)%nxf)%nxf, 2*iy  ), d2);
            ce.add(i, IF((nxf+(2*ix+2)%nxf)%nxf, 2*iy+1), d2);
            
            ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, (nyf+(2*iy-1)%nyf)%nyf), d3);
            ce.add(i, IF((nxf+(2*ix+2)%nxf)%nxf, (nyf+(2*iy-1)%nyf)%nyf), d3);
            ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, (nyf+(2*iy+2)%nyf)%nyf), d3);
            ce.add(i, IF((nxf+(2*ix+2)%nxf)%nxf, (nyf+(2*iy+2)%nyf)%nyf), d3);
            
        }
    }
    //bdry -sides
}
/*-------------------------------------------------*/
void ModelP2d::get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    double h = ug->dx(0);
    const armaicvec& n = grid->n();
    int nx(n[0]-1), ny(n[1]-1);
    assert(nx>=3);
    assert(ny>=3);
    int size = 5*nx*ny;    
    GridIndex I(n-1);
    
    int i, j;
    double d1(4.0+0.01*h*h*h), d2(-1.0);
    Construct_Elements_Matrix ce(locations, values, size);
    for (int ix = 1; ix < nx-1; ix++)
    {
        for (int iy = 1; iy < ny-1; iy++)
        {
            i = I(ix,iy);
            ce.add(i, I(ix  , iy  ), d1);
            ce.add(i, I(ix-1, iy  ), d2);
            ce.add(i, I(ix+1, iy  ), d2);
            ce.add(i, I(ix  , iy-1), d2);
            ce.add(i, I(ix  , iy+1), d2);
        }
    }
    int ix = 0;
    for (int iy = 1; iy < ny-1; iy++)
    {
        i = I(ix,iy);
        ce.add(i, I(ix  , iy  ),  d1);
        ce.add(i, I(nx-1, iy  ),  d2);
        ce.add(i, I(ix+1, iy  ),  d2);
        ce.add(i, I(ix  , iy-1),  d2);
        ce.add(i, I(ix  , iy+1),  d2);
    }
    ix = nx-1;
    for (int iy = 1; iy < ny-1; iy++)
    {
        i = I(ix,iy);
        ce.add(i, I(ix  , iy  ),  d1);
        ce.add(i, I(ix-1, iy  ),  d2);
        ce.add(i, I(0   , iy  ),  d2);
        ce.add(i, I(ix  , iy-1),  d2);
        ce.add(i, I(ix  , iy+1),  d2);
    }
    int iy = 0;
    for (int ix = 1; ix < nx-1; ix++)
    {
        i = I(ix,iy);
        ce.add(i, I(ix  , iy  ),  d1);
        ce.add(i, I(ix-1, iy  ),  d2);
        ce.add(i, I(ix+1, iy  ),  d2);
        ce.add(i, I(ix  , ny-1),  d2);
        ce.add(i, I(ix  , iy+1),  d2);
    }
    iy = ny-1;
    for (int ix = 1; ix < nx-1; ix++)
    {
        i = I(ix,iy);
        ce.add(i, I(ix  , iy  ),  d1);
        ce.add(i, I(ix-1, iy  ),  d2);
        ce.add(i, I(ix+1, iy  ),  d2);
        ce.add(i, I(ix  , iy-1),  d2);
        ce.add(i, I(ix  , 0   ),  d2);
    }
    ix = 0; iy=0;
    i = I(ix,iy);
    ce.add(i, I(ix  , iy  ),  d1);
    ce.add(i, I(nx-1, iy  ),  d2);
    ce.add(i, I(ix+1, iy  ),  d2);
    ce.add(i, I(ix  , ny-1),  d2);
    ce.add(i, I(ix  , iy+1),  d2);
    ix = 0; iy=ny-1;
    i = I(ix,iy);
    ce.add(i, I(ix  , iy  ),  d1);
    ce.add(i, I(nx-1, iy  ),  d2);
    ce.add(i, I(ix+1, iy  ),  d2);
    ce.add(i, I(ix  , iy-1),  d2);
    ce.add(i, I(ix  , 0   ),  d2);
    ix = nx-1; iy=0;
    i = I(ix,iy);
    ce.add(i, I(ix  , iy  ),  d1);
    ce.add(i, I(ix-1, iy  ),  d2);
    ce.add(i, I(0   , iy  ),  d2);
    ce.add(i, I(ix  , ny-1),  d2);
    ce.add(i, I(ix  , iy+1),  d2);
    ix = nx-1; iy=ny-1;
    i = I(ix,iy);
    ce.add(i, I(ix  , iy  ),  d1);
    ce.add(i, I(ix-1, iy  ),  d2);
    ce.add(i, I(0   , iy  ),  d2);
    ce.add(i, I(ix  , iy-1),  d2);
    ce.add(i, I(ix  , 0   ),  d2);
}

/*-------------------------------------------------*/
void ModelP3d::get_locations_values_transfer_2(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
void ModelP3d::get_locations_values_transfer_3(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
void ModelP3d::get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
// #define  C 1
// #define  D 2
// // double u_fct_p(double x, double y) {return cos(C*2*M_PI*x)*cos(D*2*M_PI*y);}
// // double u_fct_p(double x, double y) {return 0.0;}
// double u_fct_p(double x, double y) {return cos(2*M_PI*x)*cos(2*M_PI*y);}
// double rhs_fct_p(double x, double y) {return  (C*C+D*D)*4*M_PI*M_PI*cos(C*2*M_PI*x)*cos(D*2*M_PI*y);}


void ModelP2d::rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const
{
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    assert(ug);
    double vol = arma::prod(ug->dx());
    const armaicvec& n = grid->n();
    int nx(n[0]-1), ny(n[1]-1);
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            // v.at(ix,iy) = vol*rhs_fct_p(ug->xmid(ix,iy), ug->ymid(ix,iy));
            v.at(ix,iy) = vol* (*fct)(ug->xmid(ix,iy), ug->ymid(ix,iy));
        }
    }    
}
std::map<std::string,double> ModelP2d::compute_error(const GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> sol) const
{
    std::map<std::string,double> err;
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    assert(ug);
    double vol = arma::prod(ug->dx());
    const armaicvec& n = grid->n();
    int nx(n[0]-1), ny(n[1]-1);
    // GridIndex I(n-1);
    double err_glob=0.0;
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            // double err_loc = v.at(ix,iy) - u_fct_p(ug->xmid(ix,iy), ug->ymid(ix,iy));
            double err_loc = v.at(ix,iy) - (*sol)(ug->xmid(ix,iy), ug->ymid(ix,iy));
            err_glob += vol*err_loc*err_loc;
        }
    }
    err["p"] = sqrt(err_glob);
    return err;    
}

