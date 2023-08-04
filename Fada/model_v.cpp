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
#include  "model_v.hpp"


/*-------------------------------------------------*/
void ModelV::update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt)
{
    // std::cerr << "dt="<<dt << " _dt=" <<_dt<<"\n";
    // auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    // double h = ug->dx(0);
    auto spmatrix = std::dynamic_pointer_cast<ModelV::SparseMatrixDef>(matrix);
    assert(spmatrix);
    // spmatrix->add_diagonal((1/dt-1/_dt)*h*h);
    _dt = dt;
    arma::umat locations;
    armavec values;
    get_locations_values_matrix(locations, values, grid);
    spmatrix->set_elements(locations, values);
}

/*-------------------------------------------------*/
std::shared_ptr<VectorInterface> ModelV::newVector(std::shared_ptr<GridInterface const> grid) const
{
    // return std::make_shared<Vector<GridVector>>(grid->n()-1, _boundaryconditions, false);
    return std::make_shared<Vector<GridVector>>(grid->n()-1, _app->get_bc(), true);
}
/*-------------------------------------------------*/
std::shared_ptr<SmootherInterface> ModelV::newSmoother(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix) const
{
    return std::make_shared<Smoother<SmootherSimple<SparseMatrix>,Vector<armavec>>>(_smoother, matrix);
}
/*-------------------------------------------------*/
std::shared_ptr<CoarseSolverInterface> ModelV::newCoarseSolver(std::string type, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const
{
    return std::make_shared<CoarseSolver<SolverUmf, Vector<armavec>>>(matrix,_coarsesolver);
    // return std::shared_ptr<CoarseSolverInterface>(new CoarseSolver<SolverUmf, Vector<armavec>>(matrix,_coarsesolver));
}

/*-------------------------------------------------*/
std::shared_ptr<MatrixInterface> ModelV::newMatrix(std::shared_ptr<GridInterface const> grid) const
{
    arma::umat locations;
    armavec values;
    get_locations_values_matrix(locations, values, grid);
    auto p = std::make_shared<ModelV::SparseMatrixDef>(locations, values);
    return p;
}

/*-------------------------------------------------*/
std::shared_ptr<TransferInterface> ModelV::newTransfer(std::shared_ptr<GridInterface const> grid, int ref_factor) const
{
    arma::umat locations;
    armavec values;
    get_locations_values_transfer(locations, values, grid);        
    auto p = std::make_shared<Transfer<TransferByMatrix,Vector<GridVector>>>(locations, values);
    return p;    
}

/*-------------------------------------------------*/
std::shared_ptr <MatrixInterface> ModelV::newMatrixDivergence(std::shared_ptr <GridInterface const>grid) const
{
    arma::umat locations;
    armavec values;
    get_locations_values_matrix_divergence(locations, values, grid);
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
PointDataMap ModelV2d::to_point_data(const GridVector& v, std::shared_ptr<GridInterface const> grid) const
{
    const armaicvec& n = grid->n();
    // armavec p_intp(arma::prod(n), arma::fill::zeros);
    PointDataMap pdmap;
    std::stringstream ss;
    ss << "v" << _idir;
    auto p_intp = std::make_shared<armavec>(arma::prod(n), arma::fill::zeros);
    pdmap[ss.str()] = p_intp;
     
    if(_idir==0)
    {
        int nx(n[0]), ny(n[1]);
        int nxc(n[0]-1), nyc(n[1]-1);
        GridIndex IC(n-1), IN(n);
        for (int ix = 0; ix < nx-1; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                int in = IN(ix,iy);
                int iyc1 = iy%nyc;
                int iyc2 = ((iy-1)%nyc+nyc)%nyc;
            
                int i1 = IC(ix, iyc1);
                int i2 = IC(ix, iyc2);
            
                p_intp->at(in) = 0.5 * (v[i1]+v[i2]);
            }
        }
        //periodic
        for (int iy = 0; iy < ny; iy++)
        {
            int in = IN(nx-1,iy);
            int in2 = IN(0,iy);        
            p_intp->at(in) = p_intp->at(in2);
        }
    }
    else if(_idir==1)
    {
        int nx(n[0]), ny(n[1]);
        int nxc(n[0]-1), nyc(n[1]-1);
        GridIndex IC(n-1), IN(n);
        for (int ix = 0; ix < nx; ix++)
        {
            int ixc1 = ix%nxc;
            int ixc2 = ((ix-1)%nxc+nxc)%nxc;
            for (int iy = 0; iy < ny-1; iy++)
            {
                int in = IN(ix,iy);
            
                int i1 = IC(ixc1, iy);
                int i2 = IC(ixc1, iy);
            
                p_intp->at(in) = 0.5 * (v[i1]+v[i2]);
            }
        }
        for (int ix = 0; ix < nx; ix++)
        {
            int in = IN(ix,ny-1);
            int in2 = IN(ix,0);        
            p_intp->at(in) = p_intp->at(in2);
        }
    }
    else
    {
        _not_written_();
    }

    return pdmap;
}


/*-------------------------------------------------*/
void ModelV2d::get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    // std::cerr << "grid " << grid->toString() << "\n";
    const armaicvec& n = grid->n();
    int nx(n[0]-1), ny(n[1]-1);
    int nxf(2*n[0]-2), nyf(2*n[1]-2);
    // std::cerr << "nx " << nx << " ny"<< ny << "\n";
    int size = 12*nx*ny;    
    GridIndex IC(n-1), IF(2*n-2);
    
    int i, j;
    Construct_Elements_Matrix ce(locations, values, size);
    double d1(3./4.), d2(1./4.), d3(3./8.), d4(1./8.);
    if(_idir==0)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                i = IC(ix,iy);
                
                ce.add(i, IF(2*ix  , 2*iy  ), d1);
                ce.add(i, IF(2*ix  , 2*iy+1), d1);
                ce.add(i, IF(2*ix  , (nyf+(2*iy-1)%nyf)%nyf), d2);
                ce.add(i, IF(2*ix  , (nyf+(2*iy+2)%nyf)%nyf), d2);
                
                ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, 2*iy  ), d3);
                ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, 2*iy+1), d3);
                ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, (nyf+(2*iy-1)%nyf)%nyf), d4);
                ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, (nyf+(2*iy+2)%nyf)%nyf), d4);
                
                ce.add(i, IF((nxf+(2*ix+1)%nxf)%nxf, 2*iy  ), d3);
                ce.add(i, IF((nxf+(2*ix+1)%nxf)%nxf, 2*iy+1), d3);
                ce.add(i, IF((nxf+(2*ix+1)%nxf)%nxf, (nyf+(2*iy-1)%nyf)%nyf), d4);
                ce.add(i, IF((nxf+(2*ix+1)%nxf)%nxf, (nyf+(2*iy+2)%nyf)%nyf), d4);
            
            }
        }
    }
    else if(_idir==1)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                i = IC(ix,iy);
                
                ce.add(i, IF(2*ix  , 2*iy), d1);
                ce.add(i, IF(2*ix+1, 2*iy), d1);
                ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, 2*iy), d2);
                ce.add(i, IF((nxf+(2*ix+2)%nxf)%nxf, 2*iy), d2);
                
                ce.add(i, IF(2*ix                  , (nyf+(2*iy-1)%nyf)%nyf), d3);
                ce.add(i, IF(2*ix+1                , (nyf+(2*iy-1)%nyf)%nyf), d3);
                ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, (nyf+(2*iy-1)%nyf)%nyf), d4);
                ce.add(i, IF((nxf+(2*ix+2)%nxf)%nxf, (nyf+(2*iy-1)%nyf)%nyf), d4);
                
                ce.add(i, IF(2*ix                  , (nyf+(2*iy+1)%nyf)%nyf), d3);
                ce.add(i, IF(2*ix+1                , (nyf+(2*iy+1)%nyf)%nyf), d3);
                ce.add(i, IF((nxf+(2*ix-1)%nxf)%nxf, (nyf+(2*iy+1)%nyf)%nyf), d4);
                ce.add(i, IF((nxf+(2*ix+2)%nxf)%nxf, (nyf+(2*iy+1)%nyf)%nyf), d4);
                
            }
        }
    }
    else
    {
        _not_written_();
    }
}

/*-------------------------------------------------*/
void ModelV2d::get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
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
    double d1(4.0), d2(-1.0);
    if(_dt)
    {
        d1 += h*h/_dt;
    }
    else
    {
        d1 += 0.0001*h*h*h*h;
    }
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
void ModelV2d::get_locations_values_matrix_divergence(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    double h = ug->dx(0);
    const armaicvec& n = grid->n();
    int nx(n[0]-1), ny(n[1]-1);
    assert(nx>=3);
    assert(ny>=3);
    int size = 2*nx*ny;    
    GridIndex I(n-1);
    
    int i, j;
    Construct_Elements_Matrix ce(locations, values, size);
    if(_idir==0)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                i = I(ix,iy);
                ce.add(i, I(ix       , iy), -h);
                ce.add(i, I((ix+1)%nx, iy),  h);
            }
        }        
    }
    else if(_idir==1)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                i = I(ix,iy);
                ce.add(i, I(ix, iy       ), -h);
                ce.add(i, I(ix, (iy+1)%ny),  h);
            }
        }                
    }
    else
    {
        _not_written_();
    }
}

/*-------------------------------------------------*/
void ModelV3d::get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
void ModelV3d::get_locations_values_matrix(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
void ModelV3d::get_locations_values_matrix_divergence(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
// #define  A 1
// #define  B 2
// #define  C 1
// #define  D 2
// double u_fct_v_0(double x, double y) {return -2*B*M_PI*sin(A*2*M_PI*x)*cos(B*2*M_PI*y);}
// double u_fct_v_1(double x, double y) {return  2*A*M_PI*cos(A*2*M_PI*x)*sin(B*2*M_PI*y);}
//
// // double rhs_fct_v_0(double x, double y) {return -2*B*M_PI*(A*A+B*B)*4*M_PI*M_PI*sin(A*2*M_PI*x)*cos(B*2*M_PI*y) - 2*C*M_PI*sin(C*2*M_PI*x)*cos(D*2*M_PI*y);}
// // double rhs_fct_v_1(double x, double y) {return  2*A*M_PI*(A*A+B*B)*4*M_PI*M_PI*cos(A*2*M_PI*x)*sin(B*2*M_PI*y) - 2*D*M_PI*cos(C*2*M_PI*x)*sin(D*2*M_PI*y);}
// // double rhs_fct_v_0(double x, double y) {return -2*B*M_PI*(A*A+B*B)*4*M_PI*M_PI*sin(A*2*M_PI*x)*cos(B*2*M_PI*y);}
// double rhs_fct_v_0(double x, double y) {return -2*B*M_PI*(A*A+B*B)*4*M_PI*M_PI*sin(A*2*M_PI*x)*cos(B*2*M_PI*y)-2*M_PI*sin(2*M_PI*x)*cos(2*M_PI*y);}
// double rhs_fct_v_1(double x, double y) {return  2*A*M_PI*(A*A+B*B)*4*M_PI*M_PI*cos(A*2*M_PI*x)*sin(B*2*M_PI*y)-2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y);}

void ModelV2d::rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const
{
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    assert(ug);
    double vol = arma::prod(ug->dx());
    const armaicvec& n = grid->n();
    int nx(n[0]-1), ny(n[1]-1);
    if(_idir==0)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                v.at(ix,iy) = vol*(*fct)(ug->x(ix,iy), ug->ymid(ix,iy));
                // v.at(ix,iy) = vol*(*fct)(ug->x(ix,iy), ug->y(ix,iy));
            }
        }    
    }
    else if(_idir==1)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                v.at(ix,iy) = vol*(*fct)(ug->xmid(ix,iy), ug->y(ix,iy));
                // v.at(ix,iy) = vol*(*fct)(ug->x(ix,iy), ug->y(ix,iy));
            }
        }
    }
}
std::map<std::string,double> ModelV2d::compute_error(const GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> sol) const
{
    std::map<std::string,double> err;
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    assert(ug);
    double vol = arma::prod(ug->dx());
    const armaicvec& n = grid->n();
    int nx(n[0]-1), ny(n[1]-1);
    // GridIndex I(n-1);
    double err_glob=0.0;
    if(_idir==0)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                double err_loc = v.at(ix,iy) - (*sol)(ug->x(ix,iy), ug->ymid(ix,iy));
                // double err_loc = v.at(ix,iy) - (*sol)(ug->x(ix,iy), ug->y(ix,iy));
                err_glob += vol*err_loc*err_loc;
            }
        }
    }
    else if(_idir==1)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                double err_loc = v.at(ix,iy) - (*sol)(ug->xmid(ix,iy), ug->y(ix,iy));
                // double err_loc = v.at(ix,iy) - (*sol)(ug->x(ix,iy), ug->y(ix,iy));
                err_glob += vol*err_loc*err_loc;
            }
        }
    }
    std::stringstream ss;
    ss << "v" << _idir;
    err[ss.str()] = sqrt(err_glob);    
    return err;
}
