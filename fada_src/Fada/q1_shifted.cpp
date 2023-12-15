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
#include  "q1_shifted.hpp"


/*-------------------------------------------------*/
Q1shifted::Q1shifted(std::string varname, const ParameterMap& parameters, std::shared_ptr<ApplicationInterface const> app) : ModelBase(varname, parameters, app)
{
    // _idir = std::stoi(parameters.at("direction"));
    // _dt = std::stod(parameters.at("dt"));
    _idir = parameters.get<int>("direction");
    _dt = parameters.get<double>("dt");
    _periodic = app->get_bc()->all("per");
    // std::cerr << _varname << " " << "_periodic = "  << _periodic << "\n";
    // std::cerr << "_idir = "  << _idir << " _dt = "  << _dt << "\n";
}

/*-------------------------------------------------*/
void Q1shifted::update_coefficients(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix, double dt)
{
    // std::cerr << "dt="<<dt << " _dt=" <<_dt<<"\n";
    // auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    // double h = ug->dx(0);
    auto spmatrix = std::dynamic_pointer_cast<Q1shifted::SparseMatrixDef>(matrix);
    assert(spmatrix);
    // spmatrix->add_diagonal((1/dt-1/_dt)*h*h);
    _dt = dt;
    arma::umat locations;
    armavec values;
    if(_periodic)
    {
        get_locations_values_matrix_per(locations, values, grid);        
    }
    else
    {
        get_locations_values_matrix_dir(locations, values, grid);        
    }
    spmatrix->set_elements(locations, values);
}

/*-------------------------------------------------*/
std::shared_ptr<VectorInterface> Q1shifted::newVector(std::shared_ptr<GridInterface const> grid) const
{
    armaicvec n_for_index = grid->n()-1;
    if(_periodic)
    {
        auto bc = std::make_shared<BoundaryConditionsBool>(_app->get_bc());
        return std::make_shared<Vector<GridVector>>(n_for_index, bc, true);
    }
    else
    {
        n_for_index[_idir] += 1;
        // std::shared_ptr<BoundaryConditions>& get_bc() {return _boundaryconditions;}
        auto bc = std::make_shared<BoundaryConditionsBool>(_app->get_bc());
        for(int i=0;i<grid->dim();i++)
        {
            if(i==_idir) continue;
            (*bc)[i][0] = false;
            (*bc)[i][1] = false;
        }
        return std::make_shared<Vector<GridVector>>(n_for_index, bc);        
    }
}
/*-------------------------------------------------*/
std::shared_ptr<SmootherInterface> Q1shifted::newSmoother(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface> matrix) const
{
    return std::make_shared<Smoother<SmootherSimple<SparseMatrix>,Vector<armavec>>>(_parameters.get<std::string>("smoother"), matrix);
}
/*-------------------------------------------------*/
std::shared_ptr<CoarseSolverInterface> Q1shifted::newCoarseSolver(std::string type, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const
{
    return std::make_shared<CoarseSolver<SolverUmf, Vector<armavec>>>(matrix,_parameters.get<std::string>("coarsesolver"));
}

/*-------------------------------------------------*/
std::shared_ptr<MatrixInterface> Q1shifted::newMatrix(std::shared_ptr<GridInterface const> grid) const
{
    arma::umat locations;
    armavec values;
    if(_periodic)
    {
        get_locations_values_matrix_per(locations, values, grid);        
    }
    else
    {
        get_locations_values_matrix_dir(locations, values, grid);        
    }
    auto p = std::make_shared<Q1shifted::SparseMatrixDef>(locations, values);
    // p->save(std::cerr);
    return p;
}

/*-------------------------------------------------*/
std::shared_ptr<TransferInterface> Q1shifted::newTransfer(std::shared_ptr<GridInterface const> grid, int ref_factor) const
{
    arma::umat locations;
    armavec values;
    if(_periodic)
    {
        get_locations_values_transfer_per(locations, values, grid);        
    }
    else
    {
        get_locations_values_transfer_dir(locations, values, grid);        
    }
    auto p = std::make_shared<Transfer<TransferByMatrix,Vector<GridVector>>>(locations, values);
    return p;    
}

/*-------------------------------------------------*/
std::shared_ptr <MatrixInterface> Q1shifted::newMatrixDivergence(std::shared_ptr <GridInterface const>grid) const
{
    arma::umat locations;
    armavec values;
    if(_periodic)
    {
        get_locations_values_matrix_divergence_per(locations, values, grid);        
    }
    else
    {
        get_locations_values_matrix_divergence_dir(locations, values, grid);        
    }
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
PointDataMap Q1shifted::to_point_data(std::shared_ptr<GridVector const> v, std::shared_ptr<GridInterface const> grid) const
{
    if(_periodic)
    {
        return to_point_data_per(*v, grid);
    }
    return to_point_data_dir(*v, grid);
}

/*-------------------------------------------------*/
PointDataMap Q1shifted2d::to_point_data_dir(const GridVector& v, std::shared_ptr<GridInterface const> grid) const
{
    // std::shared_ptr<BoundaryConditions const> bc = _application->get_bc();
    const std::vector<std::vector<FunctionMap>>& bf = _app->get_bc()->get_bf();
    const armaicvec& n = grid->n();
    armaicvec n_for_index = grid->n()-1;
    n_for_index[_idir] += 1;
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    assert(ug);
    PointDataMap pdmap;
    auto p_intp = std::make_shared<armavec>(arma::prod(n), arma::fill::zeros);
    pdmap[_varname] = p_intp;
     
    if(_idir==0)
    {
        int nx(n[0]), ny(n[1]);
        GridIndex IC(n_for_index), IN(n);
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 1; iy < ny-1; iy++)
            {
                int in = IN(ix,iy);
                int i1 = IC(ix, iy);
                int i2 = IC(ix, iy-1);            
                p_intp->at(in) = 0.5 * (v[i1]+v[i2]);
            }
            int in = IN(ix, 0);
            // int i1 = IC(ix, 0);
            // int i2 = IC(ix, 1);
            // int i3 = IC(ix, 2);
            // p_intp->at(in) = (15./8.)*v[i1]-(5./4.)*v[i2]+(3./8.)*v[i3];
            p_intp->at(in) = (*bf[1][0].at(_varname))(ug->x(ix,0), ug->y(ix,0));
            in = IN(ix, ny-1);
            // i1 = IC(ix, ny-2);
            // i2 = IC(ix, ny-3);
            // i3 = IC(ix, ny-4);
            p_intp->at(in) = (*bf[1][1].at(_varname))(ug->x(ix,ny-1), ug->y(ix,ny-1));
            // p_intp->at(in) = (15./8.)*v[i1]-(5./4.)*v[i2]+(3./8.)*v[i3];
        }            
    }
    else if(_idir==1)
    {
        int nx(n[0]), ny(n[1]);
        GridIndex IC(n_for_index), IN(n);
        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 1; ix < nx-1; ix++)
            {
                int in = IN(ix, iy);
                int i1 = IC(ix, iy);
                int i2 = IC(ix-1, iy);            
                p_intp->at(in) = 0.5 * (v[i1]+v[i2]);
            }
            int in = IN(0, iy);
            int i1 = IC(0, iy);
            int i2 = IC(1, iy);            
            int i3 = IC(2, iy);            
            p_intp->at(in) = (15./8.)*v[i1]-(5./4.)*v[i2]+(3./8.)*v[i3];
            in = IN(nx-1, iy);
            i1 = IC(nx-2, iy);
            i2 = IC(nx-3, iy);            
            i3 = IC(nx-4, iy);            
            p_intp->at(in) = (15./8.)*v[i1]-(5./4.)*v[i2]+(3./8.)*v[i3];
        }                        
    }
    else
    {
        _not_written_();
    }

    return pdmap;
}


/*-------------------------------------------------*/
PointDataMap Q1shifted2d::to_point_data_per(const GridVector& v, std::shared_ptr<GridInterface const> grid) const
{
    const armaicvec& n = grid->n();
    armaicvec n_for_index = grid->n()-1;
    PointDataMap pdmap;
    auto p_intp = std::make_shared<armavec>(arma::prod(n), arma::fill::zeros);
    pdmap[_varname] = p_intp;
     
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
void Q1shifted2d::get_locations_values_transfer_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    const armaicvec& n = grid->n();
    armaicvec n_for_index = grid->n()-1;
    n_for_index[_idir] += 1;
    int nx(n_for_index[0]), ny(n_for_index[1]);
    
    int size = (nx-2)*( 12*(ny-2) + 18) + 16*(ny-2) + 24;    
    if(_idir==1)
    {
        size = (ny-2)*( 12*(nx-2) + 18) + 16*(nx-2) + 24; 
    }

    armaicvec n_fine = 2*(grid->n()-1);
    if(not _periodic) n_fine[_idir] += 1;

    GridIndex IC(n_for_index), IF(n_fine);
    
    int i, j;
    Construct_Elements_Matrix ce(locations, values, size);
    double d1(3./4.), d2(1./4.), d3(3./8.), d4(1./8.);
    if(_idir==0)
    {
        for (int ix = 1; ix < nx-1; ix++)
        {
            for (int iy = 1; iy < ny-1; iy++)
            {
                i = IC(ix,iy);
                
                ce.add(i, IF(2*ix  , 2*iy  ), d1);
                ce.add(i, IF(2*ix  , 2*iy+1), d1);
                ce.add(i, IF(2*ix  , 2*iy-1), d2);
                ce.add(i, IF(2*ix  , 2*iy+2), d2);
                
                ce.add(i, IF(2*ix-1, 2*iy  ), d3);
                ce.add(i, IF(2*ix-1, 2*iy+1), d3);
                ce.add(i, IF(2*ix-1, 2*iy-1), d4);
                ce.add(i, IF(2*ix-1, 2*iy+2), d4);
                
                ce.add(i, IF(2*ix+1, 2*iy  ), d3);
                ce.add(i, IF(2*ix+1, 2*iy+1), d3);
                ce.add(i, IF(2*ix+1, 2*iy-1), d4);
                ce.add(i, IF(2*ix+1, 2*iy+2), d4);            
            }
            int iy = 0;
            i = IC(ix,iy);
            
            ce.add(i, IF(2*ix  , 2*iy  ), d1);
            ce.add(i, IF(2*ix  , 2*iy+1), d1);
            ce.add(i, IF(2*ix  , 2*iy+2), d2);
            
            ce.add(i, IF(2*ix-1, 2*iy  ), d3);
            ce.add(i, IF(2*ix-1, 2*iy+1), d3);
            ce.add(i, IF(2*ix-1, 2*iy+2), d4);
            
            ce.add(i, IF(2*ix+1, 2*iy  ), d3);
            ce.add(i, IF(2*ix+1, 2*iy+1), d3);
            ce.add(i, IF(2*ix+1, 2*iy+2), d4);            
            
            iy = ny-1;
            i = IC(ix,iy);
            
            ce.add(i, IF(2*ix  , 2*iy  ), d1);
            ce.add(i, IF(2*ix  , 2*iy+1), d1);
            ce.add(i, IF(2*ix  , 2*iy-1), d2);
            
            ce.add(i, IF(2*ix-1, 2*iy  ), d3);
            ce.add(i, IF(2*ix-1, 2*iy+1), d3);
            ce.add(i, IF(2*ix-1, 2*iy-1), d4);
            
            ce.add(i, IF(2*ix+1, 2*iy  ), d3);
            ce.add(i, IF(2*ix+1, 2*iy+1), d3);
            ce.add(i, IF(2*ix+1, 2*iy-1), d4);
        }
        int ix=0;
        for (int iy = 1; iy < ny-1; iy++)
        {
            i = IC(ix,iy);
            
            ce.add(i, IF(2*ix  , 2*iy  ), d1);
            ce.add(i, IF(2*ix  , 2*iy+1), d1);
            ce.add(i, IF(2*ix  , 2*iy-1), d2);
            ce.add(i, IF(2*ix  , 2*iy+2), d2);
            
            ce.add(i, IF(2*ix+1, 2*iy  ), d3);
            ce.add(i, IF(2*ix+1, 2*iy+1), d3);
            ce.add(i, IF(2*ix+1, 2*iy-1), d4);
            ce.add(i, IF(2*ix+1, 2*iy+2), d4);            
        }
        int iy = 0;
        i = IC(ix,iy);
        
        ce.add(i, IF(2*ix  , 2*iy  ), d1);
        ce.add(i, IF(2*ix  , 2*iy+1), d1);
        ce.add(i, IF(2*ix  , 2*iy+2), d2);
        
        ce.add(i, IF(2*ix+1, 2*iy  ), d3);
        ce.add(i, IF(2*ix+1, 2*iy+1), d3);
        ce.add(i, IF(2*ix+1, 2*iy+2), d4);            
        
        iy = ny-1;
        i = IC(ix,iy);
        
        ce.add(i, IF(2*ix  , 2*iy  ), d1);
        ce.add(i, IF(2*ix  , 2*iy+1), d1);
        ce.add(i, IF(2*ix  , 2*iy-1), d2);
        
        ce.add(i, IF(2*ix+1, 2*iy  ), d3);
        ce.add(i, IF(2*ix+1, 2*iy+1), d3);
        ce.add(i, IF(2*ix+1, 2*iy-1), d4);
        
        ix=nx-1;
        for (int iy = 1; iy < ny-1; iy++)
        {
            i = IC(ix,iy);
            
            ce.add(i, IF(2*ix  , 2*iy  ), d1);
            ce.add(i, IF(2*ix  , 2*iy+1), d1);
            ce.add(i, IF(2*ix  , 2*iy-1), d2);
            ce.add(i, IF(2*ix  , 2*iy+2), d2);
            
            ce.add(i, IF(2*ix-1, 2*iy  ), d3);
            ce.add(i, IF(2*ix-1, 2*iy+1), d3);
            ce.add(i, IF(2*ix-1, 2*iy-1), d4);
            ce.add(i, IF(2*ix-1, 2*iy+2), d4);
        }
        iy = 0;
        i = IC(ix,iy);
        
        ce.add(i, IF(2*ix  , 2*iy  ), d1);
        ce.add(i, IF(2*ix  , 2*iy+1), d1);
        ce.add(i, IF(2*ix  , 2*iy+2), d2);
        
        ce.add(i, IF(2*ix-1, 2*iy  ), d3);
        ce.add(i, IF(2*ix-1, 2*iy+1), d3);
        ce.add(i, IF(2*ix-1, 2*iy+2), d4);
        
        iy = ny-1;
        i = IC(ix,iy);
        
        ce.add(i, IF(2*ix  , 2*iy  ), d1);
        ce.add(i, IF(2*ix  , 2*iy+1), d1);
        ce.add(i, IF(2*ix  , 2*iy-1), d2);
        
        ce.add(i, IF(2*ix-1, 2*iy  ), d3);
        ce.add(i, IF(2*ix-1, 2*iy+1), d3);
        ce.add(i, IF(2*ix-1, 2*iy-1), d4);
    }
    else if(_idir==1)
    {
        for (int ix = 1; ix < nx-1; ix++)
        {
            for (int iy = 1; iy < ny-1; iy++)
            {
                i = IC(ix,iy);
                
                ce.add(i, IF(2*ix  , 2*iy), d1);
                ce.add(i, IF(2*ix+1, 2*iy), d1);
                ce.add(i, IF(2*ix-1, 2*iy), d2);
                ce.add(i, IF(2*ix+2, 2*iy), d2);
                
                ce.add(i, IF(2*ix  , 2*iy-1), d3);
                ce.add(i, IF(2*ix+1, 2*iy-1), d3);
                ce.add(i, IF(2*ix-1, 2*iy-1), d4);
                ce.add(i, IF(2*ix+2, 2*iy-1), d4);
                
                ce.add(i, IF(2*ix  , 2*iy+1), d3);
                ce.add(i, IF(2*ix+1, 2*iy+1), d3);
                ce.add(i, IF(2*ix-1, 2*iy+1), d4);
                ce.add(i, IF(2*ix+2, 2*iy+1), d4);                
            }
            int iy=0;
            i = IC(ix,iy);
            
            ce.add(i, IF(2*ix  , 2*iy), d1);
            ce.add(i, IF(2*ix+1, 2*iy), d1);
            ce.add(i, IF(2*ix-1, 2*iy), d2);
            ce.add(i, IF(2*ix+2, 2*iy), d2);
            
            ce.add(i, IF(2*ix  , 2*iy+1), d3);
            ce.add(i, IF(2*ix+1, 2*iy+1), d3);
            ce.add(i, IF(2*ix-1, 2*iy+1), d4);
            ce.add(i, IF(2*ix+2, 2*iy+1), d4);                

            iy=ny-1;
            i = IC(ix,iy);
            
            ce.add(i, IF(2*ix  , 2*iy), d1);
            ce.add(i, IF(2*ix+1, 2*iy), d1);
            ce.add(i, IF(2*ix-1, 2*iy), d2);
            ce.add(i, IF(2*ix+2, 2*iy), d2);
            
            ce.add(i, IF(2*ix  , 2*iy-1), d3);
            ce.add(i, IF(2*ix+1, 2*iy-1), d3);
            ce.add(i, IF(2*ix-1, 2*iy-1), d4);
            ce.add(i, IF(2*ix+2, 2*iy-1), d4);
        }
        int ix=0;
        for (int iy = 1; iy < ny-1; iy++)
        {
            i = IC(ix,iy);
            
            ce.add(i, IF(2*ix  , 2*iy), d1);
            ce.add(i, IF(2*ix+1, 2*iy), d1);
            ce.add(i, IF(2*ix+2, 2*iy), d2);
            
            ce.add(i, IF(2*ix  , 2*iy-1), d3);
            ce.add(i, IF(2*ix+1, 2*iy-1), d3);
            ce.add(i, IF(2*ix+2, 2*iy-1), d4);
            
            ce.add(i, IF(2*ix  , 2*iy+1), d3);
            ce.add(i, IF(2*ix+1, 2*iy+1), d3);
            ce.add(i, IF(2*ix+2, 2*iy+1), d4);                
        }
        int iy=0;
        i = IC(ix,iy);
        
        ce.add(i, IF(2*ix  , 2*iy), d1);
        ce.add(i, IF(2*ix+1, 2*iy), d1);
        ce.add(i, IF(2*ix+2, 2*iy), d2);
        
        ce.add(i, IF(2*ix  , 2*iy+1), d3);
        ce.add(i, IF(2*ix+1, 2*iy+1), d3);
        ce.add(i, IF(2*ix+2, 2*iy+1), d4);                

        iy=ny-1;
        i = IC(ix,iy);
        
        ce.add(i, IF(2*ix  , 2*iy), d1);
        ce.add(i, IF(2*ix+1, 2*iy), d1);
        ce.add(i, IF(2*ix+2, 2*iy), d2);
        
        ce.add(i, IF(2*ix  , 2*iy-1), d3);
        ce.add(i, IF(2*ix+1, 2*iy-1), d3);
        ce.add(i, IF(2*ix+2, 2*iy-1), d4);
        
        ix = nx-1;
        for (int iy = 1; iy < ny-1; iy++)
        {
            i = IC(ix,iy);
            
            ce.add(i, IF(2*ix  , 2*iy), d1);
            ce.add(i, IF(2*ix+1, 2*iy), d1);
            ce.add(i, IF(2*ix-1, 2*iy), d2);
            
            ce.add(i, IF(2*ix  , 2*iy-1), d3);
            ce.add(i, IF(2*ix+1, 2*iy-1), d3);
            ce.add(i, IF(2*ix-1, 2*iy-1), d4);
            
            ce.add(i, IF(2*ix  , 2*iy+1), d3);
            ce.add(i, IF(2*ix+1, 2*iy+1), d3);
            ce.add(i, IF(2*ix-1, 2*iy+1), d4);
        }
        iy=0;
        i = IC(ix,iy);
        
        ce.add(i, IF(2*ix  , 2*iy), d1);
        ce.add(i, IF(2*ix+1, 2*iy), d1);
        ce.add(i, IF(2*ix-1, 2*iy), d2);
        
        ce.add(i, IF(2*ix  , 2*iy+1), d3);
        ce.add(i, IF(2*ix+1, 2*iy+1), d3);
        ce.add(i, IF(2*ix-1, 2*iy+1), d4);

        iy=ny-1;
        i = IC(ix,iy);
        
        ce.add(i, IF(2*ix  , 2*iy), d1);
        ce.add(i, IF(2*ix+1, 2*iy), d1);
        ce.add(i, IF(2*ix-1, 2*iy), d2);
        
        ce.add(i, IF(2*ix  , 2*iy-1), d3);
        ce.add(i, IF(2*ix+1, 2*iy-1), d3);
        ce.add(i, IF(2*ix-1, 2*iy-1), d4);
    }
}

/*-------------------------------------------------*/
void Q1shifted2d::boundary(GridVector& u, std::shared_ptr<GridInterface const> grid, std::shared_ptr<BoundaryConditions const> bc) const
{
    if(_periodic) return;
    
    const std::vector<std::vector<FunctionMap>>& bf = bc->get_bf();
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    assert(ug);
    armaicvec n_for_index = grid->n()-1;
    n_for_index[_idir] += 1;
    int nx(n_for_index(0)), ny(n_for_index(1));
    if(_idir==0)
    {
        if((*bc)[0][0]=="dir")
        {
            for(int iy=0;iy<ny;iy++)
            {
                u.at(0,iy) = (*bf[0][0].at(_varname))(ug->x(0,iy), ug->ymid(0,iy));            
            }                
        }            
        if((*bc)[0][1]=="dir")
        {
            for(int iy=0;iy<ny;iy++)
            {
                u.at(nx-1,iy) = (*bf[0][1].at(_varname))(ug->x(nx-1,iy), ug->ymid(nx-1,iy));            
            }                
        }            
        if((*bc)[1][0]=="dir")
        {
            for(int ix=1;ix<nx-1;ix++)
            {
                // std::cerr << "ix="<< ix << " iy=0 " << (*bf[1][0].at(_varname))(ug->x(ix,0), ug->y(ix,0)) <<"\n";
                u.at(ix,0) += 2.0*(*bf[1][0].at(_varname))(ug->x(ix,0), ug->y(ix,0));            
            }                
        }            
        if((*bc)[1][1]=="dir")
        {
            for(int ix=1;ix<nx-1;ix++)
            {
                // std::cerr << "ix="<< ix << " iy=ny-1 " << (*bf[1][1].at(_varname))(ug->x(ix,0), ug->y(ix,ny)) <<"\n";
                // std::cerr << "x="<< ug->x(ix,0) << " y= " << ug->y(ix,ny-1) <<"\n";
                // un de plus en iy !!!!!!!!!!!!!!!!!!!!!
                u.at(ix,ny-1) += 2.0*(*bf[1][1].at(_varname))(ug->x(ix,ny), ug->y(ix,ny));            
            }                
        }                    
    }
    else
    {
        if((*bc)[0][0]=="dir")
        {
            for(int iy=1;iy<ny-1;iy++)
            {
                u.at(0,iy) += 2*(*bf[0][0].at(_varname))(ug->x(0,iy), ug->y(0,iy));            
            }                
        }            
        if((*bc)[0][1]=="dir")
        {
            for(int iy=1;iy<ny-1;iy++)
            {
                u.at(nx-1,iy) += 2*(*bf[0][1].at(_varname))(ug->x(nx,iy), ug->y(nx,iy));            
            }                
        }            
        if((*bc)[1][0]=="dir")
        {
            for(int ix=0;ix<nx;ix++)
            {
                u.at(ix,0) = (*bf[1][0].at(_varname))(ug->xmid(ix,0), ug->y(ix,0));            
            }                
        }            
        if((*bc)[1][1]=="dir")
        {
            for(int ix=0;ix<nx;ix++)
            {
                u.at(ix,ny-1) = (*bf[1][1].at(_varname))(ug->xmid(ix,ny-1), ug->y(ix,ny-1));            
            }                
        }                    
    }
}

/*-------------------------------------------------*/
void Q1shifted2d::get_locations_values_matrix_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    const armaicvec& n = grid->n();
    armaicvec n_for_index = grid->n()-1;
    n_for_index[_idir] += 1;
    
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    double h = ug->dx(0);

    int nx(n_for_index(0)), ny(n_for_index(1));

    int size = 5*(nx-2)*(ny-2)+8*(nx-2)+2*ny;    
    if(_idir==1)
    {
        size = 5*(nx-2)*(ny-2)+8*(ny-2)+2*nx; 
    }
    // GridIndex I(n-1);
    GridIndex I(n_for_index);
    
    int i, j;
    double d1(4.0), d2(-1.0), d3(5.0);
    if(_dt)
    {
        d1 += h*h/_dt;
    }
    else
    {
        d1 += 0.0001*h*h*h*h;
    }
    Construct_Elements_Matrix ce(locations, values, size);
    if(_idir==0)
    {
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
            int iy=0;
            i = I(ix,iy);
            ce.add(i, I(ix  , iy  ), d3);
            ce.add(i, I(ix-1, iy  ), d2);
            ce.add(i, I(ix+1, iy  ), d2);
            ce.add(i, I(ix  , iy+1), d2);
            iy=ny-1;
            i = I(ix,iy);
            ce.add(i, I(ix  , iy  ), d3);
            ce.add(i, I(ix-1, iy  ), d2);
            ce.add(i, I(ix+1, iy  ), d2);
            ce.add(i, I(ix  , iy-1), d2);
        }
        int ix=0;
        for (int iy = 0; iy < ny; iy++)
        {
            i = I(ix,iy);
            ce.add(i, I(ix  , iy  ), 1.0);    
        }
        ix=nx-1;
        for (int iy = 0; iy < ny; iy++)
        {
            i = I(ix,iy);
            ce.add(i, I(ix  , iy  ), 1.0);    
        }
    }
    else
    {
        for (int iy = 1; iy < ny-1; iy++)
        {
            for (int ix = 1; ix < nx-1; ix++)
            {
                i = I(ix,iy);
                ce.add(i, I(ix  , iy  ), d1);
                ce.add(i, I(ix-1, iy  ), d2);
                ce.add(i, I(ix+1, iy  ), d2);
                ce.add(i, I(ix  , iy-1), d2);
                ce.add(i, I(ix  , iy+1), d2);
            }
            int ix=0;
            i = I(ix,iy);
            ce.add(i, I(ix  , iy  ), d3);
            ce.add(i, I(ix+1, iy  ), d2);
            ce.add(i, I(ix  , iy-1), d2);
            ce.add(i, I(ix  , iy+1), d2);
            ix=nx-1;
            i = I(ix,iy);
            ce.add(i, I(ix  , iy  ), d3);
            ce.add(i, I(ix-1, iy  ), d2);
            ce.add(i, I(ix  , iy-1), d2);
            ce.add(i, I(ix  , iy+1), d2);
        }
        int iy=0;
        for (int ix = 0; ix < nx; ix++)
        {
            i = I(ix,iy);
            ce.add(i, I(ix  , iy  ), 1.0);    
        }
        iy=ny-1;
        for (int ix = 0; ix < nx; ix++)
        {
            i = I(ix,iy);
            ce.add(i, I(ix  , iy  ), 1.0);    
        }
    }
}

/*-------------------------------------------------*/
void Q1shifted2d::get_locations_values_matrix_divergence_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    const armaicvec& n = grid->n();
    armaicvec n_for_index = grid->n()-1;
    if(not _periodic) n_for_index[_idir] += 1;
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    double h = ug->dx(0);
    
    int nx(n[0]), ny(n[1]);
    int size = 2*(nx-1)*(ny-1);    
    GridIndex I(n-1), J(n_for_index);
    
    int i, j;
    Construct_Elements_Matrix ce(locations, values, size);
    if(_idir==0)
    {
        for (int ix = 0; ix < nx-1; ix++)
        {
            for (int iy = 0; iy < ny-1; iy++)
            {
                i = I(ix,iy);
                ce.add(i, J(ix  , iy), -h);
                ce.add(i, J(ix+1, iy),  h);
            }
        }        
    }
    else if(_idir==1)
    {
        for (int ix = 0; ix < nx-1; ix++)
        {
            for (int iy = 0; iy < ny-1; iy++)
            {
                i = I(ix,iy);
                ce.add(i, J(ix, iy  ), -h);
                ce.add(i, J(ix, iy+1),  h);
            }
        }                
    }
    else
    {
        _not_written_();
    }
}


/*-------------------------------------------------*/
void Q1shifted2d::get_locations_values_transfer_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
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
void Q1shifted2d::get_locations_values_matrix_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
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
void Q1shifted2d::get_locations_values_matrix_divergence_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
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
void Q1shifted3d::get_locations_values_transfer_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
void Q1shifted3d::get_locations_values_matrix_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
void Q1shifted3d::get_locations_values_matrix_divergence_per(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
void Q1shifted3d::get_locations_values_transfer_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
void Q1shifted3d::get_locations_values_matrix_dir(arma::umat& Q1shiftedtions, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
void Q1shifted3d::get_locations_values_matrix_divergence_dir(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    _not_written_();
}

/*-------------------------------------------------*/
void Q1shifted2d::rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const
{
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    assert(ug);
    double vol = arma::prod(ug->dx());

    auto pc = std::dynamic_pointer_cast<ConstantFunction const>(fct);
    if(pc)
    {
        v.fill(vol*pc->get_constant());
        return;
    } 


    const armaicvec& n = grid->n();
    armaicvec n_for_index = grid->n()-1;
    if(not _periodic) n_for_index[_idir] += 1;
    
    int nx(n_for_index[0]), ny(n_for_index[1]);
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
/*-------------------------------------------------*/
std::map<std::string,double> Q1shifted2d::compute_error(const GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> sol) const
{
    std::map<std::string,double> err;
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    assert(ug);
    double vol = arma::prod(ug->dx());
    const armaicvec& n = grid->n();
    
    armaicvec n_for_index = grid->n()-1;
    if(not _periodic) n_for_index[_idir] += 1;
    int nx(n_for_index[0]), ny(n_for_index[1]);
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
    err[_varname] = sqrt(err_glob);    
    return err;
}
