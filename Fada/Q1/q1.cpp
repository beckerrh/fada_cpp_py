//
//  q1.cpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <cassert>
#include  <armadillo>
#include  "../coarsesolverinterface.hpp"
#include  "../coarsesolver_arma.hpp"
#include  "../construct_elements_matrix.hpp"
#include  "../feminterface.hpp"
#include  "../gridindex.hpp"
#include  "../gridvector.hpp"
#include  "../matrixinterface.hpp"
#include  "../uniformgrid.hpp"
#include  "../smoothersimple.hpp"
#include  "../smootherumf.hpp"
#include  "../sparsematrix.hpp"
#include  "../sparsematrix_arma.hpp"
#include  "../transferbymatrix.hpp"
#include  "../transferinterface.hpp"
#include  "q1.hpp"
#include  "transferq1.hpp"
#include  "stencil.hpp"

double lin2d(double x, double y) {return 3.0*x+2.0*y;}
double lin3d(double x, double y, double z) {return 3.0*x+2.0*y+z;}


/*-------------------------------------------------*/
void Q1::set_grid(std::shared_ptr<GridInterface const> grid)
{
    auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    assert(ug);
    assert(ug->dx().n_elem == ug->dim());
    assert(ug->n().n_elem == ug->dim());
    _nx = ug->nx();
    _ny = ug->ny();
    _nz = -1;
    if(ug->dim()==3) _nz = ug->nz();
    _vol = arma::prod(ug->dx());
    _ug = ug;
    if (_boundaryconditions == nullptr)
    {
        _boundaryconditions = std::make_shared<BoundaryConditions>(ug->dim());
    }
}
/*-------------------------------------------------*/
std::shared_ptr<VectorInterface> Q1::newVector(std::shared_ptr<GridInterface const> grid) const
{
    auto p = std::make_shared<Vector<GridVector>>(grid->n(), _boundaryconditions);
    p->boundary_zero();
    return p;
}
/*-------------------------------------------------*/
std::shared_ptr<SmootherInterface const> Q1::newSmoother(std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const
{
    auto stencil = std::dynamic_pointer_cast<FemAndMatrixAndSmootherInterface const>(matrix);
    if (_smoothertype=="stencil")
    {
        assert(stencil);
        return stencil;
    }
    std::shared_ptr<MatrixInterface const> matrixforumf = matrix;
    if(stencil)
    {
        arma::umat locations;
        armavec values;
        stencil->get_locations_values(locations, values);
        typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;
        matrixforumf = std::make_shared<SparseMatrixDef>(locations, values);
    }
    if(_matrixtype=="arma")
    {
        return std::make_shared<Smoother<SmootherSimple<SparseMatrix_arma>,Vector<armavec>>>(_smoother, matrixforumf);
    }
    return std::make_shared<Smoother<SmootherSimple<SparseMatrix>,Vector<armavec>>>(_smoother, matrixforumf);
}
/*-------------------------------------------------*/
std::shared_ptr<CoarseSolverInterface const> Q1::newCoarseSolver(std::string type, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> matrix) const
{
    // return std::shared_ptr<SmootherInterface>(new CoarseSolver<SmootherSimple,Vector<armavec>>("GS", matrix));
    auto armamat = std::dynamic_pointer_cast<Matrix<SparseMatrix_arma,Vector<armavec>> const>(matrix);
    if(armamat)
    {
        return std::make_shared<CoarseSolver<CoarseSolver_arma, Vector<armavec>>>(matrix,_coarsesolver);
    }
    std::shared_ptr<MatrixInterface const> matrixforumf;
    auto stencil = std::dynamic_pointer_cast<FemAndMatrixAndSmootherInterface const>(matrix);
    if(stencil)
    {
        arma::umat locations;
        armavec values;
        stencil->get_locations_values(locations, values);
        // std::cerr << "locations\n" << locations << "\n";
        // std::cerr << "values\n" << values << "\n";
        typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;
        matrixforumf = std::make_shared<SparseMatrixDef>(locations, values);
    }
    else
    {
        matrixforumf = matrix;
    }
    // return std::shared_ptr<SmootherInterface>(new Smoother<arma::sp_mat, Vector<armavec>>(matrix));
    return std::shared_ptr<CoarseSolverInterface>(new CoarseSolver<SmootherUmf, Vector<armavec>>(matrixforumf,_coarsesolver));
}

/*-------------------------------------------------*/
std::shared_ptr<MatrixInterface const> Q1::newMatrix(std::shared_ptr<GridInterface const> grid) const
{
    std::shared_ptr<FemAndMatrixAndSmootherInterface> p = newStencil(grid);
    if(_matrixtype=="stencil")
    {
        return p;
    }
    else
    {
        arma::umat locations;
        armavec values;
        p->get_locations_values(locations, values);
        if(_matrixtype=="arma")
        {
            return std::make_shared<Matrix<SparseMatrix_arma,Vector<armavec>>>(locations, values);
        }
        else
        {
            typedef Matrix<SparseMatrix,Vector<armavec>> SparseMatrixDef;
            return std::make_shared<SparseMatrixDef>(locations, values);
            // return std::shared_ptr<MatrixInterface>(new SparseMatrixDef(locations, values));
        }
    }
}

/*-------------------------------------------------*/
std::shared_ptr<TransferInterface const> Q1::newTransfer(std::shared_ptr<GridInterface const> grid) const
{
    if(_transfertype=="matrix")
    {
        arma::umat locations;
        armavec values;
        get_locations_values_transfer(locations, values, grid);
        // std::cerr << "locations=" << locations.t() <<"\n";
        // std::cerr << "values=" << values.t() <<"\n";
        return std::make_shared<Transfer<TransferByMatrix,Vector<GridVector>>>(locations, values);
    }
    else if(_transfertype=="stencil")
    {
        std::shared_ptr<UniformGrid const> ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
        assert(ug);
        const armaicvec& n = ug->n();
        const armavec& dx = ug->dx();
        if(ug->dim()==2)
        {
            return std::make_shared<Transfer<TransferQ12d,Vector<GridVector>>>(n, dx);
            // return std::shared_ptr<TransferInterface>(new Transfer<TransferQ12d,Vector<GridVector>>(n, dx, _boundaryconditions));
        }
        return std::make_shared<Transfer<TransferQ13d,Vector<GridVector>>>(n, dx);
        // return std::shared_ptr<TransferInterface>(new Transfer<TransferQ13d,Vector<GridVector>>(n, dx, _boundaryconditions));
    }
    else
    {
        _not_written_();
    }
    return nullptr;
}

/*-------------------------------------------------*/
void Q1::rhs_one(GridVector& v) const
{
    v.fill(_vol*12.0);
}
/*-------------------------------------------------*/
void Q1::rhs_random(GridVector& v) const
{
    arma::arma_rng::set_seed_random();
    // v.data().randu();
    v.randu();
    v *= 100*_vol;
}
/*-------------------------------------------------*/
void Q12d::get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    // std::cerr << "grid " << grid->toString() << "\n";
    const armaicvec& n = grid->n();
    int nx(n[0]), ny(n[1]);
    // std::cerr << "nx " << nx << " ny"<< ny << "\n";
    int size = 9*(nx-2)*(ny-2) + 6*2*( nx-2 + ny-2) + 4*4;    
    GridIndex IC(n), IF(2*n-1);
    
    int i, j;
    Construct_Elements_Matrix ce(locations, values, size);
    for (int ix = 1; ix < nx - 1; ix++)
    {
        for (int iy = 1; iy < ny - 1; iy++)
        {
            i = IC(ix,iy);
            ce.add(i, IF(2*ix  , 2*iy  ), 1);
            ce.add(i, IF(2*ix-1, 2*iy  ), 0.5);
            ce.add(i, IF(2*ix+1, 2*iy  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy-1), 0.5);
            ce.add(i, IF(2*ix  , 2*iy+1), 0.5);
            ce.add(i, IF(2*ix-1, 2*iy-1), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy+1), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy-1), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy+1), 0.25);
        }
    }
    //bdry -sides
    for (int ix = 1; ix < nx-1; ix++)
    {
        // iy=0
        i = IC(ix, 0);
        ce.add(i, IF(2*ix  , 0  ), 1);
        ce.add(i, IF(2*ix-1, 0  ), 0.5);
        ce.add(i, IF(2*ix+1, 0  ), 0.5);
        ce.add(i, IF(2*ix  , 1  ), 0.5);
        ce.add(i, IF(2*ix-1, 1  ), 0.25);
        ce.add(i, IF(2*ix+1, 1  ), 0.25);
      
        // iy=ny-1
        i = IC(ix, ny-1);
        ce.add(i, IF(2*ix  , 2*(ny-1)  ), 1);
        ce.add(i, IF(2*ix-1, 2*(ny-1)  ), 0.5);
        ce.add(i, IF(2*ix+1, 2*(ny-1)  ), 0.5);
        ce.add(i, IF(2*ix  , 2*(ny-1)-1), 0.5);
        ce.add(i, IF(2*ix-1, 2*(ny-1)-1), 0.25);
        ce.add(i, IF(2*ix+1, 2*(ny-1)-1), 0.25);
    }
    for (int iy = 1; iy < ny - 1; iy++)
    {
        // ix=0
        i = IC(0, iy);
        ce.add(i, IF(0, 2*iy  ), 1);
        ce.add(i, IF(0, 2*iy-1), 0.5);
        ce.add(i, IF(0, 2*iy+1), 0.5);
        ce.add(i, IF(1, 2*iy  ), 0.5);
        ce.add(i, IF(1, 2*iy-1), 0.25);
        ce.add(i, IF(1, 2*iy+1), 0.25);

        // ix=nx-1
        i = IC(nx-1, iy);
        ce.add(i, IF(2*(nx-1)  , 2*iy  ), 1);
        ce.add(i, IF(2*(nx-1)  , 2*iy-1), 0.5);
        ce.add(i, IF(2*(nx-1)  , 2*iy+1), 0.5);
        ce.add(i, IF(2*(nx-1)-1, 2*iy  ), 0.5);
        ce.add(i, IF(2*(nx-1)-1, 2*iy-1), 0.25);
        ce.add(i, IF(2*(nx-1)-1, 2*iy+1), 0.25);
    }
    //bdry - corners
    // ix=0 iy=0
    i = IC(0, 0);
    ce.add(i, IF(0, 0), 1);
    ce.add(i, IF(0, 1), 0.5);
    ce.add(i, IF(1, 0), 0.5);
    ce.add(i, IF(1, 1), 0.25);
    // ix=0 iy=ny-1
    i = IC(0, ny-1);
    ce.add(i, IF(0, 2*(ny-1)  ), 1);
    ce.add(i, IF(0, 2*(ny-1)-1), 0.5);
    ce.add(i, IF(1, 2*(ny-1)  ), 0.5);
    ce.add(i, IF(1, 2*(ny-1)-1), 0.25);
    // ix=nx-1 iy=0
    i = IC(nx-1, 0);
    ce.add(i, IF(2*(nx-1)  , 0), 1);
    ce.add(i, IF(2*(nx-1)  , 1), 0.5);
    ce.add(i, IF(2*(nx-1)-1, 0), 0.5);
    ce.add(i, IF(2*(nx-1)-1, 1), 0.25);
    // ix=nx-1 iy=ny-1
    i = IC(nx-1, ny-1);
    ce.add(i, IF(2*(nx-1)  , 2*(ny-1)  ), 1);
    ce.add(i, IF(2*(nx-1)  , 2*(ny-1)-1), 0.5);
    ce.add(i, IF(2*(nx-1)-1, 2*(ny-1)  ), 0.5);
    ce.add(i, IF(2*(nx-1)-1, 2*(ny-1)-1), 0.25);
}

/*-------------------------------------------------*/
void Q13d::get_locations_values_transfer(arma::umat& locations, armavec& values, std::shared_ptr <GridInterface const>grid) const
{
    // std::cerr << "grid " << grid->toString() << "\n";
    const armaicvec& n = grid->n();
    int nx(n[0]), ny(n[1]), nz(n[2]);
    // std::cerr << "nx " << nx << " ny"<< ny << "\n";
    int size = 27*(nx-2)*(ny-2)*(nz-2) + 18*2*( (ny-2)*(nz-2) + (nx-2)*(nz-2) + (ny-2)*(nz-2) ) + 4*12*((nx-2)+(ny-2)+(nz-2)) + 8*8;
    GridIndex IC(n), IF(2*n-1);
    int i, j;
    Construct_Elements_Matrix ce(locations, values, size);
    for (int ix = 1; ix < nx - 1; ix++)
    {
        for (int iy = 1; iy < ny - 1; iy++)
        {
            for (int iz = 1; iz < nz - 1; iz++)
            {
                i = IC(ix, iy, iz);
                ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
                ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
                ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
                ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
                ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
                ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
                ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
                ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
                ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
                ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
                ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
                ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
                ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
                ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
                ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
                ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
                ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
                ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
                ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
                
                ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz-1), 0.125);
                ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz-1), 0.125);
                ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz-1), 0.125);
                ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz-1), 0.125);
                ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz+1), 0.125);
                ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz+1), 0.125);
                ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz+1), 0.125);
                ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz+1), 0.125);
            }
        }
    }
    int ix = 0;
    for (int iy = 1; iy < ny - 1; iy++)
    {
        for (int iz = 1; iz < nz - 1; iz++)
        {
            i = IC(ix, iy, iz);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
            
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
            
            ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
            
            ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz+1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz+1), 0.125);
        }
    }
    ix = nx-1;
    for (int iy = 1; iy < ny - 1; iy++)
    {
        for (int iz = 1; iz < nz - 1; iz++)
        {
            i = IC(ix, iy, iz);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
            
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
            
            ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
            
            ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz+1), 0.125);
            ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz+1), 0.125);
        }
    }
    int iy=0;
    for (int ix = 1; ix < nx - 1; ix++)
    {
        for (int iz = 1; iz < nz - 1; iz++)
        {
            i = IC(ix, iy, iz);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
            ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
                
            ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz+1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz+1), 0.125);
        }
    }
    iy=ny-1;
    for (int ix = 1; ix < nx - 1; ix++)
    {
        for (int iz = 1; iz < nz - 1; iz++)
        {
            i = IC(ix, iy, iz);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
            ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
                
            ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz+1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz+1), 0.125);
        }
    }
    int iz=0;
    for (int ix = 1; ix < nx - 1; ix++)
    {
        for (int iy = 1; iy < ny - 1; iy++)
        {
            i = IC(ix, iy, iz);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
            ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
                
            ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz+1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz+1), 0.125);
            ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz+1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz+1), 0.125);
        }
    }
    iz=nz-1;
    for (int ix = 1; ix < nx - 1; ix++)
    {
        for (int iy = 1; iy < ny - 1; iy++)
        {
            i = IC(ix, iy, iz);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
            ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
                
            ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
            ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
            ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
            ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
                
            ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz-1), 0.125);
            ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz-1), 0.125);
        }
    }
    ix=0, iy=0;
    for (int iz = 1; iz < nz - 1; iz++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
                
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz-1), 0.125);
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz+1), 0.125);
    }
    ix=0, iy=ny-1;
    for (int iz = 1; iz < nz - 1; iz++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
                
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz-1), 0.125);
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz+1), 0.125);
    }
    ix=nx-1, iy=0;
    for (int iz = 1; iz < nz - 1; iz++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
                
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz-1), 0.125);
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz+1), 0.125);
    }
    ix=nx-1, iy=ny-1;
    for (int iz = 1; iz < nz - 1; iz++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
                
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz-1), 0.125);
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz+1), 0.125);
    }
    ix=0, iz=0;
    for (int iy = 1; iy < ny - 1; iy++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
                
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz+1), 0.125);
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz+1), 0.125);
    }
    ix=nx-1, iz=0;
    for (int iy = 1; iy < ny - 1; iy++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
                
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz+1), 0.125);
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz+1), 0.125);
    }
    ix=0, iz=nz-1;
    for (int iy = 1; iy < ny - 1; iy++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
                
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
                
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz-1), 0.125);
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz-1), 0.125);
    }
    ix=nx-1, iz=nz-1;
    for (int iy = 1; iy < ny - 1; iy++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
                
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
                
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz-1), 0.125);
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz-1), 0.125);
    }
    iy=0, iz=0;
    for (int ix = 1; ix < nx - 1; ix++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
                
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz+1), 0.125);
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz+1), 0.125);
    }
    iy=ny-1, iz=0;
    for (int ix = 1; ix < nx - 1; ix++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
                
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
                
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz+1), 0.125);
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz+1), 0.125);
    }
    iy=0, iz=nz-1;
    for (int ix = 1; ix < nx - 1; ix++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
                
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
                
        ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz-1), 0.125);
        ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz-1), 0.125);
    }
    iy=ny-1, iz=nz-1;
    for (int ix = 1; ix < nx - 1; ix++)
    {
        i = IC(ix, iy, iz);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
                
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
        ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
                
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
        ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
        ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
        ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
                
        ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz-1), 0.125);
        ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz-1), 0.125);
    }

    ix = 0, iy=0, iz=0;
    i = IC(ix, iy, iz);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
            
    ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
            
    ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
    ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
    ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
            
    ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz+1), 0.125);

    ix = nx-1, iy=0, iz=0;
    i = IC(ix, iy, iz);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
            
    ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
            
    ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
    ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
    ce.add(i, IF(2*ix  , 2*iy+1, 2*iz+1), 0.25);
            
    ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz+1), 0.125);

    ix = 0, iy=ny-1, iz=0;
    i = IC(ix, iy, iz);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
            
    ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
            
    ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
    ce.add(i, IF(2*ix+1, 2*iy  , 2*iz+1), 0.25);
    ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
            
    ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz+1), 0.125);

    ix = nx-1, iy=ny-1, iz=0;
    i = IC(ix, iy, iz);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
            
    ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz+1), 0.5);
            
    ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
    ce.add(i, IF(2*ix-1, 2*iy  , 2*iz+1), 0.25);
    ce.add(i, IF(2*ix  , 2*iy-1, 2*iz+1), 0.25);
            
    ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz+1), 0.125);

    ix = 0, iy=0, iz=nz-1;
    i = IC(ix, iy, iz);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
            
    ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
            
    ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz  ), 0.25);
    ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
    ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
            
    ce.add(i, IF(2*ix+1, 2*iy+1, 2*iz-1), 0.125);

    ix = nx-1, iy=0, iz=nz-1;
    i = IC(ix, iy, iz);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
            
    ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy+1, 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
            
    ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz  ), 0.25);
    ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
    ce.add(i, IF(2*ix  , 2*iy+1, 2*iz-1), 0.25);
            
    ce.add(i, IF(2*ix-1, 2*iy+1, 2*iz-1), 0.125);

    ix = 0, iy=ny-1, iz=nz-1;
    i = IC(ix, iy, iz);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
            
    ce.add(i, IF(2*ix+1, 2*iy  , 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
            
    ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz  ), 0.25);
    ce.add(i, IF(2*ix+1, 2*iy  , 2*iz-1), 0.25);
    ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
            
    ce.add(i, IF(2*ix+1, 2*iy-1, 2*iz-1), 0.125);

    ix = nx-1, iy=ny-1, iz=nz-1;
    i = IC(ix, iy, iz);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz  ), 1);
            
    ce.add(i, IF(2*ix-1, 2*iy  , 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy-1, 2*iz  ), 0.5);
    ce.add(i, IF(2*ix  , 2*iy  , 2*iz-1), 0.5);
            
    ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz  ), 0.25);
    ce.add(i, IF(2*ix-1, 2*iy  , 2*iz-1), 0.25);
    ce.add(i, IF(2*ix  , 2*iy-1, 2*iz-1), 0.25);
            
    ce.add(i, IF(2*ix-1, 2*iy-1, 2*iz-1), 0.125);
    
}


/*-------------------------------------------------*/
std::shared_ptr<FemAndMatrixAndSmootherInterface> Q12d::newStencil(std::shared_ptr<GridInterface const> grid) const
{
    std::shared_ptr<UniformGrid const> ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    assert(ug);
    const armaicvec& n = ug->n();
    const armavec& dx = ug->dx();
    assert(dx[0]==dx[1]);
    if(_stenciltype=="Full")
    {
        arma::vec::fixed<9> coef;
        coef[0] = -1.0/3.0;
        coef[1] = -1.0/3.0;
        coef[2] = -1.0/3.0;
        coef[3] = -1.0/3.0;
        coef[4] =  8.0/3.0;
        coef[5] = -1.0/3.0;
        coef[6] = -1.0/3.0;
        coef[7] = -1.0/3.0;
        coef[8] = -1.0/3.0;
        // return std::shared_ptr<MatrixInterface>(new Matrix<Stencil2d9,Vector<GridVector>>(n, coef));
        return std::shared_ptr<FemAndMatrixAndSmootherInterface>(new FemAndMatrixAndSmoother<Stencil2d9,Vector<GridVector>>(n, coef, _smoother));
    }
    else if(_stenciltype=="Trapez")
    {
        arma::vec::fixed<5> coef;
        coef[0] = -1.0;
        coef[1] = -1.0;
        coef[2] =  4.0;
        coef[3] = -1.0;
        coef[4] = -1.0;
        // return std::shared_ptr<MatrixInterface>(new Matrix<Stencil2d5,Vector<GridVector>>(n, coef));
        return std::shared_ptr<FemAndMatrixAndSmootherInterface>(new FemAndMatrixAndSmoother<Stencil2d5,Vector<GridVector>>(n, coef, _smoother));
    }
    else
    {
        std::cerr << "unknown matrix '" << _stenciltype<<"'\n";
        std::exit(1);
    }
}

/*-------------------------------------------------*/
std::shared_ptr<FemAndMatrixAndSmootherInterface> Q13d::newStencil(std::shared_ptr<GridInterface const> grid) const
{
    std::shared_ptr<UniformGrid const> ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
    // const UniformGrid* ug = dynamic_cast<const UniformGrid*>(grid);
    assert(ug);
    const armaicvec& n = ug->n();
    const armavec& dx = ug->dx();
    assert(dx[0]==dx[1]);
    assert(dx[0]==dx[2]);
    if(_stenciltype=="Full")
    {
        assert(dx[0]==dx[1]);
        double e = dx[0];
        double d0 = 8.0/3.0 * e;
        double d1 = -0.0/3.0 * e;
        double d2 = -1.0/6.0 * e;
        double d3 = -1.0/12.0 * e;
        arma::vec::fixed<27> coef;
        coef[13] = d0;
        coef[ 4] = coef[10] = coef[12] = coef[14] = coef[16] = coef[22] = d1;
        coef[ 1] = coef[ 3] = coef[ 5] = coef[ 7] = coef[ 9] = coef[11] = d2;
        coef[15] = coef[17] = coef[19] = coef[21] = coef[23] = coef[25] = d2;
        coef[ 0] = coef[ 2] = coef[ 6] = coef[ 8] = coef[18] = coef[20] = coef[24] = coef[26] = d3;
        // return std::shared_ptr<MatrixInterface>(new Matrix<Stencil3d27,Vector<GridVector>>(n, coef));
        return std::shared_ptr<FemAndMatrixAndSmootherInterface>(new FemAndMatrixAndSmoother<Stencil3d27,Vector<GridVector>>(n, coef, _smoother));
    }
    else if(_stenciltype=="Trapez")
    {
        assert(dx[0]==dx[1]);
        assert(dx[0]==dx[2]);
        arma::vec::fixed<7> coef;
        coef[0] = -1.0;
        coef[1] = -1.0;
        coef[2] = -1.0;
        coef[3] = 6.0;
        coef[4] = -1.0;
        coef[5] = -1.0;
        coef[6] = -1.0;
        coef * dx[0];
        // return std::shared_ptr<MatrixInterface>(new Matrix<Stencil3d7,Vector<GridVector>>(n, coef));
        return std::shared_ptr<FemAndMatrixAndSmootherInterface>(new FemAndMatrixAndSmoother<Stencil3d7,Vector<GridVector>>(n, coef, _smoother));
    }
    else
    {
        std::cerr << "unknown matrix '" << _stenciltype<<"'\n";
        std::exit(1);
    }
}

/*-------------------------------------------------*/
void Q12d::boundary_zero(GridVector& u) const
{
    for(int ix=0;ix<_nx;ix++)
    {
        u.at(ix,0)    = 0;
        u.at(ix,_ny-1) = 0;
    }
    for(int iy=0;iy<_ny;iy++)
    {
        u.at(0,iy)    = 0;
        u.at(_nx-1,iy) = 0;
    }
}


/*-------------------------------------------------*/
void Q13d::boundary_zero(GridVector& u) const
{
    for(int ix=0;ix<_nx;ix++)
    {
        for(int iy=0;iy<_ny;iy++)
        {
            u.at(ix,iy,0)     = 0;
            u.at(ix,iy,_nz-1) = 0;
        }
    }
    for(int ix=0;ix<_nx;ix++)
    {
        for(int iz=0;iz<_nz;iz++)
        {
            u.at(ix,0,   iz)   = 0;
            u.at(ix,_ny-1,iz)  = 0;
        }
    }
    for(int iy=0;iy<_ny;iy++)
    {
        for(int iz=0;iz<_nz;iz++)
        {
            u.at(0,   iy,iz)  = 0;
            u.at(_nx-1,iy,iz) = 0;
        }
    }
}

/*-------------------------------------------------*/
void Q12d::boundary_linear(GridVector& u) const
{
    for(int ix=0;ix<_nx;ix++)
    {
        u.at(ix,0)    = lin2d(_ug->x(ix,0), _ug->y(ix,0));
        u.at(ix,_ny-1) = lin2d(_ug->x(ix,_ny-1), _ug->y(ix,_ny-1));
    }
    for(int iy=0;iy<_ny;iy++)
    {
        u.at(0,iy)    = lin2d(_ug->x(0,iy), _ug->y(0,iy));
        u.at(_nx-1,iy) = lin2d(_ug->x(_nx-1,iy), _ug->y(_nx-1,iy));
    }
}


/*-------------------------------------------------*/
void Q13d::boundary_linear(GridVector& u) const
{
    for(int ix=0;ix<_nx;ix++)
    {
        for(int iy=0;iy<_ny;iy++)
        {
            u.at(ix,iy,0)     = lin3d(_ug->x(ix,iy,0), _ug->y(ix,iy,0), _ug->z(ix,iy,0));
            u.at(ix,iy,_nz-1) = lin3d(_ug->x(ix,iy,_nz-1), _ug->y(ix,iy,_nz-1), _ug->z(ix,iy,_nz-1));
        }
    }
    for(int ix=0;ix<_nx;ix++)
    {
        for(int iz=0;iz<_nz;iz++)
        {
            u.at(ix,0,   iz)   = lin3d(_ug->x(ix,0,   iz), _ug->y(ix,0,   iz), _ug->z(ix,0,   iz));
            u.at(ix,_ny-1,iz)  = lin3d(_ug->x(ix,_ny-1,iz), _ug->y(ix,_ny-1,iz), _ug->z(ix,_ny-1,iz));
        }
    }
    for(int iy=0;iy<_ny;iy++)
    {
        for(int iz=0;iz<_nz;iz++)
        {
            u.at(0,   iy,iz)  = lin3d(_ug->x(0,   iy,iz), _ug->y(0,   iy,iz), _ug->z(0,   iy,iz));
            u.at(_nx-1,iy,iz) = lin3d(_ug->x(_nx-1,iy,iz), _ug->y(_nx-1,iy,iz), _ug->z(_nx-1,iy,iz));
        }
    }
}
