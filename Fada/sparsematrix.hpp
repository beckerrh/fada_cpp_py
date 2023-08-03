//
//  sparsematrix.hpp
//  Fada
//
//  Created by Roland Becker on 20/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef sparsematrix_hpp
#define sparsematrix_hpp

#include  "typedefs.hpp"

class SystemMatrix;

/*-------------------------------------------------*/
class SparseMatrix
{
protected:
    arma::uword  _n_cols;
#ifdef _LONG_LONG
    arma::uvec _rows, _cols, _diags;
#else
    armaicvec _rows, _cols, _diags;
#endif
    armavec _values;
    void compute_diag();

public:
    SparseMatrix(): _n_cols(0), _rows(1), _cols(), _diags() {}
    SparseMatrix(const SparseMatrix& matrix): _rows(matrix._rows), _cols(matrix._cols), _diags(matrix._diags), _n_cols(matrix._n_cols) {}
    SparseMatrix(const arma::umat& locations, const armavec& values, bool compute_diag_flag=true) {set_elements(locations, values,compute_diag_flag);}
    SparseMatrix(std::shared_ptr<SystemMatrix const> sm) {set_elements(sm);}

    void set_elements(const arma::umat& locations, const armavec& values, bool compute_diag=true);
    void set_elements(std::shared_ptr<SystemMatrix const> sm);

#ifdef _LONG_LONG
    const arma::uvec& rows() const {return _rows;}
    const arma::uvec& cols() const {return _cols;}
#else
    const armaicvec& rows() const {return _rows;}
    const armaicvec& cols() const {return _cols;}
#endif

    const armavec& values() const {return _values;}
    int nrows() const {return _rows.n_elem-1;}
    int ncols() const {return _n_cols;}
    int nelem() const {return _values.n_elem;}
    void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const;
    std::shared_ptr<armavec> getDiag() const;
    void dot(armavec& x, const armavec& b, double d=1) const;
    void Tdot(armavec& x, const armavec& b, double d=1) const;

    void jacobi       (armavec& out, const armavec& in) const;
    void gauss_seidel1(armavec& out, const armavec& in) const;
    void gauss_seidel2(armavec& out, const armavec& in) const;
    
    void addBBT(const SparseMatrix& B, double d=1.0 , std::shared_ptr<armavec const> D=nullptr);
    void add_diagonal(double d);
    std::shared_ptr<SparseMatrix> getT() const;
};

class ConstantDiagonalSparseMatrix: public SparseMatrix
{
public:
    ConstantDiagonalSparseMatrix(const arma::umat& locations, const armavec& values, bool compute_diag=true) : SparseMatrix(locations, values) {}
    ConstantDiagonalSparseMatrix(int n, double d) : SparseMatrix()
    {
        _n_cols=n;
        _values.resize(n);
        _values.fill(d);
#ifdef _LONG_LONG
        _rows = arma::linspace<uvec>(0,n,n+1); 
        _cols = arma::linspace<uvec>(0,n-1,n); 
#else
        _rows = arma::linspace<armaicvec>(0,n,n+1); 
        _cols = arma::linspace<armaicvec>(0,n-1,n); 
#endif
    }
};


#endif /* sparsematrix_hpp */
