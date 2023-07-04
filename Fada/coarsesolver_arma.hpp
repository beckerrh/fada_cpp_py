//
//  smootherinterface.hpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef coarsesolver_arma_h
#define coarsesolver_arma_h

#include  "typedefs.hpp"
#include  "sparsematrix_arma.hpp"
#include  "matrixinterface.hpp"
#include  "vectorinterface.hpp"

/*-------------------------------------------------*/
class CoarseSolver_arma
{
protected:
  std::string _type;
  std::shared_ptr<arma::sp_mat const> _matrix;

public:
  CoarseSolver_arma() {}
  CoarseSolver_arma(const CoarseSolver_arma& matrix) : _matrix(matrix._matrix), _type(matrix._type) {}
  CoarseSolver_arma(std::shared_ptr<MatrixInterface const> matrix, std::string type) : _type(type)
  {
    auto armamat = std::dynamic_pointer_cast<Matrix<SparseMatrix_arma,Vector<armavec>> const>(matrix);
    assert(armamat);
    const arma::sp_mat&  sp_arma = static_cast<const arma::sp_mat&>(armamat->get());
    _matrix = std::make_shared<arma::sp_mat const>(sp_arma);
    assert(_matrix);
  }

  void solve(armavec& out, const armavec& in) const
  {
    arma::spsolve(out, *_matrix, in);
  }
};
#endif
