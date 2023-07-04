//
//  sparsematrix.hpp
//  Fada
//
//  Created by Roland Becker on 20/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef sparsematrix_arma_hpp
#define sparsematrix_arma_hpp

#include  "typedefs.hpp"

/*-------------------------------------------------*/
class SparseMatrix_arma : public arma::sp_mat
{
protected:
  armavec _diag;

public:
  SparseMatrix_arma() : arma::sp_mat() {}
  SparseMatrix_arma(const SparseMatrix_arma& matrix) : arma::sp_mat(matrix) {}
  SparseMatrix_arma(const arma::umat& locations, const armavec& values) : arma::sp_mat(locations, values)
  {
    const arma::sp_mat& tarma = static_cast<const arma::sp_mat&>(*this);
    _diag = tarma.diag();
  }

  void dot(armavec& out, const armavec& in, double d=1) const
  {
    const arma::sp_mat& tarma = static_cast<const arma::sp_mat&>(*this);
    // out = d*arma::dot(tarma, in);
    out += d* tarma*in;
    // std::cerr << "in=\n" << in.t() << "\n";
    // std::cerr << "tarma=\n" << tarma << "\n";
    // std::cerr << "out=\n" << out.t() << "\n";
  }
  void jacobi       (armavec& out, const armavec& in) const;
  void gauss_seidel1(armavec& out, const armavec& in) const;
  void gauss_seidel2(armavec& out, const armavec& in) const;
};


#endif /* sparsematrix_hpp */
