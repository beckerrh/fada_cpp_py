//
//  sparsematrix.cpp
//  Fada
//
//  Created by Roland Becker on 20/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include "sparsematrix_arma.hpp"

/*-------------------------------------------------*/
void SparseMatrix_arma::jacobi(armavec& out, const armavec& in) const
{
  // std::cerr << "_diag=\n" << _diag.t() << "\n";
  out = in/_diag;
}
/*-------------------------------------------------*/
void SparseMatrix_arma::gauss_seidel1(armavec& out, const armavec& in) const
{
  assert(0); std::cerr << "only Jacobi for arma::sp_mat\n"; exit(1);
  const arma::sp_mat& tarma = static_cast<const arma::sp_mat&>(*this);
  for(arma::uword c=0; c < tarma.n_cols; ++c)
  {
  const arma::uword index_start = tarma.col_ptrs[c    ];
  const arma::uword index_end   = tarma.col_ptrs[c + 1];

  for(arma::uword i=index_start; i < index_end; ++i)
    {
    const arma::uword r = tarma.row_indices[i];

    // eT& val = tarma.values[i];
    //
    // const eT result = val * B.at(r,c);
    //
    // val = result;

    }
  }

  // for(int i=0;i<tarma.n_rows;i++)
  // {
  //   double val = in[i];
  //   for(int pos=_rows[i];pos<_rows[i+1];pos++)
  //   {
  //     int j = _cols[pos];
  //     if(j>=i) break;
  //     val -= _values[pos] * out[j];
  //   }
  //   out[i] = val/_values[_diags[i]];
  // }
}
/*-------------------------------------------------*/
void SparseMatrix_arma::gauss_seidel2(armavec& out, const armavec& in) const
{
  assert(0); std::cerr << "only Jacobi for arma::sp_mat\n"; exit(1);
  // for(int i=nrows()-1;i>=0;i--)
  // {
  //   double val = in[i];
  //   for(int pos=_rows[i+1]-1;pos>=_rows[i];pos--)
  //   {
  //     int j = _cols[pos];
  //     if(j<=i) break;
  //     val -= _values[pos] * out[j];
  //   }
  //   out[i] = val/_values[_diags[i]];
  // }
}
