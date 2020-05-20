//
//  sparsematrix.cpp
//  Fada
//
//  Created by Roland Becker on 20/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include "sparsematrix.hpp"

/*-------------------------------------------------*/
SparseMatrix::SparseMatrix(arma::umat& locations, armavec& values)
{
  set_elements(locations, values);
}
/*-------------------------------------------------*/
void SparseMatrix::save(std::ostream& out, arma::file_type datatype) const
{
  for(int i=0;i<nrows();i++)
  {
    for(int pos=_rows[i];pos<_rows[i+1];pos++)
    {
      out << "(" << i << " " << _cols[pos] << ") " << _values[pos] << "\n";
    }
  }
}
/*-------------------------------------------------*/
void SparseMatrix::dot(armavec& x, const armavec& b, double d) const
{
//  x.fill(0);
  for(int i=0;i<nrows();i++)
  {
    for(int pos=_rows[i];pos<_rows[i+1];pos++)
    {
      x[i] += d*_values[pos] * b[_cols[pos]];
    }
  }
}

/*-------------------------------------------------*/
void SparseMatrix::set_elements(arma::umat& locations, armavec& values)
{
  assert(locations.n_rows==2);
  assert(locations.n_cols==values.n_elem);
  arma::uword n = locations.n_cols;
  _cols.set_size(n);
  _values.set_size(n);
  arma::uvec ind = arma::sort_index(arma::max(locations.row(1))*locations.row(0)+locations.row(1));
  arma::uword nrows = locations(0,ind[n-1])+1;
  _rows.set_size(nrows+1);
  int count = 0, prec=-1;
  for(int i=0;i<n;i++)
  {
    arma::uword index = ind[i];
    arma::uword row = locations(0,index);
    if(row!=prec)
    {
      _rows[count] = i;
      prec = row;
      count++;
    }
    _values[i] = values[index];
    _cols[i] = locations(1,index);
  }
  _rows[nrows] = n;
//  save(std::cerr);
//  assert(0);
}
