//
//  sparsematrix.cpp
//  Fada
//
//  Created by Roland Becker on 20/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include "sparsematrix.hpp"

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
  for(int i=0;i<nrows();i++)
  {
    for(int pos=_rows[i];pos<_rows[i+1];pos++)
    {
      x[i] += d*_values[pos] * b[_cols[pos]];
    }
  }
}
/*-------------------------------------------------*/
void SparseMatrix::Tdot(armavec& x, const armavec& b, double d) const
{
  for(int i=0;i<nrows();i++)
  {
    for(int pos=_rows[i];pos<_rows[i+1];pos++)
    {
      x[_cols[pos]] += d*_values[pos] * b[i];
    }
  }
}

/*-------------------------------------------------*/
void SparseMatrix::set_elements(const arma::umat& locations, const armavec& values)
{
  assert(locations.n_rows==2);
  assert(locations.n_cols==values.n_elem);
  arma::uword n = locations.n_cols;
  _cols.set_size(n);
  _values.set_size(n);
  // std::cerr << "locations\n" << locations.row(0) << "\n" << locations.row(1) << "\n";
  // sort the indices !!
  arma::uvec ind = arma::sort_index(arma::max(locations.row(1))*locations.row(0)+locations.row(1));
  arma::uword nrows = locations(0,ind[n-1])+1;
  _rows.set_size(nrows+1);
  int count = 0, prec=-1;
  for(int i=0;i<n;i++)
  {
    arma::uword index = ind[i];
    arma::uword row = locations.at(0,index);
    if(row!=prec)
    {
      _rows[count] = i;
      prec = row;
      count++;
    }
    _values[i] = values[index];
    _cols[i] = locations.at(1,index);
  }
  _rows[nrows] = n;
  // store diagonal positions
  _diags.set_size(nrows);
  for(int i=0;i<nrows;i++)
  {
    for(int pos=_rows[i];pos<_rows[i+1];pos++)
    {
      if(_cols[pos]==i)
      {
        _diags[i] = pos;
        break;
      }
    }
  }
  // for(int i=0;i<nrows;i++)
  // {
  //   std::cerr << "\nrow = " << i << " " << _cols[_diags[i]] << "\n";
  //   for(int pos=_rows[i];pos<_rows[i+1];pos++)
  //   {
  //     std::cerr << " " << _cols[pos];
  //   }
  // }
 // save(std::cerr);
//  assert(0);
}
/*-------------------------------------------------*/
void SparseMatrix::jacobi(armavec& out, const armavec& in) const
{
  for(int i=0;i<nrows();i++)
  {
    out[i] = in[i]/_values[_diags[i]];
  }
}
/*-------------------------------------------------*/
void SparseMatrix::gauss_seidel1(armavec& out, const armavec& in) const
{
  for(int i=0;i<nrows();i++)
  {
    double val = in[i];
    for(int pos=_rows[i];pos<_rows[i+1];pos++)
    {
      int j = _cols[pos];
      if(j>=i) break;
      val -= _values[pos] * out[j];
    }
    out[i] = val/_values[_diags[i]];
  }
}
/*-------------------------------------------------*/
void SparseMatrix::gauss_seidel2(armavec& out, const armavec& in) const
{
  for(int i=nrows()-1;i>=0;i--)
  {
    double val = in[i];
    for(int pos=_rows[i+1]-1;pos>=_rows[i];pos--)
    {
      int j = _cols[pos];
      if(j<=i) break;
      val -= _values[pos] * out[j];
    }
    out[i] = val/_values[_diags[i]];
  }
}
