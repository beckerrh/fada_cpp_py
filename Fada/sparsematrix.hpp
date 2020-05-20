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

/*-------------------------------------------------*/
class SparseMatrix
{
protected:
//  armaicvec _rows, _cols;
  arma::uvec _rows, _cols;
  armavec _values;
  
public:
  SparseMatrix() {}
  SparseMatrix(const SparseMatrix& matrix) {}
  SparseMatrix(arma::umat& locations, armavec& values);
  void set_elements(arma::umat& locations, armavec& values);

//  const armaicvec& rows() const {return _rows;}
//  const armaicvec& cols() const {return _cols;}
  const arma::uvec& rows() const {return _rows;}
  const arma::uvec& cols() const {return _cols;}
  const armavec& values() const {return _values;}
  arma::uword nrows() const {return _rows.n_elem-1;}
  void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const;
  void dot(armavec& x, const armavec& b, double d=1) const;
};

#endif /* sparsematrix_hpp */
