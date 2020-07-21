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
class Construct_Elements
{
protected:
  int _count;
  arma::umat _locations;
  armavec _values;
public:
  Construct_Elements() {}
  Construct_Elements(int size) {start(size);}
  const arma::umat& locations() const {return _locations;}
  const armavec& values() const {return _values;}
  void start(int size)
  {
    _locations.resize(2, size);
    _values.resize(size);
    _count=0;
  }
  void add(int i, int j, double value)
  {
    _locations(0, _count) = i;
    _locations(1, _count) = j;
    _values[_count] = value;
    _count++;
  }
};

/*-------------------------------------------------*/
class SparseMatrix
{
protected:
#ifdef _LONG_LONG
  arma::uvec _rows, _cols;
#else
  armaicvec _rows, _cols;
#endif
  armavec _values;
  
public:
  SparseMatrix() {}
  SparseMatrix(const SparseMatrix& matrix) {}
  SparseMatrix(const arma::umat& locations, const armavec& values);
  void set_elements(const arma::umat& locations, const armavec& values);

  #ifdef _LONG_LONG
    const arma::uvec& rows() const {return _rows;}
    const arma::uvec& cols() const {return _cols;}
  #else
      const armaicvec& rows() const {return _rows;}
      const armaicvec& cols() const {return _cols;}
  #endif

  const armavec& values() const {return _values;}
  arma::uword nrows() const {return _rows.n_elem-1;}
  void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const;
  void dot(armavec& x, const armavec& b, double d=1) const;
};

#endif /* sparsematrix_hpp */
