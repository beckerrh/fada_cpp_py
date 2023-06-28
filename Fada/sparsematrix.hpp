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
#include  "matrixinterface.hpp"

/*-------------------------------------------------*/
class Construct_Elements
{
protected:
  int _count;
  arma::umat& _locations;
  armavec& _values;
  void start(int size)
  {
    _locations.resize(2, size);
    _values.resize(size);
    _count=0;
  }
public:
  // Construct_Elements() {}
  Construct_Elements(arma::umat& locations, armavec& values, int size) : _locations(locations), _values(values) {start(size);}
  const arma::umat& locations() const {return _locations;}
  const armavec& values() const {return _values;}
  void add(int i, int j, double value)
  {
    _locations.at(0, _count) = i;
    _locations.at(1, _count) = j;
    _values[_count] = value;
    _count++;
  }
};

/*-------------------------------------------------*/
class SparseMatrix
{
protected:
#ifdef _LONG_LONG
  arma::uvec _rows, _cols, _diags;
#else
  armaicvec _rows, _cols, _diags;
#endif
  armavec _values;
  void set_elements(const arma::umat& locations, const armavec& values);

public:
  SparseMatrix() {}
  SparseMatrix(const SparseMatrix& matrix) {}
  SparseMatrix(const arma::umat& locations, const armavec& values) {set_elements(locations, values);}

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

  void set_grid(const armaicvec& n, const armavec& coef){assert(0);exit(0);}
  void get_locations_values(arma::umat& locations, armavec& values) const{assert(0);exit(0);}
  void jacobi       (armavec& out, const armavec& in) const;
  void gauss_seidel1(armavec& out, const armavec& in) const;
  void gauss_seidel2(armavec& out, const armavec& in) const;
};


#endif /* sparsematrix_hpp */
