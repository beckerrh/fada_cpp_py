//
//  umfmatrix.hpp
//  Fada
//
//  Created by Roland Becker on 19/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef umfmatrix_hpp
#define umfmatrix_hpp

#include  <iostream>
#include  "typedefs.hpp"
#include  "sparsematrix.hpp"

/*-------------------------------------------------*/
class UmfMatrix
{
protected:
  double* Control;
  double* Info;
  void* Symbolic, * Numeric;
  SparseMatrix _sp;
  arma::uword n;
//  const int* sb;
//  const int* cb;
  const long long* sb;
  const long long* cb;
  const double* mb;

  void init();

public:
  ~UmfMatrix();
  UmfMatrix();
  UmfMatrix( const UmfMatrix& umfmatrix);
  
  const SparseMatrix& getSparseMatrix() const {return _sp;}
  SparseMatrix& getSparseMatrix() {return _sp;}
  void save(std::ostream& os, arma::file_type datatype = arma::arma_binary) const;
  void computeLu();
  void solve(armavec& x, const armavec& b) const;
};

#endif /* umfmatrix_hpp */
