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
  arma::uword n;
  #ifdef _LONG_LONG
    const long long* sb;
    const long long* cb;
  #else
    const int* sb;
    const int* cb;
  #endif
  const double* mb;

  std::shared_ptr<MatrixInterface const> _matrix;

public:
  ~UmfMatrix();
  UmfMatrix();
  UmfMatrix( const UmfMatrix& umfmatrix);

  void save(std::ostream& os, arma::file_type datatype = arma::arma_binary) const;
  void init(std::shared_ptr<MatrixInterface const> matrix);
  void computeLu();
  void solve(armavec& x, const armavec& b) const;
};

#endif /* umfmatrix_hpp */
