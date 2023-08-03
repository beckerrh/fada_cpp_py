//
//  transferq1.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef transferbymatrix_h
#define transferbymatrix_h

#include  "sparsematrix.hpp"
#include  "gridvector.hpp"


/*-------------------------------------------------*/
class TransferByMatrix
{
protected:
  SparseMatrix _matrix;
  
public:
  ~TransferByMatrix()
  {
  }

  TransferByMatrix(const TransferByMatrix& transfer) : _matrix(transfer._matrix)
  {
  }

  TransferByMatrix(const armaicvec& n, const armavec& dx, std::shared_ptr <BoundaryConditions const> boundaryconditions)
  {
      _not_written_();
  }
  TransferByMatrix(const arma::umat& locations, const armavec& values) : _matrix(locations, values, false) {}

  const SparseMatrix& getMatrix() const {return _matrix;}
  void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const {_matrix.save(out, datatype);}
  void restrict (GridVector & out, const GridVector& in) const{out.fill(0); _matrix.dot(out,in);}
  // void restrict (GridVector & out, const GridVector& in) const{out.fill(0); _matrix.dot(out,in); out.boundary_zero();}
  void prolongate(GridVector& out, const GridVector& in) const{out.fill(0); _matrix.Tdot(out,in);}
};
#endif
