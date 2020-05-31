
//
//  matrixinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef matrixinterface_h
#define matrixinterface_h

#include  "typedefs.hpp"

class GridInterface;
class VectorInterface;
class SparseMatrix;
/*-------------------------------------------------*/
class MatrixInterface
{
public:
  virtual ~MatrixInterface() {}
  MatrixInterface() {}
  MatrixInterface(const MatrixInterface& updater) {}
  
  virtual void set_grid(const armaicvec& n, const armavec& dx)=0;
  virtual void dot(VectorInterface& out, const VectorInterface& in, double d) const=0;
  virtual void get_sparse_matrix(SparseMatrix& sp) const=0;
  
  virtual void jacobi       (VectorInterface& out, const VectorInterface& in) const=0;
  virtual void gauss_seidel1(VectorInterface& out, const VectorInterface& in) const=0;
  virtual void gauss_seidel2(VectorInterface& out, const VectorInterface& in) const=0;
};

/*-------------------------------------------------*/
template<typename MATRIX, class VECTOR>
class Matrix : public MATRIX, public MatrixInterface
{
protected:
  MATRIX& get() { return static_cast<MATRIX&>(*this); }
  MATRIX const& get() const { return static_cast<MATRIX const&>(*this); }
  const VECTOR& getVector(const VectorInterface& u) const {return static_cast<const VECTOR&>(u);}
  VECTOR& getVector(VectorInterface& u) const{return static_cast<VECTOR&>(u);}
public:
  Matrix<MATRIX, VECTOR>() : MATRIX(), MatrixInterface() {}
  Matrix<MATRIX, VECTOR>(const armaicvec& n, const armavec& dx) : MATRIX(n, dx), MatrixInterface() {}
  
  void set_grid(const armaicvec& n, const armavec& dx){get().set_grid(n, dx);}
  void dot(VectorInterface& out, const VectorInterface& in, double d) const {get().dot(getVector(out),getVector(in), d);}
  void get_sparse_matrix(SparseMatrix& sp) const{get().get_sparse_matrix(sp);}
  void jacobi(VectorInterface& out, const VectorInterface& in) const{get().jacobi(getVector(out),getVector(in));}
  void gauss_seidel1(VectorInterface& out, const VectorInterface& in) const{get().gauss_seidel1(getVector(out),getVector(in));}
  void gauss_seidel2(VectorInterface& out, const VectorInterface& in) const{get().gauss_seidel2(getVector(out),getVector(in));}
};


#endif /* matrixinterface_h */
