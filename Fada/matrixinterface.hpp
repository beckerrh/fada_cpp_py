
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
/*-------------------------------------------------*/
class MatrixInterface
{
public:
  virtual ~MatrixInterface() {}
  MatrixInterface() {}
  MatrixInterface(const MatrixInterface& matrix) {}

  virtual void dot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const=0;
  // virtual void Tdot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const=0;

  virtual void jacobi       (std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{assert(0);}
  virtual void gauss_seidel1(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{assert(0);}
  virtual void gauss_seidel2(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{assert(0);}
  virtual void set_elements(const arma::umat& locations, const armavec& values) {assert(0);exit(1);}
  virtual void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const{assert(0);exit(1);}
};

/*-------------------------------------------------*/
template<typename MATRIX, class VECTOR>
class Matrix : public  virtual MATRIX, public  virtual MatrixInterface
{
protected:
  const VECTOR& getVector(std::shared_ptr<VectorInterface const> u) const {return static_cast<const VECTOR&>(*u);}
  VECTOR& getVector(std::shared_ptr<VectorInterface> u) const{return static_cast<VECTOR&>(*u);}

public:
  Matrix<MATRIX, VECTOR>() : MATRIX(), MatrixInterface() {}
  Matrix<MATRIX, VECTOR>(const arma::umat& locations, const armavec& values): MATRIX(locations,values), MatrixInterface(){}

  MATRIX& get() { return static_cast<MATRIX&>(*this); }
  MATRIX const& get() const { return static_cast<MATRIX const&>(*this); }
  void dot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const {get().dot(getVector(out),getVector(in), d);}
  // void Tdot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const {get().Tdot(getVector(out),getVector(in), d);}
  void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const{get().save(out, datatype);}
};


#endif /* matrixinterface_h */
