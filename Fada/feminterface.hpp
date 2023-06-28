
//
//  feminterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef feminterface_h
#define feminterface_h

#include  "typedefs.hpp"
#include  "matrixinterface.hpp"

class GridInterface;
class VectorInterface;
/*-------------------------------------------------*/
class FemInterface
{
public:
  virtual ~FemInterface() {}
  FemInterface() {}
  FemInterface(const FemInterface& fem) {}

  virtual void get_locations_values(arma::umat& locations, armavec& values) const=0;
  virtual void set_grid(const armaicvec& n, const armavec& dx){assert(0);exit(1);}
};

/*-------------------------------------------------*/
template<typename FEM>
class Fem : public virtual FEM, public virtual FemInterface
{
protected:
  FEM& get() { return static_cast<FEM&>(*this); }
  FEM const& get() const { return static_cast<FEM const&>(*this); }

public:
  Fem<FEM>() : FEM(), FemInterface() {}
  Fem<FEM>(const Fem<FEM>& fem) : FEM(), FemInterface(fem) {}
  Fem<FEM>(const armaicvec& n, const armavec& coef) : FEM(n, coef), FemInterface() {}

  void get_locations_values(arma::umat& locations, armavec& values) const{get().get_locations_values(locations, values);}
  void set_grid(const armaicvec& n, const armavec& dx){get().set_grid(n, dx);}
};

/*-------------------------------------------------*/
class FemAndMatrixInterface : public virtual FemInterface, public virtual MatrixInterface
{
public:
  virtual ~FemAndMatrixInterface() {}
  FemAndMatrixInterface() : FemInterface(), MatrixInterface() {}
  FemAndMatrixInterface(const FemAndMatrixInterface& fem) : FemInterface(fem), MatrixInterface(fem) {}
};

/*-------------------------------------------------*/
template<typename FEM, class VECTOR>
class FemAndMatrix : public virtual FemAndMatrixInterface, public FEM
{
protected:
  FEM& get() { return static_cast<FEM&>(*this); }
  FEM const& get() const { return static_cast<FEM const&>(*this); }
  const VECTOR& getVector(std::shared_ptr<VectorInterface const> u) const {return static_cast<const VECTOR&>(*u);}
  VECTOR& getVector(std::shared_ptr<VectorInterface> u) const{return static_cast<VECTOR&>(*u);}

public:
  FemAndMatrix<FEM, VECTOR>() : FemAndMatrixInterface(), FEM() {}
  FemAndMatrix<FEM, VECTOR>(const armaicvec& n, const armavec& dx) : FemAndMatrixInterface(), FEM(n,dx) {}

  void get_locations_values(arma::umat& locations, armavec& values) const{get().get_locations_values(locations, values);}
  void set_grid(const armaicvec& n, const armavec& dx){get().set_grid(n, dx);}
  void dot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const {get().dot(getVector(out),getVector(in), d);}
  void jacobi(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{get().jacobi(getVector(out),getVector(in));}
  void gauss_seidel1(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{get().gauss_seidel1(getVector(out),getVector(in));}
  void gauss_seidel2(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{get().gauss_seidel2(getVector(out),getVector(in));}
  void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const{get().save(out, datatype);}
};


#endif /* feminterface_h */
