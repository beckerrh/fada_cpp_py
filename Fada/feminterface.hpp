
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
#include  "smootherinterface.hpp"

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
class FemAndMatrixAndSmootherInterface : public virtual FemInterface, public virtual MatrixInterface, public virtual SmootherInterface
{
public:
  virtual ~FemAndMatrixAndSmootherInterface() {}
  FemAndMatrixAndSmootherInterface() : FemInterface(), MatrixInterface(), SmootherInterface() {}
  FemAndMatrixAndSmootherInterface(const FemAndMatrixAndSmootherInterface& fem) : FemInterface(fem), MatrixInterface(fem), SmootherInterface(fem) {}
};

/*-------------------------------------------------*/
template<typename FEM, class VECTOR>
class FemAndMatrixAndSmoother : public virtual FemAndMatrixAndSmootherInterface, public FEM
{
protected:
  FEM& get() { return static_cast<FEM&>(*this); }
  FEM const& get() const { return static_cast<FEM const&>(*this); }
  const VECTOR& getVector(std::shared_ptr<VectorInterface const> u) const {return static_cast<const VECTOR&>(*u);}
  VECTOR& getVector(std::shared_ptr<VectorInterface> u) const{return static_cast<VECTOR&>(*u);}

public:
  FemAndMatrixAndSmoother<FEM, VECTOR>() : FemAndMatrixAndSmootherInterface(), FEM() {}
  FemAndMatrixAndSmoother<FEM, VECTOR>(const armaicvec& n, const armavec& coef, std::string smoother) : FemAndMatrixAndSmootherInterface(), FEM(n, coef, smoother) {}

  void get_locations_values(arma::umat& locations, armavec& values) const{get().get_locations_values(locations, values);}
  void set_grid(const armaicvec& n, const armavec& dx){get().set_grid(n, dx);}
  void dot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const {get().dot(getVector(out),getVector(in), d);}
  void Tdot(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in, double d=1) const {get().Tdot(getVector(out),getVector(in), d);}
  void save(std::ostream& out, arma::file_type datatype = arma::arma_ascii) const{get().save(out, datatype);}
  void presmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{get().presmooth(getVector(out),getVector(in));}
  void postsmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const{get().postsmooth(getVector(out),getVector(in));}
  void update(std::shared_ptr<MatrixInterface const> matrix) {get().update(matrix);}
  int nrows() const{_not_written_(); return 0;}
  int ncols() const{_not_written_(); return 0;}
  int nelem() const{_not_written_(); return 0;}
};


#endif /* feminterface_h */
