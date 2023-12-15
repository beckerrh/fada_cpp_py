//
//  vectorinterface.hpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef vectorinterface_h
#define vectorinterface_h

#include  "typedefs.hpp"

class BoundaryConditionsBool;

/*-------------------------------------------------*/
class VectorInterface
{
public:
  virtual ~VectorInterface() {}
  VectorInterface() {}
  VectorInterface(const VectorInterface& vector) {}

  virtual void set_size(const armaicvec& n)=0;
  virtual int get_size()const=0;
  virtual void after_prolongate()=0;
  virtual void after_restrict()=0;
  virtual void after_smooth()=0;
  virtual void after_residual()=0;
  virtual void fill(double d=0)=0;
  virtual double dot(std::shared_ptr<VectorInterface const> v)const=0;
  virtual double norm(double p=2)const=0;
  virtual void equal(std::shared_ptr<VectorInterface const> v)=0;
  virtual void add(double d, std::shared_ptr<VectorInterface const> v)=0;
  virtual void scale(double d)=0;
  virtual void save(std::ostream& os, arma::file_type ftype=arma::arma_binary) const=0;
  virtual void savehdf5(const std::string& filename) const=0;
};
/*-------------------------------------------------*/
template<typename VECTOR>
class Vector : public VECTOR, public VectorInterface
{
protected:
  VECTOR const& getVector(std::shared_ptr<VectorInterface const> v) const
  {
    auto p = std::dynamic_pointer_cast<Vector<VECTOR> const>(v);
    // const Vector<VECTOR>* p = dynamic_cast< const Vector<VECTOR>* >(v);
    return p->get();
  }
public:
  Vector<VECTOR>() : VECTOR(), VectorInterface() {}
  Vector<VECTOR>(const armaicvec& n, std::shared_ptr<BoundaryConditionsBool const> bc=nullptr, bool mean=false) : VECTOR(n,bc, mean), VectorInterface() {}

  VECTOR& get() { return static_cast<VECTOR&>(*this); }
  VECTOR const& get() const { return static_cast<VECTOR const&>(*this); }
  void set_size(const armaicvec& n) {get().set_size(n);}
  int get_size() const{return get().get_size();}
  void after_prolongate() {get().after_prolongate();}
  void after_restrict() {get().after_restrict();}
  void after_smooth() {get().after_smooth();}
  void after_residual() {get().after_residual();}
  void fill(double d=0) {get().fill(d);}
  double norm(double p=2)const {return get().norm(p);}
  void scale(double d) {get().scale(d);}
  double dot(std::shared_ptr<VectorInterface const> v)const {return get().dot(getVector(v));}
  void equal(std::shared_ptr<VectorInterface const> v) {get().equal(getVector(v));}
  void add(double d, std::shared_ptr<VectorInterface const> v) {get().add(d, getVector(v));}
  void save(std::ostream& os, arma::file_type ftype=arma::arma_binary) const{get().save(os,ftype);}
  void savehdf5(const std::string& filename) const{get().savehdf5(filename);}
};


#endif
