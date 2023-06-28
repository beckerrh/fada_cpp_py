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

/*-------------------------------------------------*/
class VectorInterface
{
public:
  virtual ~VectorInterface() {}
  VectorInterface() {}
  VectorInterface(const VectorInterface& vector) {}

  virtual void set_size(const armaicvec& n)=0;
  virtual void fill_bdry(double d=0)=0;
  virtual void fill(double d=0)=0;
  virtual double dot(const VectorInterface& v)const=0;
  virtual double norm(double p=2)const=0;
  virtual void equal(const VectorInterface& v)=0;
  virtual void add(double d, const VectorInterface& v)=0;
  virtual void scale(double d)=0;
  // virtual armavec& data() =0;
  // virtual const armavec& data() const=0;
  virtual void save(std::ostream& os, arma::file_type ftype=arma::arma_ascii) const=0;
};
/*-------------------------------------------------*/
template<typename VECTOR>
class Vector : public VECTOR, public VectorInterface
{
protected:
  VECTOR const& getVector(const VectorInterface& v) const
  {
    const Vector<VECTOR>* p = dynamic_cast< const Vector<VECTOR>* >(&v);
    return p->get();
    // return std::dynamic_pointer_cast<VECTOR const&>(&v);
    // return static_cast<VECTOR const&>(v);
  }
public:
  VECTOR& get() { return static_cast<VECTOR&>(*this); }
  VECTOR const& get() const { return static_cast<VECTOR const&>(*this); }
  Vector<VECTOR>() : VECTOR(), VectorInterface() {}
  Vector<VECTOR>(const armaicvec& n) : VECTOR(n), VectorInterface() {}

  // armavec& data()  {return get().data();}
  // const armavec& data() const {return get().data();}
  void set_size(const armaicvec& n) {get().set_size(n);}
  void fill_bdry(double d=0) {get().fill_bdry(d);}
  void fill(double d=0) {get().fill(d);}
  double norm(double p=2)const {return get().norm(p);}
  void scale(double d) {get().scale(d);}
  double dot(const VectorInterface& v)const {return get().dot(getVector(v));}
  void equal(const VectorInterface& v) {get().equal(getVector(v));}
  void add(double d, const VectorInterface& v) {get().add(d, getVector(v));}
   void save(std::ostream& os, arma::file_type ftype=arma::arma_ascii) const{get().save(os,ftype);}
};


#endif /* vectorinterface_h */
