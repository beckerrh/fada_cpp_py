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

class GridInterface;
class Vector;
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
  virtual armavec& data() =0;
  virtual const armavec& data() const=0;

};


#endif /* vectorinterface_h */
