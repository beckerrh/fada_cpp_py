
//
//  transferinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef transferinterface_h
#define transferinterface_h

#include  <memory>
#include  "typedefs.hpp"

class GridInterface;
class VectorInterface;
/*-------------------------------------------------*/
class TransferInterface
{
public:
  virtual ~TransferInterface() {}
  TransferInterface() {}
  TransferInterface(const TransferInterface& transfer) {}

  virtual void set_grid(const armaicvec& n, const armavec& dx)=0;
  virtual void restrict(VectorInterface& out, const VectorInterface& in) const=0;
  virtual void prolongate(VectorInterface& out, const VectorInterface& in) const=0;
};

/*-------------------------------------------------*/
template<typename TRANSFER, class VECTOR>
class Transfer : public TRANSFER, public TransferInterface
{
protected:
TRANSFER& get() { return static_cast<TRANSFER&>(*this); }
TRANSFER const& get() const { return static_cast<TRANSFER const&>(*this); }
const VECTOR& getVector(const VectorInterface& u) const
{
  const VECTOR* uV = dynamic_cast<const VECTOR*>(&u);
  assert(uV);
  return *uV;
}
VECTOR& getVector(VectorInterface& u) const
{
  VECTOR* uV = dynamic_cast<VECTOR*>(&u);
  assert(uV);
  return *uV;
}
public:
  Transfer<TRANSFER,VECTOR>(): TRANSFER(), TransferInterface() {}
  Transfer<TRANSFER,VECTOR>(const armaicvec& n, const armavec& dx): TRANSFER(n, dx), TransferInterface() {}

  void set_grid(const armaicvec& n, const armavec& dx) {get().set_grid(n, dx);}
  void restrict(VectorInterface& out, const VectorInterface& in) const {get().restrict(getVector(out),getVector(in));}
  void prolongate(VectorInterface& out, const VectorInterface& in) const{get().prolongate(getVector(out),getVector(in));}
};

#endif /* transferinterface_h */
