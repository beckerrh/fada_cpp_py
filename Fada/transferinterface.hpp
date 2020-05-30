
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
class Vector;
/*-------------------------------------------------*/
class TransferInterface
{
public:
  virtual ~TransferInterface() {}
  TransferInterface() {}
  TransferInterface(const TransferInterface& transfer) {}

  virtual void set_grid(const armaicvec& n, const armavec& dx)=0;
  virtual void restrict(Vector& out, const Vector& in) const=0;
  virtual void prolongate(Vector& out, const Vector& in) const=0;
};


#endif /* transferinterface_h */
