
//
//  transferinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef transferinterface_h
#define transferinterface_h

class GridInterface;
class Vector;
/*-------------------------------------------------*/
class TransferInterface
{
public:
  virtual ~TransferInterface() {}
  TransferInterface() {}
  TransferInterface(const TransferInterface& updater) {}

  virtual void set_grid(const GridInterface& grid)=0;
  virtual void restrict(Vector& out, const Vector& in) const=0;
  virtual void prolongate(Vector& out, const Vector& in) const=0;
};


#endif /* transferinterface_h */
