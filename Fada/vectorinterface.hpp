//
//  vectorinterface.hpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef vectorinterface_h
#define vectorinterface_h

class GridInterface;
class Vector;
/*-------------------------------------------------*/
class VectorInterface
{
public:
  virtual ~VectorInterface() {}
  VectorInterface() {}
  VectorInterface(const VectorInterface& vector) {}

};


#endif /* vectorinterface_h */
