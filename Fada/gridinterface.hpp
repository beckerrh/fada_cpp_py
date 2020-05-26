
//
//  gridinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef gridinterface_h
#define gridinterface_h

#include  <memory>

/*-------------------------------------------------*/
class GridInterface
{
public:
  virtual ~GridInterface() {}
  GridInterface() {}
  GridInterface(const GridInterface& updater) {}
};


#endif /* gridinterface_h */
