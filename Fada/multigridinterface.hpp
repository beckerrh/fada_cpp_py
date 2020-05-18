
//
//  multigridinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef multigridinterface_h
#define multigridinterface_h

/*-------------------------------------------------*/
class MultiGridInterface
{
public:
  virtual ~MultiGridInterface() {}
  MultiGridInterface() {}
  MultiGridInterface(const MultiGridInterface& updater) {}
};


#endif /* multigridinterface_h */
