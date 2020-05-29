
//
//  multigridinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef multigridinterface_h
#define multigridinterface_h

#include  <memory>
#include  "gridinterface.hpp"

/*-------------------------------------------------*/
class MultiGridInterface
{
public:
  virtual ~MultiGridInterface() {}
  MultiGridInterface() {}
  MultiGridInterface(const MultiGridInterface& multigrid) {}
  
  virtual std::string toString() const=0;
  virtual size_t nlevels() const=0;
  virtual size_t dim() const=0;
  virtual std::shared_ptr<GridInterface> get(int l) =0;
  virtual std::shared_ptr<const GridInterface> get(int l) const=0;
};


#endif /* multigridinterface_h */
