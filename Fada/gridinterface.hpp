
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
#include  "typedefs.hpp"

/*-------------------------------------------------*/
class GridInterface
{
public:
  virtual ~GridInterface() {}
  GridInterface() {}
  GridInterface(const GridInterface& updater) {}

  virtual std::string toString() const=0;
  virtual size_t dim() const=0;
  virtual const armaicvec& n() const=0;
  virtual size_t n(int i) const=0;
  virtual size_t nx() const=0;
  virtual size_t ny() const=0;
  virtual size_t nz() const=0;
  virtual int n_gridpoints() const=0;
  virtual void savehdf5(const std::string& filename) const=0;
  virtual void loadhdf5(const std::string& filename)=0;
};


#endif /* gridinterface_h */
