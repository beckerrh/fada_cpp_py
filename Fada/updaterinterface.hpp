//
//  updaterinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef updaterinterface_h
#define updaterinterface_h

class Operator;
class Vector;
/*-------------------------------------------------*/
class UpdaterInterface
{
public:
  virtual ~UpdaterInterface() {}
  UpdaterInterface() {}
  UpdaterInterface(const UpdaterInterface& updater) {}

  virtual void setParameters(int level, const Operator* op, int nvectors, const std::string& type="cyc", const std::string& solutiontype="gal")=0;  
  virtual void set_size(const armaicvec& n)=0;
  virtual void addUpdate(const Vector& w, Vector& u, Vector& r, bool print=false)=0;
};


#endif /* updaterinterface_h */
