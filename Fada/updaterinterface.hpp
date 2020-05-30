//
//  updaterinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef updaterinterface_h
#define updaterinterface_h

class FiniteElementInterface;
class GridInterface;
class MatrixInterface;
class VectorInterface;
/*-------------------------------------------------*/
class UpdaterInterface
{
public:
  virtual ~UpdaterInterface() {}
  UpdaterInterface() {}
  UpdaterInterface(const UpdaterInterface& updater) {}

  virtual void setParameters(const FiniteElementInterface& fem, const GridInterface& grid, std::shared_ptr<MatrixInterface> mat, int nvectors, const std::string& type="cyc", const std::string& solutiontype="gal")=0;
  virtual void addUpdate(const VectorInterface& w, VectorInterface& u, VectorInterface& r, bool print=false)=0;
};


#endif /* updaterinterface_h */
