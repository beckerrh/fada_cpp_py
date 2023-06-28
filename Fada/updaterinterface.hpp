//
//  updaterinterface.hpp
//  Fada
//
//  Created by Roland Becker on 08/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef updaterinterface_h
#define updaterinterface_h

class ModelInterface;
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

  virtual void setParameters(const ModelInterface& model, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> mat, int nvectors, const std::string& type="cyc", const std::string& solutiontype="ls")=0;
  virtual void addUpdate(std::shared_ptr<VectorInterface const> w, std::shared_ptr<VectorInterface> u, std::shared_ptr<VectorInterface> r, bool print=false)=0;
};


#endif /* updaterinterface_h */
