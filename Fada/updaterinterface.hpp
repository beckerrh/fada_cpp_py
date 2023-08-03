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
protected:
    std::string _name;
public:
    virtual ~UpdaterInterface() {}
    UpdaterInterface(std::string name) : _name(name) {}
    UpdaterInterface(const UpdaterInterface& updater) : _name(updater._name) {}

    std::string toString() const {return _name;}
    virtual void setParameters(const ModelInterface& model, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> mat)=0;
    virtual void addUpdate(std::shared_ptr<VectorInterface const> w, std::shared_ptr<VectorInterface> u, std::shared_ptr<VectorInterface> r, bool print=false)=0;
};


#endif /* updaterinterface_h */
