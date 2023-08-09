//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef modelbase_hpp
#define modelbase_hpp

#include  <map>
#include  "applicationinterface.hpp"

class AnalyticalFunctionInterface;
class GridInterface;
class GridVector;
class UniformGrid;

/*-------------------------------------------------*/
class ModelBase
{
protected:
    std::string _varname;
    std::string _stenciltype, _matrixtype, _smoothertype, _smoother, _coarsesolver, _transfertype;
    std::shared_ptr <ApplicationInterface const> _app;

public:
    ~ModelBase()
    {
    }

    ModelBase(std::string varname, const std::map <std::string, std::string>& parameters, std::shared_ptr <ApplicationInterface const> app = nullptr);
    ModelBase(const ModelBase& model) : _varname(model._varname), _stenciltype(model._stenciltype), _app(model._app)
    {
    }
    void rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const;
};

#endif
