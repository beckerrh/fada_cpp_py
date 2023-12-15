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
#include  "parametermap.hpp"

class AnalyticalFunctionInterface;
class GridInterface;
class GridVector;
class UniformGrid;

/*-------------------------------------------------*/
class ModelBase
{
protected:
    std::string _varname;
    ParameterMap _parameters;
    // std::string _stenciltype, _matrixtype, _smoothertype, _smoother, _coarsesolver, _transfertype;
    std::shared_ptr <ApplicationInterface const> _app;

public:
    ~ModelBase()
    {
    }

    // ModelBase(std::string varname, const ParameterMap& parameters, std::shared_ptr <ApplicationInterface const> app = nullptr);
    ModelBase(std::string varname, const ParameterMap& parameters, std::shared_ptr <ApplicationInterface const> app = nullptr);
    ModelBase(const ModelBase& model) : _varname(model._varname), _parameters(model._parameters), _app(model._app)
    {
    }
    void rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const;
};

#endif
