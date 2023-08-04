//
//  solverlaplace.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef  applicationinterface_hpp
#define  applicationinterface_hpp

#include  <map>
#include  <string>
#include  "boundary_conditions.hpp"
#include  "analyticalfunctioninterface.hpp"


/*-------------------------------------------------*/
class ApplicationInterface
{
protected:
    std::string _name;
    std::map<std::string,std::shared_ptr<AnalyticalFunctionInterface>> _sol, _rhs;
    std::shared_ptr<BoundaryConditions> _boundaryconditions;

public:
    virtual ~ApplicationInterface() {}
    ApplicationInterface(const ApplicationInterface& application){}      
    ApplicationInterface(std::string name, int dim, int ncomp=1) : _name(name) 
    {
        _boundaryconditions = std::make_shared<BoundaryConditions>(dim);            
    }
    bool has_solution() const {return _sol.size();}
    std::string toString() const {return _name;}
    int get_dim()const {return _boundaryconditions->size();}
    const std::shared_ptr<BoundaryConditions>& get_bc() const {return _boundaryconditions;}
    std::shared_ptr<BoundaryConditions>& get_bc() {return _boundaryconditions;}
    std::shared_ptr<AnalyticalFunctionInterface const> solution(std::string name) const {return _sol.at(name);}
    std::shared_ptr<AnalyticalFunctionInterface const> rhs(std::string name) const {return _rhs.at(name);}
};


#endif
