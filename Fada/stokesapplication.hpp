//
//  solverlaplace.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef  stokesapplication_hpp
#define  stokesapplication_hpp

#include  <sstream>
#include  <string>
#include  "applicationinterface.hpp"


/*-------------------------------------------------*/
class Sinus2DSolutionP: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Sinus2DRhsP: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Sinus2DSolutionV0: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Sinus2DRhsV0: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Sinus2DSolutionV1: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};
class Sinus2DRhsV1: public AnalyticalFunctionInterface
{
public:
    double operator()(double x, double y) const;
};

/*-------------------------------------------------*/
class StokesApplication : public ApplicationInterface
{
protected:
    // std::string _name;
    // std::shared_ptr<SystemAnalyticalFunction> _sol, _rhs;
    // std::shared_ptr<BoundaryConditions> _boundaryconditions;

    void set_application(int dim, std::string name);

public:
    StokesApplication(std::string name, int dim) : ApplicationInterface(name, dim, dim+1) 
    {
        set_application(dim, name);
    }
    // StokesApplicationInterface(const StokesApplication& application) : _name(application._name), _sol(application._sol), _Q1shiftedapplication._rhs){}
    // std::string toString() const;
    // virtual int get_dim()const {return _boundaryconditions->size();}
    // std::shared_ptr<BoundaryConditions const> get_bc() const {return _boundaryconditions;}
    // std::shared_ptr<SystemAnalyticalFunction const> solution() const
    // {
    //     return _sol;
    // }
    // std::shared_ptr<SystemAnalyticalFunction const> Q1shifted) const
    // {
    //     return _rhs;
    // }
    // std::shared_ptr<AnalyticalFunctionInterface const> solution_p() const
    // {
    //     return _sol.at("p");
    // }
    // std::shared_ptr<AnalyticalFunctionInterface const> solution_v(int i) const
    // {
    //     std::stringstream ss;
    //     ss << "v" << _idir;
    //     return _sol->at(ss.str());
    // }
    // std::shared_ptr<AnalyticalFunctionInterface const> rhs_p() const
    // {
    //     return _rhs->back();
    // }
    // std::shared_ptr<AnalyticalFunctionInterface const> rhs_v(int i) const
    // {
    //     assert(_rhs->get(i));
    //     return _rhs->get(i);
    // }
};


#endif
