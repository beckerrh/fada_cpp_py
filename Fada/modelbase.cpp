//
//  q1.cpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#include  <armadillo>
#include  <cassert>
#include  "analyticalfunctioninterface.hpp"
#include  "modelbase.hpp"
#include  "uniformgrid.hpp"
#include  "gridvector.hpp"

/*-------------------------------------------------*/
ModelBase::ModelBase(std::string varname, const std::map <std::string, std::string>& parameters, std::shared_ptr<ApplicationInterface const> app) : _app(app)
{
    _varname = varname;
    _stenciltype  = "Trapez";
    _matrixtype   = "stencil";
    _smoothertype   = "matrix";
    _smoother     = "GS";
    _coarsesolver = "direct";
    _transfertype = "matrix";
    for (std::map <std::string, std::string>::const_iterator p = parameters.begin(); p != parameters.end(); p++)
    {
        if (p->first == "stenciltype")
        {
            _stenciltype = p->second;
        }
        else if (p->first == "matrixtype")
        {
            _matrixtype = p->second;
        }
        else if (parameters.find("smoother") != parameters.end())
        {
            _smoother = p->second;
        }
        else if (parameters.find("smoothertype") != parameters.end())
        {
            _smoothertype = p->second;
        }
        else if (parameters.find("coarsesolver") != parameters.end())
        {
            _coarsesolver = p->second;
        }
        else if (parameters.find("transfertype") != parameters.end())
        {
            _transfertype = p->second;
        }
    }
    if(_matrixtype!="stencil")
    {
        assert(_smoothertype!="stencil");
    }
}

/*-------------------------------------------------*/
void ModelBase::rhs(GridVector& v, std::shared_ptr<GridInterface const> grid, std::shared_ptr<AnalyticalFunctionInterface const> fct) const
{
    auto pc = std::dynamic_pointer_cast<ConstantFunction const>(fct);
    if(pc)
    {
        auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
        assert(ug);
        double vol = arma::prod(ug->dx());
        v.fill(vol*pc->get_constant());
        return;
    } 
    auto pr = std::dynamic_pointer_cast<RandomFunction const>(fct);
    if(pr)
    {
        auto ug = std::dynamic_pointer_cast<UniformGrid const>(grid);
        assert(ug);
        double vol = arma::prod(ug->dx());
        arma::arma_rng::set_seed_random();
        if(pr->get_type()=="uniform")
        {
            v.randu();                    
        }
        else if(pr->get_type()=="normal")
        {
            v.randn();                    
        }
        else
        {
            _not_written_();
        }
        v *= 100*vol;
        return;
    } 
    _not_written_();
}
