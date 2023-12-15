//
//  stencil2d.hpp
//  Fada
//
//  Created by Roland Becker on 01/06/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef boundaryconditions_hpp
#define boundaryconditions_hpp

#include  <array>
#include  <vector>
#include  <string>
#include  <ostream>
#include  <memory>
#include  "analyticalfunctioninterface.hpp"

/*-------------------------------------------------*/
class BoundaryConditions : public std::vector<std::array<std::string,2>>
{
protected:
    std::vector<std::vector<FunctionMap>> _bdry_fct;
    
public:
  BoundaryConditions(const BoundaryConditions& bc) : std::vector<std::array<std::string,2>>(bc), _bdry_fct(bc._bdry_fct){}
  BoundaryConditions(int dim, std::string cond="dir") : std::vector<std::array<std::string,2>>(dim), _bdry_fct(dim)
  {
    for(int i=0;i<dim;i++)
    {
      (*this)[i][0] = cond;
      (*this)[i][1] = cond;
    }
    for(int i=0;i<dim;i++)
    {
        _bdry_fct[i].resize(2);
    }
  }
  FunctionMap& get_bf(int i, int j) {return _bdry_fct[i][j];}
  const std::vector<std::vector<FunctionMap>>& get_bf() const {return _bdry_fct;}
  bool all(std::string cond) const
  {
      for(auto p: *this)
      {
          if(p[0]!=cond) return false;
          if(p[1]!=cond) return false;
      }
      return true;
  }
};
static std::ostream& operator<<(std::ostream& os, const BoundaryConditions& v)
{
    for(auto p:v)
    {
        os << p[0] << " " << p[1] << std::endl;        
    }
    for(auto p:v.get_bf())
    {
        os << p[0] << " " << p[1] << std::endl;        
    }
    return os;
}

/*-------------------------------------------------*/
class BoundaryConditionsBool : public std::vector<std::array<bool,2>>
{
protected:

public:
    BoundaryConditionsBool(const BoundaryConditionsBool& bc) : std::vector<std::array<bool,2>>(bc){}
    BoundaryConditionsBool(int dim)
    {
        set_false(dim);
    }
    BoundaryConditionsBool(std::shared_ptr<BoundaryConditions const> bc)
    {
        assert(bc);
        std::vector<std::array<bool,2>>::resize(bc->size());
        for(int i=0;i<bc->size();i++)
        {
            (*this)[i][0] = ((*bc)[i][0]=="dir");
            (*this)[i][1] = ((*bc)[i][1]=="dir");
        }          
    }
    void set_false(int dim)
    {
        std::vector<std::array<bool,2>>::resize(dim);
        for(int i=0;i<dim;i++)
        {
            (*this)[i][0] = false;
            (*this)[i][1] = false;
        }              
    }
};
static std::ostream& operator<<(std::ostream& os, const BoundaryConditionsBool& v)
{
    for(auto p:v)
    {
        os << p[0] << " " << p[1] << std::endl;        
    }
    return os;
}

#endif
