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

/*-------------------------------------------------*/
class BoundaryConditions : public std::vector<std::array<std::string,2>>
{
protected:

public:
  BoundaryConditions(const BoundaryConditions& bc) : std::vector<std::array<std::string,2>>(bc){}
  BoundaryConditions(int dim) : std::vector<std::array<std::string,2>>(dim)
  {
    for(int i=0;i<dim;i++)
    {
      (*this)[i][0] = "dir";
      (*this)[i][1] = "dir";
    }
  }
};

/*-------------------------------------------------*/
class BoundaryConditionsBool : public std::vector<std::array<bool,2>>
{
protected:

public:
  BoundaryConditionsBool(const BoundaryConditionsBool& bc) : std::vector<std::array<bool,2>>(bc){}
  BoundaryConditionsBool(std::shared_ptr<BoundaryConditions const> bc) : std::vector<std::array<bool,2>>(bc->size())
  {
    for(int i=0;i<bc->size();i++)
    {
        (*this)[i][0] = ((*bc)[i][0]=="dir");
        (*this)[i][1] = ((*bc)[i][1]=="dir");
    }
  }
};

#endif
