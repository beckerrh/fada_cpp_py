//
//  uniformmultigrid.hpp
//  Fada
//
//  Created by Roland Becker on 15/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef uniformmultigrid_hpp
#define uniformmultigrid_hpp

#include  "typedefs.hpp"
#include  "multigrid.hpp"

class UniformGrid;
/*-------------------------------------------------*/
class UniformMultiGrid : public MultiGrid
{
protected:
    int _ref_factor;

public:
  UniformMultiGrid() : MultiGrid() {}
  UniformMultiGrid(const UniformMultiGrid& uniformmultigrid) : MultiGrid(uniformmultigrid), _ref_factor(uniformmultigrid._ref_factor)  {}
  UniformMultiGrid(int nlevels, const armaicvec& n0, std::shared_ptr<armamat> bounds=nullptr, int ref_factor=2) : MultiGrid(), _ref_factor(ref_factor)  {set_size(nlevels, n0, bounds, ref_factor);}
  UniformMultiGrid& operator=(const UniformMultiGrid& umg);

  std::string toString() const;

  void set_size(int nlevels, const armaicvec& n0, std::shared_ptr<armamat> bounds=nullptr, int ref_factor=2);

  int get_ref_factor() const {return _ref_factor;}
  size_t nx(size_t l) const {return _grids[l]->nx();}
  size_t ny(size_t l) const {return _grids[l]->ny();}
  size_t nz(size_t l) const {return _grids[l]->nz();}

};
std::ostream& operator<<(std::ostream& os, const UniformMultiGrid& mg);



#endif /* uniformmultigrid_hpp */
