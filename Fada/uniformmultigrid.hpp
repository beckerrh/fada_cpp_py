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

public:
  UniformMultiGrid() : MultiGrid() {}
  UniformMultiGrid(const UniformMultiGrid& uniformmultigrid) : MultiGrid(uniformmultigrid)  {}
  UniformMultiGrid(int nlevels, const armaicvec& n0, std::shared_ptr<armamat> bounds=nullptr) : MultiGrid()  {set_size(nlevels, n0, bounds);}
  UniformMultiGrid& operator=(const UniformMultiGrid& umg);

  std::string toString() const;

  void set_size(int nlevels, const armaicvec& n0, std::shared_ptr<armamat> bounds=nullptr);

  size_t nx(size_t l) const {return _grids[l]->nx();}
  size_t ny(size_t l) const {return _grids[l]->ny();}
  size_t nz(size_t l) const {return _grids[l]->nz();}

};
std::ostream& operator<<(std::ostream& os, const UniformMultiGrid& mg);



#endif /* uniformmultigrid_hpp */
