//
//  multigrid.hpp
//  Fada
//
//  Created by Roland Becker on 28/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef multigrid_hpp
#define multigrid_hpp

#include  <vector>
#include  <memory>
#include  "multigridinterface.hpp"

/*-------------------------------------------------*/
class MultiGrid : public MultiGridInterface
{
protected:
  std::vector<std::shared_ptr<GridInterface> > _grids;
  
public:
  size_t nlevels() const {return _grids.size();}
  std::shared_ptr<GridInterface> get(int l) {return _grids[l];}
  std::shared_ptr<const GridInterface> get(int l) const {return _grids[l];}
  size_t dim() const {return _grids[0]->dim();}
  int nall() const { return _grids[0]->nall();}
};

#endif /* multigrid_hpp */
