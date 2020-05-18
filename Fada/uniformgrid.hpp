//
//  uniformgrid.hpp
//  Fada
//
//  Created by Roland Becker on 15/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef uniformgrid_hpp
#define uniformgrid_hpp

#include  "typedefs.hpp"
#include  "gridinterface.hpp"

/*-------------------------------------------------*/
class UniformGrid : public GridInterface
{
protected:
  armaicvec _n;
  int _nall;
  armavec _dx;
  armamat _bounds;

public:
  UniformGrid() : GridInterface() {}
  UniformGrid(const UniformGrid& uniformgrid) : GridInterface(uniformgrid), _n(uniformgrid._n), _nall(uniformgrid._nall), _dx(uniformgrid._dx), _bounds(uniformgrid._bounds) {}
  UniformGrid(const armaicvec& n, const armamat* bounds=nullptr) {set_size(n, bounds);}

  void set_size(const armaicvec& n, const armamat* bounds=NULL);
  int dim() const {return _n.n_elem;}
  int nall()     const   { return _nall;}
  int nx() const {return _n(0);}
  int ny() const {return _n(1);}
  int nz() const {return _n(2);}
  int n(int i) const {return _n(i);}
  const armavec& dx() const {return _dx;}
  double dx(int i) const {return _dx[i];}
  double x(int ix, int iy=0, int iz=0) const {return ix*_dx(0);}
  double y(int ix, int iy, int iz=0) const {return iy*_dx(1);}
  double z(int ix, int iy, int iz) const {return iz*_dx(2);}
};


#endif /* uniformgrid_hpp */
