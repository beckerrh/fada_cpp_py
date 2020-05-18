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
#include  "multigridinterface.hpp"

class UniformGrid;
/*-------------------------------------------------*/
class UniformMultiGrid : public MultiGridInterface
{
protected:
  armaimat _n;
  armairvec _nall;
  armamat _dx;
  armamat _bounds;

public:
  UniformMultiGrid() : MultiGridInterface() {}
  UniformMultiGrid(const UniformMultiGrid& uniformmultigrid) : MultiGridInterface(), _n(uniformmultigrid._n), _nall(uniformmultigrid._nall), _dx(uniformmultigrid._dx), _bounds(uniformmultigrid._bounds)  {}
  
  void set_size(int nlevels, const armaicvec& n0, const armamat* bounds=nullptr);
  int dim() const {return _n.n_rows;}
  int minlevel() const   { return 0;}
  int maxlevel() const   { return nlevels()-1;}
  int nlevels()  const   { return _n.n_cols;}
  arma::subview_col<int>  n(int l) const { return _n.col(l);}
  arma::subview_col<int>  n() const { return _n.col(nlevels()-1);}
  int nall()     const   { return _nall(nlevels()-1);}
  int nx(int l) const {return _n(0,l);}
  int ny(int l) const {return _n(1,l);}
  int nz(int l) const {return _n(2,l);}
  int nx() const {return _n(0,_n.n_cols-1);}
  int ny() const {return _n(1,_n.n_cols-1);}
  int nz() const {return _n(2,_n.n_cols-1);}
  double x(int ix, int iy=0, int iz=0) const {return ix*_dx(0,_n.n_cols-1);}
  double y(int ix, int iy, int iz=0) const {return iy*_dx(1,_n.n_cols-1);}
  double z(int ix, int iy, int iz) const {return iz*_dx(2,_n.n_cols-1);}
  int nmax(int i) const {return _n(i,_n.n_cols-1);}
  UniformGrid grid(int l=-1) const;
};



#endif /* uniformmultigrid_hpp */
