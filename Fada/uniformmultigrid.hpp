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
  void _set_size(std::shared_ptr<armamat> bounds);

public:
  UniformMultiGrid() : MultiGridInterface() {}
  UniformMultiGrid(const UniformMultiGrid& uniformmultigrid) : MultiGridInterface(), _n(uniformmultigrid._n), _nall(uniformmultigrid._nall), _dx(uniformmultigrid._dx), _bounds(uniformmultigrid._bounds)  {}
  UniformMultiGrid(const armaimat&  n,  std::shared_ptr<armamat> bounds=nullptr) : MultiGridInterface()
  {
    set_size(n, bounds);
  }
  UniformMultiGrid& operator=(const UniformMultiGrid& umg);

  std::string toString() const;

  void set_size(int nlevelmax, int nlevel, const armaicvec& n0, std::shared_ptr<armamat> bounds=nullptr);
  void set_size(const armaimat&  n, std::shared_ptr<armamat> bounds=nullptr);
  int dim() const {return (int) _n.n_rows;}

  int nlevels()  const   { return (int) _n.n_cols;}
  const armaimat&  n() const { return _n;}
  const armamat&  bounds() const { return _bounds;}
  const armamat&  dx() const { return _dx;}
  armaimat&  n()  { return _n;}
  armamat&  bounds()  { return _bounds;}
  armamat&  dx()  { return _dx;}
  arma::subview_col<int>  n(int l) const { return _n.col(l);}
  int nall()     const   { return _nall(0);}
  int nx(int l) const {return _n(0,l);}
  int ny(int l) const {return _n(1,l);}
  int nz(int l) const {return _n(2,l);}
  double x(int ix, int iy=0, int iz=0) const {return ix*_dx(0,0);}
  double y(int ix, int iy, int iz=0) const {return iy*_dx(1,0);}
  double z(int ix, int iy, int iz) const {return iz*_dx(2,0);}
  int nmax(int i) const {return _n(i,0);}
  std::shared_ptr<UniformGrid> grid(int l=0) const;
};
std::ostream& operator<<(std::ostream& os, const UniformMultiGrid& mg);



#endif /* uniformmultigrid_hpp */
