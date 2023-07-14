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
  armavec _dx;
  armamat _bounds;

public:
  UniformGrid() : GridInterface() {}
  UniformGrid(const UniformGrid& uniformgrid) : GridInterface(uniformgrid), _n(uniformgrid._n), _dx(uniformgrid._dx), _bounds(uniformgrid._bounds) {}
  UniformGrid(const armaicvec& n, std::shared_ptr<armamat> bounds=nullptr) {set_size(n, bounds);}

  std::string toString() const;

  void set_size(const armaicvec& n, std::shared_ptr<armamat> bounds=nullptr);
  size_t dim() const {return _n.n_elem;}
  int n_fine() const { return arma::prod(_n);}
  size_t n(int i) const {return _n[i];}
  size_t nx() const {return _n[0];}
  size_t ny() const {return _n[1];}
  size_t nz() const {return _n[2];}
  const armaicvec& n() const {return _n;}
  const armavec& dx() const {return _dx;}
  double dx(int i) const {return _dx[i];}
  double x(int ix, int iy=0, int iz=0) const {return ix*_dx[0];}
  double y(int ix, int iy, int iz=0) const {return iy*_dx[1];}
  double z(int ix, int iy, int iz) const {return iz*_dx[2];}
  void savehdf5(const std::string& filename) const;
  void loadhdf5(const std::string& filename);
};


#endif /* uniformgrid_hpp */
