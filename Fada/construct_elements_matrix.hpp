//
//  sparsematrix.hpp
//  Fada
//
//  Created by Roland Becker on 20/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef construct_elements_matrix_hpp
#define construct_elements_matrix_hpp

#include  "typedefs.hpp"

/*-------------------------------------------------*/
class Construct_Elements_Matrix
{
protected:
  int _count;
  arma::umat& _locations;
  armavec& _values;
  void start(int size)
  {
    _locations.resize(2, size);
    _values.resize(size);
    _count=0;
  }
public:
  Construct_Elements_Matrix(arma::umat& locations, armavec& values, int size) : _locations(locations), _values(values) {start(size);}
  const arma::umat& locations() const {return _locations;}
  const armavec& values() const {return _values;}
  void add(int i, int j, double value)
  {
    _locations.at(0, _count) = i;
    _locations.at(1, _count) = j;
    _values[_count] = value;
    _count++;
  }
};
#endif
