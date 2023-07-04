//
//  smoothersimple.hpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef smoothersimple_hpp
#define smoothersimple_hpp

#include  "smootherinterface.hpp"
#include  "typedefs.hpp"
// #include  "sparsematrix.hpp"
#include  <string>

/*-------------------------------------------------*/
template<typename T>
class SmootherSimple
{
protected:
  std::string _type;
  // std::shared_ptr<MatrixInterface const> _matrix;
  // std::shared_ptr<SparseMatrix const> _matrix;
  std::shared_ptr<T const> _matrix;
  void set_matrix(std::shared_ptr<MatrixInterface const> matrix);

public:
  SmootherSimple<T>(std::string type) : _type(type) {}
  SmootherSimple<T>(const SmootherSimple<T>& smoother): _type(), _matrix(smoother._matrix){}
  SmootherSimple<T>(std::string type, std::shared_ptr<MatrixInterface const> matrix) : _type(type) {set_matrix(matrix);}

  // void solve(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const {presmooth(out,in);}
  // void presmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const;
  // void postsmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const;
  void presmooth(armavec& out, const armavec& in) const;
  void postsmooth(armavec& out, const armavec& in) const;
};

#endif /* smoothersimple_hpp */
