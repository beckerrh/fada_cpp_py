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
#include  <string>

/*-------------------------------------------------*/
class SmootherSimple : public SmootherInterface
{
protected:
  std::string _type;
  std::shared_ptr<MatrixInterface const> _matrix;
  void set_matrix(std::shared_ptr<MatrixInterface const> matrix);

public:
  SmootherSimple(std::string type) : SmootherInterface(), _type(type) {}
  SmootherSimple(const SmootherSimple& smoother) : SmootherInterface() {}
  SmootherSimple(std::string type, std::shared_ptr<MatrixInterface const> matrix) : SmootherInterface(), _type(type) {set_matrix(matrix);}

  void solve(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const {pre(out,in);}
  void pre(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const;
  void post(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const;
};

#endif /* smoothersimple_hpp */
