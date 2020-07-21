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
  std::shared_ptr<MatrixInterface> _matrix;
  
public:
  SmootherSimple(std::string type) : _type(type) {}
  SmootherSimple(const SmootherSimple& smoother) {}

  void set_matrix(std::shared_ptr<MatrixInterface> matrix);
  void solve(VectorInterface& out, const VectorInterface& in) const{pre(out,in);}
  void pre(VectorInterface& out, const VectorInterface& in) const;
  void post(VectorInterface& out, const VectorInterface& in) const;
};

#endif /* smoothersimple_hpp */
