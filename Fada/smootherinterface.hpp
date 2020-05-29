//
//  smootherinterface.hpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef smootherinterface_h
#define smootherinterface_h

#include  <memory>

class MatrixInterface;
class Vector;
/*-------------------------------------------------*/
class SmootherInterface
{
public:
  virtual ~SmootherInterface() {}
  SmootherInterface() {}
  SmootherInterface(const SmootherInterface& smoother) {}

  virtual void set_matrix(std::shared_ptr<MatrixInterface> matrix) = 0;
  virtual void solve(Vector& out, const Vector& in) const = 0;
  virtual void pre(Vector& out, const Vector& in) const {solve(out, in);}
  virtual void post(Vector& out, const Vector& in) const {solve(out, in);}
};


#endif /* smootherinterface_h */
