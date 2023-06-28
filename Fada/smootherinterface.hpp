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
class VectorInterface;
/*-------------------------------------------------*/
class SmootherInterface
{
public:
  virtual ~SmootherInterface() {}
  SmootherInterface() {}
  SmootherInterface(const SmootherInterface& smoother) {}

  // virtual void set_matrix(std::shared_ptr<MatrixInterface const> matrix) = 0;
  virtual void solve(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const = 0;
  virtual void pre(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const {solve(out, in);}
  virtual void post(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const {solve(out, in);}
};

/*-------------------------------------------------*/
// template<typename SMOOTHER, typename MATRIX, class VECTOR>
template<typename SMOOTHER, class VECTOR>
class Smoother : public SMOOTHER, public SmootherInterface
{
protected:
  SMOOTHER& get() { return static_cast<SMOOTHER&>(*this); }
  SMOOTHER const& get() const { return static_cast<SMOOTHER const&>(*this); }
  // MATRIX const& getMatrix(std::shared_ptr<MatrixInterface> A) const { return static_cast<MATRIX const&>(*A); }
  const VECTOR& getVector(std::shared_ptr<VectorInterface const> u) const {return static_cast<const VECTOR&>(*u);}
  VECTOR& getVector(std::shared_ptr<VectorInterface> u) const{return static_cast<VECTOR&>(*u);}
public:
  // Smoother<SMOOTHER,MATRIX, VECTOR>() : SMOOTHER(), SmootherInterface() {}
  Smoother<SMOOTHER, VECTOR>(std::shared_ptr<MatrixInterface const> matrix) : SMOOTHER(matrix), SmootherInterface() {}
  void solve(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const {get().solve(getVector(out), getVector(in));}

  // void set_matrix(std::shared_ptr<MatrixInterface const> A) {get().set_matrix(getMatrix(A));}
};


#endif /* smootherinterface_h */
