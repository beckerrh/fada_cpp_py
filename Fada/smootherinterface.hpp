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
#include  <string>

class MatrixInterface;
class VectorInterface;
/*-------------------------------------------------*/
class SmootherInterface
{
public:
  virtual ~SmootherInterface() {}
  SmootherInterface() {}
  SmootherInterface(const SmootherInterface& smoother) {}

  virtual void presmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const=0;
  virtual void postsmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const=0;
};

/*-------------------------------------------------*/
template<typename SMOOTHER, class VECTOR>
class Smoother : public virtual SMOOTHER, public virtual SmootherInterface
{
protected:
  SMOOTHER& get() { return static_cast<SMOOTHER&>(*this); }
  SMOOTHER const& get() const { return static_cast<SMOOTHER const&>(*this); }
  const VECTOR& getVector(std::shared_ptr<VectorInterface const> u) const {return static_cast<const VECTOR&>(*u);}
  VECTOR& getVector(std::shared_ptr<VectorInterface> u) const{return static_cast<VECTOR&>(*u);}
public:
  Smoother<SMOOTHER, VECTOR>(std::shared_ptr<MatrixInterface const> matrix) : SMOOTHER(matrix), SmootherInterface() {}
  Smoother<SMOOTHER, VECTOR>(std::string type, std::shared_ptr<MatrixInterface const> matrix) : SMOOTHER(type, matrix), SmootherInterface() {}
  void presmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const {get().presmooth(getVector(out), getVector(in));}
  void postsmooth(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const {get().postsmooth(getVector(out), getVector(in));}
};


#endif /* smootherinterface_h */
