//
//  smootherinterface.hpp
//  Fada
//
//  Created by Roland Becker on 27/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef coarsesolverinterface_h
#define coarsesolverinterface_h

#include  <string>
#include  <memory>

class MatrixInterface;
class VectorInterface;
/*-------------------------------------------------*/
class CoarseSolverInterface
{
public:
  virtual ~CoarseSolverInterface() {}
  CoarseSolverInterface() {}
  CoarseSolverInterface(const CoarseSolverInterface& smoother) {}

  virtual void solve(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const = 0;
  virtual void update(std::shared_ptr<MatrixInterface const> matrix)=0;
};

/*-------------------------------------------------*/
template<typename COARSESOLVER, class VECTOR>
class CoarseSolver : public COARSESOLVER, public CoarseSolverInterface
{
protected:
  COARSESOLVER& get() { return static_cast<COARSESOLVER&>(*this); }
  COARSESOLVER const& get() const { return static_cast<COARSESOLVER const&>(*this); }
  const VECTOR& getVector(std::shared_ptr<VectorInterface const> u) const {return static_cast<const VECTOR&>(*u);}
  VECTOR& getVector(std::shared_ptr<VectorInterface> u) const{return static_cast<VECTOR&>(*u);}

public:
  // CoarseSolver<COARSESOLVER, VECTOR>(std::shared_ptr<MatrixInterface const> matrix) : COARSESOLVER(matrix), CoarseSolverInterface() {}
  CoarseSolver<COARSESOLVER, VECTOR>(std::shared_ptr<MatrixInterface const> matrix, std::string type="") : COARSESOLVER(matrix, type), CoarseSolverInterface() {}
  void solve(std::shared_ptr<VectorInterface> out, std::shared_ptr<VectorInterface const> in) const {get().solve(getVector(out), getVector(in));}
  void update(std::shared_ptr<MatrixInterface const> matrix) {get().update(matrix);}
};


#endif /* smootherinterface_h */
