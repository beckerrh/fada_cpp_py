//
//  q1.hpp
//  Fada
//
//  Created by Roland Becker on 17/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef q1_hpp
#define q1_hpp

#include  "../feminterface.hpp"
#include  "../modelinterface.hpp"
#include  "nodevector.hpp"

class UniformGrid;

/*-------------------------------------------------*/
class Q1
{
protected:
  std::shared_ptr <UniformGrid const> _ug;
  std::string _stenciltype, _matrixtype, _smoothertype, _smoother, _coarsesolver;
  size_t _nx, _ny, _nz;
  double _vol;

public:
  ~Q1()
  {
  }

  // Q1() : _ug(nullptr), _stenciltype(), _matrixtype() {}
  Q1(const std::map <std::string, std::string>& parameters) : _ug(nullptr)
  {
    _stenciltype  = "Trapez";
    _matrixtype   = "stencil";
    _smoothertype   = "matrix";
    _smoother     = "GS";
    _coarsesolver = "direct";
    for (std::map <std::string, std::string>::const_iterator p = parameters.begin(); p != parameters.end(); p++)
    {
      if (p->first == "stenciltype")
      {
        _stenciltype = p->second;
      }
      else if (p->first == "matrixtype")
      {
        _matrixtype = p->second;
      }
      else if (parameters.find("smoother") != parameters.end())
      {
        _smoother = p->second;
      }
      else if (parameters.find("smoothertype") != parameters.end())
      {
        _smoothertype = p->second;
      }
      else if (parameters.find("coarsesolver") != parameters.end())
      {
        _coarsesolver = p->second;
      }
    }
    if(_matrixtype!="stencil")
    {
      assert(_smoothertype!="stencil");
    }
  }

  Q1(const Q1& model) : _ug(model._ug), _stenciltype(model._stenciltype)
  {
  }

  void set_grid(std::shared_ptr <GridInterface const> grid);
  std::shared_ptr <VectorInterface>   newVector(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <SmootherInterface const> newSmoother(std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface const> matrix) const;
  std::shared_ptr <CoarseSolverInterface const> newCoarseSolver(std::string type, std::shared_ptr <GridInterface const>grid, std::shared_ptr <MatrixInterface const> matrix) const;
  std::shared_ptr <MatrixInterface const>   newMatrix(std::shared_ptr <GridInterface const>grid) const;

  // virtual std::shared_ptr<MatrixInterface> newStencil(std::shared_ptr<GridInterface const>grid) const=0;
  virtual std::shared_ptr <FemAndMatrixAndSmootherInterface> newStencil(std::shared_ptr <GridInterface const>grid) const = 0;
  void rhs_one(NodeVector& v) const;
  void rhs_random(NodeVector& v) const;
};

/*-------------------------------------------------*/
class Q12d : public Q1
{
public:
  ~Q12d()
  {
  }

  // Q12d() : Q1() {}
  Q12d(const std::map <std::string, std::string>& parameters) : Q1(parameters)
  {
  }

  Q12d(const Q12d& model) : Q1(model)
  {
  }

  std::string toString() const
  {
    return("Q12d");
  }

  void boundary(NodeVector& v) const;

  // std::shared_ptr<MatrixInterface> newStencil(std::shared_ptr<GridInterface const>grid) const;
  std::shared_ptr <FemAndMatrixAndSmootherInterface> newStencil(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <TransferInterface const> newTransfer(std::shared_ptr <GridInterface const>grid) const;
};

/*-------------------------------------------------*/
class Q13d : public Q1
{
public:
  ~Q13d()
  {
  }

  // Q13d() : Q1() {}
  Q13d(const std::map <std::string, std::string>& parameters) : Q1(parameters)
  {
  }

  Q13d(const Q13d& model) : Q1(model)
  {
  }

  std::string toString() const
  {
    return("Q13d");
  }

  void boundary(NodeVector& v) const;

  // std::shared_ptr<MatrixInterface> newStencil(std::shared_ptr<GridInterface const> grid) const;
  std::shared_ptr <FemAndMatrixAndSmootherInterface> newStencil(std::shared_ptr <GridInterface const>grid) const;
  std::shared_ptr <TransferInterface const> newTransfer(std::shared_ptr <GridInterface const>grid) const;
};

#endif /* q1_hpp */
