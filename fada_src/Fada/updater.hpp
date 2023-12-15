//
//  updater.h
//  Fada
//
//  Created by Roland Becker on 07/05/2020.
//  Copyright © 2020 Roland Becker. All rights reserved.
//

#ifndef updater_h
#define updater_h

#include  "typedefs.hpp"
#include  "updaterinterface.hpp"
#include  "vectorinterface.hpp"


/*-------------------------------------------------*/
class UpdaterConstant : public UpdaterInterface
{
protected:
  double _omega;
  std::shared_ptr<VectorInterface> _mem;
  std::shared_ptr<MatrixInterface const> _mat;

public:
  UpdaterConstant(std::string name, double omega=1, bool print=false) : UpdaterInterface(name, print) {_omega=omega;}
  void addUpdate(std::shared_ptr<VectorInterface const> w, std::shared_ptr<VectorInterface> u, std::shared_ptr<VectorInterface> r);
  void setParameters(const ModelInterface& model, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> mat);
};

/*-------------------------------------------------*/
class UpdaterSimple : public UpdaterConstant
{
protected:
    std::string _type;
public:
    UpdaterSimple(std::string name, std::string type, bool print=false) : UpdaterConstant(name, print), _type(type) {}
  void addUpdate(std::shared_ptr<VectorInterface const> w, std::shared_ptr<VectorInterface> u, std::shared_ptr<VectorInterface> r);
};

/*-------------------------------------------------*/
class Updater : public UpdaterInterface
{
protected:
  mutable std::vector<std::shared_ptr<VectorInterface>> _mem;
  std::shared_ptr<MatrixInterface const> _mat;
  std::string _droptype, _solutiontype, _status;
  bool _scale;
  mutable double _rnorm, _condition, _conditionmax, _conditionmean;
  int _nvectors, _nvars, _nshift;
  int _nextupdate, _nextproduct;
  mutable int _nmemory, _nextmemory, _niterafterrestar;
  mutable armamat _H;
  mutable armavec _b, _x;
  std::shared_ptr<VectorInterface> getV(int i);
  std::shared_ptr<VectorInterface> getAV(int i);
  void _computeSmallSystem(int index, int nmemory, std::shared_ptr<VectorInterface const> r);
  void _matrixVectorProduct(int index) const;
  void restart();

public:
  Updater(std::string name, const std::string& droptype, const std::string& solutiontype, int nvectors, bool print=false);
  Updater& operator=( const Updater& updater);
  void setParameters(const ModelInterface& model, std::shared_ptr<GridInterface const> grid, std::shared_ptr<MatrixInterface const> mat);
  void addUpdate(std::shared_ptr<VectorInterface const> w, std::shared_ptr<VectorInterface> u, std::shared_ptr<VectorInterface> r);
};

#endif /* updater_h */
