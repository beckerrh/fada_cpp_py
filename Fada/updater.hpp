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
class UpdaterSimple : public UpdaterInterface
{
protected:
  std::shared_ptr<VectorInterface> _mem;
  std::shared_ptr<MatrixInterface> _mat;
//  int _level;

public:
  void addUpdate(const VectorInterface& w, VectorInterface& u, VectorInterface& r, bool print=false);
  void setParameters(const FiniteElementInterface& fem, const GridInterface& grid, std::shared_ptr<MatrixInterface> mat, int nvectors, const std::string& type="cyc", const std::string& solutiontype="ls");
};

/*-------------------------------------------------*/
class Updater : public UpdaterInterface
{
protected:
  mutable std::vector<std::shared_ptr<VectorInterface>> _mem;
  std::shared_ptr<MatrixInterface> _mat;
  std::string _type, _solutiontype, _status;
  bool _scale;
  mutable double _rnorm, _condition, _conditionmax, _conditionmean;
  int _nvectors, _nvars, _nshift;
  int _nextupdate, _nextproduct;
  mutable int _nmemory, _nextmemory, _niterafterrestar;
  mutable armamat _H;
  mutable armavec _b, _x;
  VectorInterface& getV(int i);
  VectorInterface& getAV(int i);
  void _computeSmallSystem(int index, int nmemory, const VectorInterface& r);
  void _matrixVectorProduct(int index) const;
  void restart();

public:
  ~Updater();
  Updater();
  Updater(const Updater& updater);
  Updater& operator=( const Updater& updater);
  void setParameters(const FiniteElementInterface& fem, const GridInterface& grid, std::shared_ptr<MatrixInterface> mat, int nvectors, const std::string& type="cyc", const std::string& solutiontype="ls");
  void addUpdate(const VectorInterface& w, VectorInterface& u, VectorInterface& r, bool print=false);
};

#endif /* updater_h */
