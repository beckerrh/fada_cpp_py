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
#include  "array.hpp"
#include  "updaterinterface.hpp"
#include  "vector.hpp"


/*-------------------------------------------------*/
class UpdaterSimple : public UpdaterInterface
{
protected:
  Vector _mem;
  const Operator* _op;
  int _level;
  
public:
  void set_size(const armaicvec& n);
  void addUpdate(const Vector& w, Vector& u, Vector& r, bool print=false);
  void setParameters(int level, const Operator* op, int nvectors, const std::string& type="cyc", const std::string& solutiontype="gal");
};

/*-------------------------------------------------*/
class Updater : public UpdaterInterface
{
protected:
  mutable Array<Vector> _mem;
  const Operator* _op;
  std::string _type, _solutiontype, _status;
  bool _scale;
  mutable double _rnorm, _condition, _conditionmax, _conditionmean;
  int _nvectors, _nvars, _nshift, _level;
  int _nextupdate, _nextproduct;
  mutable int _nmemory, _nextmemory, _niterafterrestar;
  mutable armamat _H;
  mutable armavec _b, _x;
  Vector& getV(int i);
  Vector& getAV(int i);
  void _computeSmallSystem(int index, int nmemory, const Vector& r);
  void _matrixVectorProduct(int index) const;
  void restart();
  
public:
  ~Updater();
  Updater();
  Updater(const Updater& updater);
  Updater& operator=( const Updater& updater);
  void setParameters(int level, const Operator* op, int nvectors, const std::string& type="cyc", const std::string& solutiontype="gal");
  void set_size(const armaicvec& n);
  void addUpdate(const Vector& w, Vector& u, Vector& r, bool print=false);
};

#endif /* updater_h */
